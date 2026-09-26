function [genomes, phenotypes, diagnostics] = sampleGenotypesConditional(genomeTheta, simParams, targetPhenotypes, varargin)
% sampleGenotypesConditional - Pick a genotype at random from all the genotypes
%   that have the required starting fitness AND the required trait ratio,
%   instead of building one genotype greedily.
%
% WHY THIS EXISTS
%   initializeGenomeTheta builds the initial genotype by a single greedy
%   forward pass: it walks the loci in order and sets locus ell to allele 1
%   whenever that reduces the distance to the target phenotype. The genotype
%   it returns therefore reaches the target with the FEWEST possible 1s, each of
%   them nearly aligned with the target direction.
%
%   That is not a typical genotype at that phenotype. The diagnosis is that
%   the greedy genotype is atypical in exactly the way that matters: it starts
%   the population with an unrepresentative set of available mutations. Measured
%   on the published pleiotropic GPFM at the six paper initial conditions, the
%   greedy initializer produces 39-52 1s where a genotype drawn
%   properly from the conditional distribution carries 244-288. The starting
%   genotype is a factor of ~5 sparser than it should be, so almost every locus
%   is available to flip 0 -> 1, and the mutational neighbourhood is dominated
%   by directions that were never selected for.
%
%   This routine replaces that construction with a sample from
%
%       P(g | |log W(g) - log W_0| <= tolF, |angle(x(g)) - angle(x_0)| <= tolA)
%
%   spread evenly over all genotypes that satisfy both.
%
% ALGORITHM
%   Phase 1 - the random-order walk. Start from the all-zero genome and flip loci 0 -> 1 in
%       random order until the fitness first drops to W_0. A walk that jumps
%       straight past W_0 in one step is thrown away and restarted; without that,
%       the sample is the first crossing, which is biased toward genotypes
%       with too few 1s. Repeated independently, this also gives a
%       stationarity diagnostic: the running mean of the number of 1s.
%
%   Phase 2 - getting the trait ratio right. The phase-1 genotype has the right
%       fitness but whatever trait ratio it happened to land on. A short annealing
%       run walks it to a genotype that also has the required ratio, and a
%       Metropolis-Hastings chain then samples evenly among all genotypes that
%       meet both requirements. Both proposals are symmetric - flip one random
%       locus, or swap a 1 with a 0 - and a proposal
%       is accepted only if it still meets both requirements, which is what makes
%       the sample even rather than biased.
%
%       Requiring both AT ONCE is the point. Generating a lot of genotypes at the
%       right fitness, sorting them by trait ratio and keeping the closest matches
%       does not work: almost none of them have the ratios at the extremes of our
%       six initial conditions, so those conditions end up represented by a couple
%       of near-duplicate genotypes.
%
% Inputs:
%   genomeTheta      - [1 x L] pleiotropic angles, as produced by initializeGenomeTheta
%   simParams        - needs .deltaTrait, .ellipseParams, .landscapeStdDev
%   targetPhenotypes - [nTargets x 2] target phenotypes (the W_0 contour points)
%
% Name-value options:
%   'tolLogFitness' - how close to W_0 the fitness must be, on log W. (default 0.01 * |log W_0|)
%   'tolLogRatio'   - half-width of the ratio window, on log(x2/x1).    (default 0.05)
%   'burnIn'        - MH steps discarded before the sample is taken.    (default 200 * L)
%   'numSteps'      - MH steps after burn-in.                           (default 200 * L)
%   'annealSteps'   - annealing steps used to find a genotype meeting both.  (default 400 * L)
%   'walkSamples'   - independent phase-1 walks, for the diagnostic.    (default 200)
%   'seed'          - RNG seed.                                         (default [], unset)
%   'verbose'       - print a per-target summary.                       (default true)
%
% Outputs:
%   genomes     - [nTargets x L] sampled allele states
%   phenotypes  - [nTargets x 2] phenotypes they encode
%   diagnostics - [1 x nTargets] struct with fields
%       .numOnes             number of 1s in the returned genotype
%       .numOnesGreedy       number the greedy initializer would have produced
%       .walkNumOnes         [1 x walkSamples] phase-1 landing counts
%       .walkRunningMean     running mean of walkNumOnes (stationarity check)
%       .walkStationary      true if the running mean settles within 2% by the
%                            halfway point of the sample
%       .chainNumOnes        number of 1s along the MH chain (thinned)
%       .acceptRate          MH acceptance rate
%       .feasible            true if a genotype meeting both requirements was found
%       .logFitness, .logRatio, .logFitnessTarget, .logRatioTarget
%
% Reference:
%   Kim, M., Ardell, S. M., & Kryazhimskiy, S. (2025).
%   "Module-Selection Balance in the Evolution of Modular Organisms."
%
% See also: initializeGenomeThetaSampled, initializeGenomeTheta
%
% Copyright (c) 2025 Minkyu Kim, Cornell University
% Licensed under MIT License

    ip = inputParser;
    ip.addParameter('tolLogFitness', []);
    ip.addParameter('tolLogRatio', 0.05);
    ip.addParameter('burnIn', []);
    ip.addParameter('numSteps', []);
    ip.addParameter('annealSteps', []);
    ip.addParameter('walkSamples', 200);
    ip.addParameter('seed', []);
    ip.addParameter('verbose', true);
    ip.parse(varargin{:});
    opt = ip.Results;

    if ~isempty(opt.seed), rng(opt.seed); end

    L      = numel(genomeTheta);
    delta  = simParams.deltaTrait;
    a      = simParams.ellipseParams;
    sigW   = simParams.landscapeStdDev;

    % Allele 1 at locus ell displaces the phenotype by -delta*(cos, sin).
    % This is the sign convention of initializeGenomeTheta and of Eq. (modular
    % GPM); Eq. (pleiotropic GPM) in the manuscript writes +delta.
    dxc = -delta * cos(genomeTheta(:))';
    dxs = -delta * sin(genomeTheta(:))';

    logWof = @(x) -((x(1)/a(1))^2 + (x(2)/a(2))^2) / (2 * sigW^2);

    nTargets   = size(targetPhenotypes, 1);
    genomes    = zeros(nTargets, L);
    phenotypes = zeros(nTargets, 2);
    diagCell   = cell(1, nTargets);

    if isempty(opt.burnIn),      opt.burnIn      = 200 * L; end
    if isempty(opt.numSteps),    opt.numSteps    = 200 * L; end
    if isempty(opt.annealSteps), opt.annealSteps = 400 * L; end

    for i = 1:nTargets
        xT      = targetPhenotypes(i, :);
        logW0   = logWof(xT);
        angT    = atan2(xT(2), xT(1));
        logR0   = log(xT(2) / xT(1));

        tolF = opt.tolLogFitness;
        if isempty(tolF), tolF = 0.01 * abs(logW0); end

        % A window of half-width tolLogRatio on log(x2/x1) corresponds, to first
        % order at the target, to a window of this half-width on the polar angle.
        tolA = opt.tolLogRatio * abs(sin(2 * angT)) / 2;

        %% ---- Phase 1: the random-order walk, with rejection on overshoot -------------
        walkN = zeros(1, opt.walkSamples);
        gStart = []; xStart = [];
        for k = 1:opt.walkSamples
            [gk, xk, nk] = walkToStartingFitness(dxc, dxs, logWof, logW0, tolF, L);
            walkN(k) = nk;
            if isempty(gStart) || rand() < 1/k     % keep a uniform one of them
                gStart = gk; xStart = xk;
            end
        end
        runMean = cumsum(walkN) ./ (1:opt.walkSamples);
        half    = max(2, floor(opt.walkSamples/2));
        stationary = abs(runMean(end) - runMean(half)) <= 0.02 * abs(runMean(end));

        %% ---- Phase 2a: anneal into the set of genotypes meeting both requirements ---------------
        [g, x, entered] = annealToFeasible(gStart, xStart, dxc, dxs, logWof, ...
                                           logW0, tolF, angT, tolA, opt.annealSteps, L);

        %% ---- Phase 2b: uniform MH inside the set of genotypes meeting both requirements ---------------
        if entered
            [g, x, chainN, accRate] = uniformMH(g, x, dxc, dxs, logWof, ...
                                                logW0, tolF, angT, tolA, ...
                                                opt.burnIn, opt.numSteps, L);
        else
            chainN = []; accRate = NaN;
            warning('sampleGenotypesConditional:Infeasible', ...
                ['Target %d (log R = %.3f) was not reached. The W_0 contour may ', ...
                 'not contain that trait ratio under this genotype-phenotype map. ', ...
                 'Returning the closest genotype found.'], i, logR0);
        end

        %% ---- greedy comparison, for the record --------------------------
        nGreedy = greedyNumOnes(dxc, dxs, xT, L);

        genomes(i, :)    = g;
        phenotypes(i, :) = x;

        d = struct();
        d.numOnes          = sum(g);
        d.numOnesGreedy    = nGreedy;
        d.walkNumOnes      = walkN;
        d.walkRunningMean  = runMean;
        d.walkStationary   = stationary;
        d.chainNumOnes     = chainN;
        d.acceptRate       = accRate;
        d.feasible         = entered;
        d.logFitness       = logWof(x);
        d.logFitnessTarget = logW0;
        d.logRatio         = log(x(2) / x(1));
        d.logRatioTarget   = logR0;
        diagCell{i}        = d;

        if opt.verbose
            fprintf(['  target %d: log R %+.3f -> %+.3f | log W %+.4f -> %+.4f | ', ...
                     'ones %d (greedy %d) | walk mean %.1f%s | acc %.2f%s\n'], ...
                i, logR0, d.logRatio, logW0, d.logFitness, d.numOnes, nGreedy, ...
                runMean(end), ternary(stationary, '', ' [NOT STATIONARY]'), ...
                accRate, ternary(entered, '', ' [INFEASIBLE]'));
        end
    end

    diagnostics = [diagCell{:}];
end

% ---------------------------------------------------------------------------
function [g, x, n] = walkToStartingFitness(dxc, dxs, logWof, logW0, tolF, L)
% the random-order walk: flip loci 0 -> 1 in random order until log W first lands inside
% [logW0 - tolF, logW0 + tolF]. A walk that steps past the W_0 contour without
% landing in it is discarded and restarted. Without that rejection the landing
% point is the first crossing, which is biased toward small n.
    for attempt = 1:200
        order = randperm(L);
        g = zeros(1, L);
        x = [0, 0];
        for k = 1:L
            ell = order(k);
            g(ell) = 1;
            x = x + [dxc(ell), dxs(ell)];
            lw = logWof(x);
            if abs(lw - logW0) <= tolF
                n = k;
                return;
            elseif lw < logW0 - tolF
                break;                  % overshot; discard this walk
            end
        end
    end
    % Fall back on the closest point of the last walk rather than failing.
    n = sum(g);
end

% ---------------------------------------------------------------------------
function [g, x, ok] = annealToFeasible(g, x, dxc, dxs, logWof, logW0, tolF, angT, tolA, nSteps, L)
% Simulated annealing on the constraint violation, purely to enter the feasible
% set. Nothing about the sample distribution depends on this stage: the MH chain
% that follows is what defines the distribution.
    viol = @(xx) max(0, abs(logWof(xx) - logW0) - tolF) / tolF + ...
                 max(0, abs(wrapPi(atan2(xx(2), xx(1)) - angT)) - tolA) / tolA;
    v = viol(x);
    T0 = 1; T1 = 1e-3;
    for s = 1:nSteps
        if v <= 0, break; end
        T = T0 * (T1/T0)^((s-1)/max(1, nSteps-1));
        [xp, ell, sgn] = proposeFlip(g, x, dxc, dxs, L);
        vp = viol(xp);
        if vp <= v || rand() < exp(-(vp - v)/T)
            g(ell) = g(ell) + sgn;
            x = xp;
            v = vp;
        end
    end
    ok = (v <= 0);
end

% ---------------------------------------------------------------------------
function [g, x, chainN, accRate] = uniformMH(g, x, dxc, dxs, logWof, logW0, tolF, angT, tolA, burnIn, numSteps, L)
% Uniform distribution on the set of genotypes meeting both requirements. Two symmetric proposals, mixed 50/50:
%   (a) flip one uniformly chosen locus;
%   (b) swap one 1 with one 0 (keeps the number of
%       ones fixed, and moves the phenotype much further than a single flip,
%       which is what keeps the chain mixing inside a narrow contour).
% Accept iff the proposal is feasible. Symmetric proposal + uniform target
% means the acceptance ratio is 1 on the set of genotypes meeting both requirements and 0 outside it.
    feasible = @(xx) abs(logWof(xx) - logW0) <= tolF && ...
                     abs(wrapPi(atan2(xx(2), xx(1)) - angT)) <= tolA;

    total = burnIn + numSteps;
    thin  = max(1, floor(numSteps / 2000));
    chainN = zeros(1, floor(numSteps/thin));
    nAcc = 0; idx = 0;

    for s = 1:total
        if rand() < 0.5
            [xp, ell, sgn] = proposeFlip(g, x, dxc, dxs, L);
            if feasible(xp)
                g(ell) = g(ell) + sgn; x = xp; nAcc = nAcc + 1;
            end
        else
            ones_  = find(g == 1);
            zeros_ = find(g == 0);
            if ~isempty(ones_) && ~isempty(zeros_)
                e1 = ones_(randi(numel(ones_)));
                e0 = zeros_(randi(numel(zeros_)));
                xp = x - [dxc(e1), dxs(e1)] + [dxc(e0), dxs(e0)];
                if feasible(xp)
                    g(e1) = 0; g(e0) = 1; x = xp; nAcc = nAcc + 1;
                end
            end
        end
        if s > burnIn && mod(s - burnIn, thin) == 0
            idx = idx + 1;
            if idx <= numel(chainN), chainN(idx) = sum(g); end
        end
    end
    chainN  = chainN(1:idx);
    accRate = nAcc / total;
end

% ---------------------------------------------------------------------------
function [xp, ell, sgn] = proposeFlip(g, x, dxc, dxs, L)
    ell = randi(L);
    if g(ell) == 0
        sgn = 1;  xp = x + [dxc(ell), dxs(ell)];
    else
        sgn = -1; xp = x - [dxc(ell), dxs(ell)];
    end
end

% ---------------------------------------------------------------------------
function n = greedyNumOnes(dxc, dxs, xT, L)
% What initializeGenomeTheta would have produced, for the comparison printed
% above. Same single forward pass, same accept-if-closer rule.
    x = [0, 0]; n = 0;
    for ell = 1:L
        xn = x + [dxc(ell), dxs(ell)];
        if norm(xn - xT) < norm(x - xT)
            x = xn; n = n + 1;
        end
    end
end

% ---------------------------------------------------------------------------
function d = wrapPi(d)
    d = mod(d + pi, 2*pi) - pi;
end

function s = ternary(c, a, b)
    if c, s = a; else, s = b; end
end
