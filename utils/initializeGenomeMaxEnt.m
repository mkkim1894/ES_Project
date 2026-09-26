function [genomeParams] = initializeGenomeMaxEnt(L, simParams, seed, thetaRange, opts)
% initializeGenomeMaxEnt - Initial genotypes drawn from the maximum-entropy
%   (microcanonical) ensemble at a given phenotype. Drop-in for
%   initializeGenomeSampled and initializeGenomeThetaRestricted: same first four
%   arguments, same output fields.
%
% NOTATION AGAINST THE DERIVATION
%   L here is the TOTAL number of loci, which is 2L in the derivation: his genotype runs
%   g_1 ... g_{2L}. The production runners set L = 200 loci per module and pass
%   2*L, the same convention initializeGenomeSampled uses. The derivation draws theta on
%   [pi, 3*pi/2) and sets e_l = -u_l; this code draws theta on [0, pi/2) and
%   uses e_l = (cos theta_l, sin theta_l) directly. The two agree exactly:
%   theta_SK = theta_here + pi, so -u_SK = e_here. (One typo in the note:
%   phi_l = -theta_l should read theta_l - pi; -theta_l lands in the second
%   quadrant. Nothing downstream uses phi, so the algebra is unaffected.)
%
% THE ENSEMBLE
%   With the sign convention of this codebase a locus at allele 1 displaces the
%   phenotype by -delta*(cos theta, sin theta), so writing
%
%       e_l = (cos theta_l, sin theta_l),      y(g) = -x(g)/delta = sum_l g_l e_l
%
%   the genotypes with a given phenotype are exactly the g with a given y. The
%   tilted ensemble
%
%       P_lambda(g) = exp(<lambda, y(g)>) / Z(lambda)
%
%   assigns equal probability to every genotype with the same y, whatever
%   lambda, so conditioned on y it IS the uniform (maximum-entropy) distribution
%   over genotypes at that phenotype. Because Z factorises,
%
%       Z(lambda) = prod_l (1 + exp(<lambda, e_l>)),
%
%   the loci are independent Bernoulli under P_lambda, with
%
%       p_l = 1 / (1 + exp(-<lambda, e_l>)).
%
%   Choosing lambda so that E[y] = y_0 makes the target phenotype typical rather
%   than astronomically rare, which is what lets simple rejection work. lambda
%   solves the 2-by-2 system
%
%       sum_l p_l(lambda) e_l = y_0,
%
%   whose Jacobian sum_l p_l(1-p_l) e_l e_l' is positive definite, so the
%   solution is unique and Newton converges.
%
% WHY THE TOLERANCE WINDOW NEEDS A REWEIGHT
%   This is the one deliberate addition to the maximum-entropy construction, which stops at "rejecting
%   genotypes whose phenotype y is outside of the tolerance window of y_0". Set
%   opts.reweight = false to follow the note literally.
%   y is a sum of L Bernoulli terms, so the set of genotypes hitting y_0 exactly
%   is generally empty and one has to accept a window |y - y_0| <= tol/delta.
%   Inside that window P_lambda is not flat: it varies as exp(<lambda, y - y_0>).
%   Accepted draws are therefore importance weighted by exp(-<lambda, y - y_0>)
%   before one is selected, which restores the uniform ensemble over the window.
%   The diagnostic weightRange reports how far from 1 that correction ran; with
%   the default tolerance it is small, but it is not negligible at loose
%   tolerances (it reaches ~2e4 at tol = 0.06 for a target near the cone edge).
%
% COST
%   Solving for lambda is a handful of Newton steps. Sampling is then plain
%   Bernoulli draws with an acceptance of roughly 0.04-1.4% at the default
%   tolerance, so a few hundred thousand draws per initial condition - seconds,
%   and about three orders of magnitude cheaper than the 9e5 swap moves that
%   initializeGenomeSampled uses.
%
% Inputs:
%   L          - number of loci
%   simParams  - needs .initialPhenotypes, .deltaTrait, .ellipseParams,
%                .landscapeStdDev
%   seed       - RNG seed. The angles theta_l are drawn immediately after
%                rng(seed), exactly as initializeGenomeSampled does, so the same
%                seed gives the same set of pleiotropic directions in both and
%                the two initializers can be compared locus by locus.
%   thetaRange - optional [thetaMin, thetaMax]; default [0, 2*pi) (unrestricted)
%   opts       - optional struct:
%                  .tol       phenotype tolerance, trait units (default 0.02)
%                  .nWant     accepted genotypes to pool before choosing one
%                             (default 200)
%                  .batch     draws per batch (default 5000)
%                  .maxDraws  give up after this many draws (default 4e6)
%                  .reweight  apply the importance correction (default true)
%                  .drawSeed  reseed the Bernoulli draws only, leaving the
%                             pleiotropic directions fixed by `seed`; use it to
%                             draw several genotypes from one ensemble
%
% Outputs:
%   genomeParams - .genomeTheta, .initialGenomes, .currentPhenotypes,
%                  .thetaRange, .initDiagnostics (nOnes, nOnesGreedy, residual,
%                  directionSpreadDeg, logR0Requested/Realized, and the
%                  max-entropy specific lambda, acceptRate, nDraws, weightRange,
%                  tol)
%
% See also: initializeGenomeSampled, solveMaxEntTilt
%
% Reference: Maximum-entropy (microcanonical) sampling; see the section "Sampling from the
%   microcanonical ensemble".
%
% Copyright (c) 2025 Minkyu Kim, Cornell University
% Licensed under MIT License

    if nargin < 3 || isempty(seed), seed = 1; end
    if nargin < 4 || isempty(thetaRange), thetaRange = [0, 2*pi]; end
    if nargin < 5, opts = struct(); end
    if ~isfield(opts, 'tol'),      opts.tol      = 0.02;  end
    if ~isfield(opts, 'nWant'),    opts.nWant    = 200;   end
    if ~isfield(opts, 'batch'),    opts.batch    = 5000;  end
    if ~isfield(opts, 'maxDraws'), opts.maxDraws = 4e6;   end
    if ~isfield(opts, 'reweight'), opts.reweight = true;  end
    rng(seed);

    targets = simParams.initialPhenotypes;
    delta   = simParams.deltaTrait;
    a       = simParams.ellipseParams;
    sig     = simParams.landscapeStdDev;
    nT      = size(targets, 1);

    genomeTheta = thetaRange(1) + diff(thetaRange) * rand(1, L);

    % The pleiotropic directions are fixed by `seed`. opts.drawSeed, when given,
    % reseeds only the Bernoulli sampling that follows, so repeated calls with
    % the same `seed` and different drawSeed return different genotypes drawn
    % from the same ensemble over the same loci. That is what lets one sample
    % the ensemble rather than a single representative of it.
    if isfield(opts, 'drawSeed') && ~isempty(opts.drawSeed)
        rng(opts.drawSeed);
    end

    e  = [cos(genomeTheta(:)), sin(genomeTheta(:))];     % L x 2
    ux = -delta * cos(genomeTheta);                      % allele 1 displacement
    uy = -delta * sin(genomeTheta);

    initialGenomes    = zeros(nT, L);
    currentPhenotypes = zeros(nT, 2);
    nOnes   = zeros(nT, 1);   nGreedy = zeros(nT, 1);
    resid   = zeros(nT, 1);   spread  = zeros(nT, 1);
    lamAll  = zeros(nT, 2);   accRate = zeros(nT, 1);
    nDrawn  = zeros(nT, 1);   wRange  = zeros(nT, 1);

    fprintf('initializeGenomeMaxEnt: theta ~ U[%.3f, %.3f], L = %d, tol = %g\n', ...
            thetaRange(1), thetaRange(2), L, opts.tol);

    for i = 1:nT
        xT = targets(i, :);
        y0 = -xT(:) / delta;                             % 2 x 1

        [lam, p] = solveMaxEntTilt(e, y0, sprintf('condition %d', i));
        lamAll(i, :) = lam';

        [g, d, nd, ar, wr] = drawOne(p, e, y0, delta, opts, lam);
        nDrawn(i) = nd; accRate(i) = ar; wRange(i) = wr;

        x = [sum(g .* ux), sum(g .* uy)];
        initialGenomes(i, :)    = g;
        currentPhenotypes(i, :) = x;
        nOnes(i) = sum(g);
        resid(i) = d;
        nGreedy(i) = sum(greedyGenome(ux, uy, xT, L));

        sgn  = 2*g - 1;
        dx1  = -sgn .* ux;  dx2 = -sgn .* uy;
        ben  = (logWv(x(1)+dx1, x(2)+dx2, a, sig) - logWv(x(1), x(2), a, sig)) > 0;
        spread(i) = circSpread(atan2(dx2(ben), dx1(ben)));

        fprintf(['  cond %d: 1s greedy %3d -> max-ent %3d | residual %.4f | ' ...
                 'direction spread %4.0f deg | accept %.3f%% | log R %+.3f ' ...
                 '(target %+.3f)\n'], ...
                i, nGreedy(i), nOnes(i), d, spread(i), 100*ar, ...
                log(x(2)/x(1)), log(xT(2)/xT(1)));
    end

    if any(resid > opts.tol * 1.0001)
        warning('initializeGenomeMaxEnt:Residual', ...
            'Conditions %s exceeded the tolerance %g.', ...
            mat2str(find(resid > opts.tol)'), opts.tol);
    end

    wantR = log(targets(:,2) ./ targets(:,1));
    gotR  = log(currentPhenotypes(:,2) ./ currentPhenotypes(:,1));

    genomeParams.genomeTheta       = genomeTheta;
    genomeParams.initialGenomes    = initialGenomes;
    genomeParams.currentPhenotypes = currentPhenotypes;
    genomeParams.thetaRange        = thetaRange;
    genomeParams.initDiagnostics   = struct( ...
        'nOnes',              nOnes, ...
        'nOnesGreedy',        nGreedy, ...
        'residual',           resid, ...
        'directionSpreadDeg', spread, ...
        'logR0Requested',     wantR, ...
        'logR0Realized',      gotR, ...
        'logR0SpreadRequested', max(wantR) - min(wantR), ...
        'logR0SpreadRealized',  max(gotR)  - min(gotR), ...
        'lambda',             lamAll, ...
        'acceptRate',         accRate, ...
        'nDraws',             nDrawn, ...
        'weightRange',        wRange, ...
        'tol',                opts.tol);

    fprintf('  log R_0 spread requested %.2f, realized %.2f\n\n', ...
            max(wantR) - min(wantR), max(gotR) - min(gotR));
end

% ---------------------------------------------------------------------------
function [gPick, dPick, nDraws, acceptRate, weightRange] = drawOne(p, e, y0, delta, opts, lam)
% Bernoulli draws from the tilted ensemble, keep those inside the tolerance
% window, importance weight them back to uniform, then take one.
    tolY  = opts.tol / delta;
    L     = numel(p);
    keepG = false(0, L);
    keepD = zeros(0, 1);
    keepW = zeros(0, 1);          % log importance weight, -<lambda, y - y0>
    nDraws = 0;

    while size(keepG, 1) < opts.nWant && nDraws < opts.maxDraws
        B  = min(opts.batch, opts.maxDraws - nDraws);
        G  = rand(B, L) < p';                     % B x L logical
        dy = double(G) * e - y0';                 % B x 2
        d  = hypot(dy(:,1), dy(:,2));
        hit = d <= tolY;
        nDraws = nDraws + B;
        if any(hit)
            keepG = [keepG; G(hit, :)];                    %#ok<AGROW>
            keepD = [keepD; d(hit) * delta];               %#ok<AGROW>
            keepW = [keepW; -(dy(hit, :) * lam)];          %#ok<AGROW>
        end
    end

    got = size(keepG, 1);
    if got == 0
        error('initializeGenomeMaxEnt:NoAccept', ...
            ['No draw landed within %g of the target in %d attempts. Loosen ' ...
             'opts.tol or raise opts.maxDraws.'], opts.tol, nDraws);
    end
    acceptRate = got / nDraws;

    if opts.reweight
        w = exp(keepW - max(keepW));
        weightRange = exp(max(keepW) - min(keepW));
    else
        w = ones(got, 1);
        weightRange = 1;
    end
    pick  = find(cumsum(w) >= rand * sum(w), 1, 'first');
    gPick = double(keepG(pick, :));
    dPick = keepD(pick);
end

% ---------------------------------------------------------------------------
function g = greedyGenome(ux, uy, xT, L)
% The greedy genotype, kept only as a diagnostic baseline for comparison with
% initializeGenomeSampled's nOnesGreedy. Never returned as the genotype.
    g = zeros(1, L); x = [0, 0];
    for sweep = 1:50
        improved = false;
        for ell = 1:L
            sgn = 1 - 2*g(ell);
            xn  = x + sgn * [ux(ell), uy(ell)];
            if norm(xn - xT) < norm(x - xT) - 1e-12
                g(ell) = 1 - g(ell); x = xn; improved = true;
            end
        end
        if ~improved, break; end
    end
end

% ---------------------------------------------------------------------------
function v = logWv(x1, x2, a, sig)
    v = -((x1./a(1)).^2 + (x2./a(2)).^2) ./ (2*sig^2);
end

function s = circSpread(ang)
    R = abs(mean(exp(1i*ang(:))));
    R = min(max(R, eps), 1);
    s = rad2deg(sqrt(-2*log(R)));
end
