function [genomeParams] = initializeGenomeSampled(L, simParams, seed, thetaRange)
% initializeGenomeSampled - Initial genotypes drawn by the walk-and-swap procedure instead of
%   built greedily. Drop-in for initializeGenomeTheta and for
%   initializeGenomeThetaRestricted: same arguments, same output fields.
%
% WHY
%   The greedy initializers flip a locus to 1 only when that moves the phenotype
%   closer to the target, and stop on arrival. Two things go wrong:
%
%     1. Too few 1s. The genotype reaches the target with the fewest loci it can
%        - 41-49 of 400 at our initial conditions, where a genotype of that
%        fitness normally carries ~200. Almost every locus is then free to
%        mutate at generation 1.
%     2. The wrong loci. It picks the loci whose angles point most directly at
%        the target, so the mutations available at the start all point nearly the
%        same way - a spread of 8 degrees where a typical genotype gives 31.
%
% WHAT THIS DOES
%   Stage 1 - the random-order walk. From the all-zero genotype, flip loci 0 -> 1 in random
%       order. Every genotype met whose fitness is inside a tolerance window of
%       F0 is a candidate; a walk that steps past the window is abandoned and
%       restarted from the origin. This sets the NUMBER of 1s.
%   Stage 2 - land on the target. Of the candidates nearest the target phenotype,
%       one is taken at random (not the nearest - the nearest is systematically
%       the sparsest), then loci carrying 1 are swapped for loci carrying 0 until
%       the phenotype sits on the target. Swaps never change the number of 1s, so
%       stage 1's count survives.
%   Stage 3 - decorrelate. More swaps, accepted only while the phenotype stays on
%       the target. This randomises WHICH loci carry the 1s at fixed phenotype
%       and fixed count, removing the alignment bias.
%
% WHEN THE CONE IS NARROW (thetaRange much smaller than 2*pi)
%   Stage 1 becomes inert, by construction rather than by failure. Inside a cone
%   of width w every step vector points much the same way, so steps add
%   COHERENTLY and the distance from the optimum grows LINEARLY in the number of
%   1s: |x| ~ delta * n * 2*sin(w/2)/w. The count is then fixed by the geometry.
%   At our initial conditions with w = pi/2 that is n ~ 27-37, and a genotype
%   carrying 200 1s would sit at |x| ~ 18, far past any target on the F0 contour.
%   The "~200 typical" figure quoted above belongs to the UNRESTRICTED model,
%   where steps cancel, |x| ~ delta*sqrt(n), and a high count at modest |x| is
%   possible. So under a cone, nOnes close to nOnesGreedy is the correct answer,
%   not evidence that stage 1 is broken; what stages 2 and 3 still provide is the
%   direction spread, which is the half of the bias that a cone does not fix.
%   test_restrictedTheta checks nOnes against this geometric expectation.
%
% Inputs:
%   L          - number of loci
%   simParams  - needs .initialPhenotypes, .deltaTrait, .ellipseParams,
%                .landscapeStdDev
%   seed       - RNG seed
%   thetaRange - optional [thetaMin, thetaMax]; default [0, 2*pi) (unrestricted)
%
% Outputs:
%   genomeParams - .genomeTheta, .initialGenomes, .currentPhenotypes,
%                  .thetaRange, .initDiagnostics (nOnes, nOnesGreedy, residual,
%                  directionSpreadDeg, logR0 requested/realized)
%
% See also: initializeGenomeTheta, initializeGenomeThetaRestricted,
%           sampleGenotypesConditional
%
% Copyright (c) 2025 Minkyu Kim, Cornell University
% Licensed under MIT License

    if nargin < 3 || isempty(seed), seed = 1; end
    if nargin < 4 || isempty(thetaRange), thetaRange = [0, 2*pi]; end
    rng(seed);

    targets = simParams.initialPhenotypes;
    delta   = simParams.deltaTrait;
    a       = simParams.ellipseParams;
    sig     = simParams.landscapeStdDev;
    nT      = size(targets, 1);

    genomeTheta = thetaRange(1) + diff(thetaRange) * rand(1, L);
    ux = -delta * cos(genomeTheta);
    uy = -delta * sin(genomeTheta);

    initialGenomes    = zeros(nT, L);
    currentPhenotypes = zeros(nT, 2);
    nOnes   = zeros(nT, 1);
    nGreedy = zeros(nT, 1);
    resid   = zeros(nT, 1);
    spread  = zeros(nT, 1);

    fprintf('initializeGenomeSampled: theta ~ U[%.3f, %.3f], L = %d\n', ...
            thetaRange(1), thetaRange(2), L);

    for i = 1:nT
        xT = targets(i, :);

        gG = greedyGenome(ux, uy, xT, L);
        nGreedy(i) = sum(gG);

        g = sampleAtFitness(ux, uy, xT, L, a, sig, 20000);
        [g, d] = landAndDecorrelate(g, ux, uy, xT, 6e5, 0.06);

        x = [sum(g .* ux), sum(g .* uy)];
        initialGenomes(i, :)    = g;
        currentPhenotypes(i, :) = x;
        nOnes(i) = sum(g);
        resid(i) = d;

        sgn  = 2*g - 1;
        dx1  = -sgn .* ux;  dx2 = -sgn .* uy;
        ben  = (logWv(x(1)+dx1, x(2)+dx2, a, sig) - logWv(x(1), x(2), a, sig)) > 0;
        spread(i) = circSpread(atan2(dx2(ben), dx1(ben)));

        fprintf(['  cond %d: 1s greedy %3d -> sampled %3d | residual %.3f | ' ...
                 'direction spread %4.0f deg | log R %+.3f (target %+.3f)\n'], ...
                i, nGreedy(i), nOnes(i), d, spread(i), ...
                log(x(2)/x(1)), log(xT(2)/xT(1)));
    end

    if any(resid > 0.15)
        warning('initializeGenomeSampled:Residual', ...
            ['Conditions %s did not land within 0.15 of their target phenotype. ' ...
             'Their starting ratios are not the intended ones.'], ...
            mat2str(find(resid > 0.15)'));
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
        'logR0SpreadRealized',  max(gotR)  - min(gotR));

    fprintf('  log R_0 spread requested %.2f, realized %.2f\n\n', ...
            max(wantR) - min(wantR), max(gotR) - min(gotR));
end

% ---------------------------------------------------------------------------
function g = greedyGenome(ux, uy, xT, L)
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
function gPick = sampleAtFitness(ux, uy, xT, L, a, sig, nCollect)
    logW0 = -((xT(1)/a(1))^2 + (xT(2)/a(2))^2) / (2*sig^2);
    tolF  = 0.01 * abs(logW0);

    POOL  = 600;
    poolG = zeros(POOL, L);
    poolD = inf(POOL, 1);
    worst = 1;

    got = 0; tries = 0;
    while got < nCollect && tries < 400*nCollect
        tries = tries + 1;
        ord = randperm(L);
        X1 = cumsum(ux(ord));  X2 = cumsum(uy(ord));
        lw = -((X1./a(1)).^2 + (X2./a(2)).^2) ./ (2*sig^2);

        inWin = find(abs(lw - logW0) <= tolF);
        if isempty(inWin), continue; end
        over = find(lw < logW0 - tolF, 1, 'first');
        if ~isempty(over), inWin = inWin(inWin <= over); end
        if isempty(inWin), continue; end

        for k = inWin(:)'
            got = got + 1;
            d = hypot(X1(k) - xT(1), X2(k) - xT(2));
            if d < poolD(worst)
                gk = zeros(1, L); gk(ord(1:k)) = 1;
                poolG(worst, :) = gk;
                poolD(worst)    = d;
                [~, worst] = max(poolD);
            end
        end
    end

    keep = find(isfinite(poolD));
    if isempty(keep)
        error('initializeGenomeSampled:NoCandidates', ...
              ['No genotype reached the tolerance window of F0. With L = %d and ' ...
               'delta = %g a random walk may be too short to reach that fitness.'], ...
              L, abs(ux(1))/cos(0));
    end
    [~, ordD] = sort(poolD(keep));
    nNear = max(1, round(0.1 * numel(keep)));
    gPick = poolG(keep(ordD(randi(nNear))), :);
end

% ---------------------------------------------------------------------------
function [g, d] = landAndDecorrelate(g, ux, uy, xT, nSwaps, tol)
    L = numel(g);
    x = [sum(g .* ux), sum(g .* uy)];
    d = norm(x - xT);

    half = floor(nSwaps/2);
    for k = 1:half
        if d <= tol, break; end
        i1 = randi(L); if g(i1) ~= 1, continue; end
        i0 = randi(L); if g(i0) ~= 0, continue; end
        xn = x - [ux(i1), uy(i1)] + [ux(i0), uy(i0)];
        dn = norm(xn - xT);
        if dn < d, g(i1) = 0; g(i0) = 1; x = xn; d = dn; end
    end

    for k = 1:nSwaps
        i1 = randi(L); if g(i1) ~= 1, continue; end
        i0 = randi(L); if g(i0) ~= 0, continue; end
        xn = x - [ux(i1), uy(i1)] + [ux(i0), uy(i0)];
        dn = norm(xn - xT);
        if dn <= max(tol, d), g(i1) = 0; g(i0) = 1; x = xn; d = dn; end
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
