function [genomeParams] = initializeGenomeThetaRestricted(L, simParams, seed, thetaRange)
% initializeGenomeThetaRestricted - Pleiotropic genome init with a restricted
%   range of mutational angles.
%
% Description:
%   Identical to initializeGenomeTheta except that the per-locus pleiotropic
%   angles theta_ell are drawn uniformly from thetaRange = [thetaMin, thetaMax]
%   instead of from [0, 2*pi). Everything downstream - the greedy construction of
%   the initial genome, and both simulators - is unchanged, so any difference in
%   the resulting dynamics is attributable to the angle distribution alone.
%
% WHY A RESTRICTED CONE LOWERS THE GENOTYPIC REDUNDANCY OF THE OPTIMUM
%   The allele-1 loci of a genotype sitting at phenotype x satisfy
%       sum_{ell: g_ell = 1} delta*cos(theta_ell) = |x_1|,
%       sum_{ell: g_ell = 1} delta*sin(theta_ell) = |x_2|.
%
%   Under isotropic theta these sums contain terms of both signs, so cancellation
%   is possible: a near-optimal phenotype can be reached by many large opposing
%   contributions. The number of genotypes mapping to a near-optimal phenotype is
%   therefore large, the count of available beneficial mutations does not track
%   |x_i|, and the supply never declines.
%
%   Confine theta to a single quadrant and every term is non-negative.
%   Cancellation becomes impossible, the number of allele-1 loci is tightly
%   bounded by |x|, and as x_2 approaches its optimum the surviving allele-1 loci
%   are forced toward theta ~ 0, i.e. toward being x_1-improving. This is
%   assumption (ii) of the Necessary Conditions section - the optimum realized by
%   few genetic sequences - imposed on a map in which every mutation still
%   affects both traits.
%
%   Note that theta_ell is drawn ONCE per locus and fixed. The restriction is a
%   permanent property of the genotype-phenotype map, not a property of the
%   starting genotype: every mutation available at every generation comes from
%   the same restricted set of directions.
%
% SIGN CONVENTION - READ THIS BEFORE CHOOSING thetaRange
%   Because allele 1 displaces the phenotype by -delta*(cos theta, sin theta),
%   the reachable phenotypes from the all-zero genome lie in the cone spanned by
%   -(cos theta, sin theta) over theta in thetaRange. Simulations in this project
%   start at NEGATIVE trait values, so thetaRange must be chosen such that that
%   cone covers the third quadrant:
%
%       thetaRange = [0, pi/2]        allele 1 moves (x_1, x_2) by (-, -)   OK
%       thetaRange = [pi, 3*pi/2]     allele 1 moves (x_1, x_2) by (+, +)   FAILS
%
%   The second choice is the one that reads naturally if the mutation vector is
%   written as +delta*(cos theta, sin theta), which is the manuscript's
%   convention for the displacement of a 0 -> 1 flip. Under this code's
%   convention it is the wrong quadrant: the greedy initializer below can never
%   reduce the distance to a negative target, so it sets no loci to 1 and every
%   initial phenotype collapses to the origin. That failure is silent unless it
%   is checked for, so it is checked for here and raised as an error.
%
% Inputs:
%   L          - Number of loci in the genome
%   simParams  - Structure containing simulation parameters:
%       .initialPhenotypes - Target phenotype coordinates [nAngles x 2]
%       .deltaTrait        - Mutational step size (delta)
%   seed       - Random seed for reproducibility
%   thetaRange - [thetaMin, thetaMax], the interval the angles are drawn from.
%                Defaults to [0, 2*pi), which reproduces initializeGenomeTheta.
%
% Outputs:
%   genomeParams - Structure containing:
%       .genomeTheta       - [1 x L] pleiotropic angles
%       .initialGenomes    - [nAngles x L] initial allele states
%       .currentPhenotypes - [nAngles x 2] phenotypes encoded by initialGenomes
%       .thetaRange        - the interval used, carried through for provenance
%       .initDiagnostics   - struct with per-condition initialization quality:
%           .targetPhenotype   [nAngles x 2]
%           .realizedPhenotype [nAngles x 2]
%           .residual          [nAngles x 1] Euclidean target-realized distance
%           .nOnes             [nAngles x 1] number of loci set to allele 1
%           .meanThetaOnes     [nAngles x 1] mean angle among those loci
%
% Reference:
%   Kim, M., Ardell, S. M., & Kryazhimskiy, S. (2025).
%   "Module-Selection Balance in the Evolution of Modular Organisms."
%
% See also: initializeGenomeTheta, diagnoseRestrictedTheta,
%           Run_pleiotropicRestrictedTheta
%
% Copyright (c) 2025 Minkyu Kim, Cornell University
% Licensed under MIT License

    if nargin < 4 || isempty(thetaRange)
        thetaRange = [0, 2*pi];
    end
    validateattributes(thetaRange, {'numeric'}, {'vector', 'numel', 2, 'finite'});
    if thetaRange(2) <= thetaRange(1)
        error('initializeGenomeThetaRestricted:BadRange', ...
              'thetaRange must be increasing; got [%g, %g].', thetaRange(1), thetaRange(2));
    end

    if nargin > 2 && ~isempty(seed)
        rng(seed);
    end

    targetPhenotypes = simParams.initialPhenotypes;
    delta            = simParams.deltaTrait;
    numPhenotypes    = size(targetPhenotypes, 1);

    % Angles drawn once and shared across all initial conditions and replicates,
    % exactly as in initializeGenomeTheta.
    genomeTheta = thetaRange(1) + (thetaRange(2) - thetaRange(1)) * rand(1, L);

    initialGenomes    = zeros(numPhenotypes, L);
    currentPhenotypes = zeros(numPhenotypes, 2);

    % Greedy construction, iterated to convergence.
    %
    % initializeGenomeTheta makes a SINGLE forward pass and only ever flips
    % 0 -> 1. With angles spread over the full circle that is adequate: there is
    % always some direction that reduces the distance, so one pass lands close to
    % the target. Inside a cone it is not. Precisely because cancellation is
    % unavailable - the property that makes this model useful - an early flip
    % that helps overall can overshoot one coordinate, and a single forward pass
    % has no way to undo it. Measured on the six paper initial conditions with the cone at
    % [0, pi/2], one pass leaves residuals up to 1.17 and, worse, compresses the
    % spread of realized initial ratios R_0 from the intended 0.16-5.0 down to
    % about 0.5-1.7. Trajectories launched from those points would appear to
    % converge partly because they started converged.
    %
    % Iterating passes, and allowing flips in BOTH directions, fixes this: it is
    % ordinary coordinate descent on the squared distance to the target over the
    % hypercube of genotypes, and it terminates when no single flip helps.
    maxPasses = 50;

    for i = 1:numPhenotypes
        for pass = 1:maxPasses
            changed = false;
            for ell = 1:L
                % Flipping locus ell reverses the sign of its contribution.
                signFlip = 1 - 2*initialGenomes(i, ell);   % 0 -> +1 (turn on), 1 -> -1 (turn off)
                effect   = signFlip * [-delta * cos(genomeTheta(ell)), ...
                                       -delta * sin(genomeTheta(ell))];
                newPhenotype = currentPhenotypes(i, :) + effect;

                distanceBefore = norm(currentPhenotypes(i, :) - targetPhenotypes(i, :));
                distanceAfter  = norm(newPhenotype - targetPhenotypes(i, :));

                if distanceAfter < distanceBefore - 1e-12
                    initialGenomes(i, ell)  = 1 - initialGenomes(i, ell);
                    currentPhenotypes(i, :) = newPhenotype;
                    changed = true;
                end
            end
            if ~changed
                break;
            end
        end
    end

    % ------------------------- Diagnostics --------------------------------
    residual      = sqrt(sum((currentPhenotypes - targetPhenotypes).^2, 2));
    nOnes         = sum(initialGenomes, 2);
    meanThetaOnes = NaN(numPhenotypes, 1);
    for i = 1:numPhenotypes
        sel = initialGenomes(i, :) == 1;
        if any(sel)
            meanThetaOnes(i) = mean(genomeTheta(sel));
        end
    end

    genomeParams.genomeTheta       = genomeTheta;
    genomeParams.initialGenomes    = initialGenomes;
    genomeParams.currentPhenotypes = currentPhenotypes;
    genomeParams.thetaRange        = thetaRange;
    genomeParams.initDiagnostics   = struct( ...
        'targetPhenotype',   targetPhenotypes, ...
        'realizedPhenotype', currentPhenotypes, ...
        'residual',          residual, ...
        'nOnes',             nOnes, ...
        'meanThetaOnes',     meanThetaOnes);

    % --------------------- Fail loudly, not silently -----------------------
    % A cone pointing the wrong way produces genomes that are all zeros and
    % initial phenotypes at the origin. Refusing to continue here is the whole
    % point: a run that starts at the optimum will still complete and still
    % write a results file, and the problem would only surface as an
    % inexplicable figure much later.
    targetNorm   = sqrt(sum(targetPhenotypes.^2, 2));
    relResidual  = residual ./ max(targetNorm, eps);
    badInit      = relResidual > 0.10;

    % The initial RATIOS matter more than the initial positions, because the
    % whole experiment is about whether trajectories that start at different
    % R_0 converge. If the cone cannot place the population at the intended
    % spread of R_0, apparent convergence is partly built in from the start.
    ok = all(currentPhenotypes < 0, 2) & all(targetPhenotypes < 0, 2);
    if nnz(ok) >= 2
        wantSpread = max(log(abs(targetPhenotypes(ok,2)  ./ targetPhenotypes(ok,1)))) - ...
                     min(log(abs(targetPhenotypes(ok,2)  ./ targetPhenotypes(ok,1))));
        gotSpread  = max(log(abs(currentPhenotypes(ok,2) ./ currentPhenotypes(ok,1)))) - ...
                     min(log(abs(currentPhenotypes(ok,2) ./ currentPhenotypes(ok,1))));
        genomeParams.initDiagnostics.logR0SpreadRequested = wantSpread;
        genomeParams.initDiagnostics.logR0SpreadRealized  = gotSpread;
        if gotSpread < 0.7 * wantSpread
            warning('initializeGenomeThetaRestricted:CompressedInitialSpread', ...
                    ['The realized spread of initial log R_0 is %.2f against the ' ...
                     'requested %.2f. The restricted cone cannot reach the intended ' ...
                     'initial conditions, so the populations start closer together ' ...
                     'than intended and any convergence will be partly an artifact ' ...
                     'of initialization. Widen the cone or reduce the range of ' ...
                     'initial angles before interpreting this run.'], ...
                    gotSpread, wantSpread);
        end
    end

    if all(nOnes == 0)
        error('initializeGenomeThetaRestricted:ConePointsWrongWay', ...
              ['No locus improved the distance to ANY target phenotype, so every ' ...
               'initial genome is all zeros and every initial phenotype is the ' ...
               'optimum.\nthetaRange = [%g, %g] puts the allele-1 displacement ' ...
               '-delta*(cos theta, sin theta) outside the third quadrant.\n' ...
               'For targets with x_1 < 0 and x_2 < 0, use thetaRange = [0, pi/2]. ' ...
               'See the SIGN CONVENTION note in this file.'], ...
              thetaRange(1), thetaRange(2));
    end

    if any(badInit)
        warning('initializeGenomeThetaRestricted:PoorInitialization', ...
                ['%d of %d initial phenotypes are more than half their own norm ' ...
                 'away from the requested target. The restricted cone may not be ' ...
                 'able to reach them. Inspect genomeParams.initDiagnostics before ' ...
                 'interpreting the run.'], nnz(badInit), numPhenotypes);
    end
end
