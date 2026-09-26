function report = diagnoseRestrictedTheta(varargin)
% diagnoseRestrictedTheta - Pre-flight check for the restricted-cone pleiotropic
%   GPFM. Runs in seconds and does NOT simulate anything.
%
% Purpose
%   Two things are worth establishing before committing to a long run.
%
%   1. Sign convention. A cone pointing the wrong way yields all-zero genomes and
%      six populations sitting at the optimum. The simulation would still run to
%      completion and still write a results file. Check first.
%
%   2. Does the cone actually lower the genotypic redundancy of the optimum? The
%      claim is that confining theta removes the cancellation that lets many
%      genotypes map to a near-optimal phenotype, so that the supply of
%      trait-i-improving mutations declines as trait i improves - without giving
%      up universal pleiotropy. That is a statement about the genotype-phenotype
%      map, not about the dynamics, so it can be settled by walking a genome down
%      a trajectory and counting what remains available. This function does
%      exactly that, cheaply. It also reports the realized spread of initial
%      ratios, since a cone that cannot reach the intended initial conditions
%      would build convergence in before generation 1.
%
% Usage
%   diagnoseRestrictedTheta                          % default cone [0, pi/2]
%   diagnoseRestrictedTheta('coneHalfWidth', pi/8)   % widened cone
%   diagnoseRestrictedTheta('thetaRange', [pi, 3*pi/2])   % demonstrates the failure
%   r = diagnoseRestrictedTheta('verbose', false);
%
% Name-value pairs
%   'thetaRange'    [thetaMin, thetaMax]. Overrides coneHalfWidth when given.
%   'coneHalfWidth' w >= 0; cone is [-w, pi/2 + w]. Default 0, i.e. the first
%                   quadrant. Widening the cone restores cancellation and so
%                   raises redundancy: it is the quantitative control on the
%                   factor this model is testing. At w = pi/4 the cone is a
%                   half-circle and antagonistic mutations exist.
%   'L'             Number of loci. Default 400.
%   'initialAngles' Angles defining the initial conditions to check. Defaults to
%                   the six used at paper scale. The driver passes whatever the
%                   run will actually use, so a test-mode run is diagnosed
%                   against its own two conditions rather than against six it
%                   will never visit.
%   'seed'          RNG seed. Default 1.
%   'verbose'       Print the report. Default true.
%
% Output
%   report - struct with fields .thetaRange, .init (per-condition initialization
%            diagnostics) and .supply (supply composition along a descent path).
%
% Reference:
%   Kim, M., Ardell, S. M., & Kryazhimskiy, S. (2025).
%   "Module-Selection Balance in the Evolution of Modular Organisms."
%
% See also: initializeGenomeThetaRestricted, Run_pleiotropicRestrictedTheta
%
% Copyright (c) 2025 Minkyu Kim, Cornell University
% Licensed under MIT License

p = inputParser;
addParameter(p, 'thetaRange',    [],   @(v) isempty(v) || (isnumeric(v) && numel(v) == 2));
addParameter(p, 'coneHalfWidth', 0,    @(v) isnumeric(v) && isscalar(v) && v >= 0);
addParameter(p, 'L',             400,  @(v) isnumeric(v) && isscalar(v) && v > 0);
addParameter(p, 'seed',          1,    @isnumeric);
addParameter(p, 'initialAngles', [],   @(v) isempty(v) || isnumeric(v));
addParameter(p, 'verbose',       true, @islogical);
parse(p, varargin{:});

thetaRange = p.Results.thetaRange;
if isempty(thetaRange)
    w = p.Results.coneHalfWidth;
    thetaRange = [-w, pi/2 + w];
end
L       = p.Results.L;
verbose = p.Results.verbose;

% ------------------------------- Paths ---------------------------------
thisDir  = fileparts(mfilename('fullpath'));
if isempty(thisDir), thisDir = pwd; end
addpath(thisDir);
projRoot = fileparts(thisDir);
for d = {'simulation_scripts', 'analysis_scripts', 'utils', 'figure_scripts'}
    candidate = fullfile(projRoot, d{1});
    if isfolder(candidate)
        addpath(candidate);
    end
end

% --------------------------- Parameters --------------------------------
% Matched to the paper-scale pleiotropic runs.
initialAngles = p.Results.initialAngles;
if isempty(initialAngles)
    initialAngles = [atan2(1,6.4), atan2(1,3.2), atan2(1,1.6), ...
                     atan2(1,0.8), atan2(1,0.4), atan2(1,0.2)];
end

simParams = initializeSimParams( ...
    'numIteration',    1, ...
    'initialAngles',   initialAngles, ...
    'popSize',         1e4, ...
    'ellipseRatio',    sqrt(2), ...
    'deltaTrait',      0.1, ...
    'landscapeStdDev', 2, ...
    'mutationRate',    2e-6, ...
    'omitParams',      {'geneticTargetSize'});

if verbose
    fprintf('\n=== Restricted-cone pre-flight ===\n');
    fprintf('theta ~ U[%.4f, %.4f]  (%.1f deg wide)\n', ...
            thetaRange(1), thetaRange(2), rad2deg(diff(thetaRange)));
    fprintf('L = %d loci, delta = %g\n\n', L, simParams.deltaTrait);
end

% ------------------- 1. Initialization / sign check --------------------
try
    genomeParams = initializeGenomeThetaRestricted(L, simParams, p.Results.seed, thetaRange);
catch err
    if verbose
        fprintf(2, 'INITIALIZATION FAILED\n%s\n', err.message);
    end
    report = struct('thetaRange', thetaRange, 'init', [], 'supply', [], 'error', err.message);
    return;
end

d = genomeParams.initDiagnostics;
report.thetaRange = thetaRange;
report.init       = d;

if verbose
    fprintf('--- 1. Initialization ---\n');
    fprintf('%4s %20s %20s %10s %7s %9s\n', ...
            'cond', 'target (x1,x2)', 'realized (x1,x2)', 'residual', 'n(g=1)', 'mean th');
    for i = 1:size(d.targetPhenotype, 1)
        fprintf('%4d %10.3f %9.3f %10.3f %9.3f %10.4f %7d %9.3f\n', i, ...
                d.targetPhenotype(i,1),   d.targetPhenotype(i,2), ...
                d.realizedPhenotype(i,1), d.realizedPhenotype(i,2), ...
                d.residual(i), d.nOnes(i), d.meanThetaOnes(i));
    end
    targetNorm = sqrt(sum(d.targetPhenotype.^2, 2));
    worstRel   = max(d.residual ./ max(targetNorm, eps));
    if worstRel < 0.05
        fprintf('  -> initialization OK (worst residual %.1f%% of target norm)\n', 100*worstRel);
    else
        fprintf(2, '  -> initialization POOR (worst residual %.1f%% of target norm)\n', 100*worstRel);
    end

    % The spread of initial ratios is the quantity that actually matters: the
    % experiment asks whether trajectories starting at different R_0 converge.
    sel = all(d.targetPhenotype < 0, 2) & all(d.realizedPhenotype < 0, 2);
    if nnz(sel) >= 2
        want = log(abs(d.targetPhenotype(sel,2)   ./ d.targetPhenotype(sel,1)));
        got  = log(abs(d.realizedPhenotype(sel,2) ./ d.realizedPhenotype(sel,1)));
        fprintf('  log R_0 requested: %.2f to %.2f  (spread %.2f)\n', ...
                min(want), max(want), max(want)-min(want));
        fprintf('  log R_0 realized : %.2f to %.2f  (spread %.2f)\n', ...
                min(got),  max(got),  max(got)-min(got));
        if (max(got)-min(got)) < 0.7*(max(want)-min(want))
            fprintf(2, ['  -> INITIAL SPREAD COMPRESSED. The populations would start\n' ...
                        '     closer together than intended, so apparent convergence\n' ...
                        '     would be partly built in. Widen the cone.\n']);
        end
    end
    fprintf('\n');
end

% -------------- 2. Supply composition along a descent path -------------
% Walk one genome from its initial state toward the optimum by repeatedly
% fixing the single most beneficial available mutation, which is the
% successive-mutations limit with the stochastic part removed. At intervals,
% record what the remaining supply looks like. This is not a simulation of the
% dynamics; it is an inspection of how the AVAILABLE mutations change as the
% traits improve, which is the mechanism the model is supposed to exhibit.
condIdx = size(d.targetPhenotype, 1);   % the start with x_2 furthest behind
genome  = genomeParams.initialGenomes(condIdx, :);
x       = genomeParams.currentPhenotypes(condIdx, :);
theta   = genomeParams.genomeTheta;
delta = simParams.deltaTrait;
a     = simParams.ellipseParams;
sig   = simParams.landscapeStdDev;

logW = @(v) -((v(1)/a(1))^2 + (v(2)/a(2))^2) / (2*sig^2);

snapshots = struct('x', {}, 'nBeneficial', {}, 'fracImprovingX1', {}, ...
                   'fracImprovingX2', {}, 'meanDirX1', {}, 'meanDirX2', {});

% Snapshot cadence is set from the supply actually available, not from maxSteps.
% With a narrow cone only the allele-1 loci can be beneficial, so the descent may
% be only a few dozen steps long and a cadence based on maxSteps would record a
% single row.
maxSteps  = 20000;
nOnesHere = nnz(genome);
snapEvery = max(1, floor(max(nOnesHere, 12) / 12));

for step = 1:maxSteps
    signs = 2*genome - 1;                      % +1 if allele 1 (can flip to 0)
    cand1 = x(1) + signs .* delta .* cos(theta);
    cand2 = x(2) + signs .* delta .* sin(theta);
    sAll  = -((cand1./a(1)).^2 + (cand2./a(2)).^2) / (2*sig^2) - logW(x);

    ben = find(sAll > 0);
    if isempty(ben)
        break;
    end

    if mod(step-1, snapEvery) == 0
        % Direction of the phenotypic displacement each beneficial mutation
        % would cause. Positive component = improves that trait.
        dx1 = signs(ben) .* delta .* cos(theta(ben));
        dx2 = signs(ben) .* delta .* sin(theta(ben));
        k = numel(snapshots) + 1;
        snapshots(k).x               = x;
        snapshots(k).nBeneficial     = numel(ben);
        snapshots(k).fracImprovingX1 = mean(dx1 > 0);
        snapshots(k).fracImprovingX2 = mean(dx2 > 0);
        snapshots(k).meanDirX1       = mean(dx1);
        snapshots(k).meanDirX2       = mean(dx2);
    end

    [~, best] = max(sAll(ben));
    pick = ben(best);
    x(1) = x(1) + signs(pick) * delta * cos(theta(pick));
    x(2) = x(2) + signs(pick) * delta * sin(theta(pick));
    genome(pick) = 1 - genome(pick);

    if exp(logW(x)) >= 0.99
        break;
    end
end

report.supply = snapshots;

if verbose
    fprintf('--- 2. Supply composition along a greedy descent (condition %d) ---\n', condIdx);
    fprintf('    The question: does the supply of x_i-improving mutations shrink as\n');
    fprintf('    |x_i| shrinks? If it does, the cone reproduces the modular model''s\n');
    fprintf('    mechanism without giving up pleiotropy.\n\n');
    fprintf('%9s %9s %8s %12s %12s %11s %11s\n', ...
            'x1', 'x2', 'n_ben', 'frac imp x1', 'frac imp x2', 'E[dx1]', 'E[dx2]');
    for k = 1:numel(snapshots)
        s = snapshots(k);
        fprintf('%9.3f %9.3f %8d %12.3f %12.3f %11.4f %11.4f\n', ...
                s.x(1), s.x(2), s.nBeneficial, ...
                s.fracImprovingX1, s.fracImprovingX2, s.meanDirX1, s.meanDirX2);
    end
    fprintf('\n');
    if ~isempty(snapshots)
        fprintf('Read it this way: E[dx2] falling toward zero while E[dx1] stays\n');
        fprintf('positive means the available supply has shifted toward x1-improving\n');
        fprintf('mutations as x2 approached its optimum -- the intended bias.\n');
        fprintf('If both stay flat, the cone is not producing it and the run will\n');
        fprintf('most likely behave like the unrestricted pleiotropic GPFM.\n\n');
    end
end
end
