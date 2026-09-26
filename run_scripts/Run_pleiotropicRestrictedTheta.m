function Run_pleiotropicRestrictedTheta(varargin)
% Run_pleiotropicRestrictedTheta - Driver for the pleiotropic GPFM with a
%   restricted range of mutational angles.
%
% WHAT THIS MODEL IS FOR
%   It has been argued that module-selection balance may be a property of a
%   changing mutational supply rather than of variational modularity. Answering
%   that means crossing the two factors rather than arguing about which matters.
%
%   The supply in the modular GPFM declines because of assumption (ii) of the
%   Necessary Conditions section: the optimal trait value is realized by a single
%   genetic sequence. Low genotypic redundancy near the optimum is what makes
%   b_i = |x_i|/delta shrink as module i improves.
%
%                     | low redundancy          | high redundancy
%     ----------------+-------------------------+--------------------------
%     modular GPM     | main model -> MSB       | constant supply -> MSB only
%                     |                         | under complete linkage
%     pleiotropic GPM | THIS MODEL -> ?         | standard -> no MSB
%
%   The constant-supply model relaxes assumption (ii) while keeping modularity.
%   This model tightens assumption (ii) while keeping universal pleiotropy.
%
% HOW THE CONE LOWERS REDUNDANCY
%   The allele-1 loci of a genotype at phenotype x satisfy
%       sum_{g=1} delta*cos(theta) = |x_1|,
%       sum_{g=1} delta*sin(theta) = |x_2|.
%   Under isotropic theta these sums contain terms of both signs, so cancellation
%   is possible and a near-optimal phenotype can be reached by many large
%   opposing contributions. Redundancy near the optimum stays high and the supply
%   never declines - which is why the standard pleiotropic GPFM behaves like
%   gradient ascent. Confine theta to one quadrant and every term is
%   non-negative, cancellation is impossible, few genotypes map to near-optimal
%   phenotypes, and as x_2 approaches its optimum the surviving allele-1 loci are
%   forced toward theta ~ 0, i.e. toward being x_1-improving.
%
%   Every mutation still moves BOTH traits: a quadrant is symmetric about
%   theta = pi/4, so nothing is axis-aligned and nothing affects a single trait.
%   The cone constrains the sign pattern, not the pleiotropy; what it removes is
%   antagonistic mutations that improve one trait while degrading the other.
%
% A PROPERTY OF THE MAP, NOT OF THE STARTING POINT
%   theta_ell is drawn once per locus and fixed. There is no separate mutation
%   process to restrict - the direction a mutation moves the population is
%   whatever angle its locus was assigned. The restriction is instantiated when
%   the map is built and persists for the whole run: every mutation available at
%   every generation comes from the same restricted set of directions. Worth
%   saying explicitly in the Methods, since "restricted at initialization" can be
%   misread as a transient starting-point effect.
%
% INTERPRETING THE OUTCOME
%   If this model does NOT produce a balance, reducing redundancy alone is not
%   sufficient, modularity is doing the work, and the claim stands as written.
%   If it DOES, the claim narrows from "variational modularity is necessary" to
%   "low genotypic redundancy of the trait optimum is necessary, and variational
%   modularity is its biologically documented realization". That is a sharpening
%   rather than a retraction, and it fits the rest of the section: the nested FGM
%   keeps the balance under a non-linear decline in supply, and the
%   constant-supply model loses it except under complete linkage.
%
%   'coneHalfWidth' is the quantitative control on redundancy rather than a mere
%   robustness knob: widening the cone restores cancellation. At w = pi/4 the
%   cone is a half-circle, antagonistic mutations exist, and the pre-flight shows
%   the supply decline largely disappearing (114 -> 67 available mutations,
%   against 24 -> 4 at w = 0).
%
% Usage
%   Run_pleiotropicRestrictedTheta                                  % diagnose + SSWM
%   Run_pleiotropicRestrictedTheta('mode', 'test')                  % fast sanity run
%   Run_pleiotropicRestrictedTheta('coneHalfWidth', pi/8)           % widened cone
%   Run_pleiotropicRestrictedTheta('stages', {'cm_asexual','cm_sexual','figures'})
%   Run_pleiotropicRestrictedTheta('stages', {'figures'})           % replot only
%
% Name-value pairs
%   'mode'          'reproduce' (default, paper scale) or 'test'.
%   'thetaRange'    [thetaMin, thetaMax] the angles are drawn from. Overrides
%                   coneHalfWidth when supplied. See the SIGN CONVENTION note in
%                   initializeGenomeThetaRestricted before setting this by hand:
%                   in this codebase the first quadrant [0, pi/2] is the cone
%                   that reaches negative initial phenotypes, NOT [pi, 3*pi/2].
%   'coneHalfWidth' w >= 0; cone is [-w, pi/2 + w]. Default 0.
%   'stages'        Any of 'diagnose', 'sswm', 'cm_asexual', 'cm_sexual',
%                   'figures'. Default {'diagnose', 'sswm'}.
%
%                   THE DEFAULT IS A FIRST LOOK, NOT THE FULL EXPERIMENT. All
%                   three regimes are needed, for the same reason the
%                   constant-supply model needed them: that model's answer turned
%                   out to be regime-dependent - no balance under successive
%                   mutations or free reassortment, but a balance under complete
%                   linkage, where clonal interference substitutes for the
%                   declining supply. The comparison between the two stress tests
%                   is therefore only meaningful regime by regime.
%
%                   SSWM is run first because it is minutes rather than hours,
%                   and because without clonal interference it isolates the
%                   redundancy mechanism from the selection-coupling mechanism
%                   that the constant-supply model already showed can produce a
%                   balance on its own. A balance appearing here is the
%                   conclusion-changing outcome and is worth knowing before
%                   committing to the concurrent-mutations runs.
%   'seed'          RNG seed for the angle draw and genome construction.
%                   Default 1.
%
% Outputs
%   .mat files under <this folder>/results/RestrictedTheta/<regime>/
%   (results_test/ in test mode), one per regime, with the cone encoded in the
%   filename. Each carries simParams, genomeParams (including the angle range
%   and the initialization diagnostics), the result struct, the averaged
%   trajectory, and referenceTrajectories.
%
% NOTE ON referenceTrajectories
%   The curves saved under that name are predictPleiotropicSSWM, i.e. the
%   gradient-ascent path of the UNRESTRICTED, isotropic FGM. They are not a
%   prediction for this model - restricting the cone breaks the isotropy that
%   derivation assumes. They are saved and plotted as the reference the
%   restricted model should be compared against: a population that still follows
%   them has not changed its behaviour, and a departure from them is the result.
%
% Requirements
%   The rest of the ES_Project tree (simulation_scripts, analysis_scripts,
%   utils), located automatically as the parent of this folder.
%   Also required for the SSWM stage: mybinornd, called by
%   simulatePleiotropicSSWM but not part of the repository. Add the folder
%   containing it to the path. See README.md.
%
% Reference:
%   Kim, M., Ardell, S. M., & Kryazhimskiy, S. (2025).
%   "Module-Selection Balance in the Evolution of Modular Organisms."
%
% See also: initializeGenomeThetaRestricted, diagnoseRestrictedTheta,
%           makeFigure_RestrictedTheta, Run_pleiotropicFGM
%
% Copyright (c) 2025 Minkyu Kim, Cornell University
% Licensed under MIT License

% ----------------------------- Parse inputs -----------------------------
p = inputParser;
addParameter(p, 'mode',          'reproduce', @ischar);
addParameter(p, 'thetaRange',    [],  @(v) isempty(v) || (isnumeric(v) && numel(v) == 2));
addParameter(p, 'coneHalfWidth', 0,   @(v) isnumeric(v) && isscalar(v) && v >= 0);
addParameter(p, 'stages',        {'diagnose', 'sswm'}, @iscellstr);
addParameter(p, 'seed',          1,   @isnumeric);
addParameter(p, 'init',          'greedy', @ischar);   % 'greedy'|'sampled'|'maxent'
addParameter(p, 'outRoot',       '', @ischar);          % base folder for results
parse(p, varargin{:});

initMode = lower(p.Results.init);
switch initMode
    case 'greedy'
        initFcn = @(LL, sp, sd, tr) initializeGenomeThetaRestricted(LL, sp, sd, tr);
        modelTag = 'RestrictedThetaFGM';
    case 'sampled'
        % the random-order walk sets the number of 1s; swaps land the phenotype and
        % randomise which loci carry them. See initializeGenomeSampled.
        initFcn = @(LL, sp, sd, tr) initializeGenomeSampled(LL, sp, sd, tr);
        modelTag = 'RestrictedThetaFGM-sampledinit';
    case 'maxent'
        % Uniform over all genotypes at the target phenotype, sampled by
        % exponential tilting. See initializeGenomeMaxEnt.
        initFcn = @(LL, sp, sd, tr) initializeGenomeMaxEnt(LL, sp, sd, tr);
        modelTag = 'RestrictedThetaFGM-maxentinit';
    otherwise
        error('Run_pleiotropicRestrictedTheta:BadInit', ...
              'init must be ''greedy'', ''sampled'' or ''maxent''.');
end

mode   = lower(p.Results.mode);
stages = lower(p.Results.stages);
seed   = p.Results.seed;

thetaRange = p.Results.thetaRange;
if isempty(thetaRange)
    w = p.Results.coneHalfWidth;
    thetaRange = [-w, pi/2 + w];
end

validStages = {'diagnose', 'sswm', 'cm_asexual', 'cm_sexual', 'figures'};
unknown = setdiff(stages, validStages);
if ~isempty(unknown)
    error('Run_pleiotropicRestrictedTheta:UnknownStage', ...
          'Unknown stage: %s. Valid stages are %s.', ...
          strjoin(unknown, ', '), strjoin(validStages, ', '));
end

switch mode
    case 'test',                    isTest = true;
    case {'reproduce', 'full'},     isTest = false;
    otherwise
        error('Run_pleiotropicRestrictedTheta:UnknownMode', ...
              'Unknown mode ''%s''. Use ''test'' or ''reproduce''.', mode);
end

% ------------------------------ Paths -----------------------------------
thisDir = fileparts(mfilename('fullpath'));
if isempty(thisDir), thisDir = pwd; end
addpath(thisDir);

projRoot = fileparts(thisDir);
required = {'simulation_scripts', 'analysis_scripts', 'utils'};
for d = required
    candidate = fullfile(projRoot, d{1});
    if ~isfolder(candidate)
        error('Run_pleiotropicRestrictedTheta:MissingTree', ...
              ['Expected to find %s next to this folder. This extension is ' ...
               'self-contained only in the sense that it adds no files to the ' ...
               'main tree; it still needs the simulators and utilities from it.'], ...
              candidate);
    end
    addpath(candidate);
end
figScripts = fullfile(projRoot, 'figure_scripts');
if isfolder(figScripts), addpath(figScripts); end

fprintf('Project root: %s\n', projRoot);

% --------------------------- Configuration ------------------------------
% Identical to Run_pleiotropicFGM so that the only difference from the
% published pleiotropic runs is the angle distribution.
K = 10;    % reference genetic target size; per-locus rate mu = U_ref / (2*K)
L = 200;   % loci per module, total 2*L

testCfg.numIteration      = 8;
testCfg.numTimeStamp      = 20;
testCfg.initialAngles     = [atan2(1,6.4), atan2(1,0.2)];
testCfg.mutationRateSlow  = 1e-7 * (L / K);
testCfg.mutationRateFast  = 2e-4 * (L / K);
testCfg.popSize           = 1e4;
testCfg.ellipseRatio      = sqrt(2);
testCfg.deltaTrait        = 0.1;
testCfg.landscapeStdDev   = 2;
testCfg.L_total           = 2*L;

reproduceCfg.numIteration      = 250;
reproduceCfg.numTimeStamp      = 20;
reproduceCfg.initialAngles     = [atan2(1,6.4), atan2(1,3.2), atan2(1,1.6), ...
                                  atan2(1,0.8), atan2(1,0.4), atan2(1,0.2)];
reproduceCfg.mutationRateSlow  = 1e-7 * (L / K);
reproduceCfg.mutationRateFast  = 2e-4 * (L / K);
reproduceCfg.popSize           = 1e4;
reproduceCfg.ellipseRatio      = sqrt(2);
reproduceCfg.deltaTrait        = 0.1;
reproduceCfg.landscapeStdDev   = 2;
reproduceCfg.L_total           = 2*L;

C = tern(isTest, testCfg, reproduceCfg);

% outRoot lets a caller put the results somewhere other than this folder, so
% that a run with a different initializer keeps its own tree.
outRoot = p.Results.outRoot;
if isempty(outRoot), outRoot = projRoot; end
resultsRoot = fullfile(outRoot, tern(isTest, 'results_test', 'results'));
% One subtree per initializer. The three initializers produce genuinely
% different runs of the same model, and keeping them apart means a figure
% script never has to disambiguate two files that match the same glob.
rtRoot      = fullfile(resultsRoot, ['RestrictedTheta_' lower(initMode) 'init']);

fprintf('\n=== Pleiotropic GPFM, restricted mutational cone (%s mode) ===\n', mode);
fprintf('theta ~ U[%.4f, %.4f]  (%.1f deg wide; unrestricted is 360 deg)\n', ...
        thetaRange(1), thetaRange(2), rad2deg(diff(thetaRange)));
fprintf('Stages: %s   |   initialization: %s\n\n', strjoin(stages, ', '), initMode);

% ------------------------------ Diagnose --------------------------------
if ismember('diagnose', stages)
    diagnoseRestrictedTheta('thetaRange', thetaRange, 'L', C.L_total, ...
                            'seed', seed, 'initialAngles', C.initialAngles);
end

% -------------------------------- SSWM ----------------------------------
if ismember('sswm', stages)
    if isempty(which('mybinornd'))
        error('Run_pleiotropicRestrictedTheta:MissingMybinornd', ...
              ['mybinornd is not on the MATLAB path. simulatePleiotropicSSWM ' ...
               'calls it, but it is not part of the ES_Project repository - add ' ...
               'the folder that contains it to the path. See README.md.']);
    end

    fprintf('--- SSWM ---\n'); t0 = tic;
    simParams    = baseParams(C, C.mutationRateSlow, 0);
    genomeParams = initFcn(C.L_total, simParams, seed, thetaRange);

    resultRestrictedSSWM = simulatePleiotropicSSWM(simParams, genomeParams);
    ave = computeAverageTrajectory(C.numTimeStamp, simParams, resultRestrictedSSWM.resultTable);
    fprintf('  completed in %.1f s\n', toc(t0));
    reportTermination('SSWM', resultRestrictedSSWM.terminationStatus);

    referenceTrajectories = predictPleiotropicSSWM(simParams, ave, genomeParams); %#ok<NASGU>

    fname = fullfile(ensureDir(rtRoot, 'SSWM'), ...
                     buildFilename(modelTag, 'SSWM', simParams, C.L_total, thetaRange));
    save(fname, 'simParams', 'genomeParams', 'resultRestrictedSSWM', 'ave', ...
                'referenceTrajectories', 'thetaRange');
    fprintf('  saved %s\n\n', fname);
end

% ---------------------------- CM, linked --------------------------------
if ismember('cm_asexual', stages)
    fprintf('--- Concurrent mutations, complete linkage ---\n'); t0 = tic;
    simParams    = baseParams(C, C.mutationRateFast, 0);
    genomeParams = initFcn(C.L_total, simParams, seed, thetaRange);

    resultRestrictedCM = simulatePleiotropicCM(simParams, genomeParams);
    ave = computeAverageTrajectory(C.numTimeStamp, simParams, resultRestrictedCM.resultTable);
    fprintf('  completed in %.1f s\n', toc(t0));
    reportTermination('CM linked', resultRestrictedCM.terminationStatus);

    referenceTrajectories = predictPleiotropicSSWM(simParams, ave, genomeParams); %#ok<NASGU>

    fname = fullfile(ensureDir(rtRoot, 'CM_Asexual'), ...
                     buildFilename(modelTag, 'CM_Asexual', simParams, C.L_total, thetaRange));
    save(fname, 'simParams', 'genomeParams', 'resultRestrictedCM', 'ave', ...
                'referenceTrajectories', 'thetaRange');
    fprintf('  saved %s\n\n', fname);
end

% --------------------------- CM, unlinked -------------------------------
if ismember('cm_sexual', stages)
    fprintf('--- Concurrent mutations, free reassortment ---\n'); t0 = tic;
    simParams    = baseParams(C, C.mutationRateFast, 1);
    genomeParams = initFcn(C.L_total, simParams, seed, thetaRange);

    resultRestrictedCM = simulatePleiotropicCM(simParams, genomeParams);
    ave = computeAverageTrajectory(C.numTimeStamp, simParams, resultRestrictedCM.resultTable);
    fprintf('  completed in %.1f s\n', toc(t0));
    reportTermination('CM unlinked', resultRestrictedCM.terminationStatus);

    referenceTrajectories = predictPleiotropicSSWM(simParams, ave, genomeParams); %#ok<NASGU>

    fname = fullfile(ensureDir(rtRoot, 'CM_Sexual'), ...
                     buildFilename(modelTag, 'CM_Sexual', simParams, C.L_total, thetaRange));
    save(fname, 'simParams', 'genomeParams', 'resultRestrictedCM', 'ave', ...
                'referenceTrajectories', 'thetaRange');
    fprintf('  saved %s\n\n', fname);
end

% ------------------------------ Figures ---------------------------------
if ismember('figures', stages)
    fprintf('--- Figure ---\n');
    makeFigure_RestrictedTheta('resultsRoot', resultsRoot, 'thetaRange', thetaRange, ...
                               'init', initMode);
end

fprintf('All done. Mode: %s\n', mode);
end

% ============================== Helpers ==================================

function simParams = baseParams(C, U, rho)
    args = {'numIteration',    C.numIteration, ...
            'initialAngles',   C.initialAngles, ...
            'popSize',         C.popSize, ...
            'ellipseRatio',    C.ellipseRatio, ...
            'deltaTrait',      C.deltaTrait, ...
            'landscapeStdDev', C.landscapeStdDev, ...
            'mutationRate',    U, ...
            'omitParams',      {'geneticTargetSize'}};
    if rho > 0
        args = [args, {'recombinationRate', rho}];
    end
    simParams = initializeSimParams(args{:});
end

function name = buildFilename(model, regime, sp, L, thetaRange)
    name = sprintf('%s_%s_N%.0e_L%d_M%.1e_d%.2f_eR%.2f_s%.2f', ...
        model, regime, sp.popSize, L, sp.mutationRate, sp.deltaTrait, ...
        sp.ellipseParams(1)/sp.ellipseParams(2), sp.landscapeStdDev);
    if isfield(sp, 'recombinationRate') && sp.recombinationRate > 0
        name = sprintf('%s_R%.4f', name, sp.recombinationRate);
    end
    % The cone is the whole point of this extension, so it goes in the name.
    name = sprintf('%s_th%+.3f%+.3f.mat', name, thetaRange(1), thetaRange(2));
end

function reportTermination(label, statusMatrix)
    ts = statusMatrix(:);
    fprintf('  %s termination: %d/%d reached the fitness threshold, %d ran out of beneficial mutations\n', ...
            label, sum(ts == 1), numel(ts), sum(ts == 0));
    if sum(ts == 0) > 0.1 * numel(ts)
        fprintf(2, ['  NOTE: more than 10%% of replicates stopped for want of beneficial\n' ...
                    '  mutations rather than by reaching the fitness threshold. With a\n' ...
                    '  restricted cone that is expected near the optimum, but check the\n' ...
                    '  endpoints before reading the log-ratio panel.\n']);
    end
end

function out = ensureDir(base, rel)
    out = fullfile(base, rel);
    if ~isfolder(out), mkdir(out); end
end

function y = tern(cond, a, b)
    if cond, y = a; else, y = b; end
end
