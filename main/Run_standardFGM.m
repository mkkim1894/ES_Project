%% Run_standardFGM.m
% Main driver for the Standard Genotype–Phenotype–Fitness Model (FGM).
%
% Description:
%   Runs the Standard FGM under SSWM and CM (asexual) regimes,
%   reproducing Figure for comparison with the Pleiotropic FGM.
%
% Modes:
%   Run_standardFGM          % full run (paper-scale)
%   Run_standardFGM('demo')  % fast sanity run
%
% Outputs:
%   .mat files in ./results_StandardFGM/... (or ./results_demo_StandardFGM/...)
%   and figure files via makeFigure_StandardFGM.
%
% Reference:
%   Kim, M., Ardell, S. M., & Kryazhimskiy, S. (2025).
%   "Module-Selection Balance in the Evolution of Modular Organisms."

function Run_standardFGM(mode)
    if nargin < 1, mode = 'demo'; end
    isDemo = strcmpi(mode, 'demo');

    %% Environment setup
    thisDir  = fileparts(mfilename('fullpath'));   % .../main
    projRoot = fileparts(thisDir);                 % project root
    addpath(fullfile(projRoot, 'simulation'));
    addpath(fullfile(projRoot, 'analysis'));
    addpath(fullfile(projRoot, 'utils'));
    addpath(fullfile(projRoot, 'figures'));
    fprintf('Project root: %s\n', projRoot);

    % Route outputs cleanly by mode
    resultsRoot = tern(isDemo, fullfile(projRoot, 'results_demo_StandardFGM'), fullfile(projRoot, 'results_StandardFGM'));

    % ------------------ Config knobs ------------------
    demoCfg.numIteration      = 8;
    demoCfg.numTimeStamp      = 8;
    demoCfg.mutationRateSlow  = 1e-7;
    demoCfg.mutationRateFast  = 1e-3;

    fullCfg.numIteration      = 200;
    fullCfg.numTimeStamp      = 20;
    fullCfg.mutationRateSlow  = 1e-7;
    fullCfg.mutationRateFast  = 1e-3;

    C = tern(isDemo, demoCfg, fullCfg);

    fprintf('\n=== Standard FGM: Base Regimes ===\n');

    % ======================= A1) SSWM =========================
    sectionTimer('StandardFGM SSWM');
    simParams = initializeSimParams();
    resultStandardSSWM = simulateStandardSSWM(simParams);
    ave = computeAverageTrajectory(C.numTimeStamp, simParams, resultStandardSSWM.resultTable);

    outDir = ensureDir(resultsRoot, 'SSWM');
    fname  = sprintf('%s/StandardFGM_SSWM_M%.1e.mat', outDir, simParams.mutationRate);
    save(fname, 'simParams', 'resultStandardSSWM', 'ave');

    analyticalTrajectories = predictPleiotropicSSWM(simParams, ave);
    save(fname, 'analyticalTrajectories', '-append');

    % ==================== A2) CM asexual =======================
    sectionTimer('StandardFGM CM Asexual');
    simParams = initializeSimParams('numIteration', C.numIteration, ...
                                    'mutationRate', C.mutationRateFast);
    resultStandardCM = simulateStandardCM(simParams);
    ave = computeAverageTrajectory(C.numTimeStamp, simParams, resultStandardCM.resultTable);

    outDir = ensureDir(resultsRoot, 'CM_Asexual');
    fname  = sprintf('%s/StandardFGM_CM_Asexual_M%.1e.mat', outDir, simParams.mutationRate);
    save(fname, 'simParams', 'resultStandardCM', 'ave');

    analyticalTrajectories = predictPleiotropicSSWM(simParams, ave);
    save(fname, 'analyticalTrajectories', '-append');

    % --------------------------- B) Figure ------------------------------
    fprintf('\n=== Reproducing Figure (StandardFGM) ===\n');
    makeFigure_StandardFGM(C.mutationRateSlow, C.mutationRateFast, mode);

    fprintf('\nAll done. Mode: %s\n', mode);
end

% ============================== Helpers ==================================

function out = ensureDir(base, rel)
    out = fullfile(base, rel);
    if ~exist(out, 'dir'), mkdir(out); end
end

function sectionTimer(label)
    persistent lastT
    if isempty(lastT), lastT = tic; fprintf('\n'); return; end
    dt = toc(lastT);
    fprintf('[%s] finished in %.2f s\n', label, dt);
    lastT = tic;
end

function y = tern(cond, a, b)
    if cond, y = a; else, y = b; end
end
