%% Run_standardFGM.m
% Main driver for the Standard Genotype–Phenotype–Fitness Model (FGM).
%
% Description:
%   Runs the Standard FGM under SSWM and CM (asexual/sexual) regimes,
%   reproducing Figure 2 equivalents for comparison with the Pleiotropic FGM.
%
% Modes:
%   Run_standardFGM          % full run (paper-scale)
%   Run_standardFGM('demo')  % fast sanity run, outputs -> ./results_demo_StandardFGM/...
%
% Outputs:
%   .mat files in ./results_StandardFGM/... (or ./results_demo_StandardFGM/...)
%   and figure files via makeFigure2.
%
% Reference:
%   Kim, M., Ardell, S. M., & Kryazhimskiy, S. (2025).
%   "Module-Selection Balance in the Evolution of Modular Organisms."
%
% -------------------------------------------------------------------------

function Run_standardFGM(mode)
    if nargin < 1, mode = 'demo'; end
    isDemo = strcmpi(mode, 'demo');

    % Route outputs cleanly by mode
    resultsRoot = tern(isDemo, '../results_demo_StandardFGM', '../results_StandardFGM');

    %% Environment setup
    thisDir  = fileparts(mfilename('fullpath'));   % .../StandardGPFM/main
    projRoot = fileparts(thisDir);                 % .../StandardGPFM
    addpath(genpath(projRoot));
    rootDisp(projRoot);

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
    % ---------------------------------------------------------------------

    fprintf('\n=== Standard FGM: Base Regimes ===\n');

    % ======================= A1) SSWM =========================
    sectionTimer('StandardFGM SSWM');
    simParams = initializeSimParams();
    resultStandardSSWM = simulateStandardSSWM(simParams);
    ave = computeAverageTrajectory(C.numTimeStamp, simParams, resultStandardSSWM.resultTable);

    outDir = ensureDir(thisDir, [resultsRoot '/SSWM']);
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

    outDir = ensureDir(thisDir, [resultsRoot '/CM_Asexual']);
    fname  = sprintf('%s/StandardFGM_CM_Asexual_M%.1e.mat', outDir, simParams.mutationRate);
    save(fname, 'simParams', 'resultStandardCM', 'ave');

    analyticalTrajectories = predictPleiotropicSSWM(simParams, ave);
    save(fname, 'analyticalTrajectories', '-append');

    % --------------------------- B) Figure 2 ------------------------------
    fprintf('\n=== Reproducing Figure 2 (StandardFGM) ===\n');
    makeFigure_StandardFGM(C.mutationRateSlow, C.mutationRateFast);

    fprintf('\nAll done. Mode: %s\n', mode);
end

% ============================== Helpers ==================================

function out = ensureDir(here, rel)
    out = fullfile(here, rel);
    if ~exist(out, 'dir'), mkdir(out); end
end

function sectionTimer(label)
    persistent lastT
    if isempty(lastT), lastT = tic; fprintf('\n'); return; end
    dt = toc(lastT);
    fprintf('[%s] finished in %.2f s\n', label, dt);
    lastT = tic;
end

function rootDisp(root)
    fprintf('Added all paths under: %s\n', root);
end

function y = tern(cond, a, b)
    if cond, y = a; else, y = b; end
end
