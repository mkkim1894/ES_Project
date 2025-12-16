%% Run_pleiotropicFGM.m
% Main driver for the Pleiotropic Genotype–Phenotype–Fitness Model (GPFM).
%
% Description:
%   Runs pleiotropic FGM under SSWM and CM (asexual/sexual), and reproduces Figure 2.
%
% Modes:
%   Run_pleiotropicFGM          % full run (paper-scale)
%   Run_pleiotropicFGM('demo')  % fast sanity run
%
% Outputs:
%   .mat files in ./results/... (or ./results_demo/...) and figure files via makeFigure2.
%
% Reference:
%   Kim, M., Ardell, S. M., & Kryazhimskiy, S. (2025).
%   "Module-Selection Balance in the Evolution of Modular Organisms."

function Run_pleiotropicFGM(mode)
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
    resultsRoot = tern(isDemo, fullfile(projRoot, 'results_demo'), fullfile(projRoot, 'results'));

    % ------------------ Config knobs ------------------
    demoCfg.numIteration      = 8;
    demoCfg.numTimeStamp      = 8;
    demoCfg.recombList        = [0, 1];
    demoCfg.initialAngles     = [atan2(1,6.4), atan2(1,0.2)];
    demoCfg.mutationRateSlow  = 1e-7;
    demoCfg.mutationRateFast  = 2e-4;
    demoCfg.L_SSWM            = 4000;
    demoCfg.L_CM              = 400;

    fullCfg.numIteration      = 200;
    fullCfg.numTimeStamp      = 20;
    fullCfg.recombList        = [0, 1];
    fullCfg.initialAngles     = [atan2(1,6.4), atan2(1,3.2), atan2(1,1.6), ...
                                 atan2(1,0.8), atan2(1,0.4), atan2(1,0.2)];
    fullCfg.mutationRateSlow  = 1e-7;
    fullCfg.mutationRateFast  = 2e-4;
    fullCfg.L_SSWM            = 4000;
    fullCfg.L_CM              = 400;

    C = tern(isDemo, demoCfg, fullCfg);

    fprintf('\n=== Pleiotropic FGM: Base Regimes ===\n');

    % ======================= A1) SSWM (L = 4000) =========================
    sectionTimer('PleiotropicFGM SSWM');
    simParams = initializeSimParams('omitParams', {'geneticTargetSize'});
    genomeParams = initializeGenomeTheta(C.L_SSWM, simParams, 1);

    resultPleiotropicSSWM = simulatePleiotropicSSWM(simParams, genomeParams);
    ave = computeAverageTrajectory(C.numTimeStamp, simParams, resultPleiotropicSSWM.resultTable);

    outDir = ensureDir(resultsRoot, 'SSWM');
    fname  = sprintf('%s/PleiotropicFGM_SSWM_L%d_M%.1e.mat', outDir, ...
                     length(genomeParams.genomeTheta), simParams.mutationRate);
    save(fname, 'simParams', 'genomeParams', 'resultPleiotropicSSWM', 'ave');
    analyticalTrajectories = predictPleiotropicSSWM(simParams, ave);
    save(fname, 'analyticalTrajectories', '-append');

    % ==================== A2) CM asexual (L = 400) =======================
    sectionTimer('PleiotropicFGM CM Asexual');
    simParams = initializeSimParams('numIteration', C.numIteration, ...
                                    'mutationRate', C.mutationRateFast, ...
                                    'omitParams', {'geneticTargetSize'});
    genomeParams = initializeGenomeTheta(C.L_CM, simParams, 1);

    resultPleiotropicCM = simulatePleiotropicCM(simParams, genomeParams);
    ave = computeAverageTrajectory(C.numTimeStamp, simParams, resultPleiotropicCM.resultTable);

    outDir = ensureDir(resultsRoot, 'CM_Asexual');
    fname  = sprintf('%s/PleiotropicFGM_CM_Asexual_L%d_M%.1e.mat', outDir, ...
                     length(genomeParams.genomeTheta), simParams.mutationRate);
    save(fname, 'simParams', 'genomeParams', 'resultPleiotropicCM', 'ave');
    analyticalTrajectories = predictPleiotropicSSWM(simParams, ave);
    save(fname, 'analyticalTrajectories', '-append');

    % ==================== A3) CM sexual (L = 400, r=1) ===================
    sectionTimer('PleiotropicFGM CM Sexual (full recombination)');
    simParams = initializeSimParams('numIteration', C.numIteration, ...
                                    'mutationRate', C.mutationRateFast, ...
                                    'recombinationRate', 1, ...
                                    'omitParams', {'geneticTargetSize'});
    genomeParams = initializeGenomeTheta(C.L_CM, simParams, 1);

    resultPleiotropicCM = simulatePleiotropicCM(simParams, genomeParams);
    ave = computeAverageTrajectory(C.numTimeStamp, simParams, resultPleiotropicCM.resultTable);

    outDir = ensureDir(resultsRoot, 'CM_Sexual');
    fname  = sprintf('%s/PleiotropicFGM_CM_Sexual_L%d_M%.1e_R%.4f.mat', outDir, ...
                     length(genomeParams.genomeTheta), simParams.mutationRate, simParams.recombinationRate);
    save(fname, 'simParams', 'genomeParams', 'resultPleiotropicCM', 'ave');
    analyticalTrajectories = predictPleiotropicSSWM(simParams, ave);
    save(fname, 'analyticalTrajectories', '-append');

    % --------------------------- B) Figure 2 ------------------------------
    fprintf('\n=== Reproducing Figure 2 ===\n');
    makeFigure2('PleiotropicFGM', C.mutationRateSlow, C.mutationRateFast, mode);

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
