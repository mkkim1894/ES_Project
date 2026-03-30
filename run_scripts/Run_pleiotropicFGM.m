%% Run_pleiotropicFGM.m
% Standalone driver for Pleiotropic FGM simulations.
%
% Usage:
%   Run_pleiotropicFGM              % reproduce mode (default)
%   Run_pleiotropicFGM('test')      % fast sanity run
%   Run_pleiotropicFGM('reproduce') % full paper-scale run
%
% Outputs:
%   .mat files in ./results/... (or ./results_test/...)
%
% Reference:
%   Kim, M., Ardell, S. M., & Kryazhimskiy, S. (2025).
%   "Module-Selection Balance in the Evolution of Modular Organisms."

function Run_pleiotropicFGM(mode)
    if nargin < 1, mode = 'reproduce'; end
    isTest = strcmpi(mode, 'test');

    %% Environment setup
    thisDir  = fileparts(mfilename('fullpath'));
    projRoot = fileparts(thisDir);
    addpath(fullfile(projRoot, 'simulation_scripts'));
    addpath(fullfile(projRoot, 'analysis_scripts'));
    addpath(fullfile(projRoot, 'utils'));
    addpath(fullfile(projRoot, 'figure_scripts'));
    fprintf('Project root: %s\n', projRoot);

    resultsRoot = tern(isTest, fullfile(projRoot, 'results_test'), fullfile(projRoot, 'results'));

    % ------------------ Config knobs ------------------
    K = 10;   % reference genetic target size; defines per-locus rate mu = U_ref / (2*K)
    L = 200;  % number of loci per module (total loci = 2*L)
    % Per-locus mutation rate mu = U_ref / (2*K), independent of L:
    %   SSWM: U_ref = 1e-7 (at K=10), mu = 1e-7 / (2*10) = 5e-9
    %   CM:   U_ref = 2e-4 (at K=10), mu = 2e-4 / (2*10) = 1e-5
    % U is scaled as U = U_ref * (L/K) to keep mu fixed as L changes.

    testCfg.numIteration      = 8;
    testCfg.numTimeStamp      = 20;
    testCfg.initialAngles     = [atan2(1,6.4), atan2(1,0.2)];
    testCfg.mutationRateSlow  = 1e-7 * (L / K);
    testCfg.mutationRateFast  = 2e-4 * (L / K);
    testCfg.L_SSWM            = 2*L;
    testCfg.L_CM              = 2*L;
    testCfg.popSize           = 10^4;
    testCfg.ellipseRatio      = sqrt(2);
    testCfg.deltaTrait        = 0.1;
    testCfg.landscapeStdDev   = 2;

    reproduceCfg.numIteration      = 250;
    reproduceCfg.numTimeStamp      = 20;
    reproduceCfg.initialAngles     = [atan2(1,6.4), atan2(1,3.2), atan2(1,1.6), ...
                                atan2(1,0.8), atan2(1,0.4), atan2(1,0.2)];
    reproduceCfg.mutationRateSlow  = 1e-7 * (L / K);
    reproduceCfg.mutationRateFast  = 2e-4 * (L / K);
    reproduceCfg.L_SSWM            = 2*L;
    reproduceCfg.L_CM              = 2*L;
    reproduceCfg.popSize           = 10^4;
    reproduceCfg.ellipseRatio      = sqrt(2);
    reproduceCfg.deltaTrait        = 0.1;
    reproduceCfg.landscapeStdDev   = 2;

    C = tern(isTest, testCfg, reproduceCfg);

    % Proximity cutoff for log-ratio plots (panels D-F).
    % Set to 0 to disable cutoff.
    proximityCutoff = 0.1;

    fprintf('\n=== Pleiotropic FGM (%s mode) ===\n', mode);

    % ======================= SSWM (L = 400) =========================
    fprintf('Running PleiotropicFGM SSWM...\n');
    tic;
    simParams = initializeSimParams('numIteration', C.numIteration, ...
                                'initialAngles', C.initialAngles, ...
                                'popSize', C.popSize, ...
                                'ellipseRatio', C.ellipseRatio, ...
                                'deltaTrait', C.deltaTrait, ...
                                'landscapeStdDev', C.landscapeStdDev, ...
                                'mutationRate', C.mutationRateSlow, ...
                                'omitParams', {'geneticTargetSize'});
    genomeParams = initializeGenomeTheta(C.L_SSWM, simParams, 1);

    resultPleiotropicSSWM = simulatePleiotropicSSWM(simParams, genomeParams);
    ave = computeAverageTrajectory(C.numTimeStamp, simParams, resultPleiotropicSSWM.resultTable);
    fprintf('  SSWM completed in %.2f s\n', toc);
    reportTermination('SSWM', resultPleiotropicSSWM.terminationStatus);

    outDir = ensureDir(resultsRoot, 'SSWM');
    fname = fullfile(outDir, buildFilename('PleiotropicFGM', 'SSWM', simParams, C.L_SSWM));
    save(fname, 'simParams', 'genomeParams', 'resultPleiotropicSSWM', 'ave');
    analyticalTrajectories = predictPleiotropicSSWM(simParams, ave, genomeParams);
    save(fname, 'analyticalTrajectories', '-append');

    % ==================== CM asexual (L = 400) =======================
    fprintf('Running PleiotropicFGM CM Asexual...\n');
    tic;
    simParams = initializeSimParams('numIteration', C.numIteration, ...
                                    'initialAngles', C.initialAngles, ...
                                    'popSize', C.popSize, ...
                                    'ellipseRatio', C.ellipseRatio, ...
                                    'deltaTrait', C.deltaTrait, ...
                                    'landscapeStdDev', C.landscapeStdDev, ...
                                    'mutationRate', C.mutationRateFast, ...
                                    'omitParams', {'geneticTargetSize'});
    genomeParams = initializeGenomeTheta(C.L_CM, simParams, 1);

    resultPleiotropicCM = simulatePleiotropicCM(simParams, genomeParams);
    ave = computeAverageTrajectory(C.numTimeStamp, simParams, resultPleiotropicCM.resultTable);
    fprintf('  CM Asexual completed in %.2f s\n', toc);
    reportTermination('CM Asexual', resultPleiotropicCM.terminationStatus);

    outDir = ensureDir(resultsRoot, 'CM_Asexual');
    fname = fullfile(outDir, buildFilename('PleiotropicFGM', 'CM_Asexual', simParams, C.L_CM));
    save(fname, 'simParams', 'genomeParams', 'resultPleiotropicCM', 'ave');
    % SSWM prediction shown in all panels (A-C) to illustrate little
    % variations across the regimes
    analyticalTrajectories = predictPleiotropicSSWM(simParams, ave, genomeParams);
    save(fname, 'analyticalTrajectories', '-append');

    % ==================== CM sexual (L = 400, r=1) ===================
    fprintf('Running PleiotropicFGM CM Sexual (full recombination)...\n');
    tic;
    simParams = initializeSimParams('numIteration', C.numIteration, ...
                                    'initialAngles', C.initialAngles, ...
                                    'popSize', C.popSize, ...
                                    'ellipseRatio', C.ellipseRatio, ...
                                    'deltaTrait', C.deltaTrait, ...
                                    'landscapeStdDev', C.landscapeStdDev, ...
                                    'mutationRate', C.mutationRateFast, ...
                                    'recombinationRate', 1, ...
                                    'omitParams', {'geneticTargetSize'});
    genomeParams = initializeGenomeTheta(C.L_CM, simParams, 1);

    resultPleiotropicCM = simulatePleiotropicCM(simParams, genomeParams);
    ave = computeAverageTrajectory(C.numTimeStamp, simParams, resultPleiotropicCM.resultTable);
    fprintf('  CM Sexual completed in %.2f s\n', toc);
    reportTermination('CM Sexual', resultPleiotropicCM.terminationStatus);

    outDir = ensureDir(resultsRoot, 'CM_Sexual');
    fname = fullfile(outDir, buildFilename('PleiotropicFGM', 'CM_Sexual', simParams, C.L_CM));
    save(fname, 'simParams', 'genomeParams', 'resultPleiotropicCM', 'ave');
    % SSWM prediction shown in all panels (A-C) to illustrate little
    % variations across the regimes
    analyticalTrajectories = predictPleiotropicSSWM(simParams, ave, genomeParams);
    save(fname, 'analyticalTrajectories', '-append');

    fprintf('\nAll done. Mode: %s\n', mode);
end

% ============================== Helpers ==================================

function name = buildFilename(model, regime, sp, L)
    % Build a parameter-encoded filename
    % Format: Model_Regime_N1e4_L400_M1e-7_d0.10_eR1.41_s2.00_R0.0000.mat
    name = sprintf('%s_%s_N%.0e_L%d_M%.1e_d%.2f_eR%.2f_s%.2f', ...
        model, regime, sp.popSize, L, sp.mutationRate, sp.deltaTrait, ...
        sp.ellipseParams(1)/sp.ellipseParams(2), sp.landscapeStdDev);
    if isfield(sp, 'recombinationRate') && sp.recombinationRate > 0
        name = sprintf('%s_R%.4f', name, sp.recombinationRate);
    end
    name = [name, '.mat'];
end

function reportTermination(label, statusMatrix)
    ts = statusMatrix(:);
    nTotal = numel(ts);
    fprintf('  %s termination: %d/%d reached fitness, %d no beneficial muts\n', ...
        label, sum(ts == 1), nTotal, sum(ts == 0));
end

function out = ensureDir(base, rel)
    out = fullfile(base, rel);
    if ~exist(out, 'dir'), mkdir(out); end
end

function y = tern(cond, a, b)
    if cond, y = a; else, y = b; end
end