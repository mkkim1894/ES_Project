%% Run_modularFGM.m
% Standalone driver for Modular FGM simulations.
%
% Usage:
%   Run_modularFGM              % reproduce mode (default)
%   Run_modularFGM('test')      % fast sanity run
%   Run_modularFGM('reproduce') % full paper-scale run
%
% Outputs:
%   .mat files in ./results/... (or ./results_test/...)
%
% Reference:
%   Kim, M., Ardell, S. M., & Kryazhimskiy, S. (2025).
%   "Module-Selection Balance in the Evolution of Modular Organisms."

function Run_modularFGM(mode)
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
    testCfg.numIteration      = 8;
    testCfg.preRunSteps       = 40;
    testCfg.numTimeStamp      = 20;
    testCfg.recombList        = [0, 0.01, 1];
    testCfg.initialAngles     = [atan2(1,6.4), atan2(1,3.2), atan2(1,1.6), ...
                                 atan2(1,0.8), atan2(1,0.4), atan2(1,0.2)];
    testCfg.mutationRateSlow  = 1e-7;
    testCfg.mutationRateFast  = 2e-4;
    testCfg.popSize           = 10^4;
    testCfg.ellipseRatio      = sqrt(2);
    testCfg.deltaTrait        = 0.1;
    testCfg.landscapeStdDev   = 2;
    testCfg.geneticTargetSize = [10, 10];

    reproduceCfg.numIteration      = 250;
    reproduceCfg.preRunSteps       = 200;
    reproduceCfg.numTimeStamp      = 20;
    reproduceCfg.recombList        = [0, 2e-4, 1e-3, 1e-2, 1e-1, 1];
    reproduceCfg.initialAngles     = [atan2(1,6.4), atan2(1,3.2), atan2(1,1.6), ...
                                      atan2(1,0.8), atan2(1,0.4), atan2(1,0.2)];
    reproduceCfg.mutationRateSlow  = 1e-7;
    reproduceCfg.mutationRateFast  = 2e-4;
    reproduceCfg.popSize           = 10^4;
    reproduceCfg.ellipseRatio      = sqrt(2);
    reproduceCfg.deltaTrait        = 0.1;
    reproduceCfg.landscapeStdDev   = 2;
    reproduceCfg.geneticTargetSize = [10, 10];

    C = tern(isTest, testCfg, reproduceCfg);

    fprintf('\n=== Modular FGM (%s mode) ===\n', mode);

    % SSWM
    fprintf('Running ModularFGM SSWM...\n');
    tic;
    simParams = initializeSimParams('numIteration', C.numIteration, ...
                                    'initialAngles', C.initialAngles, ...
                                    'popSize', C.popSize, ...
                                    'ellipseRatio', C.ellipseRatio, ...
                                    'deltaTrait', C.deltaTrait, ...
                                    'landscapeStdDev', C.landscapeStdDev, ...
                                    'geneticTargetSize', C.geneticTargetSize);
    resultModularSSWM = simulateModularSSWM(simParams);
    ave = computeAverageTrajectory(C.numTimeStamp, simParams, resultModularSSWM.resultTable);
    fprintf('  SSWM completed in %.2f s\n', toc);
    reportTermination('SSWM', resultModularSSWM.terminationStatus);

    outDir = ensureDir(resultsRoot, 'SSWM');
    fname = fullfile(outDir, buildFilename('ModularFGM', 'SSWM', simParams));
    save(fname, 'simParams', 'resultModularSSWM', 'ave');
    analyticalTrajectories = predictModularSSWM(discretizeInitialPhenotypes(simParams), ave);
    save(fname, 'analyticalTrajectories', '-append');

    % CM asexual
    fprintf('Running ModularFGM CM Asexual...\n');
    tic;
    simParams = initializeSimParams('numIteration', C.numIteration, ...
                                    'initialAngles', C.initialAngles, ...
                                    'popSize', C.popSize, ...
                                    'ellipseRatio', C.ellipseRatio, ...
                                    'deltaTrait', C.deltaTrait, ...
                                    'landscapeStdDev', C.landscapeStdDev, ...
                                    'geneticTargetSize', C.geneticTargetSize, ...
                                    'mutationRate', C.mutationRateFast);
    resultModularCM = simulateModularCM(simParams);
    ave = computeAverageTrajectory(C.numTimeStamp, simParams, resultModularCM.resultTable);
    fprintf('  CM Asexual completed in %.2f s\n', toc);
    reportTermination('CM Asexual', resultModularCM.terminationStatus);

    outDir = ensureDir(resultsRoot, 'CM_Asexual');
    fname = fullfile(outDir, buildFilename('ModularFGM', 'CM_Asexual', simParams));
    save(fname, 'simParams', 'resultModularCM', 'ave');
    analyticalTrajectories = predictModularCM(discretizeInitialPhenotypes(simParams));
    save(fname, 'analyticalTrajectories', '-append');

    % CM sexual (full recombination)
    fprintf('Running ModularFGM CM Sexual (full recombination)...\n');
    tic;
    simParams = initializeSimParams('numIteration', C.numIteration, ...
                                    'initialAngles', C.initialAngles, ...
                                    'popSize', C.popSize, ...
                                    'ellipseRatio', C.ellipseRatio, ...
                                    'deltaTrait', C.deltaTrait, ...
                                    'landscapeStdDev', C.landscapeStdDev, ...
                                    'geneticTargetSize', C.geneticTargetSize, ...
                                    'mutationRate', C.mutationRateFast, ...
                                    'recombinationRate', 1);
    initialParameters = freezeParam(simParams);
    simParams.populationMatrices = preRunSimulation(initialParameters, simParams, C.preRunSteps, 'rngSeed', 1);
    resultModularCM = simulateModularCM(simParams);
    ave = computeAverageTrajectory(C.numTimeStamp, simParams, resultModularCM.resultTable);
    fprintf('  CM Sexual completed in %.2f s\n', toc);
    reportTermination('CM Sexual', resultModularCM.terminationStatus);

    outDir = ensureDir(resultsRoot, 'CM_Sexual');
    fname = fullfile(outDir, buildFilename('ModularFGM', 'CM_Sexual', simParams));
    save(fname, 'simParams', 'resultModularCM', 'ave');
    analyticalTrajectories = predictFullRecomb(discretizeInitialPhenotypes(simParams), ave);
    save(fname, 'analyticalTrajectories', '-append');

    fprintf('\nAll done. Mode: %s\n', mode);

end

% ============================== Helpers ==================================

function name = buildFilename(model, regime, sp)
    % Build parameter-encoded filename
    % Format: Model_Regime_N1e4_M1e-7_d0.10_eR1.41_s2.00[_K10-10][_n10-10][_R0.0000].mat

    name = sprintf('%s_%s_N%.0e_M%.1e_d%.2f_eR%.2f_s%.2f', ...
        model, regime, sp.popSize, sp.mutationRate, sp.deltaTrait, ...
        sp.ellipseParams(1)/sp.ellipseParams(2), sp.landscapeStdDev);

    % ModularFGM parameter
    if isfield(sp, 'geneticTargetSize')
        name = sprintf('%s_K%d-%d', name, sp.geneticTargetSize(1), sp.geneticTargetSize(2));
    end

    % NestedFGM parameter
    if isfield(sp, 'moduleDimension')
        name = sprintf('%s_n%d-%d', name, sp.moduleDimension(1), sp.moduleDimension(2));
    end

    % Recombination
    if isfield(sp, 'recombinationRate') && sp.recombinationRate > 0
        name = sprintf('%s_R%.4f', name, sp.recombinationRate);
    end

    % Discordant environment angles
    if isfield(sp, 'discordantAngles') && ~isempty(sp.discordantAngles)
        name = sprintf('%s_T%.4f-%.4f', name, sp.discordantAngles(1), sp.discordantAngles(2));
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
    if ~exist(out, 'dir')
        mkdir(out);
    end
end

function sp = discretizeInitialPhenotypes(sp)
    delta = sp.deltaTrait;
    x0    = sp.initialPhenotypes;
    x0    = -delta * round(-x0 / delta);
    x0    = min(x0, 0);
    sp.initialPhenotypes = x0;
end

function y = tern(cond, a, b)
    if cond
        y = a;
    else
        y = b;
    end
end