%% Run_nestedFGM.m
% Standalone driver for Nested FGM simulations.
%
% Usage:
%   Run_nestedFGM                                % reproduce, all regimes, [10,10]
%   Run_nestedFGM('test')                        % test mode, all regimes, [10,10]
%   Run_nestedFGM('reproduce')                   % reproduce, all regimes, [10,10]
%   Run_nestedFGM('reproduce', {'SSWM'})         % reproduce, SSWM only, [10,10]
%   Run_nestedFGM('reproduce', {}, [10,20])      % reproduce, all regimes, [10,20]
%   Run_nestedFGM('reproduce', {'CM_Asexual'}, [10,20])
%
% Regime names:
%   'SSWM', 'CM_Asexual', 'CM_Sexual'
%
% Outputs:
%   .mat files in ./results/Generalization/NestedFGM/... (or ./results_test/...)
%
% Reference:
%   Kim, M., Ardell, S. M., & Kryazhimskiy, S. (2025).
%   "Module-Selection Balance in the Evolution of Modular Organisms."

function Run_nestedFGM(mode, regimes, moduleDimension)

    if nargin < 1 || isempty(mode)
        mode = 'reproduce';
    end
    if nargin < 2 || isempty(regimes)
        regimes = {'SSWM', 'CM_Asexual', 'CM_Sexual'};
    end
    if nargin < 3 || isempty(moduleDimension)
        moduleDimension = [10, 10];
    end

    isTest = strcmpi(mode, 'test');

    %% Environment setup
    thisDir  = fileparts(mfilename('fullpath'));
    projRoot = fileparts(thisDir);
    addpath(fullfile(projRoot, 'simulation_scripts'));
    addpath(fullfile(projRoot, 'analysis_scripts'));
    addpath(fullfile(projRoot, 'utils'));
    addpath(fullfile(projRoot, 'figure_scripts'));

    fprintf('Project root: %s\n', projRoot);

    resultsRoot = tern(isTest, fullfile(projRoot, 'results_test'), ...
                               fullfile(projRoot, 'results'));

    %% Configuration — identical to master runner
    testCfg.numIteration      = 8;
    testCfg.preRunSteps       = 40;
    testCfg.numTimeStamp      = 20;
    testCfg.initialAngles     = [atan2(1,6.4), atan2(1,3.2), atan2(1,1.6), ...
                                 atan2(1,0.8), atan2(1,0.4), atan2(1,0.2)];
    testCfg.mutationRateSlow  = 1e-7;
    testCfg.mutationRateFast  = 2e-4;
    testCfg.popSize           = 10^4;
    testCfg.ellipseRatio      = sqrt(2);
    testCfg.deltaTrait        = 0.1;
    testCfg.landscapeStdDev   = 2;
    testCfg.moduleDimension   = [10, 10];

    reproduceCfg.numIteration      = 250;
    reproduceCfg.preRunSteps       = 200;
    reproduceCfg.numTimeStamp      = 20;
    reproduceCfg.initialAngles     = [atan2(1,6.4), atan2(1,3.2), atan2(1,1.6), ...
                                      atan2(1,0.8), atan2(1,0.4), atan2(1,0.2)];
    reproduceCfg.mutationRateSlow  = 1e-7;
    reproduceCfg.mutationRateFast  = 2e-4;
    reproduceCfg.popSize           = 10^4;
    reproduceCfg.ellipseRatio      = sqrt(2);
    reproduceCfg.deltaTrait        = 0.1;
    reproduceCfg.landscapeStdDev   = 2;
    reproduceCfg.moduleDimension   = [10, 10];

    C = tern(isTest, testCfg, reproduceCfg);

    % Override moduleDimension with argument if provided
    C.moduleDimension = moduleDimension;

    % Calibrate nested mutation-vector magnitude m so that mean |delta_i| = delta
    % i.e. m^2/2 = delta  =>  m = sqrt(2*delta)
    nestedM = sqrt(2 * C.deltaTrait);

    fprintf('\n=== Nested FGM only (%s mode) ===\n', mode);
    fprintf('deltaTrait = %.4f, nestedM = %.4f, moduleDimension = [%d, %d]\n', ...
        C.deltaTrait, nestedM, C.moduleDimension(1), C.moduleDimension(2));

    %% Regime selection
    runSSWM      = any(strcmpi(regimes, 'SSWM'));
    runCMAsexual = any(strcmpi(regimes, 'CM_Asexual'));
    runCMSexual  = any(strcmpi(regimes, 'CM_Sexual'));

    if ~(runSSWM || runCMAsexual || runCMSexual)
        error('No valid regimes selected. Use ''SSWM'', ''CM_Asexual'', and/or ''CM_Sexual''.');
    end

    %% Nested SSWM
    if runSSWM
        fprintf('Running NestedFGM SSWM...\n');
        tic;
        simParams = initializeSimParams('numIteration', C.numIteration, ...
                                        'initialAngles', C.initialAngles, ...
                                        'popSize', C.popSize, ...
                                        'ellipseRatio', C.ellipseRatio, ...
                                        'deltaTrait', nestedM, ...
                                        'landscapeStdDev', C.landscapeStdDev, ...
                                        'mutationRate', C.mutationRateSlow);
        simParams.moduleDimension = C.moduleDimension;

        resultNestedSSWM = simulateNestedSSWM(simParams);
        ave = computeAverageTrajectory(C.numTimeStamp, simParams, resultNestedSSWM.resultTable);
        fprintf('  NestedFGM SSWM completed in %.2f s\n', toc);
        reportTermination('NestedFGM SSWM', resultNestedSSWM.terminationStatus);

        outDir = ensureDir(resultsRoot, fullfile('Generalization', 'NestedFGM', 'SSWM'));
        fname = fullfile(outDir, buildFilename('NestedFGM', 'SSWM', simParams));
        save(fname, 'simParams', 'resultNestedSSWM', 'ave');
    end

    %% Nested CM asexual
    if runCMAsexual
        fprintf('Running NestedFGM CM Asexual...\n');
        tic;
        simParams = initializeSimParams('numIteration', C.numIteration, ...
                                        'initialAngles', C.initialAngles, ...
                                        'popSize', C.popSize, ...
                                        'ellipseRatio', C.ellipseRatio, ...
                                        'deltaTrait', nestedM, ...
                                        'landscapeStdDev', C.landscapeStdDev, ...
                                        'mutationRate', C.mutationRateFast);
        simParams.moduleDimension = C.moduleDimension;

        resultNestedCM = simulateNestedCM(simParams);
        ave = computeAverageTrajectory(C.numTimeStamp, simParams, resultNestedCM.resultTable);
        fprintf('  NestedFGM CM Asexual completed in %.2f s\n', toc);
        reportTermination('NestedFGM CM Asexual', resultNestedCM.terminationStatus);

        outDir = ensureDir(resultsRoot, fullfile('Generalization', 'NestedFGM', 'CM_Asexual'));
        fname = fullfile(outDir, buildFilename('NestedFGM', 'CM_Asexual', simParams));
        save(fname, 'simParams', 'resultNestedCM', 'ave');
    end

    %% Nested CM sexual
    if runCMSexual
        fprintf('Running NestedFGM CM Sexual...\n');
        tic;
        simParams = initializeSimParams('numIteration', C.numIteration, ...
                                        'initialAngles', C.initialAngles, ...
                                        'popSize', C.popSize, ...
                                        'ellipseRatio', C.ellipseRatio, ...
                                        'deltaTrait', nestedM, ...
                                        'landscapeStdDev', C.landscapeStdDev, ...
                                        'mutationRate', C.mutationRateFast, ...
                                        'recombinationRate', 1);
        simParams.moduleDimension = C.moduleDimension;

        % No preRunSimulation for nested FGM (not implemented)
        resultNestedCM = simulateNestedCM(simParams);
        ave = computeAverageTrajectory(C.numTimeStamp, simParams, resultNestedCM.resultTable);
        fprintf('  NestedFGM CM Sexual completed in %.2f s\n', toc);
        reportTermination('NestedFGM CM Sexual', resultNestedCM.terminationStatus);

        outDir = ensureDir(resultsRoot, fullfile('Generalization', 'NestedFGM', 'CM_Sexual'));
        fname = fullfile(outDir, buildFilename('NestedFGM', 'CM_Sexual', simParams));
        save(fname, 'simParams', 'resultNestedCM', 'ave');
    end

    fprintf('\nNested FGM run complete. Mode: %s\n', mode);

end

% ============================== Helpers ==================================

function name = buildFilename(model, regime, sp)
    name = sprintf('%s_%s_N%.0e_M%.1e_d%.2f_eR%.2f_s%.2f', ...
        model, regime, sp.popSize, sp.mutationRate, sp.deltaTrait, ...
        sp.ellipseParams(1)/sp.ellipseParams(2), sp.landscapeStdDev);

    if isfield(sp, 'geneticTargetSize')
        name = sprintf('%s_K%d-%d', name, sp.geneticTargetSize(1), sp.geneticTargetSize(2));
    end

    if isfield(sp, 'moduleDimension')
        name = sprintf('%s_n%d-%d', name, sp.moduleDimension(1), sp.moduleDimension(2));
    end

    if isfield(sp, 'recombinationRate') && sp.recombinationRate > 0
        name = sprintf('%s_R%.4f', name, sp.recombinationRate);
    end

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

function y = tern(cond, a, b)
    if cond
        y = a;
    else
        y = b;
    end
end
