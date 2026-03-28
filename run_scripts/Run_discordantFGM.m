%% Run_discordantFGM.m
% Standalone driver for Discordant Module simulations.
%
% Initial phenotypes are filtered in y-space: only angles whose preimage
% y0 = M^{-1}*x0 satisfies y0_k <= -delta for both chromosomes k are
% retained. This ensures both chromosomes have at least one occupied locus
% (b_k >= 1) at the start, consistent with the binary locus model.
%
% Usage:
%   Run_discordantFGM                    % reproduce mode, all regimes
%   Run_discordantFGM('test')            % test mode, all regimes
%   Run_discordantFGM('reproduce')       % reproduce mode, all regimes
%   Run_discordantFGM('test', {'SSWM'})  % test mode, only SSWM
%   Run_discordantFGM('reproduce', {'CM_Asexual', 'CM_Sexual'})
%
% Regime names:
%   'SSWM', 'CM_Asexual', 'CM_Sexual'
%
% Outputs:
%   .mat files in ./results/Generalization/DiscordantFGM/... (or ./results_test/...)
%
% Reference:
%   Kim, M., Ardell, S. M., & Kryazhimskiy, S. (2025).
%   "Module-Selection Balance in the Evolution of Modular Organisms."

function Run_discordantFGM(mode, regimes)

    if nargin < 1 || isempty(mode)
        mode = 'reproduce';
    end
    if nargin < 2 || isempty(regimes)
        regimes = {'SSWM', 'CM_Asexual', 'CM_Sexual'};
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

    %% Configuration
    K = 10;   % reference genetic target size; defines per-locus rate mu = U_ref / (2*K)
    L = 200;  % number of loci per module (controls deleterious mutation supply)
    % Per-locus mutation rate mu = U_ref / (2*K), independent of L:
    %   SSWM: U_ref = 1e-7 (at K=10), mu = 1e-7 / (2*10) = 5e-9
    %   CM:   U_ref = 2e-4 (at K=10), mu = 2e-4 / (2*10) = 1e-5
    % U is scaled as U = U_ref * (L/K) to keep mu fixed as L changes.

    testCfg.numIteration      = 8;
    testCfg.preRunSteps       = 40;
    testCfg.numTimeStamp      = 20;
    testCfg.initialAngles     = [atan2(1,6.4), atan2(1,3.2), atan2(1,1.6), ...
                             atan2(1,0.8), atan2(1,0.4), atan2(1,0.2)];
    testCfg.mutationRateSlow  = 1e-7 * (L / K);
    testCfg.mutationRateFast  = 2e-4 * (L / K);
    testCfg.popSize           = 10^4;
    testCfg.ellipseRatio      = sqrt(2);
    testCfg.deltaTrait        = 0.1;
    testCfg.landscapeStdDev   = 2;
    testCfg.geneticTargetSize = [L, L];
    testCfg.discordantAngles  = [pi/16, 7*pi/16];

    reproduceCfg.numIteration      = 250;
    reproduceCfg.preRunSteps       = 200;
    reproduceCfg.numTimeStamp      = 20;
    reproduceCfg.initialAngles     = [atan2(1,6.4), atan2(1,3.2), atan2(1,1.6), ...
                                  atan2(1,0.8), atan2(1,0.4), atan2(1,0.2)];
    reproduceCfg.mutationRateSlow  = 1e-7 * (L / K);
    reproduceCfg.mutationRateFast  = 2e-4 * (L / K);
    reproduceCfg.popSize           = 10^4;
    reproduceCfg.ellipseRatio      = sqrt(2);
    reproduceCfg.deltaTrait        = 0.1;
    reproduceCfg.landscapeStdDev   = 2;
    reproduceCfg.geneticTargetSize = [L, L];
    reproduceCfg.discordantAngles  = [pi/16, 7*pi/16];

    C = tern(isTest, testCfg, reproduceCfg);
    fprintf('\n=== Discordant FGM only (%s mode) ===\n', mode);

    %% Filter admissible initial angles in y-space
    % Accept angle theta only if y0 = M^{-1}*x0 satisfies y0_k <= -delta
    % for both chromosomes k=1,2. This ensures b_k >= 1 at the start.
    masterAngles    = C.initialAngles;
    thetaDiscordant = C.discordantAngles;

    M_disc = [cos(thetaDiscordant(1)), cos(thetaDiscordant(2));
              sin(thetaDiscordant(1)), sin(thetaDiscordant(2))];

    % Compute ellipseParams consistent with C.ellipseRatio
    if C.ellipseRatio >= 1
        ellipseParams = [1, 1/C.ellipseRatio];
    else
        ellipseParams = [C.ellipseRatio, 1];
    end

    allInitialPhenotypes = findInitialPhenotypes(masterAngles, ellipseParams, 0.25, C.landscapeStdDev);

    discordantMask = false(1, numel(masterAngles));
    for ii = 1:numel(masterAngles)
        y0 = M_disc \ allInitialPhenotypes(ii, :)';
        discordantMask(ii) = all(y0 <= -C.deltaTrait);
    end

    discordantInitialAngles = masterAngles(discordantMask);
    discordantAngleIdx      = find(discordantMask);

    if isempty(discordantInitialAngles)
        error('No admissible initial angles remain for the chosen discordantAngles and delta.');
    end

    fprintf('DiscordantFGM admissible initial angles: %d / %d retained (y-space filter, delta=%.2f).\n', ...
        numel(discordantInitialAngles), numel(masterAngles), C.deltaTrait);

    %% Regime selection
    runSSWM      = any(strcmpi(regimes, 'SSWM'));
    runCMAsexual = any(strcmpi(regimes, 'CM_Asexual'));
    runCMSexual  = any(strcmpi(regimes, 'CM_Sexual'));

    if ~(runSSWM || runCMAsexual || runCMSexual)
        error('No valid regimes selected. Use ''SSWM'', ''CM_Asexual'', and/or ''CM_Sexual''.');
    end

    %% Discordant SSWM
    if runSSWM
        fprintf('Running DiscordantFGM SSWM...\n');
        tic;
        simParams = initializeSimParams('numIteration', C.numIteration, ...
                                        'initialAngles', discordantInitialAngles, ...
                                        'popSize', C.popSize, ...
                                        'ellipseRatio', C.ellipseRatio, ...
                                        'deltaTrait', C.deltaTrait, ...
                                        'landscapeStdDev', C.landscapeStdDev, ...
                                        'geneticTargetSize', C.geneticTargetSize, ...
                                        'discordantAngles', C.discordantAngles);
        simParams.masterInitialAngles = masterAngles;
        simParams.initialAngleIdx     = discordantAngleIdx;

        resultDiscordantSSWM = simulateDiscordantSSWM(simParams);
        ave = computeAverageTrajectory(C.numTimeStamp, simParams, resultDiscordantSSWM.resultTable);
        fprintf('  Discordant SSWM completed in %.2f s\n', toc);
        reportTermination('Discordant SSWM', resultDiscordantSSWM.terminationStatus);

        outDir = ensureDir(resultsRoot, fullfile('Generalization', 'DiscordantFGM', 'SSWM'));
        fname  = fullfile(outDir, buildFilename('DiscordantFGM', 'SSWM', simParams));
        save(fname, 'simParams', 'resultDiscordantSSWM', 'ave');
    end

    %% Discordant CM asexual
    if runCMAsexual
        fprintf('Running DiscordantFGM CM Asexual...\n');
        tic;
        simParams = initializeSimParams('numIteration', C.numIteration, ...
                                        'initialAngles', discordantInitialAngles, ...
                                        'popSize', C.popSize, ...
                                        'ellipseRatio', C.ellipseRatio, ...
                                        'deltaTrait', C.deltaTrait, ...
                                        'landscapeStdDev', C.landscapeStdDev, ...
                                        'geneticTargetSize', C.geneticTargetSize, ...
                                        'mutationRate', C.mutationRateFast, ...
                                        'discordantAngles', C.discordantAngles);
        simParams.masterInitialAngles = masterAngles;
        simParams.initialAngleIdx     = discordantAngleIdx;

        resultDiscordantCM = simulateDiscordantCM(simParams);
        ave = computeAverageTrajectory(C.numTimeStamp, simParams, resultDiscordantCM.resultTable);
        fprintf('  Discordant CM Asexual completed in %.2f s\n', toc);
        reportTermination('Discordant CM Asexual', resultDiscordantCM.terminationStatus);

        outDir = ensureDir(resultsRoot, fullfile('Generalization', 'DiscordantFGM', 'CM_Asexual'));
        fname  = fullfile(outDir, buildFilename('DiscordantFGM', 'CM_Asexual', simParams));
        save(fname, 'simParams', 'resultDiscordantCM', 'ave');
    end

    %% Discordant CM sexual
    if runCMSexual
        fprintf('Running DiscordantFGM CM Sexual...\n');
        tic;
        simParams = initializeSimParams('numIteration', C.numIteration, ...
                                        'initialAngles', discordantInitialAngles, ...
                                        'popSize', C.popSize, ...
                                        'ellipseRatio', C.ellipseRatio, ...
                                        'deltaTrait', C.deltaTrait, ...
                                        'landscapeStdDev', C.landscapeStdDev, ...
                                        'geneticTargetSize', C.geneticTargetSize, ...
                                        'mutationRate', C.mutationRateFast, ...
                                        'recombinationRate', 1, ...
                                        'discordantAngles', C.discordantAngles);
        simParams.masterInitialAngles = masterAngles;
        simParams.initialAngleIdx     = discordantAngleIdx;

        % No preRunSimulation: discordant y-space pre-run not implemented
        resultDiscordantCM = simulateDiscordantCM(simParams);
        ave = computeAverageTrajectory(C.numTimeStamp, simParams, resultDiscordantCM.resultTable);
        fprintf('  Discordant CM Sexual completed in %.2f s\n', toc);
        reportTermination('Discordant CM Sexual', resultDiscordantCM.terminationStatus);

        outDir = ensureDir(resultsRoot, fullfile('Generalization', 'DiscordantFGM', 'CM_Sexual'));
        fname  = fullfile(outDir, buildFilename('DiscordantFGM', 'CM_Sexual', simParams));
        save(fname, 'simParams', 'resultDiscordantCM', 'ave');
    end

    fprintf('\nDiscordant-only run complete. Mode: %s\n', mode);

end

% ============================== Helpers ==================================

function name = buildFilename(model, regime, sp)
    name = sprintf('%s_%s_N%.0e_M%.1e_d%.2f_eR%.2f_s%.2f', ...
        model, regime, sp.popSize, sp.mutationRate, sp.deltaTrait, ...
        sp.ellipseParams(1)/sp.ellipseParams(2), sp.landscapeStdDev);

    if isfield(sp, 'geneticTargetSize')
        name = sprintf('%s_L%d-%d', name, sp.geneticTargetSize(1), sp.geneticTargetSize(2));
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
    ts     = statusMatrix(:);
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
    if cond; y = a; else; y = b; end
end
