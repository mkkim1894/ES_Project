%% Run_supplementary.m
% Standalone driver for supplementary simulations (steady-state CM validation).
%
% Usage:
%   Run_supplementary              % reproduce mode (default)
%   Run_supplementary('test')      % fast test run
%   Run_supplementary('reproduce') % full paper-scale run
%
% Outputs:
%   .mat files in ./results_supplementary/
%
% Reference:
%   Kim, M., Ardell, S. M., & Kryazhimskiy, S. (2025).
%   "Module-Selection Balance in the Evolution of Modular Organisms."

function Run_supplementary(mode)
    if nargin < 1 || isempty(mode)
        mode = 'reproduce';
    end

    isTest = strcmpi(mode, 'test');

    %% Setup paths
    thisDir  = fileparts(mfilename('fullpath'));
    projRoot = fileparts(thisDir);
    addpath(fullfile(projRoot, 'simulation_scripts'));
    addpath(fullfile(projRoot, 'analysis_scripts'));
    addpath(fullfile(projRoot, 'utils'));

    resultsDir = fullfile(projRoot, 'results_supplementary');
    if ~isfolder(resultsDir), mkdir(resultsDir); end

    fprintf('\n========================================================\n');
    fprintf('  SUPPLEMENTARY SIMULATIONS - Module-Selection Balance\n');
    fprintf('========================================================\n');
    fprintf('Mode: %s\n', mode);
    fprintf('Results directory: %s\n\n', resultsDir);

    %% Configuration
    testCfg.N           = 10^3;
    testCfg.generations = 5000;
    testCfg.t_burnin    = 2000;

    reproduceCfg.N           = 10^4;
    reproduceCfg.generations = 15000;
    reproduceCfg.t_burnin    = 10000;

    C = tern(isTest, testCfg, reproduceCfg);

    %% S1. Steady-state CM simulation
    fprintf('--- S1. Steady-State CM Simulation ---\n');
    fprintf('N = %d, generations = %d, burnin = %d\n', C.N, C.generations, C.t_burnin);
    tStart = tic;
    try
        simulateSteadyStateCM(C.N, ...
            'generations', C.generations, ...
            't_burnin',    C.t_burnin, ...
            'savePath',    resultsDir);
    catch ME
        fprintf('ERROR in S1 simulation: %s\n', ME.message);
    end
    fprintf('S1 completed in %.1f s\n\n', toc(tStart));

    fprintf('========================================================\n');
    fprintf('  SUPPLEMENTARY SIMULATIONS COMPLETE\n');
    fprintf('  Run Run_supplementaryFigures to generate figures.\n');
    fprintf('========================================================\n\n');

    matFiles = dir(fullfile(resultsDir, '*.mat'));
    fprintf('Generated data files:\n');
    for i = 1:length(matFiles)
        fprintf('  - %s\n', matFiles(i).name);
    end
    fprintf('\n');
end

function y = tern(cond, a, b)
    if cond; y = a; else; y = b; end
end
