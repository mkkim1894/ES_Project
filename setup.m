%% setup.m
% Run this script from the project root to add all paths and verify setup.
%
% Usage:
%   cd /path/to/project
%   setup

function setup()
    % Get the directory where setup.m lives (project root)
    projRoot = fileparts(mfilename('fullpath'));
    
    % Add all subdirectories
    addpath(fullfile(projRoot, 'main'));
    addpath(fullfile(projRoot, 'simulation'));
    addpath(fullfile(projRoot, 'analysis'));
    addpath(fullfile(projRoot, 'utils'));
    addpath(fullfile(projRoot, 'figures'));
    
    fprintf('=== FGM Project Setup ===\n');
    fprintf('Project root: %s\n\n', projRoot);
    
    % Verify required functions exist
    requiredFunctions = {
        % Main drivers
        'Run_modularFGM', 'Run_pleiotropicFGM', 'Run_standardFGM', ...
        % Simulations
        'simulateModularCM', 'simulateModularSSWM', ...
        'simulateNestedCM', 'simulateNestedSSWM', ...
        'simulatePleiotropicCM', 'simulatePleiotropicSSWM', ...
        'simulateStandardCM', 'simulateStandardSSWM', ...
        % Analysis
        'computeAverageTrajectory', 'predictFullRecomb', ...
        'predictModularCM', 'predictModularSSWM', 'predictPleiotropicSSWM', ...
        % Utils
        'findInitialPhenotypes', 'freezeParam', 'initializeGenomeTheta', ...
        'initializeSimParams', 'preRunSimulation', ...
        % Figures
        'makeFigure2', 'makeFigure3', 'makeFigure4', 'makeFigure5', 'makeFigure6', ...
        'makeFigure_StandardFGM'
    };
    
    fprintf('Checking required functions...\n');
    allFound = true;
    missingFuncs = {};
    for i = 1:length(requiredFunctions)
        fn = requiredFunctions{i};
        if exist(fn, 'file') == 2
            fprintf('  [OK] %s\n', fn);
        else
            fprintf('  [MISSING] %s\n', fn);
            allFound = false;
            missingFuncs{end+1} = fn;
        end
    end
    
    fprintf('\n');
    if allFound
        fprintf('All %d functions found. Setup complete!\n\n', length(requiredFunctions));
        fprintf('To run a demo simulation:\n');
        fprintf('  >> cd main\n');
        fprintf('  >> Run_modularFGM(''demo'')\n\n');
        fprintf('To run full simulations:\n');
        fprintf('  >> Run_modularFGM(''full'')\n');
    else
        fprintf('WARNING: %d functions missing:\n', length(missingFuncs));
        for i = 1:length(missingFuncs)
            fprintf('  - %s\n', missingFuncs{i});
        end
        fprintf('\nCheck your installation.\n');
    end
end
