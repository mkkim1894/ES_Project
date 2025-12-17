%% Run_supplementary.m
% Master driver for all supplementary analyses and figures.
%
% Description:
%   Runs all supplementary simulations, analyses, and figure generation for
%   the modular FGM paper. Organizes results into the supplementary/ directory.
%
% Usage:
%   Run_supplementary              % Full run (all analyses)
%   Run_supplementary('demo')      % Quick demo run (reduced parameters)
%   Run_supplementary('figures')   % Generate figures only (requires existing data)
%
% Outputs:
%   Data files saved to: ./results_supplementary/
%   Figures saved to:    ./results_supplementary/Figures/
%
% Supplementary Analyses:
%   S1. Steady-state CM validation (simulateSteadyStateCM)
%   S2. Steady-state recombination/variance (simulateSteadyStateRecomb)
%   S3. Threshold D sensitivity analysis (Run_ThresholdDAnalysis)
%   S4. Nested FGM distribution visualization
%
% Reference:
%   Kim, M., Ardell, S. M., & Kryazhimskiy, S. (2025).
%   "Module-Selection Balance in the Evolution of Modular Organisms."
%
% Copyright (c) 2025 Minkyu Kim, Cornell University
% Licensed under MIT License

function Run_supplementary(mode)
    if nargin < 1
        mode = 'full';
    end
    
    isDemo = strcmpi(mode, 'demo');
    figuresOnly = strcmpi(mode, 'figures');
    
    %% Setup paths
    thisDir = fileparts(mfilename('fullpath'));
    addpath(genpath(fileparts(thisDir)));
    
    resultsDir = fullfile(thisDir, '../results_supplementary');
    figuresDir = fullfile(resultsDir, 'Figures');
    
    if ~isfolder(resultsDir), mkdir(resultsDir); end
    if ~isfolder(figuresDir), mkdir(figuresDir); end
    
    fprintf('\n');
    fprintf('========================================================\n');
    fprintf('  SUPPLEMENTARY ANALYSES - Module-Selection Balance\n');
    fprintf('========================================================\n');
    fprintf('Mode: %s\n', mode);
    fprintf('Results directory: %s\n', resultsDir);
    fprintf('Figures directory: %s\n', figuresDir);
    fprintf('\n');
    
    %% Configuration
    if isDemo
        cfg.N = 10^3;                    % Reduced population size
        cfg.generations = 5000;           % Reduced generations
        cfg.t_burnin = 2000;
        cfg.numRepeat = 10;               % Reduced replicates
        cfg.D_values = [1000, 100];       % Reduced D values
    else
        cfg.N = 10^4;
        cfg.generations = 15000;
        cfg.t_burnin = 10000;
        cfg.numRepeat = 100;
        cfg.D_values = [10000, 1000, 100, 10];
    end
    
    %% ================================================================
    %  S1. STEADY-STATE CM VALIDATION
    %  ================================================================
    if ~figuresOnly
        fprintf('\n--- S1. Steady-State CM Validation ---\n');
        sectionTimer('S1 start');
        
        try
            % simulateSteadyStateCM saves the file directly with correct variable names
            % when savePath is provided, so we don't need to save again
            [~, ~, ~, ~, ~, ~] = simulateSteadyStateCM(cfg.N, ...
                'generations', cfg.generations, ...
                't_burnin', cfg.t_burnin, ...
                'savePath', resultsDir);
            
        catch ME
            fprintf('Warning: S1 failed - %s\n', ME.message);
        end
        
        sectionTimer('S1 complete');
    end
    
    % Generate S1 figures
    fprintf('\nGenerating S1 figures...\n');
    dataFile = fullfile(resultsDir, sprintf('SteadyStateCM_N_%d.mat', cfg.N));
    if isfile(dataFile)
        try
            makeFigureS_SteadyStateCM(dataFile, 'outputDir', figuresDir);
        catch ME
            fprintf('Warning: S1 figures failed - %s\n', ME.message);
        end
    else
        fprintf('Skipping S1 figures - data file not found\n');
    end
    
    %% ================================================================
    %  S2. STEADY-STATE RECOMBINATION/VARIANCE ANALYSIS
    %  ================================================================
    if ~figuresOnly
        fprintf('\n--- S2. Steady-State Recombination Analysis ---\n');
        sectionTimer('S2 start');
        
        try
            % Run three conditions
            fprintf('  Running two-trait model...\n');
            [~, ~] = simulateSteadyStateRecomb('both', ...
                'numRepeat', cfg.numRepeat, ...
                'generations', cfg.generations, ...
                't_burnin', cfg.t_burnin, ...
                'savePath', resultsDir);
            
            fprintf('  Running trait 1 only...\n');
            [~, ~] = simulateSteadyStateRecomb('trait1', ...
                'numRepeat', cfg.numRepeat, ...
                'generations', cfg.generations, ...
                't_burnin', cfg.t_burnin, ...
                'savePath', resultsDir);
            
            fprintf('  Running trait 2 only...\n');
            [~, ~] = simulateSteadyStateRecomb('trait2', ...
                'numRepeat', cfg.numRepeat, ...
                'generations', cfg.generations, ...
                't_burnin', cfg.t_burnin, ...
                'savePath', resultsDir);
            
        catch ME
            fprintf('Warning: S2 failed - %s\n', ME.message);
        end
        
        sectionTimer('S2 complete');
    end
    
    % Generate S2 figures
    fprintf('\nGenerating S2 figures...\n');
    dataFiles.both = fullfile(resultsDir, 'SteadyStateRecomb_TwoTraits.mat');
    dataFiles.trait1 = fullfile(resultsDir, 'SteadyStateRecomb_Trait1.mat');
    dataFiles.trait2 = fullfile(resultsDir, 'SteadyStateRecomb_Trait2.mat');
    
    if isfile(dataFiles.both) && isfile(dataFiles.trait1) && isfile(dataFiles.trait2)
        try
            makeFigureS_SteadyStateRecomb(dataFiles, 'outputDir', figuresDir);
        catch ME
            fprintf('Warning: S2 figures failed - %s\n', ME.message);
        end
    else
        fprintf('Skipping S2 figures - data files not found\n');
    end
    
    %% ================================================================
    %  S3. THRESHOLD D SENSITIVITY ANALYSIS
    %  ================================================================
    if ~figuresOnly
        fprintf('\n--- S3. Threshold D Sensitivity Analysis ---\n');
        sectionTimer('S3 start');
        
        try
            Run_ThresholdDAnalysis('D_values', cfg.D_values, ...
                                   'outputDir', resultsDir, ...
                                   'outputFile', 'ThresholdD_Analysis.mat');
        catch ME
            fprintf('Warning: S3 failed - %s\n', ME.message);
        end
        
        sectionTimer('S3 complete');
    end
    
    % Generate S3 figures (requires simulation data from Run_modularFGM)
    fprintf('\nGenerating S3 figures...\n');
    theoryFile = fullfile(resultsDir, 'ThresholdD_Analysis.mat');
    
    % Look for simulation data in multiple possible locations
    possibleSimFiles = {
        fullfile(resultsDir, 'DynSim_2.mat'), ...
        fullfile(thisDir, '../results/CM_Asexual/DynSim_2.mat'), ...
        fullfile(thisDir, '../results_demo/CM_Asexual/DynSim_2.mat')
    };
    
    simFile = '';
    for i = 1:length(possibleSimFiles)
        if isfile(possibleSimFiles{i})
            simFile = possibleSimFiles{i};
            break;
        end
    end
    
    if isfile(theoryFile)
        if ~isempty(simFile)
            try
                makeFigureS_ThresholdDTrajectories(simFile, theoryFile, ...
                    'D_values', cfg.D_values, 'outputDir', figuresDir);
            catch ME
                fprintf('Warning: S3 figures failed - %s\n', ME.message);
            end
        else
            fprintf('Skipping S3 trajectory figures - simulation data not found\n');
            fprintf('  To generate S3 figures, first run Run_modularFGM to create DynSim_2.mat\n');
            fprintf('  Searched locations:\n');
            for i = 1:length(possibleSimFiles)
                fprintf('    - %s\n', possibleSimFiles{i});
            end
        end
    else
        fprintf('Skipping S3 figures - theory data not found\n');
    end
    
    %% ================================================================
    %  S4. NESTED FGM DISTRIBUTION VISUALIZATION
    %  ================================================================
    fprintf('\n--- S4. Nested FGM Distributions ---\n');
    
    try
        makeFigureS_NestedFGMDistributions('outputDir', figuresDir);
    catch ME
        fprintf('Warning: S4 figures failed - %s\n', ME.message);
    end
    
    %% ================================================================
    %  SUMMARY
    %  ================================================================
    fprintf('\n');
    fprintf('========================================================\n');
    fprintf('  SUPPLEMENTARY ANALYSES COMPLETE\n');
    fprintf('========================================================\n');
    fprintf('\n');
    
    % List generated files
    fprintf('Generated data files:\n');
    dataFilesList = dir(fullfile(resultsDir, '*.mat'));
    for i = 1:length(dataFilesList)
        fprintf('  - %s\n', dataFilesList(i).name);
    end
    
    fprintf('\nGenerated figures:\n');
    figFiles = dir(fullfile(figuresDir, '*.eps'));
    for i = 1:length(figFiles)
        fprintf('  - %s\n', figFiles(i).name);
    end
    
    fprintf('\nAll supplementary analyses complete. Mode: %s\n', mode);
end

%% Helper function
function sectionTimer(label)
    persistent lastT
    if isempty(lastT)
        lastT = tic;
        return;
    end
    if contains(label, 'start')
        lastT = tic;
    else
        dt = toc(lastT);
        fprintf('[%s] completed in %.2f s\n', label, dt);
    end
end
