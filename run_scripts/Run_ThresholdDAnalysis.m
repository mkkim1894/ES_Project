%% Run_ThresholdDAnalysis.m
% Computes analytical trajectory predictions across multiple threshold D values.
% No simulations are run. Requires existing Run_modularFGM output.
%
% Usage:
%   Run_ThresholdDAnalysis              % reproduce mode (default)
%   Run_ThresholdDAnalysis('test')      % uses results_test/
%
% Outputs:
%   results/ThresholdD_Analysis.mat (or results_test/)
%
% Reference:
%   Kim, M., Ardell, S. M., & Kryazhimskiy, S. (2025).
%   "Module-Selection Balance in the Evolution of Modular Organisms."

function Run_ThresholdDAnalysis(mode)
    if nargin < 1 || isempty(mode)
        mode = 'reproduce';
    end

    isTest = strcmpi(mode, 'test');

    thisDir  = fileparts(mfilename('fullpath'));
    projRoot = fileparts(thisDir);
    addpath(fullfile(projRoot, 'simulation_scripts'));
    addpath(fullfile(projRoot, 'analysis_scripts'));
    addpath(fullfile(projRoot, 'utils'));
    addpath(fullfile(projRoot, 'figure_scripts'));

    resultsRoot = tern(isTest, ...
        fullfile(projRoot, 'results_test'), ...
        fullfile(projRoot, 'results'));

    fprintf('=== Threshold D Analysis ===\n');
    fprintf('Results root: %s\n', resultsRoot);

    % ----------------------------------------------------------------
    % Locate CM Asexual simulation file
    % ----------------------------------------------------------------
    cmDir   = fullfile(resultsRoot, 'CM_Asexual');
    pattern = fullfile(cmDir, 'ModularFGM_CM_Asexual_*.mat');
    files   = dir(pattern);

    if isempty(files)
        error('No ModularFGM CM_Asexual file found in %s.\nRun Run_modularFGM first.', cmDir);
    end
    if numel(files) > 1
        [~, newest] = max([files.datenum]);
        files = files(newest);
        warning('Multiple CM_Asexual files found. Using most recent: %s', files.name);
    end

    simFile = fullfile(files(1).folder, files(1).name);
    fprintf('Loading simulation file: %s\n', simFile);

    d        = load(simFile);
    simParams = d.simParams;

    % ----------------------------------------------------------------
    % D values to sweep (descending for figure panel order A-D)
    % ----------------------------------------------------------------
    D_values = [10000, 1000, 100, 10];

    % ----------------------------------------------------------------
    % Compute predictions for each D value
    % discretize initial phenotypes consistent with simulation start
    % ----------------------------------------------------------------
    spDiscrete = discretizeInitialPhenotypes(simParams);

    nAngles      = numel(simParams.initialAngles);
    nD           = numel(D_values);
    Theory_D_check = cell(nAngles, nD);

    for di = 1:nD
        D_val = D_values(di);
        fprintf('  Computing predictions for D = %d ...\n', D_val);
        traj = predictModularCM(spDiscrete, 'D', D_val);
        for j = 1:nAngles
            Theory_D_check{j, di} = traj{j};
        end
    end

    % ----------------------------------------------------------------
    % Save
    % ----------------------------------------------------------------
    outFile = fullfile(resultsRoot, 'ThresholdD_Analysis.mat');
    save(outFile, 'D_values', 'Theory_D_check', 'simParams');
    fprintf('Saved: %s\n', outFile);
    fprintf('\nDone. Run Run_supplementaryFigures to generate the threshold D figure.\n');
end

% ============================== Helpers ==================================

function sp = discretizeInitialPhenotypes(sp)
    delta = sp.deltaTrait;
    x0    = sp.initialPhenotypes;
    x0    = -delta * round(-x0 / delta);
    x0    = min(x0, 0);
    sp.initialPhenotypes = x0;
end

function y = tern(cond, a, b)
    if cond; y = a; else; y = b; end
end
