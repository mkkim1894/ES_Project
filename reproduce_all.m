%% reproduce_all.m
% Master script to reproduce all simulations and figures.
% Run from the project root directory.
% Expected runtime: ~10 hours on a machine with parallel computing
% toolbox (with 10 cores).
%
% Run test_all first. It runs the same stages in the same order at reduced
% parameter sets into results_test/, in minutes, and reports what breaks.
%
% Figure 1 is drawn by hand. Figure 7 and Figure S6 come from the Python
% analysis in LTEE_analysis/, which is not driven from here.
%
% Outputs:
%   results/                        main simulation data
%   results/Supplementary/          supplementary simulation data
%   results/Figures/                main figures
%   results/Supplementary/Figures/  supplementary figures
%   results/main_figures/           the manuscript set, collected

cd(fileparts(mfilename('fullpath')));
pathList = strsplit(genpath('.'), pathsep);
pathList = pathList(~contains(pathList, '_relegated') & ~contains(pathList, '_to_delete'));
addpath(strjoin(pathList, pathsep));

%% Simulations
fprintf('=== Running simulations ===\n');

% NOTE ON INITIALIZATION. 'sampled' is the initializer behind the published
% pleiotropic results, so this is what reproduces the current Figures 2-4. It
% becomes 'maxent' once the maximum-entropy migration lands; see
% docs/maximum-entropy-initialization.md. A bare Run_pleiotropicFGM defaults
% to 'greedy', which is the original initializer and reproduces nothing in the
% manuscript.
Run_pleiotropicFGM('reproduce', 'sampled')
Run_modularFGM('reproduce')
Run_nestedFGM('reproduce')
Run_nestedFGM('reproduce', {}, [10, 20])   % asymmetric dimensionality

Run_modularConstantSupply('mode', 'reproduce', ...
    'stages', {'sswm', 'cm_asexual', 'cm_sexual'})

Run_pleiotropicRestrictedTheta('mode', 'reproduce', 'init', 'maxent', ...
    'stages', {'diagnose', 'sswm', 'cm_asexual', 'cm_sexual'})

Run_supplementary('reproduce')
Run_ThresholdDAnalysis('reproduce')        % no simulation; needs Run_modularFGM
Run_nModule('stages', {'modular', 'pleiotropic'})

%% Figures
fprintf('=== Generating figures ===\n');

Run_allMainFigures('mode', 'reproduce')    % F2-F6 and the restricted cone
Run_supplementaryFigures('reproduce')      % S1-S5
Run_nModule('stages', {'figures'})         % S7, S8

fprintf('=== Done ===\n');
fprintf('F1 is drawn by hand; F7 and S6 come from LTEE_analysis/ltee_analysis.py\n');
