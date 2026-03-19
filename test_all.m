%% test_all.m
% Quick pipeline verification in test mode.
% Run from the project root directory.
% Checks that all simulation and figure scripts run without errors
% using small parameter sets. Expected runtime: ~15 minutes.
%
% Does NOT reproduce publication-quality results.
% For full reproduction, run reproduce_all.m instead.

cd(fileparts(mfilename('fullpath')));
addpath(genpath('.'));

%% Simulations (test mode)
fprintf('=== Running simulations (test mode) ===\n');

Run_pleiotropicFGM('test')
Run_modularFGM('test')
Run_discordantFGM('test')
Run_nestedFGM('test')
Run_nestedFGM('test', {}, [10, 20]) % asymmetric dimensionality
Run_supplementary('test')
Run_ThresholdDAnalysis('test')

%% Figures
fprintf('=== Generating figures (test mode) ===\n');

Run_mainFigures('test')
Run_supplementaryFigures('test')

fprintf('=== Test complete. Check results_test/ and results_supplementary/Figures/ ===\n');