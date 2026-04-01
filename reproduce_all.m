%% reproduce_all.m
% Master script to reproduce all simulations and figures.
% Run from the project root directory.
% Expected runtime: ~10 hours on a machine with parallel computing
% toolbox (with 10 cores).
%
% Outputs:
%   results/           - Main simulation data
%   results_supplementary/ - Supplementary simulation data
%   results/Figures/   - Main figures
%   results_supplementary/Figures/ - Supplementary figures

cd(fileparts(mfilename('fullpath')));
addpath(genpath('.'));

%% Simulations
fprintf('=== Running simulations ===\n');

Run_pleiotropicFGM
Run_modularFGM
Run_discordantFGM
Run_nestedFGM
Run_nestedFGM('reproduce', {}, [10, 20]) % asymmetric dimensionality
Run_supplementary
Run_ThresholdDAnalysis

%% Figures
fprintf('=== Generating figures ===\n');

Run_mainFigures('reproduce')
Run_supplementaryFigures('reproduce')

fprintf('=== Done ===\n');