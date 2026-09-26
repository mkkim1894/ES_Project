%% test_all.m
% Run every simulation and every figure script at test scale, end to end.
%
% The mirror of reproduce_all: the same stages in the same order, at reduced
% parameter sets, with all output written to results_test/ so that nothing in
% results/ is touched. Use it to confirm the pipeline runs before committing
% to the full ~10 hour reproduction.
%
% One failure does not stop the rest, so a single pass reports everything that
% is broken. Each stage is timed, and its expected output is checked for a
% modification time inside that stage, so a stage that silently wrote nothing
% is reported rather than passing on a leftover file.
%
% Run from the project root directory.
%
% Outputs (all under results_test/):
%   results_test/{SSWM,CM_Asexual,CM_Sexual}/   pleiotropic and modular
%   results_test/ConstantSupply/                constant-supply control
%   results_test/RestrictedTheta_*init/         restricted cone
%   results_test/Generalization/NestedFGM/      nested FGM
%   results_test/nModule/                       ten-module, at n = 4
%   results_test/Supplementary/                 parameter-grid data
%   results_test/Figures/                       main figures
%   results_test/Supplementary/Figures/         supplementary figures
%   results_test/main_figures/                  the manuscript set, collected
%
% Expected runtime: minutes, not hours.
%
% Reference:
%   Kim, M., Ardell, S. M., & Kryazhimskiy, S. (2025).
%   "Module-Selection Balance in the Evolution of Modular Organisms."
%
% See also: reproduce_all, Run_smokeTest

cd(fileparts(mfilename('fullpath')));
pathList = strsplit(genpath('.'), pathsep);
pathList = pathList(~contains(pathList, '_relegated') & ~contains(pathList, '_to_delete'));
addpath(strjoin(pathList, pathsep));

testResults = Run_smokeTest();

% testResults is a table of stage, status, seconds, output and message.
% Anything other than status 'ok' and output 'fresh' needs attention before
% reproduce_all is worth starting.
