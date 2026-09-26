%% Run_maxEntSimulations - the three evolutionary regimes, max-entropy initial genotypes
%
% Drives the production simulations with initial genotypes drawn from the
% maximum-entropy ensemble at the target phenotype instead of the
% random-walk-plus-swaps procedure:
%
%   SSWM         strong selection, weak mutation
%   CM asexual   concurrent mutations, complete linkage
%   CM sexual    concurrent mutations, free reassortment
%
% Nothing about the simulators changes. This calls the same
% Run_pleiotropicRestrictedTheta and Run_pleiotropicFGM that produce the
% published runs, with 'init' set to 'maxent', so any difference in the results
% is attributable to the initial genotypes alone.
%
% WHERE THE OUTPUT GOES
%   restricted      results[_test]/RestrictedTheta_maxentinit/<REGIME>/
%                   as RestrictedThetaFGM-maxentinit_*.mat
%   unrestricted    the main results tree, as PleiotropicFGM-maxentinit_*.mat
%
%   Both filenames carry the initializer, so neither can overwrite an existing
%   run. (The published greedy and sampled pleiotropic runs share one filename
%   between them, which is why the greedy ones were archived under
%   _relegated/from_results_tree/. This does not add to that problem.)
%
% USAGE
%   MODE = 'test'; Run_maxEntSimulations          % 8 replicates, minutes
%   MODE = 'reproduce'; Run_maxEntSimulations     % 250 replicates, hours
%
%   MODELS = {'restricted'};        % default is {'restricted'}
%   MODELS = {'restricted', 'unrestricted'};
%   STAGES = {'sswm'};              % default is all three regimes
%   SEED   = 1;
%
% RUN TEST FIRST. The test configuration uses 8 replicates and 2 initial
% conditions against 250 and 6, so it finishes in minutes and exercises every
% code path. Only when that comes back clean is the reproduce run worth
% starting.
%
% WHAT TO WATCH
%   Each regime prints the initializer's per-condition line before simulating,
%   so an unreachable initial condition or a poor landing shows up in seconds
%   rather than after the simulation. Check that
%     - every condition reports a residual at or below the tolerance (0.02)
%     - the number of 1s is what the geometry allows: under the cone,
%       n is near |x|/(delta * R) and cannot fall below |x|/delta
%     - termination reports few replicates stopping for want of beneficial
%       mutations
%
% Copyright (c) 2025 Minkyu Kim, Cornell University
% Licensed under MIT License

if ~exist('MODE',   'var') || isempty(MODE),   MODE   = 'test';           end
if ~exist('MODELS', 'var') || isempty(MODELS), MODELS = {'restricted'};   end
if ~exist('STAGES', 'var') || isempty(STAGES)
    STAGES = {'sswm', 'cm_asexual', 'cm_sexual'};
end
if ~exist('SEED',   'var') || isempty(SEED),   SEED   = 1;                end

HERE = fileparts(mfilename('fullpath'));
if isempty(HERE), HERE = pwd; end
ROOT = fileparts(HERE);
pp = strsplit(genpath(ROOT), pathsep);
pp = pp(~contains(pp, '_to_delete') & ~contains(pp, '_relegated'));
addpath(strjoin(pp, pathsep));
addpath(fullfile(ROOT, 'utils'), fullfile(ROOT, 'run_scripts'), ...
        fullfile(ROOT, 'figure_scripts'), fullfile(ROOT, 'analysis_scripts'));

if isempty(which('initializeGenomeMaxEnt'))
    error('Run_maxEntSimulations:MissingInit', ...
          'initializeGenomeMaxEnt is not on the path.');
end

fprintf('\n#############################################################\n');
fprintf('  Run_maxEntSimulations\n');
fprintf('  mode   : %s\n', MODE);
fprintf('  models : %s\n', strjoin(MODELS, ', '));
fprintf('  regimes: %s\n', strjoin(STAGES, ', '));
fprintf('  seed   : %d\n', SEED);
fprintf('#############################################################\n');

tAll = tic;

%% restricted pleiotropy, theta in [0, pi/2]
if ismember('restricted', MODELS)
    fprintf('\n===== restricted pleiotropy, max-entropy init =====\n');
    Run_pleiotropicRestrictedTheta( ...
        'mode',    MODE, ...
        'init',    'maxent', ...
        'stages',  STAGES, ...
        'seed',    SEED);
end

%% universal pleiotropy, theta in [0, 2*pi) - the control
if ismember('unrestricted', MODELS)
    fprintf('\n===== universal pleiotropy, max-entropy init =====\n');
    fprintf(['  NOTE: this is a control run, not a replacement for the\n' ...
             '  published pleiotropic results, which are already produced with\n' ...
             '  initializeGenomeSampled and start at n = [148 293 212 230 140\n' ...
             '  178], mean 200 - the same place the max-entropy ensemble puts\n' ...
             '  them. What max-entropy changes here is the scatter, not the\n' ...
             '  centre. (The greedy initializer gave n around 43; those runs\n' ...
             '  are archived under _relegated/from_results_tree.) Output is\n' ...
             '  tagged PleiotropicFGM-maxentinit and cannot be picked up by\n' ...
             '  makeRegimeFigure_Generations, whose glob is PleiotropicFGM_*.\n\n']);
    Run_pleiotropicFGM(MODE, 'maxent');
end

fprintf('\n#############################################################\n');
fprintf('  finished in %.1f min\n', toc(tAll)/60);
if strcmpi(MODE, 'test')
    fprintf('  Test configuration complete. Review the initializer and\n');
    fprintf('  termination summaries, then re-run with MODE = ''reproduce''.\n');
end
fprintf('#############################################################\n\n');
