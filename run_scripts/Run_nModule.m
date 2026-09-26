function Run_nModule(varargin)
% Run_nModule  Run the n-module simulations and generate the figures.
%
% Executes the modular and pleiotropic simulations and the figure script in
% sequence, reporting the wall-clock time of each stage. Simulation output is
% written to results/nModule and figures to results/Figures.
%
% Name-value pairs
%   'n'                  Number of modules. Default 10.
%   'popSize'            Population size N. Default 1e4.
%   'numGenerations'     Generations per run. Default 1e4.
%   'numReplicates'      Replicates per initial condition. Default 30.
%   'conditions'         nCond-by-2 matrix of initial conditions, each row
%                        [k, c]. Default [], giving three conditions: k = 1,
%                        c = 0.15 (concentrated), k = 3, c = 0.25
%                        (intermediate) and a uniform control. See
%                        runModularND for the two-level profile.
%   'geneticTargetSize'  Loci per module L_i. Default [], in which case L_i
%                        accommodates the largest deficit required by any
%                        initial condition.
%   'referenceModel'     Model shown in the reference figure, 'Modular'
%                        (default) or 'Pleiotropic'.
%   'referenceCond'      Condition shown in the reference figure. Default 1.
%   'stages'             Cell array selecting stages to run, any of
%                        'modular', 'pleiotropic', 'figures'. Default all
%                        three; use this to regenerate figures without
%                        repeating the simulations.
%
% Dependencies
%   runModularND, runPleiotropicND, makeFigure_nModule, and the Statistics
%   and Machine Learning Toolbox.
%
% Reference
%   Kim, M., Ardell, S. M., & Kryazhimskiy, S. (2025). Module-Selection
%   Balance in the Evolution of Modular Organisms.

% ----------------------------- Parse inputs -----------------------------
p = inputParser;
addParameter(p, 'n',                10,        @isscalar);
addParameter(p, 'popSize',          1e4,       @isscalar);
addParameter(p, 'numGenerations',   1e4,       @isscalar);
addParameter(p, 'numReplicates',    30,        @isscalar);
addParameter(p, 'conditions',       [],        @(v) isempty(v) || size(v,2) == 2);
addParameter(p, 'geneticTargetSize', [],       @(v) isempty(v) || isscalar(v));
addParameter(p, 'referenceModel',   'Modular', @ischar);
addParameter(p, 'referenceCond',    1,         @isscalar);
addParameter(p, 'stages', {'modular', 'pleiotropic', 'figures'}, @iscellstr);
addParameter(p, 'outRoot',          '',        @ischar);   % '' = <project>/results
parse(p, varargin{:});

opts   = p.Results;
stages = lower(opts.stages);

validStages = {'modular', 'pleiotropic', 'figures'};
unknown     = setdiff(stages, validStages);
if ~isempty(unknown)
    error('Run_nModule:UnknownStage', ...
          'Unknown stage: %s. Valid stages are %s.', ...
          strjoin(unknown, ', '), strjoin(validStages, ', '));
end

% --------------------------- Path and toolbox ---------------------------
thisDir = fileparts(mfilename('fullpath'));
if isempty(thisDir), thisDir = pwd; end
addpath(thisDir);

required = {'runModularND', 'runPleiotropicND', 'makeFigure_nModule'};
missing  = required(cellfun(@(f) isempty(which(f)), required));
if ~isempty(missing)
    error('Run_nModule:MissingDependencies', ...
          'Not on the MATLAB path: %s', strjoin(missing, ', '));
end

if isempty(which('poissrnd')) || isempty(which('mnrnd')) || isempty(which('datasample'))
    error('Run_nModule:MissingToolbox', ...
          'Statistics and Machine Learning Toolbox is required.');
end

% ------------------------------- Header ---------------------------------
fprintf('\n========================================================\n');
fprintf('  n-module simulations\n');
fprintf('  started %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf('  n = %d, N = %g, %g generations per run\n', ...
        opts.n, opts.popSize, opts.numGenerations);
fprintf('  %d replicates per condition\n', opts.numReplicates);
fprintf('  stages: %s\n', strjoin(stages, ', '));
fprintf('========================================================\n');

timings   = struct('stage', {}, 'seconds', {});
totalTic  = tic;

% ------------------------------- Stages ---------------------------------
if ismember('modular', stages)
    t0 = tic;
    runModularND('n', opts.n, 'popSize', opts.popSize, ...
                 'numGenerations', opts.numGenerations, ...
                 'numReplicates', opts.numReplicates, ...
                 'conditions', opts.conditions, ...
                 'geneticTargetSize', opts.geneticTargetSize, ...
                 'outRoot', opts.outRoot);
    timings(end+1) = struct('stage', 'Modular simulation', 'seconds', toc(t0));
end

if ismember('pleiotropic', stages)
    t0 = tic;
    runPleiotropicND('n', opts.n, 'popSize', opts.popSize, ...
                     'numGenerations', opts.numGenerations, ...
                     'numReplicates', opts.numReplicates, ...
                     'conditions', opts.conditions, ...
                     'geneticTargetSize', opts.geneticTargetSize, ...
                     'outRoot', opts.outRoot);
    timings(end+1) = struct('stage', 'Pleiotropic simulation', 'seconds', toc(t0));
end

if ismember('figures', stages)
    t0 = tic;
    makeFigure_nModule('n', opts.n, 'popSize', opts.popSize, ...
                       'referenceModel', opts.referenceModel, ...
                       'referenceCond', opts.referenceCond, ...
                       'resultsRoot', opts.outRoot);
    timings(end+1) = struct('stage', 'Figures', 'seconds', toc(t0));
end

% ------------------------------- Summary --------------------------------
totalSeconds = toc(totalTic);

fprintf('\n========================================================\n');
fprintf('  Summary\n');
fprintf('--------------------------------------------------------\n');
for i = 1:numel(timings)
    fprintf('  %-24s %10s\n', timings(i).stage, formatDuration(timings(i).seconds));
end
fprintf('--------------------------------------------------------\n');
fprintf('  %-24s %10s\n', 'Total', formatDuration(totalSeconds));
fprintf('  finished %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf('========================================================\n\n');
end

% ======================================================================
% Helper subfunctions
% ======================================================================

function s = formatDuration(seconds)
% Format a duration in seconds, minutes or hours, whichever is most readable.
    if seconds < 60
        s = sprintf('%.1f s', seconds);
    elseif seconds < 3600
        s = sprintf('%dm %02ds', floor(seconds/60), round(mod(seconds, 60)));
    else
        s = sprintf('%dh %02dm', floor(seconds/3600), ...
                    round(mod(seconds, 3600) / 60));
    end
end
