function T = Run_smokeTest(varargin)
% Run_smokeTest - run every simulation and every figure script at test scale.
%
%   The engine behind test_all. Exercises each code path once, in dependency
%   order, at the smallest configuration that still visits every branch, and
%   writes everything into results_test/ so that nothing in results/ is
%   touched. Nothing it produces is a result.
%
%   A failure in one stage is reported and does not stop the others, so one
%   pass tells you everything that is broken rather than the first thing.
%
%   Each stage is timed, and its expected output is checked for existence AND
%   for a modification time later than the moment that stage began. A stage
%   that silently wrote nothing reports 'none'; a stage that left an older
%   file in place reports 'stale' instead of passing on the strength of it.
%
%   The ten-module stages run at n = 4, N = 1e3 rather than the published
%   n = 10, N = 1e4, and are pointed at results_test/ by 'outRoot'.
%
% Name-value pairs
%   'only'     Cell array of stage names to run. Default all. The names are
%              the first column of the returned table.
%   'skip'     Cell array of stage names to leave out.
%   'root'     Project root. Default inferred from this file's location.
%
% Output
%   T  Table, one row per stage: stage, status, seconds, output, message.
%
% Usage
%   cd ~/ES_Project
%   test_all                                        % the usual way in
%
%   addpath(genpath(pwd));
%   Run_smokeTest('skip', {'nmodule_sims', 'fig_nmodule'})
%   Run_smokeTest('only', {'restricted_maxent', 'fig_restrictedtheta'})
%
% See also: test_all, reproduce_all
%
% Copyright (c) 2025 Minkyu Kim, Cornell University
% Licensed under MIT License

p = inputParser;
addParameter(p, 'only', {}, @iscell);
addParameter(p, 'skip', {}, @iscell);
addParameter(p, 'root', '', @ischar);
parse(p, varargin{:});
P = p.Results;

root = P.root;
if isempty(root)
    root = fileparts(fileparts(mfilename('fullpath')));
    if isempty(root) || ~isfolder(fullfile(root, 'run_scripts')), root = pwd; end
end

% Retained-but-inactive directories hold superseded copies of these functions;
% leaving them on the path lets one shadow the version that should run.
pl = strsplit(genpath(root), pathsep);
pl = pl(~contains(pl, '_relegated') & ~contains(pl, '_to_delete'));
addpath(strjoin(pl, pathsep));

rt = fullfile(root, 'results_test');

% -------------------------------------------------------------------- stages
%  name                   what it runs                                        expected output (glob, relative to root)
S = {
 'pleiotropic_sampled',  @() Run_pleiotropicFGM('test', 'sampled'),           'results_test/SSWM/PleiotropicFGM_*.mat'
 'pleiotropic_maxent',   @() Run_pleiotropicFGM('test', 'maxent'),            'results_test/SSWM/PleiotropicFGM-maxentinit_*.mat'
 'maxent_wrapper',       @() maxEntWrapper(),                                 'results_test/SSWM/PleiotropicFGM-maxentinit_*.mat'
 'modular',              @() Run_modularFGM('test'),                          'results_test/SSWM/ModularFGM_*.mat'
 'nested_symmetric',     @() Run_nestedFGM('test'),                           'results_test/Generalization/NestedFGM/SSWM/NestedFGM_SSWM_*n10-10*.mat'
 'nested_asymmetric',    @() Run_nestedFGM('test', {}, [10 20]),              'results_test/Generalization/NestedFGM/SSWM/NestedFGM_SSWM_*n10-20*.mat'
 'constantsupply_sims',  @() Run_modularConstantSupply('mode','test', ...
                              'stages', {'sswm','cm_asexual','cm_sexual'}),   'results_test/ConstantSupply/SSWM/*.mat'
 'restricted_maxent',    @() Run_pleiotropicRestrictedTheta('mode','test', ...
                              'init','maxent', 'stages', ...
                              {'diagnose','sswm','cm_asexual','cm_sexual'}),  'results_test/RestrictedTheta_maxentinit/SSWM/*.mat'
 'restricted_sampled',   @() Run_pleiotropicRestrictedTheta('mode','test', ...
                              'init','sampled', 'stages', {'sswm'}),          'results_test/RestrictedTheta_sampledinit/SSWM/*.mat'
 'supplementary_sims',   @() Run_supplementary('test'),                       'results_test/Supplementary/SteadyStateCM_N_1000.mat'
 'thresholdD',           @() Run_ThresholdDAnalysis('test'),                  'results_test/ThresholdD_Analysis.mat'
 'nmodule_sims',         @() Run_nModule('n', 4, 'popSize', 1e3, ...
                              'numGenerations', 300, 'numReplicates', 3, ...
                              'outRoot', rt, ...
                              'stages', {'modular','pleiotropic'}),           'results_test/nModule/Pleiotropic_nModule_n4_N1e+03.mat'
 % ------------------------------------------------------------------ figures
 'fig_main',             @() Run_mainFigures('test'),                         'results_test/Figures/Figure2_SSWM.pdf'
 'fig_constantsupply',   @() Run_modularConstantSupply('mode','test', ...
                              'stages', {'figures'}),                         'results_test/Figures/Figure_ModularConstantSupply_Generations.pdf'
 'fig_restrictedtheta',  @() makeFigure_RestrictedTheta('resultsRoot', rt, ...
                              'init', 'maxent'),                              'results_test/Figures/Figure_RestrictedTheta_Generations.pdf'
 'fig_collector',        @() Run_allMainFigures('mode', 'test'),              'results_test/main_figures/F2.pdf'
 'fig_supplementary',    @() Run_supplementaryFigures('test'),                'results_test/Supplementary/Figures/FigureS_CM_Main.pdf'
 'fig_nmodule',          @() Run_nModule('n', 4, 'popSize', 1e3, ...
                              'outRoot', rt, 'stages', {'figures'}),          'results_test/Figures/Figure_nModule.pdf'
};

names = S(:,1);
keep  = true(size(names));
if ~isempty(P.only), keep = keep & ismember(names, lower(P.only)); end
if ~isempty(P.skip), keep = keep & ~ismember(names, lower(P.skip)); end
S = S(keep, :);
if isempty(S)
    error('Run_smokeTest:NoStages', ...
          'No stage matched. Valid names:\n  %s', strjoin(names', sprintf('\n  ')));
end

% ----------------------------------------------------------------------- run
fprintf('\n=============================================================\n');
fprintf('  test_all / Run_smokeTest    %d stages, test scale\n', size(S,1));
fprintf('  every output goes to %s\n', rt);
fprintf('=============================================================\n');

n      = size(S,1);
status = cell(n,1);  secs = zeros(n,1);  outcome = cell(n,1);  msg = cell(n,1);
tAll   = tic;

for i = 1:n
    fprintf('\n--- [%d/%d] %s ---\n', i, n, S{i,1});
    t0 = now;  tic;
    try
        S{i,2}();
        status{i} = 'ok';
        msg{i}    = '';
    catch err
        status{i} = 'FAILED';
        msg{i}    = firstLine(err.message);
        fprintf('  FAILED: %s\n', msg{i});
    end
    secs(i)    = toc;
    outcome{i} = checkOutput(root, S{i,3}, t0);
    fprintf('  %s | %.1f s | output: %s\n', status{i}, secs(i), outcome{i});
end

T = table(S(:,1), status, secs, outcome, msg, ...
    'VariableNames', {'stage','status','seconds','output','message'});

% ------------------------------------------------------------------- summary
fprintf('\n=============================================================\n');
fprintf('  finished in %.1f min\n', toc(tAll)/60);
bad = ~strcmp(status,'ok') | ~strcmp(outcome,'fresh');
fprintf('  %d stages, %d failed, %d wrote nothing, %d left a stale file\n', ...
        n, sum(~strcmp(status,'ok')), sum(strcmp(outcome,'none')), ...
        sum(strcmp(outcome,'stale')));
if any(bad)
    fprintf('\n  needs attention:\n');
    for i = find(bad)'
        fprintf('    %-22s %-7s %-6s %s\n', S{i,1}, status{i}, outcome{i}, msg{i});
    end
else
    fprintf('\n  every stage ran and wrote a fresh output. Pipeline is sound.\n');
end
fprintf('  Nothing in results/ was touched. Delete %s when done.\n', rt);
fprintf('=============================================================\n\n');
disp(T);
end

% ===========================================================================
function maxEntWrapper()
% Run_maxEntSimulations is a script driven by workspace variables, so it is
% called from a function of its own: the four below are its input and they do
% not leak into the caller.
MODE   = 'test';            %#ok<NASGU>
MODELS = {'unrestricted'};  %#ok<NASGU>
STAGES = {'sswm'};          %#ok<NASGU>
SEED   = 1;                 %#ok<NASGU>
Run_maxEntSimulations
end

% ---------------------------------------------------------------------------
function s = checkOutput(root, pattern, t0)
% 'fresh' the expected file exists and was written during this stage
% 'stale' it exists but predates the stage, so this run did not produce it
% 'none'  nothing matches
d = dir(fullfile(root, pattern));
d = d(~[d.isdir]);
if isempty(d)
    s = 'none';
elseif max([d.datenum]) >= t0 - 1/86400
    s = 'fresh';
else
    s = 'stale';
end
end

function s = firstLine(m)
k = find(m == sprintf('\n'), 1);
if isempty(k), s = m; else, s = m(1:k-1); end
end
