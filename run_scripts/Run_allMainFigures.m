function Run_allMainFigures(varargin)
% Run_allMainFigures - Regenerate every MATLAB main figure and collect them.
%
%   Redraws the main figures from existing simulation output and copies the
%   results into one directory under the manuscript's own numbering, so that the
%   set can be uploaded without picking it out of the results tree by hand. No
%   simulations are run.
%
%   Figure 1 is drawn by hand and Figure 7 is produced by the Python analysis in
%   LTEE_analysis, so neither is generated here. Both are reported as missing in
%   the closing summary, which lists what was collected and what was not.
%
%   Figures produced
%     F2   successive mutations, pleiotropic and modular
%     F3   concurrent mutations, linked chromosomes
%     F4   concurrent mutations, unlinked chromosomes
%     F5   nested Fisher's geometric model
%     F6   constant supply of module-improving mutations
%     restricted cone, across the three regimes (manuscript number not settled;
%     collected under its descriptive name)
%
%   A failure in one figure is reported and does not stop the others, so a
%   partial results tree still yields whatever it can support.
%
% Name-value pairs
%   'mode'     'reproduce' (default) or 'test'. Selects which results tree is
%              read and where the figures are written.
%   'collect'  Copy the figures into one directory. Default true.
%   'outDir'   Where to collect them. Default results/main_figures/, or
%              results_test/main_figures/ in test mode.
%   'only'     Cell array selecting a subset of {'main', 'constantsupply',
%              'restrictedtheta'}. Default all three.
%
% Usage
%   Run_allMainFigures
%   Run_allMainFigures('mode', 'test')
%   Run_allMainFigures('only', {'restrictedtheta'})
%   Run_allMainFigures('collect', false)
%
% Reference:
%   Kim, M., Ardell, S. M., & Kryazhimskiy, S. (2025).
%   "Module-Selection Balance in the Evolution of Modular Organisms."
%
% See also: Run_mainFigures, Run_modularConstantSupply, makeFigure_RestrictedTheta
%
% Copyright (c) 2025 Minkyu Kim, Cornell University
% Licensed under MIT License

p = inputParser;
addParameter(p, 'mode',    'reproduce', @ischar);
addParameter(p, 'collect', true,  @(v) islogical(v) || isnumeric(v));
addParameter(p, 'outDir',  '',    @ischar);
addParameter(p, 'only',    {'main', 'constantsupply', 'restrictedtheta'}, @iscell);
parse(p, varargin{:});
P = p.Results;

isTest = strcmpi(P.mode, 'test');
root   = fileparts(fileparts(mfilename('fullpath')));
if isempty(root), root = pwd; end

if isempty(P.outDir)
    if isTest
        P.outDir = fullfile(root, 'results_test', 'main_figures');
    else
        P.outDir = fullfile(root, 'results', 'main_figures');
    end
end

% Retained-but-inactive directories hold superseded copies of figure functions;
% leaving them on the path lets one shadow the version that should be used.
pathList = strsplit(genpath(root), pathsep);
pathList = pathList(~contains(pathList, '_relegated') & ~contains(pathList, '_to_delete'));
addpath(strjoin(pathList, pathsep));

if isTest, tree = 'results_test'; else, tree = 'results'; end

fprintf('\n=============================================================\n');
fprintf('  Run_allMainFigures   mode: %s\n', P.mode);
fprintf('=============================================================\n');

made = {};   % {manuscript name, source path}
skipped = {};

tAll = tic;

% ---------------------------------------------------------------- F2 - F5
if ismember('main', lower(P.only))
    fprintf('\n--- Figures 2-5 ---\n');
    try
        Run_mainFigures(P.mode);
        src = fullfile(root, tree, 'Figures');
        made = [made; {
            'F2.pdf', fullfile(src, 'Figure2_SSWM.pdf');
            'F3.pdf', fullfile(src, 'Figure3_CMLinked.pdf');
            'F4.pdf', fullfile(src, 'Figure4_CMUnlinked.pdf');
            'F5.pdf', fullfile(src, 'Figure_NestedFGM_Generations.pdf')}];
    catch err
        skipped = [skipped; {'F2-F5', err.message}];
        fprintf('  FAILED: %s\n', err.message);
    end
end

% ------------------------------------------------------------------- F6
if ismember('constantsupply', lower(P.only))
    fprintf('\n--- Figure 6, constant supply ---\n');
    try
        Run_modularConstantSupply('mode', P.mode, 'stages', {'figures'});
        made = [made; {'F6.pdf', fullfile(root, tree, 'Figures', ...
            'Figure_ModularConstantSupply_Generations.pdf')}];
    catch err
        skipped = [skipped; {'F6', err.message}];
        fprintf('  FAILED: %s\n', err.message);
    end
end

% -------------------------------------------------------- restricted cone
if ismember('restrictedtheta', lower(P.only))
    fprintf('\n--- Restricted cone ---\n');
    try
        makeFigure_RestrictedTheta('resultsRoot', fullfile(root, tree), ...
                                   'init', 'maxent');
        made = [made; {'F_RestrictedTheta.pdf', fullfile(root, tree, 'Figures', ...
            'Figure_RestrictedTheta_Generations.pdf')}];
    catch err
        skipped = [skipped; {'restricted cone', err.message}];
        fprintf('  FAILED: %s\n', err.message);
    end
end

% ---------------------------------------------------------------- collect
if P.collect && ~isempty(made)
    if ~isfolder(P.outDir), mkdir(P.outDir); end
    fprintf('\n--- Collecting into %s ---\n', P.outDir);
    for i = 1:size(made, 1)
        if isfile(made{i,2})
            copyfile(made{i,2}, fullfile(P.outDir, made{i,1}));
            fprintf('  %-24s <- %s\n', made{i,1}, shorten(made{i,2}, root));
        else
            skipped = [skipped; {made{i,1}, 'figure file not found after drawing'}]; %#ok<AGROW>
            fprintf('  %-24s NOT FOUND at %s\n', made{i,1}, shorten(made{i,2}, root));
        end
    end
end

% ---------------------------------------------------------------- summary
fprintf('\n=============================================================\n');
fprintf('  finished in %.1f s\n', toc(tAll));
fprintf('  Not produced here:\n');
fprintf('    F1   drawn by hand\n');
fprintf('    F7   LTEE_analysis/ltee_analysis.py\n');
if ~isempty(skipped)
    fprintf('  Incomplete:\n');
    for i = 1:size(skipped, 1)
        fprintf('    %-18s %s\n', skipped{i,1}, firstLine(skipped{i,2}));
    end
end
fprintf('  The restricted-cone figure is collected under its descriptive name;\n');
fprintf('  rename it once its manuscript number is settled. See docs/figures.md.\n');
fprintf('=============================================================\n\n');
end

% ---------------------------------------------------------------------------
function s = shorten(pth, root)
    if strncmp(pth, root, numel(root))
        s = pth(numel(root)+2:end);
    else
        s = pth;
    end
end

function s = firstLine(msg)
    nl = find(msg == sprintf('\n'), 1);
    if isempty(nl), s = msg; else, s = msg(1:nl-1); end
end
