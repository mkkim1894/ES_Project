function T = compareLogRatioConvergence(varargin)
% compareLogRatioConvergence - Does the module performance ratio converge?
%
%   Measures how far apart the six initial conditions remain in log R = log(x_2/x_1)
%   at the end of the usable part of a run, and compares that with how far apart
%   they started. A genotype-phenotype map that admits a module-selection balance
%   collapses the six onto one value; a map that does not preserves their order and
%   roughly their separation.
%
%   The quantity reported is the SPREAD, max_j log R_j - min_j log R_j over the
%   initial conditions j, of the replicate-mean log R. Its value at the end of the
%   window divided by its value at the start is the convergence ratio: near 0 means
%   the six have collapsed, near 1 means they have not.
%
%   Each dataset is evaluated over its own usable window, which ends at the last
%   generation where every initial condition still retains at least 'retention' of
%   its replicates. Replicates stop contributing once either trait comes within
%   delta of its optimum, which is the exclusion the figures apply; without it,
%   log R is dominated by lattice granularity near the optimum.
%
% Name-value pairs
%   'datasets'  n x 3 cell array {label, resultsRoot, modelDir}. Default compares
%               the restricted cone, the modular map and universal pleiotropy.
%   'regimes'   subset of {'SSWM', 'CM_Asexual', 'CM_Sexual'}. Default all three.
%   'retention' fraction of replicates that must survive. Default 0.80.
%   'nStamps'   generations sampled along each trajectory. Default 200.
%
% Output
%   T  table, one row per dataset and regime, with the spread at the start and at
%      the end of the window, their ratio, and the window length.
%
% Reference:
%   Kim, M., Ardell, S. M., & Kryazhimskiy, S. (2025).
%   "Module-Selection Balance in the Evolution of Modular Organisms."
%
% Copyright (c) 2025 Minkyu Kim, Cornell University
% Licensed under MIT License

p = inputParser;
addParameter(p, 'datasets',  defaultDatasets(), @iscell);
addParameter(p, 'regimes',   {'SSWM', 'CM_Asexual', 'CM_Sexual'}, @iscell);
addParameter(p, 'retention', 0.80, @isnumeric);
addParameter(p, 'nStamps',   200,  @isnumeric);
parse(p, varargin{:});
P = p.Results;

rows = {};
fprintf('\n%-24s %-12s %8s %8s %8s %10s %7s\n', 'model', 'regime', ...
        'spread0', 'spreadT', 'ratio', 'window', 'nCond');
fprintf('%s\n', repmat('-', 1, 84));

for i = 1:size(P.datasets, 1)
    label = P.datasets{i,1};
    root  = P.datasets{i,2};
    mdl   = P.datasets{i,3};

    for k = 1:numel(P.regimes)
        regime = P.regimes{k};
        f = locate(root, mdl, regime);
        if isempty(f)
            fprintf('%-24s %-12s %s\n', label, regime, '(no results)');
            continue;
        end

        [s0, sT, tEnd, nCond] = spreadOverWindow(f, P.retention, P.nStamps);
        if isnan(s0)
            fprintf('%-24s %-12s %s\n', label, regime, '(no usable window)');
            continue;
        end

        fprintf('%-24s %-12s %8.3f %8.3f %8.3f %10.0f %7d\n', ...
                label, regime, s0, sT, sT/s0, tEnd, nCond);
        rows(end+1, :) = {label, regime, s0, sT, sT/s0, tEnd, nCond}; %#ok<AGROW>
    end
end

fprintf('\nratio = spread at the end of the window / spread at the start.\n');
fprintf('Small means the initial conditions converged on a common ratio.\n\n');

if isempty(rows)
    T = [];
else
    T = cell2table(rows, 'VariableNames', ...
        {'model', 'regime', 'spread0', 'spreadT', 'ratio', 'window', 'nCond'});
end
end

% ---------------------------------------------------------------------------
function d = defaultDatasets()
    here = fileparts(mfilename('fullpath'));
    if isempty(here), here = pwd; end
    root = fileparts(here);
    res = fullfile(root, 'results');
    d = { ...
      'restricted cone',  res, 'RestrictedTheta_maxentinit'; ...
      'modular',          res, 'ModularFGM';                 ...
      'universal pleio.', res, 'PleiotropicFGM'              };
end

function f = locate(root, modelDir, regime)
% The results trees are organised either as <root>/<model>/<regime>/ or as
% <root>/<regime>/ with the model in the filename. Both are handled.
    f = '';
    cands = {fullfile(root, modelDir, regime), fullfile(root, regime)};
    for i = 1:numel(cands)
        if ~isfolder(cands{i}), continue; end
        d = dir(fullfile(cands{i}, '*.mat'));
        d = d(~[d.isdir]);
        if i == 2
            d = d(contains({d.name}, modelDir));
        end
        if isscalar(d)
            f = fullfile(d(1).folder, d(1).name);
            return;
        elseif numel(d) > 1
            error('compareLogRatioConvergence:Ambiguous', ...
                  ['%d result files in %s.\nMove superseded runs into a ' ...
                   '_relegated subfolder.'], numel(d), cands{i});
        end
    end
end

function [s0, sT, tEnd, nCond] = spreadOverWindow(file, retention, nStamps)
    s0 = NaN; sT = NaN; tEnd = NaN;

    d  = load(file);
    sp = d.simParams;
    R  = resultTableOf(d);
    [nCond, nRep] = size(R);
    need = max(2, ceil(retention * nRep));

    maxGen = 0;
    for j = 1:nCond
        for r = 1:nRep
            if ~isempty(R{j,r}), maxGen = max(maxGen, max(R{j,r}(:,2))); end
        end
    end
    if maxGen <= 1, return; end

    stamps = floor(linspace(1, maxGen, nStamps));
    meanLR = NaN(nCond, nStamps);
    ok     = false(nCond, nStamps);

    for j = 1:nCond
        LR = NaN(nRep, nStamps);
        for r = 1:nRep
            traj = R{j,r};
            if isempty(traj), continue; end
            tt = traj(:,2);
            for s = 1:nStamps
                if stamps(s) > max(tt), break; end
                idx = find(tt <= stamps(s), 1, 'last');
                if isempty(idx), continue; end
                x1 = traj(idx,3); x2 = traj(idx,4);
                % Same exclusion as the figures: drop a replicate once either
                % trait is within delta of its optimum.
                if abs(x1) < sp.deltaTrait || abs(x2) < sp.deltaTrait, break; end
                if x1 < -1e-9 && x2 < -1e-9
                    LR(r,s) = log(abs(x2) / abs(x1));
                end
            end
        end
        meanLR(j,:) = mean(LR, 1, 'omitnan');
        ok(j,:)     = sum(~isnan(LR), 1) >= need;
    end

    usable = find(all(ok, 1));
    if numel(usable) < 2, return; end

    a  = usable(1);
    b  = usable(end);
    s0 = max(meanLR(:,a)) - min(meanLR(:,a));
    sT = max(meanLR(:,b)) - min(meanLR(:,b));
    tEnd = stamps(b);
end

function R = resultTableOf(d)
    fn = fieldnames(d);
    for i = 1:numel(fn)
        v = d.(fn{i});
        if isstruct(v) && isfield(v, 'resultTable')
            R = v.resultTable;
            return;
        end
    end
    error('compareLogRatioConvergence:NoResultTable', ...
          'No resultTable in the loaded file.');
end
