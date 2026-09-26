function makeFigure_nModule(varargin)
% makeFigure_nModule  Figures for the n-module simulations.
%
% Reads the output of runPleiotropicND and runModularND and writes
%
%   Figure_nModule.pdf            Panels A and B, the effective number of
%                                 lagging modules n_x against time for the
%                                 pleiotropic and modular GPFMs; panels C and
%                                 D, the effective number of module targets
%                                 n_g for the same two models.
%   Figure_nModule_Reference.pdf  Module performances and mean fitness for a
%                                 single run.
%
% Name-value pairs
%   'n'               Number of modules. Default 10.
%   'popSize'         Population size, used to match result files. Default 1e4.
%   'referenceModel'  'Modular' (default) or 'Pleiotropic'.
%   'referenceCond'   Condition shown in the reference figure. Default 1.
%   'numTimeStamp'    Generation stamps per n_x trajectory. Default 200.
%   'numWindows'      Windows spanning the longest run, for n_g. Default 6.
%   'numBoot'         Rarefaction subsamples for n_g. Default 2e3.
%   'targetDepth'     Rarefaction depth. Default 104, the depth used in the
%                     LTEE analysis. Windows that cannot support this depth
%                     are excluded.
%   'maxFitness'      Mean fitness at which n_x trajectories are truncated.
%                     Default 0.99, the criterion used to terminate the
%                     two-module simulations.
%   'rarefactionCheck'  Logical. Default false. When true, additionally
%                     writes Figure_nModule_RarefactionCheck.pdf and prints a
%                     table comparing n_g with and without rarefaction. The
%                     main figure is not affected.
%
% Statistics
%   n_g is the effective number of modules improved by selection over a
%   window of Delta t generations. With the improvement in module i written
%   as m_i = max(x_i(t + Dt) - x_i(t), 0) / delta and its share as
%   p_i = m_i / sum_j m_j,
%
%       n_g = 1 / sum_i p_i^2,
%
%   the counterpart of the effective number of gene targets computed for the
%   LTEE data. It is defined on both GPFMs. On the modular GPFM each
%   beneficial mutation changes exactly one x_i by delta, so m_i counts
%   fixations in module i. On the pleiotropic GPFM a mutation cannot be
%   attributed to a single module, but m_i remains defined as the share of
%   the population's trait-space improvement carried by trait i.
%
%   n_x is the effective number of lagging modules, the inverse Simpson index
%   of the deficit shares q_i = |x_i| / sum_j |x_j|. It describes where the
%   population lies in trait space rather than where selection is acting, and
%   generalizes the module performance ratio R of Figures 2-4.
%
%   Two properties of these estimators bear on their interpretation. The
%   simulation logs record the population mean phenotype, so m_i is a net
%   change and undercounts wherever beneficial and deleterious mutations
%   offset. Simpson's index is also compressed toward small values at small
%   sample size, by approximately n/K at depth K, so each window is rarefied
%   to a common depth before n_g is computed; at K = 104 the compression is
%   about 10 per cent for n = 10.
%
% Reference
%   Kim, M., Ardell, S. M., & Kryazhimskiy, S. (2025). Module-Selection
%   Balance in the Evolution of Modular Organisms.

% ----------------------------- Parse inputs -----------------------------
p = inputParser;
addParameter(p, 'n',              10,        @isscalar);
addParameter(p, 'popSize',        1e4,       @isscalar);
addParameter(p, 'referenceModel', 'Modular', @ischar);
addParameter(p, 'referenceCond',  1,         @isscalar);
addParameter(p, 'numTimeStamp',   200,       @isscalar);
addParameter(p, 'numWindows',     6,         @isscalar);
addParameter(p, 'numBoot',        2e3,       @isscalar);
addParameter(p, 'targetDepth',    104,       @isscalar);
addParameter(p, 'maxFitness',     0.99,      @isscalar);
addParameter(p, 'rarefactionCheck', false, @islogical);
addParameter(p, 'resultsRoot',    '',        @ischar);   % '' = <project>/results
parse(p, varargin{:});

n              = p.Results.n;
popSize        = p.Results.popSize;
referenceModel = p.Results.referenceModel;
referenceCond  = p.Results.referenceCond;
numTimeStamp   = p.Results.numTimeStamp;
numWindows     = p.Results.numWindows;
numBoot        = p.Results.numBoot;
targetDepth    = p.Results.targetDepth;
maxFitness     = p.Results.maxFitness;
rarefactionCheck = p.Results.rarefactionCheck;

% Sequential ramp; conditions are ordered by initial concentration.
condColors = [0.09 0.28 0.44; 0.15 0.46 0.66; 0.45 0.70 0.85];
markers    = {'d', '^', 'v', '>', '<', 'o'};

% Trace and reference-line widths for the main figure.
LW_TRACE = 2.6;
LW_REF   = 1.3;

% Minimum total remaining deficit, in units of delta, for which n_x is
% evaluated; below this the deficit shares are dominated by quantization.
% The threshold applies to the total sum_i |x_i| and not to individual |x_i|,
% since concentrated initial conditions begin with most modules at small or
% zero deficit.
MIN_DEFICIT_UNITS = 10;

% ---------------------------- Load result files -------------------------
resultsRoot = p.Results.resultsRoot;
if isempty(resultsRoot)
    resultsRoot = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'results');
end
dataDir     = fullfile(resultsRoot, 'nModule');
models      = {'Pleiotropic', 'Modular'};

data = struct();
for m = 1:2
    f = fullfile(dataDir, sprintf('%s_nModule_n%d_N%.0e.mat', models{m}, n, popSize));
    if ~isfile(f)
        error('makeFigure_nModule:MissingResults', ...
              'Not found: %s\nRun run%sND first.', f, models{m});
    end
    data.(models{m}) = load(f);
    fprintf('  Loading: %s\n', f);
end

cPle = data.Pleiotropic.simParams.conditions;
cMod = data.Modular.simParams.conditions;
if ~isequal(size(cPle), size(cMod)) || max(abs(cPle(:) - cMod(:))) > 1e-12
    error('makeFigure_nModule:ConditionMismatch', ...
        ['The two result files use different initial conditions. The ' ...
         'definitions of twoLevelDirection and defaultConditions must ' ...
         'agree between runModularND and runPleiotropicND.']);
end
conditions = cMod;
nCond      = size(conditions, 1);
deltaTrait = data.Modular.simParams.deltaTrait;

figDir = fullfile(resultsRoot, 'Figures');
if ~isfolder(figDir), mkdir(figDir); end

% ======================================================================
% n_g: gather pooled event counts, choose one common depth, bootstrap
%
% POOLING IS DELIBERATE. n_sel measures how broadly adaptive change is spread
% across modules, and the LTEE counterpart pools mutations across the six
% populations before computing Simpson's index. Computing it per replicate and
% then averaging measures something else - how broadly one population happened
% to move in one window - and it does NOT recover the same answer: on the
% modular map, condition 1, the pooled statistic rises to 9.8 by generation
% 5800 while the per-replicate average peaks at 6.4 and then falls back to 4.7.
% ======================================================================
maxGen = 0;
for m = 1:2
    rt = data.(models{m}).resultTable;
    for idx = 1:numel(rt)
        if ~isempty(rt{idx}), maxGen = max(maxGen, max(rt{idx}(:, 2))); end
    end
end

winWidth = maxGen / numWindows;
winShift = winWidth / 4;
winStart = 1:winShift:max(1, maxGen - winWidth);
winMid   = winStart + winWidth / 2;

counts = struct();
for m = 1:2
    counts.(models{m}) = gatherCounts(data.(models{m}).resultTable, ...
                                      deltaTrait, winStart, winWidth, n);
end

% A single rarefaction depth is applied to every window of every condition
% of both models, so that all panels are comparable. Windows that cannot
% support the depth are excluded rather than allowed to lower it. The LTEE
% analysis follows the same rule, its K = 104 being the minimum across
% windows that had already passed a count threshold.
depth = targetDepth;

printSamplingDiagnostic(counts, models, winMid, conditions, depth);

qualifying = zeros(1, 2);
for m = 1:2
    T = sum(counts.(models{m}), 3);
    qualifying(m) = nnz(T >= depth);
end

fprintf('\n  Rarefaction depth: K = %d (fixed; LTEE analysis uses 104)\n', depth);
fprintf('  Windows at or above K: Pleiotropic %d, Modular %d (of %d x %d)\n', ...
        qualifying(1), qualifying(2), nCond, numel(winStart));

if any(qualifying < 3)
    warning('makeFigure_nModule:TooFewWindows', ...
        ['Fewer than 3 windows reach the depth K = %d in at least one ' ...
         'model, which is not enough to show a trend. Scale numReplicates ' ...
         'by the factor printed in the diagnostic above.'], depth);
end

ng = struct();
for m = 1:2
    ng.(models{m}) = rarefyBootstrap(counts.(models{m}), depth, numBoot, ...
                                     depth, winMid);
end

% ======================================================================
% Main figure: n_x (A, B) over n_g (C, D)
%
% ROW ORDER. The lagging-module statistic n_x is the state of the population
% (how unevenly the remaining deficit is spread) and n_g is what selection is
% doing to it, so n_x is drawn on top and n_g below. This is the reverse of the
% earlier version, in which n_g was the top row.
% ======================================================================
figure('Units', 'centimeters', 'Position', [1, 1, 17.8, 14]);

pos = { [0.09, 0.58, 0.38, 0.32], [0.58, 0.58, 0.38, 0.32], ...
        [0.09, 0.10, 0.38, 0.32], [0.58, 0.10, 0.38, 0.32] };
lab = {'A', 'B', 'C', 'D'};

addColumnAnnotations({'Pleiotropic GPFM', 'Modular GPFM'}, [0.09, 0.58], 0.38, 0.93, 0.05);

% --- A, B: n_x, the effective number of LAGGING modules -------------
for m = 1:2
    subplot('Position', pos{m}); hold on;
    addSubplotLabel(lab{m}, pos{m});

    rt = data.(models{m}).resultTable;
    for j = 1:nCond
        col = condColors(min(j, size(condColors,1)), :);
        [gen, nxMean, nxSd] = computeNx(rt(j, :), numTimeStamp, deltaTrait, ...
                                        MIN_DEFICIT_UNITS, maxGen, maxFitness);

        % Spread across replicates, plus or minus one standard deviation, the
        % same convention used for the trait-space and log-ratio bands in the
        % main figures. Panels C and D carry a bootstrap interval over
        % rarefaction subsamples, which answers a different question (sampling
        % depth) and so is not directly comparable to this one.
        ok = ~isnan(nxMean) & ~isnan(nxSd);
        if any(ok)
            fill([gen(ok), fliplr(gen(ok))], ...
                 [nxMean(ok) + nxSd(ok), fliplr(nxMean(ok) - nxSd(ok))], ...
                 col, 'FaceAlpha', 0.18, 'EdgeColor', 'none');
        end

        plot(gen, nxMean, '-', 'Color', col, 'LineWidth', LW_TRACE);
        first = find(~isnan(nxMean), 1, 'first');
        if ~isempty(first)
            plot(gen(first), nxMean(first), 'Marker', markers{min(j, numel(markers))}, ...
                 'Color', col, 'MarkerFaceColor', col, 'MarkerSize', 7);
        end
    end

    yline(n, '--', 'LineWidth', LW_REF, 'Color', [0.4 0.4 0.4]);
    yline(1, ':',  'LineWidth', LW_REF, 'Color', [0.6 0.6 0.6]);
    xlim([1, maxGen]); ylim([0.5, n + 0.5]);
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 10);
    if m == 1
        ylabel({'Effective number of'; 'lagging modules, $n_\mathrm{lag}$'}, ...
               'Interpreter', 'latex', 'FontSize', 12);
    end
end

% --- C, D: n_g, the effective number of MODULE TARGETS --------------
for m = 1:2
    subplot('Position', pos{m+2}); hold on;
    addSubplotLabel(lab{m+2}, pos{m+2});

    st = ng.(models{m});
    for j = 1:nCond
        if isempty(st(j).centers), continue; end
        col = condColors(min(j, size(condColors,1)), :);
        fill([st(j).centers, fliplr(st(j).centers)], [st(j).lo, fliplr(st(j).hi)], ...
             col, 'FaceAlpha', 0.18, 'EdgeColor', 'none');
        plot(st(j).centers, st(j).mean, '-', 'Color', col, 'LineWidth', LW_TRACE);
        plot(st(j).centers(1), st(j).mean(1), 'Marker', markers{min(j, numel(markers))}, ...
             'Color', col, 'MarkerFaceColor', col, 'MarkerSize', 7);
    end

    yline(n, '--', 'LineWidth', LW_REF, 'Color', [0.4 0.4 0.4]);
    yline(1, ':',  'LineWidth', LW_REF, 'Color', [0.6 0.6 0.6]);
    xlim([1, maxGen]); ylim([0.5, n + 0.5]);
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 10);
    xlabel('Generations', 'FontName', 'Helvetica', 'FontSize', 12);
    if m == 1
        ylabel({'Effective number of'; 'module targets, $n_\mathrm{sel}$'}, ...
               'Interpreter', 'latex', 'FontSize', 12);
    end
end

outName = fullfile(figDir, 'Figure_nModule.pdf');
print(outName, '-dpdf', '-vector');
fprintf('Saved %s\n', outName);

% ======================================================================
% Reference figure
% ======================================================================
traj = data.(referenceModel).resultTable{referenceCond, 1};
if isempty(traj)
    error('makeFigure_nModule:EmptyReference', ...
          'Reference run %d of the %s model is empty.', referenceCond, referenceModel);
end

gen = traj(:, 2);  X = traj(:, 3:end);  W = traj(:, 1);
refCol = condColors(min(referenceCond, size(condColors,1)), :);

figure('Units', 'centimeters', 'Position', [1, 1, 17.8, 8]);
refPos = {[0.09, 0.16, 0.38, 0.68], [0.58, 0.16, 0.38, 0.68]};
addColumnAnnotations({ {'Module performance'; [referenceModel ' GPFM']}, ...
                       {'Mean fitness';       [referenceModel ' GPFM']} }, ...
                     [0.09, 0.58], 0.38, 0.85, 0.12);

subplot('Position', refPos{1});  addSubplotLabel('A', refPos{1});
plot(gen, X, '-', 'Color', refCol, 'LineWidth', 0.7);
xlim([1, max(gen)]);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 10);
xlabel('Generations', 'FontName', 'Helvetica', 'FontSize', 12);
ylabel('Module performance, $x_i$', 'Interpreter', 'latex', 'FontSize', 12);

subplot('Position', refPos{2});  addSubplotLabel('B', refPos{2});
plot(gen, W, '-', 'Color', refCol, 'LineWidth', 0.9);
xlim([1, max(gen)]);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 10);
xlabel('Generations', 'FontName', 'Helvetica', 'FontSize', 12);
ylabel('Mean fitness, $\bar{W}$', 'Interpreter', 'latex', 'FontSize', 12);

outName = fullfile(figDir, 'Figure_nModule_Reference.pdf');
print(outName, '-dpdf', '-vector');
fprintf('Saved %s\n', outName);

% ======================================================================
% Optional rarefaction check (diagnostic only; the main figure above is
% written first and is not affected by anything below).
%
% Overlays n_g computed WITHOUT rarefaction - Simpson's index taken from
% every event in the window - on the rarefied estimate the main figure uses.
% Same pooling, same windows, same qualifying rule. The only difference is
% whether the window is subsampled to K before 1/D is formed.
%
% Simpson's index is compressed at small sample size, so if window depths
% differ a lot the unrarefied curve will sit low where windows are shallow
% and rise as they deepen, producing a trend that is an artifact of sampling.
% If instead the two curves lie on top of each other, rarefaction is not
% doing any work on these data and the procedure can be simplified away.
% ======================================================================
if rarefactionCheck
    ngRaw = struct();
    for m = 1:2
        ngRaw.(models{m}) = rawSimpson(counts.(models{m}), depth, winMid);
    end

    figure('Units', 'centimeters', 'Position', [1, 1, 17.8, 8]);
    chkPos = {[0.09, 0.16, 0.38, 0.68], [0.58, 0.16, 0.38, 0.68]};
    addColumnAnnotations({'Pleiotropic GPFM', 'Modular GPFM'}, ...
                         [0.09, 0.58], 0.38, 0.85, 0.12);

    fprintf('\n  ---- rarefaction check: n_g with and without rarefaction ----\n');
    fprintf('    solid = rarefied to K = %d (main figure), dashed = all events\n', depth);
    fprintf('    %-12s %-14s %8s %10s %8s %9s %8s\n', ...
            'model', 'condition', 'gen', 'rarefied', 'raw', 'raw-rare', 'depth');

    for m = 1:2
        subplot('Position', chkPos{m}); hold on;
        addSubplotLabel(lab{m}, chkPos{m});

        st  = ng.(models{m});
        raw = ngRaw.(models{m});
        for j = 1:nCond
            if isempty(st(j).centers), continue; end
            col = condColors(min(j, size(condColors,1)), :);
            fill([st(j).centers, fliplr(st(j).centers)], ...
                 [st(j).lo, fliplr(st(j).hi)], ...
                 col, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
            plot(st(j).centers,  st(j).mean,  '-',  'Color', col, ...
                 'LineWidth', LW_TRACE);
            plot(raw(j).centers, raw(j).mean, '--', 'Color', col, ...
                 'LineWidth', 1.6);

            for q = 1:numel(raw(j).centers)
                fprintf('    %-12s k=%2d c=%.2f %8d %10.2f %8.2f %9+.2f %8d\n', ...
                        models{m}, conditions(j,1), conditions(j,2), ...
                        round(raw(j).centers(q)), st(j).mean(q), ...
                        raw(j).mean(q), raw(j).mean(q) - st(j).mean(q), ...
                        raw(j).depth(q));
            end
        end

        yline(n, '--', 'LineWidth', LW_REF, 'Color', [0.4 0.4 0.4]);
        yline(1, ':',  'LineWidth', LW_REF, 'Color', [0.6 0.6 0.6]);
        xlim([1, maxGen]); ylim([0.5, n + 0.5]);
        set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 10);
        xlabel('Generations', 'FontName', 'Helvetica', 'FontSize', 12);
        if m == 1
            ylabel({'Effective number of'; 'module targets, $n_\mathrm{sel}$'}, ...
                   'Interpreter', 'latex', 'FontSize', 12);
        end
    end

    % Size of the discrepancy, against the width of the rarefaction interval
    % that the main figure already draws. A discrepancy well inside that band
    % is not worth a procedure.
    dAll = []; ciAll = []; depAll = [];
    for m = 1:2
        st = ng.(models{m}); raw = ngRaw.(models{m});
        for j = 1:nCond
            if isempty(st(j).centers), continue; end
            dAll   = [dAll,   raw(j).mean - st(j).mean];      %#ok<AGROW>
            ciAll  = [ciAll,  (st(j).hi - st(j).lo) / 2];     %#ok<AGROW>
            depAll = [depAll, raw(j).depth];                  %#ok<AGROW>
        end
    end
    if isempty(dAll)
        warning('makeFigure_nModule:NoCheckWindows', ...
                ['No window reaches K = %d, so there is nothing to compare. ' ...
                 'Lower targetDepth or raise numReplicates.'], depth);
        return;
    end

    fprintf('\n    windows compared          : %d\n', numel(dAll));
    fprintf('    window depth range        : %d to %d events (ratio %.1fx)\n', ...
            min(depAll), max(depAll), max(depAll) / max(min(depAll), 1));
    fprintf('    mean(raw - rarefied)      : %+.3f\n', mean(dAll));
    fprintf('    max |raw - rarefied|      : %.3f\n', max(abs(dAll)));
    fprintf('    mean half-width of the rarefaction CI : %.3f\n', mean(ciAll));
    if max(abs(dAll)) < mean(ciAll)
        fprintf(['    -> every discrepancy is smaller than the interval the ' ...
                 'figure already shows.\n']);
    else
        fprintf(['    -> at least one discrepancy exceeds the interval the ' ...
                 'figure already shows.\n']);
    end
    fprintf('    correlation of (raw - rarefied) with window depth : %+.2f\n', ...
            corrDepth(dAll, depAll));

    outName = fullfile(figDir, 'Figure_nModule_RarefactionCheck.pdf');
    print(outName, '-dpdf', '-vector');
    fprintf('Saved %s\n', outName);
end

end

% ======================================================================
% Helper subfunctions
% ======================================================================

function S = rawSimpson(C, minEvents, winMid)
% n_g from every event in the window, with no subsampling. The qualifying
% rule matches rarefyBootstrap, so the two sets of windows coincide and the
% curves can be compared point by point. Diagnostic use only.
    [nCond, ~, ~] = size(C);
    S = repmat(struct('centers', [], 'mean', [], 'depth', []), nCond, 1);

    for i = 1:nCond
        tot  = sum(squeeze(C(i, :, :)), 2)';
        keep = find(tot >= minEvents);
        if isempty(keep), continue; end

        mm = zeros(1, numel(keep));
        for q = 1:numel(keep)
            share = squeeze(C(i, keep(q), :))' / tot(keep(q));
            mm(q) = 1 / sum(share.^2);
        end
        S(i).centers = winMid(keep);
        S(i).mean    = mm;
        S(i).depth   = tot(keep);
    end
end

function r = corrDepth(d, dep)
% Pearson correlation, guarded against a degenerate vector.
    if numel(d) < 3 || std(d) == 0 || std(dep) == 0
        r = NaN;
    else
        x = d(:) - mean(d(:));
        y = double(dep(:)) - mean(double(dep(:)));
        r = (x' * y) / sqrt((x' * x) * (y' * y));
    end
end

function C = gatherCounts(resultTable, deltaTrait, winStart, winWidth, n)
% Pooled improvement counts C(condition, window, module). Replicates are
% pooled within a window, as the six LTEE populations are.
    [nCond, nRep] = size(resultTable);
    C = zeros(nCond, numel(winStart), n);

    for i = 1:nCond
        for rep = 1:nRep
            traj = resultTable{i, rep};
            if isempty(traj), continue; end
            simTime = traj(:, 2);
            X       = traj(:, 3:end);

            for w = 1:numel(winStart)
                t0 = winStart(w);  t1 = t0 + winWidth;
                if t1 > max(simTime), break; end
                i0 = find(simTime <= t0, 1, 'last');
                i1 = find(simTime <= t1, 1, 'last');
                if isempty(i0) || isempty(i1) || i1 <= i0, continue; end

                % Improvement only, the counterpart of driver mutations.
                % Accumulated across replicates in real arithmetic and
                % rounded once, so that rounding does not accumulate.
                C(i, w, :) = squeeze(C(i, w, :))' + ...
                             max(X(i1, :) - X(i0, :), 0) / deltaTrait;
            end
        end
    end
    C = round(C);
end

function S = rarefyBootstrap(C, depth, numBoot, minEvents, winMid)
% Rarefy every qualifying window to the common depth and bootstrap n_g.
    [nCond, nWin, n] = size(C);
    S = repmat(struct('centers', [], 'mean', [], 'lo', [], 'hi', []), nCond, 1);

    for i = 1:nCond
        tot  = sum(squeeze(C(i, :, :)), 2)';
        keep = find(tot >= minEvents & tot >= depth);
        if isempty(keep), continue; end

        mm = zeros(1, numel(keep));  lo = mm;  hi = mm;
        for q = 1:numel(keep)
            w      = keep(q);
            labels = repelem(1:n, squeeze(C(i, w, :))');
            total  = numel(labels);

            vals = zeros(1, numBoot);
            for b = 1:numBoot
                sub   = accumarray(labels(randperm(total, depth))', 1, [n 1]);
                share = sub / depth;
                vals(b) = 1 / sum(share.^2);
            end
            % Aggregate as ltee_analysis.py does: average Simpson's index
            % across subsamples and invert once, rather than averaging the
            % inverse. The interval is unaffected, because 1/x is monotone and
            % the percentiles of 1/D are the reciprocals of the reversed
            % percentiles of D, which is exactly the band the LTEE figure draws.
            mm(q) = 1 / mean(1 ./ vals);
            lo(q) = prctile(vals, 2.5);
            hi(q) = prctile(vals, 97.5);
        end
        S(i).centers = winMid(keep);
        S(i).mean = mm;  S(i).lo = lo;  S(i).hi = hi;
    end
end

function [genStamps, nxMean, nxSd, nContrib] = computeNx(repRow, numTimeStamp, deltaTrait, ...
                                         minDeficitUnits, maxGen, maxFitness)
% Mean n_x across replicates of one condition, on a common time grid.
%
% Each replicate is truncated once its mean fitness reaches maxFitness. Past
% that point the population is at mutation-selection balance rather than
% adapting, and the deficit shares are set by that balance rather than by the
% trajectory.
    genStamps = floor(linspace(1, maxGen, numTimeStamp));
    acc = NaN(numel(repRow), numTimeStamp);

    for rep = 1:numel(repRow)
        traj = repRow{rep};
        if isempty(traj), continue; end
        simTime = traj(:, 2);
        X       = traj(:, 3:end);

        stop = find(traj(:, 1) >= maxFitness, 1, 'first');
        if ~isempty(stop)
            simTime = simTime(1:stop);  X = X(1:stop, :);
        end

        for t = 1:numTimeStamp
            if genStamps(t) > max(simTime), break; end
            idx = find(simTime <= genStamps(t), 1, 'last');
            if isempty(idx), idx = 1; end

            v     = abs(X(idx, :));
            total = sum(v);
            if total >= minDeficitUnits * deltaTrait
                share = v / total;
                acc(rep, t) = 1 / sum(share.^2);
            end
        end
    end
    nxMean   = mean(acc, 1, 'omitnan');
    nxSd     = std(acc, 0, 1, 'omitnan');
    nContrib = sum(~isnan(acc), 1);

    % Replicates leave the average at different times, because each is truncated
    % where its own mean fitness reaches maxFitness. Once only a handful remain,
    % the mean is dominated by whichever replicates happen to still be adapting.
    % The log-ratio panels of the main figures handle the same problem with an
    % 80 percent retention rule; the equivalent here is MIN_RETAINED_FRACTION.
    %
    % Set to 0.80, the same rule used in Figures 2-4, Figure 6 and the
    % restricted-theta figure. Without it the pleiotropic curves turn sharply
    % upward at the end of the run: those populations reach W = 0.99 around
    % generation 5000 and leave the average at different times, and the ones
    % still running are the ones whose traits are still far from their optima,
    % so the surviving mean n_lag rises. That upturn is an artifact of which
    % replicates remain, not a feature of the dynamics.
    MIN_RETAINED_FRACTION = 0.80;

    if MIN_RETAINED_FRACTION > 0
        minContrib = max(2, ceil(MIN_RETAINED_FRACTION * numel(repRow)));
        nxMean(nContrib < minContrib) = NaN;
        nxSd(nContrib < minContrib)   = NaN;
    end
end

function printSamplingDiagnostic(counts, models, winMid, conditions, depth)
% Per-window pooled event counts, with the factor by which the replicate
% count would have to be scaled for every window through the middle of the
% run to reach the rarefaction depth. Adaptive activity decays over a run, so
% the later windows are the sparse ones.
    fprintf('\n  ---- sampling diagnostic (pooled events per window) ----\n');
    for m = 1:2
        C = counts.(models{m});
        fprintf('  %s\n', models{m});
        for i = 1:size(C, 1)
            tot = sum(squeeze(C(i, :, :)), 2)';
            fprintf('    k=%2d c=%.2f :%s\n', conditions(i,1), conditions(i,2), ...
                    sprintf('%6d', tot));
        end
    end
    fprintf('    window mid-points   :%s\n', sprintf('%6d', round(winMid)));

    fprintf('    replicate scale-up needed to reach K = %d:\n', depth);
    for m = 1:2
        C = counts.(models{m});
        need = zeros(1, size(C, 1));
        for i = 1:size(C, 1)
            tot = sum(squeeze(C(i, :, :)), 2)';
            half = tot(1:ceil(numel(tot)/2));
            need(i) = depth / max(min(half(half > 0)), 1);
        end
        fprintf('      %-12s x%.1f  (to carry every window through mid-run)\n', ...
                models{m}, max(need));
    end
end

function addColumnAnnotations(titles, xPos, width, yPos, height)
    if nargin < 4, yPos   = 0.93; end
    if nargin < 5, height = 0.05; end
    for i = 1:numel(titles)
        annotation('textbox', [xPos(i), yPos, width, height], ...
            'String', titles{i}, 'FontSize', 14, 'FontName', 'Helvetica', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
            'EdgeColor', 'none');
    end
end

function addSubplotLabel(label, pos)
    annotation('textbox', [pos(1)-0.055, pos(2)+pos(4), 0.03, 0.03], ...
        'String', label, 'FontSize', 14, 'FontWeight', 'bold', 'EdgeColor', 'none');
end
