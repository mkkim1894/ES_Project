function makeFigure_RestrictedTheta(varargin)
% makeFigure_RestrictedTheta - Figure for the pleiotropic GPFM with a restricted
%   mutational cone.
%
% Layout:
%   A, B, C - the distribution of phenotypic effects available at three of the
%             initial phenotypes, identified by the marker that marks the start
%             of the corresponding trajectory in the panels below
%   D, E, F - trajectories in the trait space, one per evolutionary regime
%   G, H, I - log(x_2/x_1) against generations, same three regimes
%
% Set 'showDPE' false for the six-panel figure as it stood before the top row
% was added; D-I then revert to A-F and the layout is unchanged.
%
% HOW TO READ IT
%   This figure fills the fourth cell of the stress test:
%   universal pleiotropy combined with LOW genotypic redundancy of the optimum.
%
%   Black dashed: the path of steepest ascent through each initial phenotype,
%   x_2 = x_20 (x_1/x_10)^(a_1^2/a_2^2), equation "FGM trait path". This curve
%   follows from the phenotype-fitness map alone, so it is the SAME reference
%   for the restricted and the unrestricted GPM; what differs between them is
%   whether a population follows it. It is on by default ('showGradientPath').
%
%   The same curve is ALSO the SSWM prediction for the unrestricted pleiotropic
%   GPFM, and drawing it in yellow ('showReference', off by default) would say
%   that instead. Everywhere else in this project a yellow curve means "the
%   prediction for the model in this panel", so yellow here would invert the
%   convention and read as a failed theory. Black dashed says "reference
%   geometry", which is what is meant.
%
%   Orange: the equal fitness benefits line x_1/a_1^2 = x_2/a_2^2, and the
%   horizontal log(a_2^2/a_1^2) in the ratio panels - the module-selection
%   balance of the modular GPFM. Trajectories that collapse onto it show a
%   balance arising from low redundancy alone, with no modular encoding, which
%   would narrow the paper's claim. Trajectories that do not show that reducing
%   redundancy under pleiotropy is not sufficient, and modularity is doing the
%   work.
%
%   Judge convergence in the log-ratio panels, not the trajectory panels, where
%   trajectories crowd together near the optimum for reasons that have nothing
%   to do with an attractor.
%
%   The top row is a property of the phenotype, not of a simulation: at each of
%   those initial phenotypes it is the distribution over theta of the beneficial
%   mutations available to the genotypes that realise it. Under universal
%   pleiotropy that distribution is flat over the window shown; under the
%   restricted cone it is skewed, and its lobe sits on the direction to the
%   optimum rather than on the fitness gradient. See drawDPEPanel for the key.
%   It is computed by computeDPEAtPhenotype, in analysis_scripts.
%
% Usage
%   makeFigure_RestrictedTheta
%   makeFigure_RestrictedTheta('resultsRoot', '/path/to/results')
%   makeFigure_RestrictedTheta('thetaRange', [0, pi/2])
%
% Name-value pairs
%   'resultsRoot'  Directory containing RestrictedTheta/{SSWM, CM_Asexual,
%                  CM_Sexual}. Default: results/ beside this file, falling back
%                  to results_test/.
%   'thetaRange'   Select the run with this cone when several are present.
%   'yLimits'      y-limits of panels D-F. Default [-3, 2], as in Figures 2-4.
%   'showGradientPath' Overlay the path of steepest ascent in black dashed.
%                  Default TRUE - see HOW TO READ IT above.
%   'showReference' Overlay the same curve in yellow, as the unrestricted-model
%                  prediction. Default FALSE - see HOW TO READ IT above.
%   'outputFile'   Filename without extension. Default
%                  Figure_RestrictedTheta_Generations.
%
% Regimes with no results file are drawn as an empty labelled panel rather than
% raising an error, so the figure can be produced after the SSWM stage alone.
%
% Reference:
%   Kim, M., Ardell, S. M., & Kryazhimskiy, S. (2025).
%   "Module-Selection Balance in the Evolution of Modular Organisms."
%
% See also: Run_pleiotropicRestrictedTheta, makeFigure_ConstantSupply
%
% Copyright (c) 2025 Minkyu Kim, Cornell University
% Licensed under MIT License

% ----------------------------- Parse inputs -----------------------------
p = inputParser;
addParameter(p, 'resultsRoot', '',  @ischar);
addParameter(p, 'thetaRange',  [],  @(v) isempty(v) || (isnumeric(v) && numel(v) == 2));
addParameter(p, 'yLimits',     [-3, 2], @(v) isnumeric(v) && numel(v) == 2);
addParameter(p, 'outputFile',  'Figure_RestrictedTheta_Generations', @ischar);
% Off by default: see HOW TO READ IT above. Set true to overlay the unrestricted
% isotropic FGM gradient path saved as referenceTrajectories.
addParameter(p, 'showReference', false, @(v) islogical(v) || isnumeric(v));
addParameter(p, 'showGradientPath', true, @(v) islogical(v) || isnumeric(v));
% The top row. Off gives the six-panel figure exactly as it stood before.
addParameter(p, 'showDPE',     true, @(v) islogical(v) || isnumeric(v));
% WHICH initial conditions the top row shows, as indices into the initial
% conditions the run itself used - so the marker over each panel is the marker
% on that trajectory by construction, not by matching a value. Default: the
% shallowest, the middle one and the steepest.
addParameter(p, 'dpeConditions', [1, 3, 6], @isnumeric);
addParameter(p, 'dpeNGeno',    2000, @isnumeric);
% The ensemble draws take a couple of minutes, so the result is cached beside
% the figure and reused when the parameters match. Set false to force a redraw.
addParameter(p, 'dpeCache',    true, @(v) islogical(v) || isnumeric(v));
% Which initializer's runs to draw. Each has its own results subtree.
addParameter(p, 'init',        'sampled', @ischar);
parse(p, varargin{:});

yLimits       = p.Results.yLimits;
outputFile    = p.Results.outputFile;
thetaRange    = p.Results.thetaRange;
showReference    = logical(p.Results.showReference);
showGradientPath = logical(p.Results.showGradientPath);
showDPE          = logical(p.Results.showDPE);
dpeConditions    = p.Results.dpeConditions;
dpeNGeno         = p.Results.dpeNGeno;
dpeCache         = logical(p.Results.dpeCache);
rtSubtree        = ['RestrictedTheta_' lower(p.Results.init) 'init'];

thisDir = fileparts(mfilename('fullpath'));
if isempty(thisDir), thisDir = pwd; end

resultsRoot = p.Results.resultsRoot;
if isempty(resultsRoot)
    projRoot  = fileparts(thisDir);
    candidate = fullfile(projRoot, 'results');
    if isfolder(fullfile(candidate, rtSubtree))
        resultsRoot = candidate;
    else
        resultsRoot = fullfile(projRoot, 'results_test');
    end
end
fprintf('makeFigure_RestrictedTheta loading from: %s\n', resultsRoot);

markers        = {'d','^','v','>','<','o'};
if showDPE
    subplot_labels = {'A','B','C','D','E','F','G','H','I'};
else
    subplot_labels = {'A','B','C','D','E','F'};
end
dirs           = {'SSWM', 'CM_Asexual', 'CM_Sexual'};
column_titles  = {'Successive mutations', ...
                  'Concurrent mutations, linked chromosomes', ...
                  'Concurrent mutations, unlinked chromosomes'};

% --------------------------- Locate files -------------------------------
file_names = cell(1,3);
for i = 1:3
    dirPath = fullfile(resultsRoot, rtSubtree, dirs{i});
    files = [];
    if isfolder(dirPath)
        files = dir(fullfile(dirPath, sprintf('RestrictedThetaFGM*_%s_*.mat', dirs{i})));
        files = files(~[files.isdir]);
    end

    if isempty(files)
        file_names{i} = '';
        fprintf('  %-11s no results yet - panel will be left empty\n', dirs{i});
        continue;
    end

    if ~isempty(thetaRange)
        tag = sprintf('_th%+.3f%+.3f', thetaRange(1), thetaRange(2));
        hit = contains({files.name}, tag);
        if any(hit), files = files(hit); end
    end

    if ~isscalar(files)
        error('makeFigure_RestrictedTheta:AmbiguousResults', ...
              ['%d result files in %s match:\n  %s\n' ...
               'Pass ''thetaRange'' to select one, or move the superseded runs ' ...
               'into a _relegated/ subfolder.'], ...
              numel(files), dirPath, strjoin({files.name}, sprintf('\n  ')));
    end

    file_names{i} = fullfile(files(1).folder, files(1).name);
    fprintf('  %-11s %s\n', dirs{i}, files(1).name);
end

if all(cellfun(@isempty, file_names))
    error('makeFigure_RestrictedTheta:NoResults', ...
          ['No results found under %s.\nRun Run_pleiotropicRestrictedTheta first.'], ...
          fullfile(resultsRoot, rtSubtree));
end

% ------------------------- Top row, if asked for ------------------------
% Computed before the figure is opened: the ensemble draws take a couple of
% minutes, and leaving a half-drawn figure on screen for that long is
% indistinguishable from a hang.
if showDPE
    [D, dpeMarkers] = getDPE(file_names, resultsRoot, dpeConditions, dpeNGeno, ...
                             dpeCache, markers);
end

% --------------------------- Figure layout ------------------------------
if showDPE
    figH = 18.0;
else
    figH = 12;
end
hFig = figure('Units','centimeters','Position',[1,1,17.8,figH]);
% Otherwise the PDF is a figure printed onto a letter page, with a wide band of
% white below it. The margin is not decoration: several labels are placed
% outside the axes they belong to and overflow the figure box - "density" to the
% left of panel A, "Module 2 Performance" to the right of the trajectory panels
% - and on a letter page they used to land in the page margin. Crop to the
% figure box exactly and they are cut off.
mrg = [0.8, 0.5, 0.5, 0.3];                    % left, bottom, right, top, cm
set(hFig, 'PaperUnits', 'centimeters', 'PaperPositionMode', 'manual', ...
          'PaperPosition', [mrg(1), mrg(2), 17.8, figH], ...
          'PaperSize',     [17.8 + mrg(1) + mrg(3), figH + mrg(2) + mrg(4)]);
positions = defineSubplotPositions(showDPE);
addColumnAnnotations(column_titles, positions((1:3) + 3*showDPE), showDPE);

nPanel = 6 + 3*showDPE;
for i = 1:nPanel
    subplot('Position', positions{i});
    addSubplotLabel(subplot_labels{i}, positions{i});

    if showDPE && i <= 3
        drawDPEPanel(D, i, 'showKey', i == 1, 'showYLabel', i == 1, ...
                     'tickSize', 9, 'labelSize', 11);
        continue;
    end

    j   = i - 3*showDPE;              % 1-6 within the original six panels
    idx = mod(j-1, 3) + 1;
    if isempty(file_names{idx})
        drawEmptyPanel(dirs{idx});
        continue;
    end

    if j <= 3
        plotTrajectoryPanel(file_names{idx}, markers, showReference, ...
                            showGradientPath, showDPE);
    else
        plotLogRatioPanel(file_names{idx}, markers, yLimits);
    end
end

% The markers go on last. subplot('Position', ...) deletes any axes it overlaps,
% and these sit in the gap above the top row.
if showDPE
    for i = 1:3
        drawStartMarker(positions{i}, dpeMarkers{i});
    end
end

figDir = fullfile(resultsRoot, 'Figures');
if ~isfolder(figDir), mkdir(figDir); end
print(fullfile(figDir, [outputFile '.pdf']), '-dpdf', '-vector');
fprintf('Saved %s\n', fullfile(figDir, [outputFile '.pdf']));
end

% ======================================================================
% Panels
% ======================================================================

function plotTrajectoryPanel(file, markers, showReference, showGradientPath, threeRow)
    d  = load(file);
    sp = d.simParams;
    av = getAverageTrajectory(d);

    plotFitnessContours(sp); hold on;
    plotReferenceRays(sp);

    % Two passes, so the layer order is global rather than per-trajectory:
    % every reference curve first, every simulated mean on top. The reference is
    % suppressed unless asked for - it belongs to the unrestricted model, not
    % this one.
    % Path of steepest ascent through each initial phenotype, drawn BLACK
    % DASHED so it cannot be mistaken for a yellow model prediction. The curve
    % is x_2 = x_20 (x_1/x_10)^(a_1^2/a_2^2), equation "FGM trait path": it
    % follows from the phenotype-fitness map alone and is therefore the same
    % reference for the restricted and unrestricted GPMs. What differs between
    % them is whether a population follows it, which is exactly what the panel
    % is for.
    if showGradientPath && isfield(d, 'referenceTrajectories')
        for j = 1:numel(d.referenceTrajectories)
            R = d.referenceTrajectories{j,1};
            if isempty(R) || size(R,1) < 2, continue; end
            R = R(all(isfinite(R), 2), :);
            plot(R(:,1), R(:,2), '--', 'Color', [0 0 0], 'LineWidth', 1.1);
        end
    end

    if showReference && isfield(d, 'referenceTrajectories')
        for j = 1:numel(d.referenceTrajectories)
            R = d.referenceTrajectories{j,1};
            if isempty(R) || size(R,1) < 2, continue; end
            R = R(all(isfinite(R), 2), :);
            plot(R(:,1), R(:,2), '-', 'Color', '#EDB120', 'LineWidth', 2);
        end
    end

    % Replicate spread around the mean path, drawn the same way as Figures 2-5:
    % each replicate's deviation is projected onto the axis perpendicular to the
    % local direction of travel, so the band shows variation in the SHAPE of the
    % path rather than in how far along it a replicate happens to be.
    resultTable = getResultTable(d);
    for j = 1:numel(av.averageTimeStamp)
        ts = av.averageTimeStamp{j};
        bandData = computeTrajectoryPerpendicularBand(resultTable(j,:), size(ts, 2));
        fill([bandData.upperX1, fliplr(bandData.lowerX1)], ...
             [bandData.upperX2, fliplr(bandData.lowerX2)], ...
             [0.6118 0.7333 0.8431], 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    end

    for j = 1:numel(av.averageTimeStamp)
        ts = av.averageTimeStamp{j};
        plot(ts(1,:), ts(2,:), '-', 'Color', '#2776A8', 'LineWidth', 1.0);
        scatter(ts(1,1), ts(2,1), 70, markers{min(j, numel(markers))}, ...
                'MarkerEdgeColor', '#2776A8', 'MarkerFaceColor', '#2776A8');
    end

    if threeRow, tickFS = 9; else, tickFS = 10; end
    customizeAxes(1.2.*[-2.8, 0.05], 1.2.*[-2.1875, 0.05], tickFS);
    if threeRow
        % Left to itself MATLAB ticks these panels every 0.5 at this size. The
        % y tick labels are drawn to the RIGHT of the axis, which sits at the
        % right edge of the panel, so "-2.5" is half again as wide as "-2" and
        % eats the space "Module 2 Performance" needs.
        set(gca, 'XTick', -3:1:-1, 'YTick', -2:1:-1);
    end
    placeAxisTitles(threeRow);
end

function placeAxisTitles(threeRow)
% "Module 1 Performance" and "Module 2 Performance" are axis titles, but the
% axes are drawn at the origin, which in these panels is the top right corner,
% so they cannot be xlabel and ylabel and are placed by hand.
%
% Placing them at fixed DATA coordinates works only at one panel size. The tick
% labels they have to clear are sized in points and do not scale with the panel,
% and "Module 1 Performance" left-aligned from x = -3 is nearly as long as the
% panel is wide, so it also overhangs the right edge. In the nine-panel layout
% both show. There they are placed in normalized axes units instead, centred on
% the panel and offset far enough to clear a 10 pt tick label, at 11 pt.
    if ~threeRow
        % Unchanged, so the six-panel figure is pixel for pixel what it was.
        text(-3, 0.4, 'Module 1 Performance', 'FontName','Helvetica', 'FontSize',12);
        text(0.4, -0.1, 'Module 2 Performance', 'FontName','Helvetica', ...
             'FontSize',12, 'Rotation',270);
        text(0.08, 0.12, '0', 'FontName','Helvetica', 'FontSize',12);
        return;
    end
    text(0.45, 1.19, 'Module 1 Performance', 'Units','normalized', ...
         'FontName','Helvetica', 'FontSize',10, ...
         'HorizontalAlignment','center', 'VerticalAlignment','middle');
    text(1.18, 0.50, 'Module 2 Performance', 'Units','normalized', ...
         'FontName','Helvetica', 'FontSize',10, 'Rotation',270, ...
         'HorizontalAlignment','center', 'VerticalAlignment','middle');
    text(1.008, 1.03, '0', 'Units','normalized', ...
         'FontName','Helvetica', 'FontSize',10);
end

function plotLogRatioPanel(file, markers, yLimits)
    d  = load(file);
    sp = d.simParams;
    resultTable = getResultTable(d);

    numTimeStamp = 200;
    numSims      = size(resultTable, 2);
    maxGenAll    = 0;

    % Same rule as Figures 2-4 and 6: stop displaying log R once more than 20%
    % of replicates have dropped out. Replicates do not leave at random - the
    % ones that persist are those whose trait 1 is still far from its optimum -
    % so the surviving mean drifts upward and manufactures an apparent late-time
    % rise in R. Requiring 80% retention shortens the curves and removes it.
    minSimsCutoff = max(2, ceil(0.80 * numSims));

    hold on;

    for j = 1:size(resultTable, 1)
        maxGen = 0;
        for k = 1:numSims
            traj = resultTable{j, k};
            if ~isempty(traj), maxGen = max(maxGen, max(traj(:, 2))); end
        end
        if maxGen <= 1, continue; end

        genStamps   = floor(linspace(1, maxGen, numTimeStamp));
        logRatioAll = NaN(numSims, numTimeStamp);

        for k = 1:numSims
            traj = resultTable{j, k};
            if isempty(traj), continue; end
            simTime = traj(:, 2);
            for t = 1:numTimeStamp
                if genStamps(t) > max(simTime), break; end
                idx = find(simTime <= genStamps(t), 1, 'last');
                if isempty(idx), continue; end
                x1val = traj(idx, 3);
                x2val = traj(idx, 4);
                % Same rule as Figures 2-4 and 6, and as stated in Methods:
                % exclude data points once either trait is within delta of its
                % optimum, where log R is dominated by lattice granularity.
                if abs(x1val) < sp.deltaTrait || abs(x2val) < sp.deltaTrait
                    break;
                end
                % Guard against the floating-point residue that repeated lattice
                % additions leave in place of an exact zero, and against sign
                % changes once a trait overshoots its optimum.
                if x1val < -1e-9 && x2val < -1e-9
                    logRatioAll(k, t) = log(abs(x2val) / abs(x1val));
                end
            end
        end

        numContributing = sum(~isnan(logRatioAll), 1);
        valid           = numContributing >= minSimsCutoff;

        meanLogRatio = mean(logRatioAll, 1, 'omitnan');
        stdLogRatio  = std(logRatioAll, 0, 1, 'omitnan');
        meanLogRatio(~valid) = NaN;
        stdLogRatio(~valid)  = NaN;

        if any(valid), maxGenAll = max(maxGenAll, max(genStamps(valid))); end

        ok = ~isnan(meanLogRatio);
        if any(ok)
            fill([genStamps(ok), fliplr(genStamps(ok))], ...
                 [meanLogRatio(ok) + stdLogRatio(ok), fliplr(meanLogRatio(ok) - stdLogRatio(ok))], ...
                 [0.15, 0.46, 0.66], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
            plot(genStamps(ok), meanLogRatio(ok), '-', 'Color', '#2776A8', 'LineWidth', 0.9);
            first = find(ok, 1, 'first');
            plot(genStamps(first), meanLogRatio(first), ...
                 'Marker', markers{min(j, numel(markers))}, ...
                 'Color', '#2776A8', 'MarkerFaceColor', '#2776A8', 'MarkerSize', 7);
        end
    end

    if maxGenAll <= 1, maxGenAll = 100; end

    % Module-selection balance of the MODULAR GPFM. The question this figure
    % asks is whether the restricted cone reaches it without modular encoding.
    yline(log(sp.ellipseParams(2)^2 / sp.ellipseParams(1)^2), ...
          '-', 'LineWidth', 2, 'Color', [0.8, 0.3, 0]);
    yline(0, '-', 'LineWidth', 1.4, 'Color', [0.4 0.4 0.4]);

    xlim([1, maxGenAll]); ylim(yLimits);
    set(gca, 'TickLabelInterpreter','latex', 'FontSize', 10, 'Box', 'on');
    xlabel('Generations', 'FontName','Helvetica', 'FontSize', 12);
    ylabel('$\log(x_2/x_1)$', 'Interpreter', 'latex', 'FontSize', 12);
end

function drawEmptyPanel(regimeName)
    axis off;
    text(0.5, 0.5, sprintf('%s\nnot run', strrep(regimeName, '_', '\_')), ...
         'Units', 'normalized', 'HorizontalAlignment', 'center', ...
         'FontName', 'Helvetica', 'FontSize', 11, 'Color', [0.55 0.55 0.55]);
end

% ======================================================================
% Helpers
% ======================================================================

function plotFitnessContours(sp)
% Initial and final fitness contours, drawn numerically. The main-tree figure
% scripts obtain the inner level symbolically; doing it numerically here keeps
% this extension free of a Symbolic Math Toolbox dependency.
    W0 = 0.25;
    if isfield(sp, 'initialFitness') && ~isempty(sp.initialFitness)
        W0 = sp.initialFitness;
    end

    % Evaluated on an explicit grid rather than with fcontour, which keeps the
    % level list unambiguous and avoids a dependency on a newer graphics
    % function than the rest of this extension needs.
    g  = linspace(-4, 1.2, 400);
    [X1, X2] = meshgrid(g, g);
    W  = exp(-((X1./sp.ellipseParams(1)).^2 + (X2./sp.ellipseParams(2)).^2) ...
             ./ (2 * sp.landscapeStdDev^2));
    contour(X1, X2, W, sort([W0, 0.99]), 'LineColor', [0.7 0.7 0.7], 'LineWidth', 1.4);
end

function plotReferenceRays(sp)
    x1 = linspace(-3, 0.5, 1000);
    plot(x1, x1, '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 1);
    text(-2.3, -2.4, '$\mathbf{x_1 = x_2}$', 'Interpreter','latex', ...
         'FontSize',12, 'Color',[0.4 0.4 0.4]);

    R_bar = sp.ellipseParams(2)^2 / sp.ellipseParams(1)^2;
    plot(x1, R_bar * x1, '-', 'Color', [0.8, 0.3, 0], 'LineWidth', 2);
    text(-3.3, -1.8, '$\mathbf{s_1 = s_2}$', 'Interpreter','latex', ...
         'FontSize',12, 'Color',[0.8,0.3,0]);
end

function av = getAverageTrajectory(d)
    if isfield(d, 'averageTrajectory'), av = d.averageTrajectory;
    elseif isfield(d, 'ave'),           av = d.ave;
    else, error('makeFigure_RestrictedTheta:NoTrajectory', 'No trajectory data in file.');
    end
end

function resultTable = getResultTable(d)
    fn = fieldnames(d);
    hits = fn(contains(fn, 'result', 'IgnoreCase', true));
    for i = 1:numel(hits)
        if isstruct(d.(hits{i})) && isfield(d.(hits{i}), 'resultTable')
            resultTable = d.(hits{i}).resultTable;
            return;
        end
    end
    error('makeFigure_RestrictedTheta:NoResultTable', 'Could not find resultTable in file.');
end

function customizeAxes(xl, yl, fs)
    if nargin < 3, fs = 10; end
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    ax.Box       = 'off';
    ax.XColor    = 'k';
    ax.YColor    = 'k';
    ax.LineWidth = 1;
    xlim(xl); ylim(yl);
    set(ax, 'TickLabelInterpreter', 'latex', 'FontSize', fs);
end

function pos = defineSubplotPositions(showDPE)
    x = [0.06, 0.38, 0.70];
    w = 0.25;
    if ~showDPE
        % Unchanged from before the top row was added.
        s = 10/12;
        pos = { [x(1), 0.56*s, w, 0.38*s];
                [x(2), 0.56*s, w, 0.38*s];
                [x(3), 0.56*s, w, 0.38*s];
                [x(1), 0.08*s, w, 0.38*s];
                [x(2), 0.08*s, w, 0.38*s];
                [x(3), 0.08*s, w, 0.38*s] };
        return;
    end
    % Three rows in a 17.8 x 18.0 cm figure. The trajectory and log-ratio panels
    % keep roughly the physical size they had in the two-row layout; the
    % histogram row is shorter, since it carries no aspect ratio of its own.
    %
    % The gap between rows 1 and 2 is the wide one, and it has three things to
    % hold, bottom to top: the "Module 1 Performance" label, which sits ABOVE
    % the trajectory panel it belongs to rather than inside it; the column
    % titles, which wrap onto two lines; and the histogram row's x label. Each
    % needs about 0.05 of figure height, which is what sets the 0.180 gap.
    % Narrower panels than the two-row layout, to open the gaps between
    % columns: "Module 2 Performance" hangs about 0.9 cm off the right of each
    % trajectory panel and has to clear the next column.
    x    = [0.055, 0.385, 0.715];
    w    = 0.235;
    yRow = [0.760, 0.355, 0.055];
    hRow = [0.165, 0.225, 0.225];
    pos  = cell(9, 1);
    for r = 1:3
        for c = 1:3
            pos{3*(r-1) + c} = [x(c), yRow(r), w, hRow(r)];
        end
    end
end

function addColumnAnnotations(titles, rowPos, showDPE)
% Placed just above the row of trajectory panels, which is what they describe.
%
% These strings wrap onto two lines in a quarter-width box, and a textbox fills
% downward from its top edge, so they are anchored by their BOTTOM and grow
% upward. The offset clears the "Module 1 Performance" label, which sits ABOVE
% the trajectory panel it belongs to - plotReferenceRays places it at y = 0.4
% while the panel's y limit is 0.06 - and so occupies the first 0.03 of the gap.
    if ~showDPE
        for i = 1:numel(titles)
            p = rowPos{i};
            annotation('textbox', [p(1), 0.95, p(3), 0.05], 'String', titles{i}, ...
                       'FontSize', 14, 'FontName', 'Helvetica', ...
                       'HorizontalAlignment', 'center', 'EdgeColor', 'none');
        end
        return;
    end
    yTop = rowPos{1}(2) + rowPos{1}(4);
    for i = 1:numel(titles)
        p = rowPos{i};
        annotation('textbox', [p(1) - 0.03, yTop + 0.070, p(3) + 0.06, 0.085], ...
                   'String', titles{i}, ...
                   'FontSize', 12, 'FontName', 'Helvetica', ...
                   'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
                   'EdgeColor', 'none');
    end
end

function [D, mk] = getDPE(file_names, resultsRoot, sel, nGeno, useCache, markers)
% The top row, at the initial conditions the run itself used.
%
% sel indexes those initial conditions, so the marker returned for each panel is
% the marker drawn at the start of that trajectory in the panels below. The
% phenotypes and the cone come from the results file rather than from defaults,
% so the histograms describe the run in the figure and not a similar one.
    if isempty(which('computeDPEAtPhenotype'))
        % Added here so that the figures stage works from a bare path.
        sib = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'analysis_scripts');
        if isfolder(sib), addpath(sib); end
    end
    if isempty(which('computeDPEAtPhenotype'))
        error('makeFigure_RestrictedTheta:NoDPE', ...
              ['computeDPEAtPhenotype is not on the path; add ' ...
               'analysis_scripts, or pass ''showDPE'', false for the ' ...
               'six-panel figure.']);
    end

    src = file_names(~cellfun(@isempty, file_names));
    d   = load(src{1}, 'simParams', 'genomeParams', 'thetaRange');
    sp  = d.simParams;

    if ~isfield(sp, 'initialPhenotypes') || isempty(sp.initialPhenotypes)
        error('makeFigure_RestrictedTheta:NoPhenotypes', ...
              'No initialPhenotypes in %s.', src{1});
    end
    nCond = size(sp.initialPhenotypes, 1);
    if any(sel < 1) || any(sel > nCond)
        error('makeFigure_RestrictedTheta:BadConditions', ...
              ['''dpeConditions'' must index the %d initial conditions of ' ...
               'the run; got [%s].'], nCond, num2str(sel));
    end
    X0 = sp.initialPhenotypes(sel, :);
    mk = markers(min(sel, numel(markers)));

    tr = [0, pi/2];
    if isfield(d, 'thetaRange') && numel(d.thetaRange) == 2
        tr = d.thetaRange;
    end

    % The run's own loci, so the histograms describe that genome rather than
    % another draw from the same cone.
    if ~isfield(d, 'genomeParams') || ~isfield(d.genomeParams, 'genomeTheta')
        error('makeFigure_RestrictedTheta:NoGenome', ...
              'No genomeParams.genomeTheta in %s.', src{1});
    end
    theta = d.genomeParams.genomeTheta;

    args = {'theta', theta, 'delta', sp.deltaTrait, 'N', sp.popSize, ...
            'a', sp.ellipseParams, 'sigma', sp.landscapeStdDev, ...
            'thetaRange', tr, 'X0', X0, 'nGeno', nGeno};

    % The ensemble draws take a couple of minutes, so the result is kept beside
    % the figure and reused whenever the inputs above are identical.
    cacheFile = fullfile(resultsRoot, 'Figures', 'DPE_cache.mat');
    key = jsonencode({theta, sp.deltaTrait, sp.popSize, sp.ellipseParams, ...
                      sp.landscapeStdDev, tr, X0, nGeno});
    if useCache && isfile(cacheFile)
        C = load(cacheFile);
        if isfield(C, 'key') && strcmp(C.key, key)
            fprintf('  top row: reusing %s\n', cacheFile);
            D = C.D;
            return;
        end
    end

    D = computeDPEAtPhenotype(args{:});
    if useCache
        if ~isfolder(fileparts(cacheFile)), mkdir(fileparts(cacheFile)); end
        save(cacheFile, 'D', 'key');
        fprintf('  top row: cached to %s\n', cacheFile);
    end
end

function drawStartMarker(pos, mk)
% The initial condition a top-row panel belongs to, named by the marker that
% marks the start of its trajectory in the panels below rather than by its
% value of R_0.
    w  = 0.030;
    ax = axes('Position', [pos(1) + 0.5*pos(3) - w/2, pos(2) + pos(4) + 0.006, w, w]);
    plot(ax, 0, 0, mk, 'MarkerEdgeColor', '#2776A8', 'MarkerFaceColor', '#2776A8', ...
         'MarkerSize', 8);
    xlim(ax, [-1 1]); ylim(ax, [-1 1]);
    axis(ax, 'off');
end

function addSubplotLabel(label, pos)
    annotation('textbox', [pos(1)-0.055, pos(2)+pos(4)+0.012, 0.045, 0.045], ...
               'String', label, 'FontSize', 14, 'FontWeight', 'bold', 'EdgeColor', 'none');
end
