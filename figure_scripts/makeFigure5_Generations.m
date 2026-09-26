function makeFigure5_Generations(figureType, simParamsRef, varargin)
% makeFigure5_Generations — Generate Figure 5 (Nested FGM).
%
% Layout (3x3):
%   Row 1: A, B, C -- transcend regime (same FGM parameters apply to all
%          three regimes, so this row sits ABOVE the regime titles).
%          A: schematic (not drawn here, added by hand later)
%          B: fraction of mutations that are beneficial vs. distance to
%             the optimum, |x_i|, for fixed (m, n_i)
%          C: expected effect size of a beneficial mutation
%             (E[delta_i | delta_i > 0]) vs. distance to the optimum
%   Row 2: phase-plane trajectory plot (SSWM | CM linked | CM unlinked),
%          labeled by the regime titles directly above it
%   Row 3: log(x2/x1) vs. generations (SSWM | CM linked | CM unlinked)
%
% Usage:
%   makeFigure5_Generations('NestedFGM', simParams)
%   makeFigure5_Generations('NestedFGM', simParams, 'test')
%   makeFigure5_Generations('NestedFGM', simParams, 'full', 'proximityCutoff', 0.1)
%
% Inputs:
%   figureType   - 'NestedFGM'
%   simParamsRef - simParams struct used to match files when multiple exist in the results directory
%
% Optional positional input:
%   mode              - 'test', 'full', or 'auto'
%
% Optional name-value input:
%   'proximityCutoff' - Exclude ratio data when |x_i| < cutoff. Default: 0 (no cutoff).
%   'outputFile'      - Override default output filename (without extension).
%
% Notes:
%   Simulation trajectories are shown in blue. The orange line marks the
%   theoretical module-selection balance log(a2^2/a1^2).
%   If simParams.initialAngleIdx exists, markers are assigned using the
%   original indices from the master angle list so marker identity is
%   preserved across figures.
%
% Outputs:
%   Saves Figure_NestedFGM_Generations.pdf to results/Figures/
%   (or outputFile.pdf if outputFile override is specified)

% ----------------------------- Parse inputs -----------------------------
mode = 'auto';
proximityCutoff = 0;
outputFileOverride = '';   % if non-empty, overrides the default output filename
showEffectPanels = true;   % row A/B/C (schematic + mutation-effect panels)

idx = 1;
if nargin >= 3 && ischar(varargin{1}) && ismember(lower(varargin{1}), {'test','full','auto'})
    mode = lower(varargin{1});
    idx = 2;
end

while idx <= length(varargin) - 1
    key = varargin{idx};
    val = varargin{idx + 1};
    switch lower(key)
        case 'proximitycutoff'
            proximityCutoff = val;
        case 'outputfile'
            outputFileOverride = val;
        case 'showeffectpanels'
            % false -> original 2x3 layout (phase-plane + log-ratio only),
            % labelled A-F on a 17.8 x 12 cm canvas. Used for the asymmetric
            % nested-FGM supplementary figure, whose caption reads "same as
            % Figure 5 but with n1=10, n2=20" and for which the regime-agnostic
            % A/B/C row would be redundant.
            showEffectPanels = val;
        otherwise
            error('Unknown parameter: %s', key);
    end
    idx = idx + 2;
end

% ----------------------------- Settings --------------------------------
if showEffectPanels
    subplot_labels = {'A','B','C','D','E','F','G','H','I'};
else
    subplot_labels = {'A','B','C','D','E','F'};
end
markers        = {'d','^','v','>','<','o'};

% ----------------------- Resolve results root --------------------------
figureRoot  = fileparts(mfilename('fullpath'));
resultsRoot = resolveResultsRoot(figureRoot, mode);
fprintf('makeFigure5_Generations loading data from: %s\n', resultsRoot);

baseDirs = struct( ...
    'SSWM',       fullfile(resultsRoot, 'Generalization', 'NestedFGM', 'SSWM'), ...
    'CM_Asexual', fullfile(resultsRoot, 'Generalization', 'NestedFGM', 'CM_Asexual'), ...
    'CM_Sexual',  fullfile(resultsRoot, 'Generalization', 'NestedFGM', 'CM_Sexual'));

dirs          = {'SSWM', 'CM_Asexual', 'CM_Sexual'};
column_titles = {'Successive mutations', ...
                 'Concurrent mutations, linked modules', ...
                 'Concurrent mutations, unlinked modules'};
output_file   = 'Figure_NestedFGM';
regime_colors = {'#2776A8', '#2776A8', '#2776A8'};

% --------------------------- Locate files ------------------------------
file_names = cell(1,3);
for i = 1:3
    dirPath = baseDirs.(dirs{i});
    files = dir(fullfile(dirPath, sprintf('%s_%s_*.mat', figureType, dirs{i})));

    if isempty(files)
        error('No matching file for %s in %s.', figureType, dirPath);
    end

    if isscalar(files)
        file_names{i} = fullfile(files(1).folder, files(1).name);
    else
        file_names{i} = pickBestMatchingFile(files, simParamsRef);
    end

    fprintf('  Loading: %s\n', file_names{i});
end

% --------------------------- Figure layout -----------------------------
% A/B/C transcend regime (same underlying FGM parameters apply to all
% three), so they sit ABOVE the regime titles, at the very top of the
% figure. The regime titles sit directly above rows 2-3, which are the
% regime-specific panels.
tmp = load(file_names{1}, 'simParams');
sharedSimParams = tmp.simParams;   % m, moduleDimension -- shared across regimes

if showEffectPanels
    figure('Units','centimeters','Position',[1,1,17.8,17.8]);
else
    figure('Units','centimeters','Position',[1,1,17.8,12]);
end
addColumnAnnotations(column_titles, showEffectPanels);
subplot_positions = defineSubplotPositions(showEffectPanels);

for i = 1:numel(subplot_positions)
    subplot('Position', subplot_positions{i});
    addSubplotLabel(subplot_labels{i}, subplot_positions{i});

    % 1 = A/B/C (regime-agnostic), 2 = phase-plane, 3 = log-ratio.
    % Without the effect panels the grid is 2x3 and row 1 is skipped.
    row = ceil(i / 3);
    if ~showEffectPanels
        row = row + 1;
    end
    col = mod(i-1, 3) + 1;

    if row == 1
        if col == 1
            % A: reserved for schematic (added by hand later). Deliberately
            % blank -- just claim the layout space and show the panel label.
            axis off;
        elseif col == 2
            plotFractionBeneficial(sharedSimParams);
        else
            plotTypicalBeneficialEffect(sharedSimParams);
        end
    elseif row == 2
        plotFirstThreeSubplots_Nested(file_names{col}, markers, regime_colors{col});
    else
        plotLastThreeSubplots_Nested_Generations(file_names{col}, markers, proximityCutoff);
    end
end

figDir = fullfile(resultsRoot, 'Figures');
if ~isfolder(figDir)
    mkdir(figDir);
end

if ~isempty(outputFileOverride)
    saveFile = outputFileOverride;
else
    saveFile = [output_file '_Generations'];
end

print(fullfile(figDir, [saveFile '.pdf']), '-dpdf', '-vector');
fprintf('Saved %s\n', fullfile(figDir, [saveFile '.pdf']));
end

% ======================================================================
% Helper subfunctions
% ======================================================================

function resultsRoot = resolveResultsRoot(figureRoot, mode)
    switch mode
        case 'test'
            candidates = {fullfile(figureRoot,'../results_test'), fullfile(figureRoot,'../../results_test')};
        case 'full'
            candidates = {fullfile(figureRoot,'../results'), fullfile(figureRoot,'../../results')};
        otherwise
            candidates = {fullfile(figureRoot,'../results'), fullfile(figureRoot,'../../results'), ...
                          fullfile(figureRoot,'../results_test'), fullfile(figureRoot,'../../results_test')};
    end

    for i = 1:numel(candidates)
        if isfolder(candidates{i})
            resultsRoot = candidates{i};
            return;
        end
    end

    error('No results directory found for mode ''%s''.', mode);
end

function picked = pickBestMatchingFile(files, simParamsRef)
    picked = '';

    tags = {};

    if isfield(simParamsRef, 'popSize')
        tags{end+1} = sprintf('N%.0e', simParamsRef.popSize);
    end
    if isfield(simParamsRef, 'mutationRate')
        tags{end+1} = sprintf('M%.1e', simParamsRef.mutationRate);
    end
    if isfield(simParamsRef, 'deltaTrait')
        tags{end+1} = sprintf('d%.2f', simParamsRef.deltaTrait);
    end
    if isfield(simParamsRef, 'ellipseParams')
        tags{end+1} = sprintf('eR%.2f', simParamsRef.ellipseParams(1) / simParamsRef.ellipseParams(2));
    end
    if isfield(simParamsRef, 'landscapeStdDev')
        tags{end+1} = sprintf('s%.2f', simParamsRef.landscapeStdDev);
    end
    if isfield(simParamsRef, 'moduleDimension')
        tags{end+1} = sprintf('n%d-%d', simParamsRef.moduleDimension(1), simParamsRef.moduleDimension(2));
    end
    if isfield(simParamsRef, 'recombinationRate') && simParamsRef.recombinationRate > 0
        tags{end+1} = sprintf('R%.4f', simParamsRef.recombinationRate);
    end

    scores = zeros(1, numel(files));
    for f = 1:numel(files)
        name = files(f).name;
        for t = 1:numel(tags)
            if contains(name, tags{t})
                scores(f) = scores(f) + 1;
            end
        end
    end

    bestIdx = find(scores == max(scores));
    if isscalar(bestIdx)
        picked = fullfile(files(bestIdx).folder, files(bestIdx).name);
        return;
    end

    % A tie used to fall through to "most recent" with a warning, which is easy
    % to miss in a long run log. Superseded runs belong in a _relegated/
    % subfolder, which dir() does not match; anything still tied here is a
    % genuine ambiguity the caller has to resolve.
    tied = {files(bestIdx).name};
    error('makeFigure5_Generations:AmbiguousResults', ...
          ['%d result files match the requested parameters equally well:\n  %s\n' ...
           'Move the superseded ones into a _relegated/ subfolder.'], ...
          numel(tied), strjoin(tied, sprintf('\n  ')));
end

function pos = defineSubplotPositions(showEffectPanels)
    % 3x3 grid on a 17.8 x 17.8 cm canvas. Panel size (4.45 x 3.8 cm) is
    % unchanged from the original 2x3 version. Bottom-up in cm:
    %   0.0 - 0.8   bottom margin
    %   0.8 - 4.6   Row 3 (log-ratio, G/H/I)
    %   4.6 - 5.6   gap (1.0 cm)
    %   5.6 - 9.4   Row 2 (phase-plane, D/E/F)
    %   9.4 - 12.0  regime-title zone (2.0 cm gap + 0.6 cm title text) --
    %               titles label rows 2-3, which ARE regime-specific
    %  12.0 - 13.0  gap (1.0 cm)
    %  13.0 - 16.8  Row 1 (A/B/C, regime-agnostic -- sits above the titles
    %               since it applies to all three regimes at once)
    %  16.8 - 17.8  top margin (clearance for the A/B/C subplot labels)
    colW  = 0.25;    % 4.45 cm / 17.8 cm
    rowH  = 0.2135;  % 3.8 cm / 17.8 cm
    xCol1 = 0.06;
    xCol2 = 0.38;
    xCol3 = 0.70;
    yRow3 = 0.0449;  % 0.8 cm / 17.8 cm
    yRow2 = 0.3146;  % 5.6 cm / 17.8 cm
    yRow1 = 0.7303;  % 13.0 cm / 17.8 cm

    if ~showEffectPanels
        % Original 2x3 layout on a 17.8 x 12 cm canvas: phase-plane over
        % log-ratio, panels A-F, matching the pre-restructure figure.
        s = 10/12;
        pos = {
            [0.06, 0.56*s, 0.25, 0.38*s];
            [0.38, 0.56*s, 0.25, 0.38*s];
            [0.70, 0.56*s, 0.25, 0.38*s];
            [0.06, 0.08*s, 0.25, 0.38*s];
            [0.38, 0.08*s, 0.25, 0.38*s];
            [0.70, 0.08*s, 0.25, 0.38*s];
        };
        return;
    end

    pos = {
        [xCol1, yRow1, colW, rowH];   % A: blank (schematic)
        [xCol2, yRow1, colW, rowH];   % B: fraction of beneficial mutations vs. distance
        [xCol3, yRow1, colW, rowH];   % C: typical beneficial effect size vs. distance
        [xCol1, yRow2, colW, rowH];   % D: phase-plane, SSWM
        [xCol2, yRow2, colW, rowH];   % E: phase-plane, CM linked
        [xCol3, yRow2, colW, rowH];   % F: phase-plane, CM unlinked
        [xCol1, yRow3, colW, rowH];   % G: log-ratio, SSWM
        [xCol2, yRow3, colW, rowH];   % H: log-ratio, CM linked
        [xCol3, yRow3, colW, rowH];   % I: log-ratio, CM unlinked
    };
end

function addColumnAnnotations(titles, showEffectPanels)
    % Sits between Row 1 (A/B/C, regime-agnostic) and Row 2 (phase-plane,
    % regime-specific) -- these titles label rows 2-3, not row 1.
    x = [0.06, 0.38, 0.70];
    if showEffectPanels
        yTitle = 0.6404;   % between row 1 and row 2 of the 3x3 layout
    else
        yTitle = 0.95;     % top of the 2x3 layout, as before the restructure
    end
    for i = 1:numel(titles)
        annotation('textbox', [x(i), yTitle, 0.25, 0.0337], ...
            'String', titles{i}, ...
            'FontSize', 16, ...
            'FontName', 'Helvetica', ...
            'HorizontalAlignment', 'center', ...
            'EdgeColor', 'none');
    end
end

function addSubplotLabel(label, pos)
    annotation('textbox', [pos(1)-0.055, pos(2)+pos(4)+0.012, 0.045, 0.045], ...
        'String', label, ...
        'FontSize', 14, ...
        'FontWeight', 'bold', ...
        'EdgeColor', 'none');
end

function plotFirstThreeSubplots_Nested(file, markers, color)
    d = load(file);
    sp = d.simParams;
    av = getAverageTrajectory(d);
    resultTable = getResultTable(d);

    [~, lvl] = defineGaussianPDF(sp);
    plotPDFContour(sp, lvl, [0.7 0.7 0.7]);
    hold on;

    angleIdx = getInitialAngleIdx(sp, numel(av.averageTimeStamp));

    % Simulation averages, now with a perpendicular-projected error band
    % showing spread across replicate runs
    plotAverageTrajectories(av, markers, angleIdx, color, resultTable);

    % Reference lines, drawn last so they render on top of the band and
    % the simulated average line: gray x1=x2 line and the orange
    % theoretical module-selection balance line, both labeled -- matching
    % the Modular column in makeRegimeFigure_Generations.m.
    plotReferenceLines_Nested(sp);

    customizeAxes(1.2.*[-2.8, 0.05], 1.2.*[-2.1875, 0.05]);
    text(-3, 0.4, 'Module 1 Performance', 'FontName','Helvetica', 'FontSize',12);
    text(0.4, -0.1, 'Module 2 Performance', 'FontName','Helvetica', 'FontSize',12, 'Rotation',270);
    text(0.08, 0.12, '0', 'FontName','Helvetica', 'FontSize',12);
end

function plotReferenceLines_Nested(sp)
    x1 = linspace(-3, 0.5, 1000);
    plot(x1, x1, '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 1);
    text(-2.3, -2.4, '$\mathbf{x_1 = x_2}$', 'Interpreter','latex', 'FontSize',12, 'Color',[0.4 0.4 0.4]);

    Rstar = getNestedBalanceRatio(sp);
    xplot = linspace(1.2*(-2.8), 0, 200);
    yplot = Rstar .* xplot;
    plot(xplot, yplot, '-', 'Color', '#D95319', 'LineWidth', 1.5);
    text(-3.3, -1.8, '$\mathbf{s_1 = s_2}$', 'Interpreter','latex', 'FontSize',12, 'Color','#D95319');
end

function plotFractionBeneficial(sp)
    % Panel B: fraction of mutations that are beneficial, as a function of
    % distance to the optimum |x_i|, for fixed FGM parameters (m, n_i).
    %
    % Model (matches manuscript Section "Nested Fisher's geometric model"
    % and simulateNestedSSWM.m): a mutation's effect delta_i on module i is
    %   delta_i ~ Normal(mu_eff, sigma_i^2)
    %   mu_eff  = -m^2/2                        (fixed)
    %   sigma_i = m * sqrt(2*|x_i| / n_i)        (grows with distance)
    % x_i is negative distance-from-optimum, so delta_i > 0 is beneficial
    % (moves toward the optimum). This is module-agnostic here because
    % moduleDimension is symmetric (n1 == n2) for this figure's parameters
    % -- if n1 != n2 the two modules would trace different curves.
    m = sp.deltaTrait;
    n = sp.moduleDimension(1);
    mu_eff = -(m^2) / 2;

    distances = linspace(0, 1.2*2.8, 300);
    sigma = m .* sqrt(2 .* distances ./ n);
    fracBeneficial = normcdf(0, mu_eff, sigma, 'upper');

    plot(distances, fracBeneficial, '-', 'Color', 'k', 'LineWidth', 1.5);

    xlim([0, max(distances)]);
    ylim([0, 1]);
    set(gca, 'TickLabelInterpreter','latex', 'FontSize', 10, 'Box', 'off');
    xlabel('Distance to optimum, $|x_i|$', 'Interpreter','latex', 'FontSize', 12);
    ylabel('$\mathrm{Pr}(\delta_i > 0)$', 'Interpreter','latex', 'FontSize', 14);
end

function plotTypicalBeneficialEffect(sp)
    % Panel C: expected effect size of a BENEFICIAL mutation,
    % E[delta_i | delta_i > 0], as a function of distance to the
    % optimum |x_i|. This is the mean of the same Normal(mu_eff, sigma_i^2)
    % distribution as panel B, truncated to its beneficial (>0) tail.
    m = sp.deltaTrait;
    n = sp.moduleDimension(1);
    mu_eff = -(m^2) / 2;

    distances = linspace(0, 1.2*2.8, 300);
    sigma = m .* sqrt(2 .* distances ./ n);

    % Truncated-normal conditional mean above threshold 0:
    %   E[X | X>0] = mu + sigma * phi(alpha) / (1 - Phi(alpha)),  alpha = -mu/sigma
    alpha = -mu_eff ./ max(sigma, eps);
    condMean = mu_eff + sigma .* normpdf(alpha) ./ max(1 - normcdf(alpha), eps);
    condMean(sigma < 1e-9) = 0;   % limiting value as distance -> 0

    plot(distances, condMean, '-', 'Color', 'k', 'LineWidth', 1.5);

    xlim([0, max(distances)]);
    set(gca, 'TickLabelInterpreter','latex', 'FontSize', 10, 'Box', 'off');
    xlabel('Distance to optimum, $|x_i|$', 'Interpreter','latex', 'FontSize', 12);
    ylabel('$\mathrm{E}(\delta_i \mid \delta_i>0)$', 'Interpreter','latex', 'FontSize', 14);
end

function av = getAverageTrajectory(d)
    if isfield(d, 'averageTrajectory')
        av = d.averageTrajectory;
    elseif isfield(d, 'ave')
        av = d.ave;
    else
        error('No trajectory data found.');
    end
end

function angleIdx = getInitialAngleIdx(sp, nTraj)
    if isfield(sp, 'initialAngleIdx') && ~isempty(sp.initialAngleIdx)
        angleIdx = sp.initialAngleIdx(:)';
    else
        angleIdx = 1:nTraj;
    end

    if numel(angleIdx) ~= nTraj
        warning('Length of simParams.initialAngleIdx does not match trajectory count. Falling back to sequential markers.');
        angleIdx = 1:nTraj;
    end
end

function plotPDFContour(sp, lvl, lineColor)
    fc = fcontour(@(x1,x2) exp(-sqrt((x1./sp.ellipseParams(1)).^2 + (x2./sp.ellipseParams(2)).^2).^2 / ...
        (2 * sp.landscapeStdDev^2)), ...
        'LineColor', lineColor, 'LineWidth', 1.4);
    fc.LevelList = [0.99, lvl];
end

function plotAverageTrajectories(av, markers, angleIdx, color, resultTable)
    nMarkers = numel(markers);
    bandColor = hex2rgbLocal(color);

    % Two passes so the layer order is global: bands underneath, all
    % simulated means and their start markers on top.
    for j = 1:numel(av.averageTimeStamp)
        ts = av.averageTimeStamp{j};
        numTimeStamp = size(ts, 2);
        bandData = computeTrajectoryPerpendicularBand(resultTable(j,:), numTimeStamp);
        fill([bandData.upperX1, fliplr(bandData.lowerX1)], ...
             [bandData.upperX2, fliplr(bandData.lowerX2)], ...
             bandColor, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    end

    for j = 1:numel(av.averageTimeStamp)
        ts = av.averageTimeStamp{j};
        plot(ts(1,:), ts(2,:), '-', 'Color', color, 'LineWidth', 1.0);
        mkIdx = angleIdx(j);
        if mkIdx < 1 || mkIdx > nMarkers
            warning('Marker index %d is out of range. Using sequential fallback.', mkIdx);
            mkIdx = min(j, nMarkers);
        end
        scatter(ts(1,1), ts(2,1), 70, markers{mkIdx}, ...
            'MarkerEdgeColor', color, 'MarkerFaceColor', color);
    end
end

function rgb = hex2rgbLocal(hexStr)
    % fill() requires an RGB triplet, not a hex string, when passed
    % positionally.
    hexStr = strrep(hexStr, '#', '');
    rgb = [hex2dec(hexStr(1:2)), hex2dec(hexStr(3:4)), hex2dec(hexStr(5:6))] / 255;
end

function plotLastThreeSubplots_Nested_Generations(file, markers, proximityCutoff)
    d = load(file);
    sp = d.simParams;
    resultTable = getResultTable(d);

    angleIdx = getInitialAngleIdx(sp, size(resultTable, 1));

    numTimeStamp = 200;
    maxGenAll = 0;
    numSims = size(resultTable, 2);
    % Stop displaying log R once more than 20% of replicates have dropped out
    % (either finished, or excluded by proximityCutoff). Replicates do not drop
    % out at random: the ones that persist longest are those whose trait 1 is
    % still far from the optimum, which biases the surviving mean upward and
    % produces an apparent late-time increase in R. Requiring 80% retention
    % shortens the curves but removes that bias.
    %
    % NOTE: the previous rule, min(40, floor(0.5*numSims)), was capped at 40 and
    % so retained only 40/250 = 16% at paper scale, not the 50% it appears to say.
    minRetainedFraction = 0.80;
    minSimsCutoff = max(2, ceil(minRetainedFraction * numSims));

    hold on;

    for j = 1:size(resultTable, 1)
        maxGen = 0;
        for k = 1:numSims
            traj = resultTable{j, k};
            if ~isempty(traj)
                maxGen = max(maxGen, max(traj(:, 2)));
            end
        end

        if maxGen <= 1
            continue;
        end

        genStamps = floor(linspace(1, maxGen, numTimeStamp));
        logRatioAll = NaN(numSims, numTimeStamp);

        for k = 1:numSims
            traj = resultTable{j, k};
            if isempty(traj)
                continue;
            end

            simTime = traj(:, 2);
            trait1 = traj(:, 3);
            trait2 = traj(:, 4);
            simMaxGen = max(simTime);

            for t = 1:numTimeStamp
                if genStamps(t) <= simMaxGen
                    idx = find(simTime <= genStamps(t), 1, 'last');
                    if isempty(idx)
                        idx = 1;
                    end

                    x1val = trait1(idx);
                    x2val = trait2(idx);

                    if proximityCutoff > 0 && (abs(x1val) < proximityCutoff || abs(x2val) < proximityCutoff)
                        break;
                    end

                    if x1val ~= 0 && x2val ~= 0
                        logRatioAll(k, t) = log(abs(x2val) / abs(x1val));
                    end
                end
            end
        end

        numContributing = sum(~isnan(logRatioAll), 1);
        validTimeIdx = numContributing >= minSimsCutoff;

        meanLogRatio = mean(logRatioAll, 1, 'omitnan');
        stdLogRatio  = std(logRatioAll, 0, 1, 'omitnan');

        meanLogRatio(~validTimeIdx) = NaN;
        stdLogRatio(~validTimeIdx)  = NaN;

        validGens = genStamps(validTimeIdx);
        if ~isempty(validGens)
            maxGenAll = max(maxGenAll, max(validGens));
        end

        validIdx = ~isnan(meanLogRatio);
        if any(validIdx)
            xFill = [genStamps(validIdx), fliplr(genStamps(validIdx))];
            yFill = [meanLogRatio(validIdx) + stdLogRatio(validIdx), ...
                     fliplr(meanLogRatio(validIdx) - stdLogRatio(validIdx))];
            fill(xFill, yFill, [0.15, 0.46, 0.66], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        end

        plot(genStamps(validIdx), meanLogRatio(validIdx), '-', ...
            'Color', '#2776A8', 'LineWidth', 0.9);

        if any(validIdx)
            firstValid = find(validIdx, 1, 'first');

            mkIdx = angleIdx(j);
            if mkIdx < 1 || mkIdx > numel(markers)
                warning('Marker index %d is out of range. Using sequential fallback.', mkIdx);
                mkIdx = min(j, numel(markers));
            end

            plot(genStamps(firstValid), meanLogRatio(firstValid), ...
                'Marker', markers{mkIdx}, ...
                'Color', '#2776A8', ...
                'MarkerFaceColor', '#2776A8', ...
                'MarkerSize', 7);
        end
    end

    if maxGenAll <= 1
        maxGenAll = 100;
    end

    % Reference line and theoretical module-selection balance
    Rstar = getNestedBalanceRatio(sp);
    yline(log(Rstar), '-', 'Color', '#D95319', 'LineWidth', 1.5);

    xlim([1, maxGenAll]);
    ylim([-3, 2]);
    set(gca, 'TickLabelInterpreter','latex', 'FontSize', 10, 'Box', 'on');
    xlabel('Generations', 'FontName','Helvetica', 'FontSize', 12);
    ylabel('$\log(x_2/x_1)$', 'Interpreter', 'latex', 'FontSize', 12);
end

function [f, lvl] = defineGaussianPDF(sp)
    syms x1 x2
    f = exp(-sqrt((x1./sp.ellipseParams(1)).^2 + (x2./sp.ellipseParams(2)).^2).^2 / ...
        (2 * sp.landscapeStdDev^2));
    dH = det(hessian(f, [x1, x2]));
    dHf = matlabFunction(dH, 'Vars', [x1, x2]);

    xyStar = fminsearch(@(x) abs(dHf(x(1), x(2))), [0, 0]);
    lvl = double(subs(f, [x1, x2], xyStar));
end

function Rstar = getNestedBalanceRatio(sp)
    % Same balance condition as the simple modular model for the symmetric
    % nested case used in Figure 5.
    Rstar = (sp.ellipseParams(2)^2) / (sp.ellipseParams(1)^2);
end

function customizeAxes(xl, yl)
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    ax.Box = 'off';
    ax.XColor = 'k';
    ax.YColor = 'k';
    ax.LineWidth = 1;
    xlim(xl);
    ylim(yl);
    ax.XTick = [-2, -1];
    ax.YTick = [-2, -1];
    ax.TickLength = [0.015, 0.015];
    ax.FontSize = 10;
    ax.TickLabelInterpreter = 'latex';
end

function resultTable = getResultTable(d)
    fieldNames = fieldnames(d);
    resultFields = fieldNames(contains(fieldNames, 'result', 'IgnoreCase', true));

    for i = 1:length(resultFields)
        fn = resultFields{i};
        if isstruct(d.(fn)) && isfield(d.(fn), 'resultTable')
            resultTable = d.(fn).resultTable;
            return;
        end
    end

    error('Could not find resultTable in data file.');
end