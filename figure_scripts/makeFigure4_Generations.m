function makeFigure4_Generations(figureType, simParamsRef, varargin)
% makeFigure4_Generations — Generate Figure 4 (Discordant-module FGM).
%
% Usage:
%   makeFigure4_Generations('DiscordantFGM', simParams)
%   makeFigure4_Generations('DiscordantFGM', simParams, 'test')
%   makeFigure4_Generations('DiscordantFGM', simParams, 'full', 'proximityCutoff', 0.1)
%
% Inputs:
%   figureType   - 'DiscordantFGM'
%   simParamsRef - simParams struct used to match files when multiple exist in the results directory
%
% Optional positional input:
%   mode              - 'test', 'full', or 'auto'
%
% Optional name-value input:
%   'proximityCutoff' - Exclude ratio data when |x_i| < cutoff. Default: 0 (no cutoff).
%
% Notes:
%   Shows simulation trajectories with theoretical attractor overlay (orange line).
%   If simParams.initialAngleIdx exists, markers are assigned using the
%   original indices from the master angle list so marker identity is
%   preserved across figures.
%
% Outputs:
%   Saves Figure_DiscordantFGM_Generations.eps to results/Figures/

% ----------------------------- Parse inputs -----------------------------
mode = 'auto';
proximityCutoff = 0;

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
        otherwise
            error('Unknown parameter: %s', key);
    end
    idx = idx + 2;
end

% ----------------------------- Settings --------------------------------
subplot_labels = {'A','B','C','D','E','F'};
markers        = {'d','^','v','>','<','o'};

% ----------------------- Resolve results root --------------------------
figureRoot  = fileparts(mfilename('fullpath'));
resultsRoot = resolveResultsRoot(figureRoot, mode);
fprintf('makeFigure4_Generations loading data from: %s\n', resultsRoot);

baseDirs = struct( ...
    'SSWM',       fullfile(resultsRoot, 'Generalization', 'DiscordantFGM', 'SSWM'), ...
    'CM_Asexual', fullfile(resultsRoot, 'Generalization', 'DiscordantFGM', 'CM_Asexual'), ...
    'CM_Sexual',  fullfile(resultsRoot, 'Generalization', 'DiscordantFGM', 'CM_Sexual'));

dirs          = {'SSWM', 'CM_Asexual', 'CM_Sexual'};
column_titles = {'Successive mutations', ...
                 'Concurrent mutations, linked modules', ...
                 'Concurrent mutations, unlinked modules'};
output_file   = 'Figure_DiscordantFGM';
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
figure('Units','centimeters','Position',[1,1,17.8,12]);
addColumnAnnotations(column_titles);
subplot_positions = defineSubplotPositions(10/12);

for i = 1:6
    subplot('Position', subplot_positions{i});
    addSubplotLabel(subplot_labels{i}, subplot_positions{i});

    if i <= 3
        plotFirstThreeSubplots_Discordant(file_names{i}, markers, regime_colors{i});
    else
        plotLastThreeSubplots_Discordant_Generations(file_names{i-3}, markers, proximityCutoff);
    end
end

figDir = fullfile(resultsRoot, 'Figures');
if ~isfolder(figDir)
    mkdir(figDir);
end

print(fullfile(figDir, [output_file '_Generations.eps']), '-depsc', '-r300');
fprintf('Saved %s\n', fullfile(figDir, [output_file '_Generations.eps']));
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
    if isfield(simParamsRef, 'geneticTargetSize')
        tags{end+1} = sprintf('L%d-%d', simParamsRef.geneticTargetSize(1), simParamsRef.geneticTargetSize(2));
    end
    if isfield(simParamsRef, 'recombinationRate') && simParamsRef.recombinationRate > 0
        tags{end+1} = sprintf('R%.4f', simParamsRef.recombinationRate);
    end
    if isfield(simParamsRef, 'discordantAngles') && ~isempty(simParamsRef.discordantAngles)
        tags{end+1} = sprintf('T%.4f-%.4f', simParamsRef.discordantAngles(1), simParamsRef.discordantAngles(2));
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
    else
        [~, newestLocal] = max([files(bestIdx).datenum]);
        picked = fullfile(files(bestIdx(newestLocal)).folder, files(bestIdx(newestLocal)).name);
        warning('Multiple matching files found. Using most recent: %s', files(bestIdx(newestLocal)).name);
    end
end

function pos = defineSubplotPositions(s)
    pos = {
        [0.06, 0.56*s, 0.25, 0.38*s];
        [0.38, 0.56*s, 0.25, 0.38*s];
        [0.70, 0.56*s, 0.25, 0.38*s];
        [0.06, 0.08*s, 0.25, 0.38*s];
        [0.38, 0.08*s, 0.25, 0.38*s];
        [0.70, 0.08*s, 0.25, 0.38*s];
    };
end

function addColumnAnnotations(titles)
    x = [0.06, 0.38, 0.70];
    for i = 1:numel(titles)
        annotation('textbox', [x(i), 0.95, 0.25, 0.05], ...
            'String', titles{i}, ...
            'FontSize', 14, ...
            'FontName', 'Helvetica', ...
            'HorizontalAlignment', 'center', ...
            'EdgeColor', 'none');
    end
end

function addSubplotLabel(label, pos)
    annotation('textbox', [pos(1)-0.04, pos(2)+pos(4), 0.03, 0.03], ...
        'String', label, ...
        'FontSize', 12, ...
        'FontWeight', 'bold', ...
        'EdgeColor', 'none');
end

function plotFirstThreeSubplots_Discordant(file, markers, color)
    d = load(file);
    sp = d.simParams;
    av = getAverageTrajectory(d);

    [~, lvl] = defineGaussianPDF(sp);
    plotPDFContour(sp, lvl, [0.7 0.7 0.7]);
    hold on;

    % ---------------- Axis limits used in this panel ----------------
    xlimVals = 1.2 .* [-2.8, 0.05];
    ylimVals = 1.2 .* [-2.1875, 0.05];

    % ---------------- Discordant geometry ----------------
    theta1 = sp.discordantAngles(1);
    theta2 = sp.discordantAngles(2);

    R1 = tan(theta1);
    R2 = tan(theta2);
    [~, R3] = getDiscordantBalanceRatios(sp);

    xplot = linspace(xlimVals(1), 0, 200);

    % ---------------- Shade inaccessible regions ----------------
    % Accessible wedge is between x2 = R1*x1 and x2 = R2*x1.
    % Since x1 is negative in these panels, note that R2*x1 is more negative
    % than R1*x1, so:
    %   upper boundary = R1*x1
    %   lower boundary = R2*x1

    % Region above the upper boundary
    x_top = [xplot, fliplr(xplot)];
    y_top = [R1 .* xplot, repmat(ylimVals(2), size(xplot))];
    fill(x_top, y_top, [0.85 0.85 0.85], ...
        'EdgeColor', 'none', 'FaceAlpha', 0.5);

    % Region below the lower boundary
    x_bot = [xplot, fliplr(xplot)];
    y_bot = [repmat(ylimVals(1), size(xplot)), fliplr(R2 .* xplot)];
    fill(x_bot, y_bot, [0.85 0.85 0.85], ...
        'EdgeColor', 'none', 'FaceAlpha', 0.5);

    % ---------------- Plot simulation trajectories ----------------
    angleIdx = getInitialAngleIdx(sp, numel(av.averageTimeStamp));
    plotAverageTrajectories(av, markers, angleIdx, color);

    % ---------------- Plot boundaries and attractor ----------------
    plot(xplot, R3 .* xplot, '-',  'Color', '#D95319', 'LineWidth', 1.5);
    text(-3.2, R3 * (-3.2) - 0.2, '$\mathbf{x_2 = R^{*} x_1}$', ...
        'Interpreter','latex', ...
        'FontSize',10, ...
        'Color',[0.8 0.3 0]);

    customizeAxes(xlimVals, ylimVals);
    text(-3, 0.4, 'Module 1 Performance', 'FontName','Helvetica', 'FontSize',10);
    text(0.4, -0.1, 'Module 2 Performance', 'FontName','Helvetica', 'FontSize',10, 'Rotation',270);
    text(0.08, 0.12, '0', 'FontName','Helvetica', 'FontSize',10);

    % ---------------- y1/y2 chromosome direction inset ----------------
    % Draw two arrows from a common origin in the bottom-left corner
    % showing the directions of chromosomes 1 and 2 in x-space.
    ox = xlimVals(1) * 0.94;   % inset origin x
    oy = ylimVals(1) * 0.98;   % inset origin y
    arrowLen = 0.4;            % arrow length in data units

    dx1 = arrowLen * cos(theta1);  dy1 = arrowLen * sin(theta1);
    dx2 = arrowLen * cos(theta2);  dy2 = arrowLen * sin(theta2);

    quiver(ox, oy, dx1, dy1, 0, 'k', 'LineWidth', 1.2, 'MaxHeadSize', 0.5);
    quiver(ox, oy, dx2, dy2, 0, 'k', 'LineWidth', 1.2, 'MaxHeadSize', 0.5);

    text(ox + dx1*1.25, oy + dy1*1.25, '$y_1$', 'Interpreter','latex', ...
        'FontSize', 8, 'HorizontalAlignment','center');
    text(ox + dx2*1.25, oy + dy2*1.25, '$y_2$', 'Interpreter','latex', ...
        'FontSize', 8, 'HorizontalAlignment','center');
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

function plotAverageTrajectories(av, markers, angleIdx, color)
    nMarkers = numel(markers);

    for j = 1:numel(av.averageTimeStamp)
        ts = av.averageTimeStamp{j};
        plot(ts(1,:), ts(2,:), '-', 'Color', color, 'LineWidth', 1.0);

        mkIdx = angleIdx(j);
        if mkIdx < 1 || mkIdx > nMarkers
            warning('Marker index %d is out of range. Using sequential fallback.', mkIdx);
            mkIdx = min(j, nMarkers);
        end

        scatter(ts(1,1), ts(2,1), 30, markers{mkIdx}, ...
            'MarkerEdgeColor', color, 'MarkerFaceColor', color);
    end
end

function plotLastThreeSubplots_Discordant_Generations(file, markers, proximityCutoff)
    d = load(file);
    sp = d.simParams;
    resultTable = getResultTable(d);

    angleIdx = getInitialAngleIdx(sp, size(resultTable, 1));

    % Reconstruct M^{-1} for y-space delta-band filter
    theta1_disc = sp.discordantAngles(1);
    theta2_disc = sp.discordantAngles(2);
    M_disc = [cos(theta1_disc), cos(theta2_disc);
              sin(theta1_disc), sin(theta2_disc)];
    Minv   = inv(M_disc);

    numTimeStamp = 200;
    maxGenAll = 0;
    numSims = size(resultTable, 2);
    minSimsCutoff = min(40, max(2, floor(numSims * 0.5)));

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

                    % Filter in y-space: stop trajectory when either chromosome
                    % enters the delta-band near zero (b_k < 1, approximation breaks down)
                    yval = Minv * [x1val; x2val];
                    if proximityCutoff > 0 && (abs(yval(1)) < proximityCutoff || abs(yval(2)) < proximityCutoff)
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
                'MarkerSize', 5);
        end
    end

    if maxGenAll <= 1
        maxGenAll = 100;
    end

    % Discordant-module boundaries and theoretical attractor
    theta1 = sp.discordantAngles(1);
    theta2 = sp.discordantAngles(2);

    R1 = tan(theta1);
    R2 = tan(theta2);
    [~, R3] = getDiscordantBalanceRatios(sp);

    % y-axis limits — change here to update both shading and axis together
    ymin = -3;
    ymax =  2;

    % gray inaccessible zone below lower boundary
    fill([1, maxGenAll, maxGenAll, 1], ...
        [ymin, ymin, log(R1), log(R1)], ...
        [0.85 0.85 0.85], ...
        'EdgeColor', 'none', 'FaceAlpha', 0.5);

    % gray inaccessible zone above upper boundary
    fill([1, maxGenAll, maxGenAll, 1], ...
        [log(R2), log(R2), ymax, ymax], ...
        [0.85 0.85 0.85], ...
        'EdgeColor', 'none', 'FaceAlpha', 0.5);

    % theoretical attractor only
    yline(log(R3), '-', 'Color', '#D95319', 'LineWidth', 1.5);

    xlim([1, maxGenAll]);
    ylim([ymin, ymax]);
    set(gca, 'TickLabelInterpreter','latex', 'FontSize', 8, 'Box', 'on');
    xlabel('Generations', 'FontName','Helvetica', 'FontSize', 10);
    ylabel('$\log(x_2/x_1)$', 'Interpreter', 'latex', 'FontSize', 10);
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

function [R1, R3] = getDiscordantBalanceRatios(sp)
    theta1 = sp.discordantAngles(1);
    theta2 = sp.discordantAngles(2);
    a1 = sp.ellipseParams(1);
    a2 = sp.ellipseParams(2);

    R1 = tan(theta1);
    R3 = (a2^2 / a1^2) * (cos(theta1) - cos(theta2)) / (sin(theta2) - sin(theta1));
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
    ax.FontSize = 8;
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