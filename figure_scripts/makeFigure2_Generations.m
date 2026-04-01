function makeFigure2_Generations(figureType, simParamsRef, varargin)
% makeFigure2_Generations — Generate Figure 2 (Pleiotropic FGM).
%
% Usage:
%   makeFigure2_Generations('PleiotropicFGM', simParams)
%   makeFigure2_Generations('PleiotropicFGM', simParams, 'test')
%   makeFigure2_Generations('PleiotropicFGM', simParams, 'full', 'proximityCutoff', 0.1)
%
% Inputs:
%   figureType   - 'PleiotropicFGM'
%   simParamsRef - simParams struct used to match files when multiple exist in the results directory
%
% Optional:
%   mode              - 'test', 'full', or 'auto' (positional, after simParamsRef)
%   'proximityCutoff' - Exclude data from ratio plots (D-F) when |x_i| < cutoff.
%                       Set to 0 to disable. Default: 0 (no cutoff).
%
% Outputs:
%   Saves Figure_PleiotropicFGM_Generations.pdf to results/Figures/

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

% -------------------------- Resolve results root ------------------------
figureRoot  = fileparts(mfilename('fullpath'));
resultsRoot = resolveResultsRoot(figureRoot, mode);

baseDirs = struct('SSWM',       fullfile(resultsRoot,'SSWM'), ...
                  'CM_Asexual', fullfile(resultsRoot,'CM_Asexual'), ...
                  'CM_Sexual',  fullfile(resultsRoot,'CM_Sexual'));

dirs          = {'SSWM', 'CM_Asexual', 'CM_Sexual'};
column_titles = {'Successive mutations', ...
                 'Concurrent mutations, linked genome', ...
                 'Concurrent mutations, unlinked genome'};
output_file   = 'Figure_PleiotropicFGM';

% --------------------------- Locate result files ------------------------
file_names = cell(1,3);
for i = 1:3
    dirPath = baseDirs.(dirs{i});
    % Use wildcard matching on all .mat files in the directory
    files = dir(fullfile(dirPath, sprintf('%s_%s_*.mat', figureType, dirs{i})));
    if isempty(files)
        error('No matching file for %s in %s.', figureType, dirPath);
    end
    if isscalar(files)
        file_names{i} = fullfile(files(1).folder, files(1).name);
    else
        % Multiple files: pick the one matching current parameters
        matched = false;
        paramTag = sprintf('N%.0e', simParamsRef.popSize);
        for f = 1:length(files)
            if contains(files(f).name, paramTag)
                file_names{i} = fullfile(files(f).folder, files(f).name);
                matched = true;
                break;
            end
        end
        if ~matched
            % Fallback: use most recent file
            [~, newest] = max([files.datenum]);
            file_names{i} = fullfile(files(newest).folder, files(newest).name);
            warning('Multiple files found in %s. Using most recent: %s', dirPath, files(newest).name);
        end
    end
    fprintf('  Loading: %s\n', file_names{i});
end

% --------------------------- Figure layout ------------------------------
figure('Units','centimeters','Position',[1,1,17.8,12]);
addColumnAnnotations(column_titles);
subplot_positions = defineSubplotPositions(10/12);
subplot_labels = {'A','B','C','D','E','F'};
markers = {'d','^','v','>','<','o'};

for i = 1:6
    subplot('Position', subplot_positions{i});
    addSubplotLabel(subplot_labels{i}, subplot_positions{i});
    if i <= 3
        plotFirstThreeSubplots(file_names{i}, markers);
    else
        plotLastThreeSubplots_Generations(file_names{i-3}, markers, proximityCutoff);
    end
end

figDir = fullfile(resultsRoot,'Figures');
if ~isfolder(figDir), mkdir(figDir); end
print(fullfile(figDir, [output_file '_Generations.pdf']), '-dpdf', '-vector');
fprintf('Saved %s\n', fullfile(figDir, [output_file '_Generations.pdf']));
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
        annotation('textbox', [x(i), 0.95, 0.25, 0.05], 'String', titles{i}, ...
            'FontSize', 14, 'FontName', 'Helvetica', 'HorizontalAlignment', 'center', 'EdgeColor', 'none');
    end
end

function addSubplotLabel(label, pos)
    annotation('textbox', [pos(1)-0.04, pos(2)+pos(4), 0.03, 0.03], ...
        'String', label, 'FontSize', 12, 'FontWeight', 'bold', 'EdgeColor', 'none');
end

function plotFirstThreeSubplots(file, markers)
    color = '#EDB120';
    d = load(file);
    sp = d.simParams;
    av = getAverageTrajectory(d);
    [~, lvl] = defineGaussianPDF(sp);
    plotPDFContour(sp, lvl, [0.7 0.7 0.7]); hold on;
    if isfield(d, 'analyticalTrajectories')
        plotTrajectories(d.analyticalTrajectories, av, markers, color);
    else
        plotAverageTrajectories(av, markers, color);
    end
    plotReferenceLines();
    customizeAxes(1.2.*[-2.8,0.05], 1.2.*[-2.1875,0.05]);
    text(-3, 0.4, 'Module 1 Performance', 'FontName','Helvetica', 'FontSize',10);
    text(0.4, -0.1, 'Module 2 Performance', 'FontName','Helvetica', 'FontSize',10, 'Rotation',270);
    text(0.08, 0.12, '0', 'FontName','Helvetica', 'FontSize',10);
end

function av = getAverageTrajectory(d)
    if isfield(d,'averageTrajectory'), av = d.averageTrajectory;
    elseif isfield(d,'ave'), av = d.ave;
    else, error('No trajectory data found.');
    end
end

function plotPDFContour(sp, lvl, lineColor)
    fc = fcontour(@(x1,x2) exp(-sqrt((x1./sp.ellipseParams(1)).^2 + (x2./sp.ellipseParams(2)).^2).^2 / ...
        (2 * sp.landscapeStdDev^2)), 'LineColor', lineColor, 'LineWidth', 1.4);
    fc.LevelList = [0.99, lvl];
end

function plotReferenceLines()
    x1 = linspace(-3, 0.5, 1000);
    plot(x1, x1, '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 1);
    text(-2.3, -2.4, '$\mathbf{x_1 = x_2}$', 'Interpreter','latex', 'FontSize',10, 'Color',[0.4 0.4 0.4]);
end

function plotAverageTrajectories(av, markers, color)
    for j = 1:numel(av.averageTimeStamp)
        ts = av.averageTimeStamp{j};
        plot(ts(1,:), ts(2,:), '-', 'Color', color, 'LineWidth', 1);
        scatter(ts(1,1), ts(2,1), 30, markers{j}, 'MarkerEdgeColor', color, 'MarkerFaceColor', color);
    end
end

function plotTrajectories(analytic, av, markers, color)
    for j = 1:numel(analytic)
        if ~isempty(analytic{j,1}) && size(analytic{j,1}, 1) > 0
            R = analytic{j,1}';
            if size(R, 1) >= 2 && size(R, 2) > 0
                plot(R(1,:), R(2,:), '-', 'Color', color, 'LineWidth', 2);
            end
        end
        ts = av.averageTimeStamp{j};
        plot(ts(1,:), ts(2,:), '-', 'Color', '#2776A8', 'LineWidth', 1);
        scatter(ts(1,1), ts(2,1), 30, 'Marker', markers{j}, ...
            'MarkerEdgeColor', '#2776A8', 'MarkerFaceColor', '#2776A8');
    end
end

function plotLastThreeSubplots_Generations(file, markers, proximityCutoff)
    d = load(file);
    sp = d.simParams;
    resultTable = getResultTable(d);

    numTimeStamp = 200;
    maxGenAll = 0;
    numSims = size(resultTable, 2);
    minSimsCutoff = min(40, max(2, floor(numSims * 0.5)));

    for j = 1:length(sp.initialAngles)
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
            if isempty(traj), continue; end

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
        stdLogRatio = std(logRatioAll, 0, 1, 'omitnan');

        meanLogRatio(~validTimeIdx) = NaN;
        stdLogRatio(~validTimeIdx) = NaN;

        validGens = genStamps(validTimeIdx);
        if ~isempty(validGens)
            maxGenAll = max(maxGenAll, max(validGens));
        end

        validIdx = ~isnan(meanLogRatio);
        if any(validIdx)
            xFill = [genStamps(validIdx), fliplr(genStamps(validIdx))];
            yFill = [meanLogRatio(validIdx) + stdLogRatio(validIdx), ...
                     fliplr(meanLogRatio(validIdx) - stdLogRatio(validIdx))];
            fill(xFill, yFill, [0.15, 0.46, 0.66], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
        end

        plot(genStamps(validIdx), meanLogRatio(validIdx), '-', 'Color', '#2776A8', 'LineWidth', 0.9); hold on;
        if any(validIdx)
            firstValid = find(validIdx, 1, 'first');
            plot(genStamps(firstValid), meanLogRatio(firstValid), 'Marker', markers{j}, 'Color', '#2776A8', ...
                 'MarkerFaceColor', '#2776A8', 'MarkerSize', 5);
        end
    end

    if maxGenAll <= 1
        maxGenAll = 100;
    end

    yline(0, '-', 'LineWidth', 1.4, 'Color', [0.4 0.4 0.4]);
    xlim([1, maxGenAll]);
    ylim([-3, 2]);
    set(gca, 'TickLabelInterpreter','latex','FontSize',8);
    xlabel('Generations','FontName','Helvetica','FontSize',10);
    ylabel('$\log(x_2/x_1)$', 'Interpreter', 'latex', 'FontSize', 10);
end

function [f, lvl] = defineGaussianPDF(sp)
    syms x1 x2
    f = exp(-sqrt((x1./sp.ellipseParams(1)).^2 + (x2./sp.ellipseParams(2)).^2).^2 / ...
        (2 * sp.landscapeStdDev^2));
    dH = det(hessian(f, [x1, x2]));
    dHf = matlabFunction(dH, 'Vars', [x1, x2]);
    lvl = double(subs(f, [x1,x2], fminsearch(@(x) abs(dHf(x(1),x(2))), [0,0])));
end

function customizeAxes(xl, yl)
    ax = gca;
    ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin';
    ax.Box = 'off'; ax.XColor = 'k'; ax.YColor = 'k'; ax.LineWidth = 1;
    xlim(xl); ylim(yl);
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
