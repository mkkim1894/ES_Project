function makeFigure2(figureType, SSWMRate, CMRate, varargin)
% makeFigure2 — Generate publication-quality Figure 2 for all FGM variants.
%
% Usage:
%   makeFigure2('PleiotropicFGM', 1e-7, 2e-4)
%   makeFigure2('PleiotropicFGM', 1e-7, 2e-4, 'demo')
%   makeFigure2('PleiotropicFGM', 1e-7, 2e-4, 'full')

% ----------------------------- Mode parsing -----------------------------
mode = 'auto';
if nargin >= 4
    tmp = lower(varargin{1});
    if ~ischar(tmp) || ~ismember(tmp, {'demo','full','auto'})
        error('Mode must be ''demo'', ''full'', or ''auto''.');
    end
    mode = tmp;
end

% -------------------------- Resolve results root ------------------------
figureRoot  = fileparts(mfilename('fullpath'));
resultsRoot = resolveResultsRootFlat(figureRoot, mode);

baseDirs = struct('SSWM',       fullfile(resultsRoot,'SSWM'), ...
                  'CM_Asexual', fullfile(resultsRoot,'CM_Asexual'), ...
                  'CM_Sexual',  fullfile(resultsRoot,'CM_Sexual'));

[dirs, column_titles, output_file, regime_colors] = specForFigureType(figureType);

% --------------------------- Locate result files ------------------------
file_names = cell(1,3);
for i = 1:3
    isSSWM  = strcmp(dirs{i}, 'SSWM');
    dirPath = baseDirs.(dirs{i});
    rate    = tern(isSSWM, SSWMRate, CMRate);
    patterns = {
        sprintf('%s_%s_L*_M%.1e*.mat', figureType, dirs{i}, rate), ...
        sprintf('%s_%s_L*_M*.mat', figureType, dirs{i}), ...
        sprintf('%s_%s_M%.1e*.mat', figureType, dirs{i}, rate)
    };
    files = [];
    for p = 1:numel(patterns)
        files = dir(fullfile(dirPath, patterns{p}));
        if ~isempty(files), break; end
    end
    if isempty(files)
        error('No matching file for %s in %s.', figureType, dirPath);
    end
    file_names{i} = fullfile(files(1).folder, files(1).name);
end

% --------------------------- Figure layout ------------------------------
figure('Units','centimeters','Position',[1,1,17.8,12]);
addColumnAnnotations(column_titles);
subplot_positions = defineSubplotPositions(10/12);
subplot_labels = {'A','B','C','D','E','F'};
markers = {'d','^','v','>','<','o'};
gray = [0.4 0.4 0.4];

for i = 1:6
    subplot('Position', subplot_positions{i});
    addSubplotLabel(subplot_labels{i}, subplot_positions{i});
    if i <= 3
        plotFirstThreeSubplots(file_names{i}, markers, regime_colors{i}, gray, figureType);
    else
        plotLastThreeSubplots(file_names{i-3}, markers, gray, figureType);
        if strcmp(figureType,'PleiotropicFGM'), ylim([-4, 2]); end
        if strcmp(figureType,'NestedFGM'),     ylim([-2, 2]); end
    end
end

figDir = fullfile(resultsRoot,'Figures');
if ~isfolder(figDir), mkdir(figDir); end
print(fullfile(figDir, [output_file '.eps']), '-depsc', '-r300');
end

% ======================================================================
% Helper subfunctions
% ======================================================================

function resultsRoot = resolveResultsRootFlat(figureRoot, mode)
    switch mode
        case 'demo'
            candidates = {fullfile(figureRoot,'../results_demo'), fullfile(figureRoot,'../../results_demo')};
        case 'full'
            candidates = {fullfile(figureRoot,'../results'), fullfile(figureRoot,'../../results')};
        otherwise
            candidates = {fullfile(figureRoot,'../results_demo'), fullfile(figureRoot,'../../results_demo'), ...
                          fullfile(figureRoot,'../results'), fullfile(figureRoot,'../../results')};
    end
    for i = 1:numel(candidates)
        if isfolder(candidates{i}), resultsRoot = candidates{i}; return; end
    end
    error('results/ or results_demo/ not found.');
end

function [dirs, titles, outfile, colors] = specForFigureType(figureType)
    dirs = {'SSWM','CM_Asexual','CM_Sexual'};
    switch figureType
        case 'ModularFGM'
            titles = {'Successive mutations','Concurrent, linked modules','Concurrent, unlinked modules'};
            outfile = 'Figure_ModularFGM';
        case 'PleiotropicFGM'
            titles = {'Successive mutations','Concurrent mutations, linked genome','Concurrent mutations, unlinked genome'};
            outfile = 'Figure_PleiotropicFGM';
        case 'NestedFGM'
            titles = {'Successive mutations','Concurrent mutations, linked genome','Concurrent mutations, unlinked genome'};
            outfile = 'Figure_NestedFGM';
        otherwise
            error('Invalid figureType.');
    end
    colors = {'#EDB120','#EDB120','#EDB120'};
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

function plotFirstThreeSubplots(file, markers, color, gray, figureType)
    d = load(file);
    sp = d.simParams; av = getAverageTrajectory(d);
    [~, lvl] = defineGaussianPDF(sp);
    plotPDFContour(sp, lvl, [0.7 0.7 0.7]); hold on;
    if strcmp(figureType,'NestedFGM')
        plotAverageTrajectories(av, markers, color);
    elseif isfield(d,'analyticalTrajectories')
        plotTrajectories(d.analyticalTrajectories, av, markers, color);
    else
        plotAverageTrajectories(av, markers, color);
    end
    plotReferenceLines(sp, gray, figureType);
    customizeAxes(1.2.*[-2.8,0.05], 1.2.*[-2.1875,0.05]);
    text(-3, 0.2, 'Module 1 Performance', 'FontName','Helvetica', 'FontSize',10);
    text(0.2, -0.1, 'Module 2 Performance', 'FontName','Helvetica', 'FontSize',10, 'Rotation',270);
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

function plotReferenceLines(sp, color, figType)
    x1 = linspace(-3,0.5,1000);
    plot(x1, x1, '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 1);
    text(-2.3, -2.4, '$\mathbf{x_1 = x_2}$', 'Interpreter','latex', 'FontSize',10, 'Color',[0.4 0.4 0.4]);
    if any(strcmp(figType,{'ModularFGM','NestedFGM'})) && isfield(sp,'geneticTargetSize')
        if sp.geneticTargetSize(1) ~= sp.geneticTargetSize(2)
            Rb = (sp.ellipseParams(2)^2*sp.geneticTargetSize(2)) / ...
                 (sp.ellipseParams(1)^2*sp.geneticTargetSize(1));
            plot(x1, Rb*x1, '-', 'Color', color, 'LineWidth', 1.2);
        end
    end
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
        R = analytic{j,1}';
        plot(R(1,:), R(2,:), '-', 'Color', color, 'LineWidth', 1.8);
        ts = av.averageTimeStamp{j};
        plot(ts(1,:), ts(2,:), '-', 'Color', '#2776A8', 'LineWidth', 1);
        scatter(ts(1,1), ts(2,1), 30, 'Marker', markers{j}, ...
            'MarkerEdgeColor', '#2776A8', 'MarkerFaceColor', '#2776A8');
    end
end

function plotLastThreeSubplots(file, markers, ~, ~)
    d = load(file); av = getAverageTrajectory(d);
    for j = 1:numel(av.logRatio)
        lr = av.logRatio{j};
        nt = numel(lr);
        plot(1:nt, lr, '-', 'Color', '#2776A8', 'LineWidth', 0.9); hold on;
        plot(1, lr(1), 'Marker', markers{j}, 'Color', '#2776A8', ...
             'MarkerFaceColor', '#2776A8', 'MarkerSize', 5);
    end
    yline(0, '-', 'LineWidth', 1.4, 'Color', [0.4 0.4 0.4]);
    nt = numel(av.logRatio{1});
    xticks(linspace(1, nt, 5)); xticklabels(linspace(0, 1, 5));
    xlim([1, nt]); ylim([-2, 2]);
    set(gca, 'TickLabelInterpreter','latex','FontSize',8);
    xlabel('Normalized Time','FontName','Helvetica','FontSize',10);
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
    ax.XTick = []; ax.YTick = [];
end

function out = tern(c,a,b), out = a; if ~c, out = b; end, end
