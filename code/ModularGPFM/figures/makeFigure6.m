function makeFigure6(figureType, SSWMRate, CMRate, varargin)
% makeFigure6 — Publication-ready Figure 6 (trajectories only; heatmap suppressed)
%
% Usage:
%   makeFigure6('ModularFGM',    1e-7, 2e-4)          % auto-detect results root (prefers full)
%   makeFigure6('PleiotropicFGM',1e-7, 2e-4,'demo')  % force demo
%   makeFigure6('NestedFGM',     1e-7, 2e-4,'full')  % force full
%
% Notes:
% - Loads from <resultsRoot>/Generalization/<modelDir>/<regime>
%   where modelDir is modularFGM | pleiotropicFGM | nestedFGM
% - Saves to <resultsRoot>/Figures/Figure_<Type>_Generalization.eps

% ----------------------------- Parse mode ------------------------------
mode = 'auto';
if nargin >= 4
    tmp = lower(varargin{1});
    if ~ismember(tmp, {'demo','full','auto'})
        error('Invalid mode: use ''demo'', ''full'', or ''auto''.');
    end
    mode = tmp;
end

% ----------------------- Resolve results root --------------------------
figureRoot  = fileparts(mfilename('fullpath'));
resultsRoot = resolveResultsRoot(figureRoot, mode);
fprintf('makeFigure6 loading data from: %s\n', resultsRoot);

% ------------------- Regime dirs and figure metadata -------------------
subplot_labels = {'A','B','C','D','E','F'};
markers         = {'d','^','v','>','<','o'};

modelDir = [lower(figureType(1)) figureType(2:end)]; % e.g. 'NestedFGM' -> 'nestedFGM'

baseDirs = struct( ...
    'SSWM',       fullfile(resultsRoot,'Generalization',modelDir,'SSWM'), ...
    'CM_Asexual', fullfile(resultsRoot,'Generalization',modelDir,'CM_Asexual'), ...
    'CM_Sexual',  fullfile(resultsRoot,'Generalization',modelDir,'CM_Sexual'));

[dirs, column_titles, ~, regime_colors] = specForFigureType(figureType);

% --------------------------- Locate files ------------------------------
file_names = cell(1,3);
for i = 1:3
    isSSWM  = strcmp(dirs{i}, 'SSWM');
    dirPath = baseDirs.(dirs{i});
    rate    = tern(isSSWM, SSWMRate, CMRate);
    fname   = findFirstFile(dirPath, figureType, dirs{i}, rate);
    if isempty(fname)
        fprintf(2, 'No files found. Looked in: %s\n', dirPath);
        error('No files found in %s for %s', dirPath, figureType);
    end
    file_names{i} = fname;
end

% --------------------------- Figure layout -----------------------------
figureWidth = 17.8; figureHeight = 12;
scaling_factor   = 10/12;
subplot_positions = defineSubplotPositions(scaling_factor);
custom_gray = [0.5, 0.5, 0.5];

figure('Units','centimeters','Position',[1,1,figureWidth,figureHeight]);
addColumnAnnotations(column_titles);

% --------------------------- Panels A–F --------------------------------
for i = 1:6
    subplot('Position', subplot_positions{i});
    addSubplotLabel(subplot_labels{i}, subplot_positions{i});

    if i <= 3
        % Top row: phenotype-space trajectories (no heatmap)
        plotTopPanels(file_names{i}, markers, regime_colors{i}, figureType);
    else
        % Bottom row: log R over normalized time
        plotBottomPanels(file_names{i-3}, markers, custom_gray, figureType);
        if strcmp(figureType,'PleiotropicFGM'), ylim([-4, 2]); end
        if strcmp(figureType,'NestedFGM'),     ylim([-2, 2]); end
    end
end

% ----------------------------- Save figure -----------------------------
figDir = fullfile(resultsRoot, 'Figures');
if ~isfolder(figDir), mkdir(figDir); end
outPath = fullfile(figDir, sprintf('Figure_%s_Generalization.eps', figureType));
saveFigure(outPath, figureWidth, figureHeight);
fprintf('Saved %s\n', outPath);

end % ========================= end main =================================


% ======================================================================
% Helper subfunctions
% ======================================================================

function resultsRoot = resolveResultsRoot(figureRoot, mode)
    switch mode
        case 'demo', resultsRoot = fullfile(figureRoot,'../results_demo');
        case 'full', resultsRoot = fullfile(figureRoot,'../results');
        otherwise
            if isfolder(fullfile(figureRoot,'../results'))
                resultsRoot = fullfile(figureRoot,'../results');
            elseif isfolder(fullfile(figureRoot,'../results_demo'))
                resultsRoot = fullfile(figureRoot,'../results_demo');
            else
                error('Neither results/ nor results_demo/ found relative to %s', figureRoot);
            end
    end
    if ~isfolder(resultsRoot)
        resultsRoot = fullfile(figureRoot,'../../results');
    end
end

function [dirs, column_titles, output_file, regime_colors] = specForFigureType(figureType)
    commonDirs   = {'SSWM','CM_Asexual','CM_Sexual'};
    commonColors = {'#EDB120','#EDB120','#EDB120'}; 
    switch figureType
        case 'ModularFGM'
            dirs = commonDirs;
            column_titles = {'Successive mutations', ...
                             'Concurrent mutations, linked modules', ...
                             'Concurrent mutations, unlinked modules'};
            output_file   = 'Figure_ModularFGM_Generalization';
            regime_colors = commonColors;
        case 'PleiotropicFGM'
            dirs = commonDirs;
            column_titles = {'Successive mutations', ...
                             'Concurrent mutations, linked genome', ...
                             'Concurrent mutations, unlinked genome'};
            output_file   = 'Figure_PleiotropicFGM_Generalization';
            regime_colors = commonColors;
        case 'NestedFGM'
            dirs = commonDirs;
            column_titles = {'Successive mutations', ...
                             'Concurrent mutations, linked genome', ...
                             'Concurrent mutations, unlinked genome'};
            output_file   = 'Figure_NestedFGM_Generalization';
            regime_colors = commonColors;
        otherwise
            error('Invalid figureType.');
    end
end

function fname = findFirstFile(dirPath, figureType, regime, rate)
    % Try multiple patterns to match numeric formatting differences
    patterns = { ...
        sprintf('%s_%s_M%.1e*.mat', figureType, regime, rate), ...
        sprintf('%s_%s_M%.0e*.mat', figureType, regime, rate), ...
        sprintf('%s_%s_M%s*.mat',   figureType, regime, num2str(rate)) ...
    };
    fname = '';
    for k = 1:numel(patterns)
        files = dir(fullfile(dirPath, patterns{k}));
        if ~isempty(files)
            fname = fullfile(files(1).folder, files(1).name);
            return
        end
    end
end

function positions = defineSubplotPositions(s)
    positions = {
        [0.06, 0.56 * s, 0.25, 0.38 * s];
        [0.38, 0.56 * s, 0.25, 0.38 * s];
        [0.70, 0.56 * s, 0.25, 0.38 * s];
        [0.06, 0.08 * s, 0.25, 0.38 * s];
        [0.38, 0.08 * s, 0.25, 0.38 * s];
        [0.70, 0.08 * s, 0.25, 0.38 * s];
    };
end

function addColumnAnnotations(titles)
    x_pos = [0.06, 0.38, 0.70];
    for i = 1:numel(titles)
        annotation('textbox', [x_pos(i), 0.95, 0.25, 0.05], ...
            'String', titles{i}, 'FontSize', 14, 'FontName', 'Helvetica', ...
            'HorizontalAlignment', 'center', 'EdgeColor', 'none');
    end
end

function addSubplotLabel(label, position)
    annotation('textbox', [position(1) - 0.04, position(2) + position(4), 0.03, 0.03], ...
        'String', label, 'FontSize', 12, 'FontWeight', 'bold', 'EdgeColor', 'none');
end

% --------------------------- Panel helpers -----------------------------

function plotTopPanels(file_name, markers, ~, figureType)
    % --- Fixed palette (match empirical blue from makeFigure3) ---
    col_blue = '#2776A8';   % empirical/average trajectories only
    col_cnt  = [0.7, 0.7, 0.7]; % contour gray

    % Load + robust average trajectory extraction
    data      = load(file_name);
    simParams = data.simParams;
    avgTraj   = getAverageTrajectory(data);

    % Draw the Gaussian landscape contours
    [~, inflection_pdf_level] = defineGaussianPDF(simParams);
    plotPDFContour(simParams, inflection_pdf_level, col_cnt);
    hold on;

    % Always plot average trajectories in BLUE
    plotAverageTrajectories(avgTraj, markers, col_blue);

    % Reference lines (no diagonal)
    plotReferenceLines(simParams, figureType);

    % Axes and labels
    customizeAxes(1.2.*[-2.8, 0.05], 1.2.*[-2.1875, 0.05]);
    text(-3, 0.2, 'Module 1 Performance', 'FontName','Helvetica', 'FontSize',10);
    text(0.2, -0.1, 'Module 2 Performance', 'FontName','Helvetica', 'FontSize',10, 'Rotation',270);
    text(0.08, 0.12, '0', 'FontName','Helvetica', 'FontSize',10);
end


function plotBottomPanels(file_name, markers, ~, ~)
    data      = load(file_name);
    simParams = data.simParams;
    avgTraj   = getAverageTrajectory(data);

    % Trajectories (consistent blue accent)
    num_conditions = length(avgTraj.logRatio);
    for j = 1:num_conditions
        lr = avgTraj.logRatio{j};
        nt = length(lr);
        plot(1:nt, lr, '-', 'Color', '#2776A8', 'LineWidth', 0.9); hold on;
        plot(1, lr(1), 'Marker', markers{j}, 'Color', '#2776A8', ...
             'MarkerFaceColor', '#2776A8', 'MarkerSize', 5);
    end

    col_msl  = [0.8, 0.3, 0];
    y_esbl = log( (simParams.ellipseParams(2)^2) / (simParams.ellipseParams(1)^2) );
    yline(y_esbl, '-', 'LineWidth', 2, 'Color', col_msl);

    nt = length(avgTraj.logRatio{1});
    xticks(linspace(1, nt, 5)); xticklabels(linspace(0, 1, 5));
    xlim([1, nt]); ylim([-2, 2]);
    set(gca, 'TickLabelInterpreter','latex','FontSize',8);
    xlabel('Normalized Time', 'FontName','Helvetica','FontSize',10);
    ylabel('$\log(x_2/x_1)$', 'Interpreter', 'latex', 'FontSize', 10);
end

% ---------------------- Model/plot primitives --------------------------

function avgTraj = getAverageTrajectory(data)
    if isfield(data,'averageTrajectory')
        avgTraj = data.averageTrajectory;
    elseif isfield(data,'ave')
        avgTraj = data.ave;
    else
        error('Neither "averageTrajectory" nor "ave" found.');
    end
end

function [f, inflection_pdf_level] = defineGaussianPDF(simParams)
    syms x1 x2
    f = exp(-sqrt((x1./simParams.ellipseParams(1)).^2 + ...
                  (x2./simParams.ellipseParams(2)).^2).^2 / ...
            (2 * simParams.landscapeStdDev^2));
    det_H = det(hessian(f, [x1, x2]));
    det_H_func = matlabFunction(det_H, 'Vars', [x1, x2]);
    inflection_level = fminsearch(@(x) abs(det_H_func(x(1), x(2))), [0, 0]);
    inflection_pdf_level = double(subs(f, [x1, x2], inflection_level));
end

function plotPDFContour(simParams, inflection_pdf_level, lineColor)
    fc = fcontour(@(x1,x2) exp(-sqrt((x1./simParams.ellipseParams(1)).^2 + ...
                                     (x2./simParams.ellipseParams(2)).^2).^2 / ...
                            (2 * simParams.landscapeStdDev^2)), ...
                  'LineColor', lineColor, 'LineStyle','-', 'LineWidth',1.4);
    fc.LevelList = [0.99, inflection_pdf_level];
end

function plotReferenceLines(simParams, ~)
    x1 = linspace(-3, 0.5, 1000);
    col_msl  = [0.8, 0.3, 0];    
    x2_esbl = (simParams.ellipseParams(2)^2)/(simParams.ellipseParams(1)^2) * x1;
    plot(x1, x2_esbl, '-', 'Color', col_msl, 'LineWidth', 2);
    text(-3.3, -1.8, '$\mathbf{s_1 = s_2}$', 'Interpreter','latex', ...
         'FontSize',10, 'Color', col_msl);
end

function plotAverageTrajectories(averageTrajectory, markers, color)
    for j = 1:length(averageTrajectory.averageTimeStamp)
        ts = averageTrajectory.averageTimeStamp{j};
        plot(ts(1,:), ts(2,:), '-', 'Color', color, 'LineWidth', 1);
        scatter(ts(1,1), ts(2,1), 30, 'Marker', markers{j}, ...
            'MarkerEdgeColor', color, 'MarkerFaceColor', color);
    end
end

% ----------------------------- Aesthetics ------------------------------

function saveFigure(file_name, width, height)
    set(gcf,'PaperUnits','centimeters');
    set(gcf,'PaperPosition',[0,0,width,height]);
    set(gcf,'PaperSize',[width,height]);
    print(file_name,'-depsc','-r300');
end

function customizeAxes(x_limits, y_limits)
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    ax.Box = 'off';
    ax.XColor = 'k';
    ax.YColor = 'k';
    ax.LineWidth = 1;
    xlim(x_limits);
    ylim(y_limits);
    ax.XTick = [];
    ax.YTick = [];
end

function out = tern(cond, a, b)
    if cond, out = a; else, out = b; end
end
