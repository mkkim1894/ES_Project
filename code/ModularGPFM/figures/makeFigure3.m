function makeFigure3(figureType, SSWMRate, CMRate, varargin)
% makeFigure3 — Generate Figure 3 comparing evolutionary regimes
%
% Usage:
%   makeFigure3('ModularFGM', 1e-7, 2e-4)           % auto-detect full/demo (prefers full)
%   makeFigure3('ModularFGM', 1e-7, 2e-4, 'demo')   % force demo mode
%   makeFigure3('ModularFGM', 1e-7, 2e-4, 'full')   % force full mode
%
% Inputs:
%   figureType : 'ModularFGM', 'PleiotropicFGM', or 'NestedFGM'
%   SSWMRate   : Mutation rate for SSWM regime
%   CMRate     : Mutation rate for CM (asexual/sexual) regimes
%   mode       : optional ('demo', 'full', or 'auto'); default = 'auto'
%
% Works whether called from repo root or any subdirectory (main/figures/etc.).
%
% Reference:
%   Kim, M., Ardell, S. M., & Kryazhimskiy, S. (2025).
%   "Module-Selection Balance in the Evolution of Modular Organisms."

% ----------------------------- Parse mode ------------------------------
mode = 'auto';
if nargin >= 4
    tmp = lower(varargin{1});
    if ~ismember(tmp, {'demo','full','auto'}), error('Invalid mode: use ''demo'', ''full'', or ''auto''.'); end
    mode = tmp;
end

% ------------------------- Global fig settings -------------------------
subplot_labels = {'A','B','C','D','E','F'};
markers        = {'d','^','v','>','<','o'};
customcolor    = [0.5, 0.5, 0.5];

% ----------------------- Resolve results root --------------------------
figureRoot  = fileparts(mfilename('fullpath'));
resultsRoot = resolveResultsRoot(figureRoot, mode);
fprintf('makeFigure3 loading data from: %s\n', resultsRoot);

% ------------------- Regime dirs and figure metadata -------------------
baseDirs = struct('SSWM', fullfile(resultsRoot,'SSWM'), ...
                  'CM_Asexual', fullfile(resultsRoot,'CM_Asexual'), ...
                  'CM_Sexual',  fullfile(resultsRoot,'CM_Sexual'));

[dirs, column_titles, output_file, regime_colors] = specForFigureType(figureType);

% --------------------------- Locate files ------------------------------
file_names = cell(1,3);
for i = 1:3
    isSSWM  = strcmp(dirs{i}, 'SSWM');
    dirPath = baseDirs.(dirs{i});
    rate    = tern(isSSWM, SSWMRate, CMRate);
    pattern = sprintf('%s_%s_M%.1e*.mat', figureType, dirs{i}, rate);
    files   = dir(fullfile(dirPath, pattern));
    if isempty(files)
        missingRate = tern(isSSWM, SSWMRate, CMRate);
        error('No files found in %s for %s with mutation rate %.1e', dirPath, figureType, missingRate);
    end
    file_names{i} = fullfile(files(1).folder, files(1).name);
end

% --------------------------- Figure layout -----------------------------
figureWidth = 17.8; figureHeight = 12;
scaling_factor   = 10/12;
subplot_positions = defineSubplotPositions(scaling_factor);

figure('Units','centimeters','Position',[1,1,figureWidth,figureHeight]);
addColumnAnnotations(column_titles);

% --------------------------- Plot subplots -----------------------------
for i = 1:6
    subplot('Position', subplot_positions{i});
    addSubplotLabel(subplot_labels{i}, subplot_positions{i});

    if i <= 3
        plotFirstThreeSubplots(i, file_names{i}, markers, regime_colors{i}, customcolor, figureType);
    else
        plotLastThreeSubplots(i, file_names{i-3}, markers, customcolor, figureType);
        if strcmp(figureType,'PleiotropicFGM'), ylim([-4, 2]); end
        if strcmp(figureType,'NestedFGM'),     ylim([-2, 2]); end
    end
end

% ----------------------------- Save figure -----------------------------
figDir = fullfile(resultsRoot,'Figures');
if ~isfolder(figDir), mkdir(figDir); end
outPath = fullfile(figDir, [output_file '.eps']);
saveFigure(outPath, figureWidth, figureHeight);
fprintf('Saved %s\n', outPath);

end % ========================= end main =================================


% ======================================================================
% Helper subfunctions (behavior-preserving refactors only)
% ======================================================================

function resultsRoot = resolveResultsRoot(figureRoot, mode)
    switch mode
        case 'demo', resultsRoot = fullfile(figureRoot,'../results_demo');
        case 'full', resultsRoot = fullfile(figureRoot,'../results');
        otherwise    % 'auto' — prefer full for publication
            if isfolder(fullfile(figureRoot,'../results'))
                resultsRoot = fullfile(figureRoot,'../results');
            elseif isfolder(fullfile(figureRoot,'../results_demo'))
                resultsRoot = fullfile(figureRoot,'../results_demo');
            else
                error('Neither results/ nor results_demo/ found.');
            end
    end
    if ~isfolder(resultsRoot) % fallback for deeper directory structure
        resultsRoot = fullfile(figureRoot,'../../results');
    end
end

function [dirs, column_titles, output_file, regime_colors] = specForFigureType(figureType)
    % Maps figureType to labels/colors exactly as before.
    commonDirs = {'SSWM','CM_Asexual','CM_Sexual'};
    switch figureType
        case 'ModularFGM'
            dirs = commonDirs;
            column_titles = {'Successive mutations', ...
                             'Concurrent mutations, linked modules', ...
                             'Concurrent mutations, unlinked modules'};
            output_file   = 'Figure_ModularFGM';
            regime_colors = {'#EDB120','#EDB120','#EDB120'};
        case 'PleiotropicFGM'
            dirs = commonDirs;
            column_titles = {'Successive mutations', ...
                             'Concurrent mutations, linked genome', ...
                             'Concurrent mutations, unlinked genome'};
            output_file   = 'Figure_PleiotropicFGM';
            regime_colors = {'#EDB120','#EDB120','#EDB120'};
        case 'NestedFGM'
            dirs = commonDirs;
            column_titles = {'Successive mutations', ...
                             'Concurrent mutations, linked genome', ...
                             'Concurrent mutations, unlinked genome'};
            output_file   = 'Figure_NestedFGM';
            regime_colors = {'#EDB120','#EDB120','#EDB120'};
        otherwise
            error('Invalid figureType. Choose ModularFGM, PleiotropicFGM, or NestedFGM.');
    end
end

function positions = defineSubplotPositions(scaling_factor)
    positions = {
        [0.06, 0.56 * scaling_factor, 0.25, 0.38 * scaling_factor];
        [0.38, 0.56 * scaling_factor, 0.25, 0.38 * scaling_factor];
        [0.70, 0.56 * scaling_factor, 0.25, 0.38 * scaling_factor];
        [0.06, 0.08 * scaling_factor, 0.25, 0.38 * scaling_factor];
        [0.38, 0.08 * scaling_factor, 0.25, 0.38 * scaling_factor];
        [0.70, 0.08 * scaling_factor, 0.25, 0.38 * scaling_factor];
    };
end

function addColumnAnnotations(column_titles)
    x_pos = [0.06, 0.38, 0.70];
    for i = 1:numel(column_titles)
        annotation('textbox', [x_pos(i), 0.95, 0.25, 0.05], ...
            'String', column_titles{i}, 'FontSize', 14, 'FontName', 'Helvetica', ...
            'HorizontalAlignment', 'center', 'EdgeColor', 'none');
    end
end

function addSubplotLabel(label, position)
    annotation('textbox', [position(1) - 0.04, position(2) + position(4), 0.03, 0.03], ...
        'String', label, 'FontSize', 12, 'FontWeight', 'bold', 'EdgeColor', 'none');
end

function plotFirstThreeSubplots(~, file_name, markers, color, customcolor, figureType)
    data       = load(file_name);
    simParams  = data.simParams;
    avgTraj    = getAverageTrajectory(data);

    % Gaussian fitness contours (unchanged behavior)
    [~, inflection_pdf_level] = defineGaussianPDF(simParams);
    plotPDFContour(simParams, inflection_pdf_level, [0.7, 0.7, 0.7]); hold on;

    % Trajectories (unchanged behavior)
    if strcmp(figureType,'NestedFGM')
        plotAverageTrajectories(avgTraj, markers, color);
    else
        if isfield(data,'analyticalTrajectories')
            plotTrajectories(data.analyticalTrajectories, avgTraj, markers, color);
        else
            plotAverageTrajectories(avgTraj, markers, color);
        end
    end

    % Reference geometry (unchanged behavior)
    plotReferenceLines(simParams, customcolor, figureType);

    % Axes and labeling (unchanged behavior)
    customizeAxes(1.2.*[-2.8, 0.05], 1.2.*[-2.1875, 0.05]);
    text(-3, 0.2, 'Module 1 Performance', 'FontName','Helvetica', 'FontSize',10);
    text(0.2, -0.1, 'Module 2 Performance', 'FontName','Helvetica', 'FontSize',10, 'Rotation',270);
    text(0.08, 0.12, '0', 'FontName','Helvetica', 'FontSize',10);
end

function avgTraj = getAverageTrajectory(data)
    if isfield(data,'averageTrajectory')
        avgTraj = data.averageTrajectory;
    elseif isfield(data,'ave')
        avgTraj = data.ave;
    else
        error('Neither "averageTrajectory" nor "ave" found in %s', inputname(1));
    end
end

function plotPDFContour(simParams, inflection_pdf_level, lineColor)
    fc = fcontour(@(x1,x2) exp(-sqrt((x1./simParams.ellipseParams(1)).^2 + ...
                                     (x2./simParams.ellipseParams(2)).^2).^2 / ...
                            (2 * simParams.landscapeStdDev^2)), ...
                  'LineColor', lineColor, 'LineStyle','-', 'LineWidth',1.4);
    fc.LevelList = [0.99, inflection_pdf_level];
end

function plotReferenceLines(simParams, ~, figureType)
    x1 = linspace(-3, 0.5, 1000);
    x2_diag = x1;
    plot(x1, x2_diag, '-', 'Color',[0.4 0.4 0.4], 'LineWidth',1);  % gray diagonal
    text(-2.3, -2.4, '$\mathbf{x_1 = x_2}$', 'Interpreter','latex', 'FontSize',10, 'Color',[0.4 0.4 0.4]);

    if any(strcmp(figureType, {'ModularFGM','NestedFGM'})) && isfield(simParams,'geneticTargetSize')
        if simParams.geneticTargetSize(1) == simParams.geneticTargetSize(2)
            if strcmp(figureType,'NestedFGM')
                x2_ESBL = (simParams.ellipseParams(2)^2) / (simParams.ellipseParams(1)^2) * x1;
                plot(x1, x2_ESBL, '-', 'Color', [0.8, 0.3, 0], 'LineWidth', 2);
                text(-3.3, -1.8, '$\mathbf{s_1 = s_2}$', 'Interpreter','latex', 'FontSize',10, 'Color',[0.8,0.3,0]);
            else
                R_bar  = (simParams.ellipseParams(2)^2 * simParams.geneticTargetSize(2)) / ...
                         (simParams.ellipseParams(1)^2 * simParams.geneticTargetSize(1));
                x2_MSL = R_bar * x1;
                plot(x1, x2_MSL, '-', 'Color', [0.8, 0.3, 0], 'LineWidth', 2);
                text(-3.3, -1.8, '$\mathbf{s_1 = s_2}$', 'Interpreter','latex', 'FontSize',10, 'Color',[0.8,0.3,0]);
            end
        else
            R_bar  = (simParams.ellipseParams(2)^2 * simParams.geneticTargetSize(2)) / ...
                     (simParams.ellipseParams(1)^2 * simParams.geneticTargetSize(1));
            x2_MSL = R_bar * x1;
            plot(x1, x2_MSL, '-',  'Color', [0.8, 0.3, 0], 'LineWidth', 2);

            x2_ESBL = (simParams.ellipseParams(2)^2) / (simParams.ellipseParams(1)^2) * x1;
            plot(x1, x2_ESBL, '--', 'Color', [0.1 0.2 0.8], 'LineWidth', 1.2);

            text(-2.4, -2.5, '$\mathbf{s_1/K_1 = s_2/K_2}$', 'Interpreter','latex', 'FontSize',10, 'Color',[0.1 0.2 0.8]);
            text(-3.3, -1.8, '$\mathbf{s_1 = s_2}$',        'Interpreter','latex', 'FontSize',10, 'Color',[0.8, 0.3, 0]);
        end
    end
end

function plotAverageTrajectories(averageTrajectory, markers, color)
    for j = 1:length(averageTrajectory.averageTimeStamp)
        ts = averageTrajectory.averageTimeStamp{j};
        plot(ts(1,:), ts(2,:), '-', 'Color', color, 'LineWidth', 1);
        scatter(ts(1,1), ts(2,1), 30, 'Marker', markers{j}, ...
            'MarkerEdgeColor', color, 'MarkerFaceColor', color);
    end
end

function plotTrajectories(analyticalTrajectories, averageTrajectory, markers, color)
    for j = 1:length(analyticalTrajectories)
        Record = analyticalTrajectories{j,1}';
        plot(Record(1,:), Record(2,:), '-', 'Color', color, 'LineWidth', 2);
        ts = averageTrajectory.averageTimeStamp{j};
        plot(ts(1,:), ts(2,:), '-', 'Color', '#2776A8', 'LineWidth', 1);
        scatter(ts(1,1), ts(2,1), 30, 'Marker', markers{j}, ...
            'MarkerEdgeColor', '#2776A8', 'MarkerFaceColor', '#2776A8');
    end
end

function plotLastThreeSubplots(~, file_name, markers, ~, figureType)
    data      = load(file_name);
    simParams = data.simParams;
    avgTraj   = getAverageTrajectory(data);

    num_conditions = length(avgTraj.logRatio);
    for j = 1:num_conditions
        lr = avgTraj.logRatio{j};
        nt = length(lr);
        plot(1:nt, lr, '-', 'Color', '#2776A8', 'LineWidth', 0.9); hold on;
        plot(1, lr(1), 'Marker', markers{j}, 'Color', '#2776A8', ...
             'MarkerFaceColor', '#2776A8', 'MarkerSize', 5);
    end

    if any(strcmp(figureType, {'ModularFGM','NestedFGM'}))
        yline(0, '-', 'LineWidth', 1.4, 'Color', [0.4 0.4 0.4]);
        if isfield(simParams,'geneticTargetSize') && simParams.geneticTargetSize(1) ~= simParams.geneticTargetSize(2)
            yline(log( (simParams.geneticTargetSize(2)*simParams.ellipseParams(2)^2) / ...
                        (simParams.geneticTargetSize(1)*simParams.ellipseParams(1)^2) ), ...
                  '-',  'LineWidth', 2, 'Color', [0.8, 0.3, 0]); % MSL
        else
            yline(log(simParams.ellipseParams(2)^2 / simParams.ellipseParams(1)^2), ...
                  '-',  'LineWidth', 2, 'Color', [0.8, 0.3, 0]); % equal-fitness line
        end
    else
        yline(0, '-', 'LineWidth', 1.4, 'Color', [0.4 0.4 0.4]); % Pleiotropic baseline only
    end

    nt = length(avgTraj.logRatio{1});
    xticks(linspace(1, nt, 5)); xticklabels(linspace(0, 1, 5));
    xlim([1, nt]); ylim([-2, 2]);
    set(gca, 'TickLabelInterpreter','latex','FontSize',8);
    xlabel('Normalized Time', 'FontName','Helvetica','FontSize',10);
    ylabel('$\log(x_2/x_1)$', 'Interpreter', 'latex', 'FontSize', 10);
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

function saveFigure(file_name, width, height)
    set(gcf, 'PaperUnits','centimeters');
    set(gcf, 'PaperPosition', [0, 0, width, height]);
    set(gcf, 'PaperSize',     [width, height]);
    print(file_name, '-depsc', '-r300');
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
