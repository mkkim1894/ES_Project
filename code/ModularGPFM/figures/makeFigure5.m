function makeFigure5(figureType, SSWMRate, CMRate, varargin)
% makeFigure5 — Generate Figure 5 (Generalization performance across regimes)
%
% Usage:
%   makeFigure5('ModularFGM', 1e-7, 2e-4)             % auto-detect full/demo (prefers full)
%   makeFigure5('PleiotropicFGM', 1e-7, 2e-4, 'demo') % force demo mode
%   makeFigure5('NestedFGM', 1e-7, 2e-4, 'full')      % force full mode
%
% Inputs:
%   figureType : 'ModularFGM', 'PleiotropicFGM', or 'NestedFGM'
%   SSWMRate   : Mutation rate for SSWM regime
%   CMRate     : Mutation rate for CM regimes
%   mode       : optional ('full','demo','auto'); default = 'auto'
%
% Notes:
% - Loads from results/Generalization/... or results_demo/Generalization/...
% - Saves EPS to <resultsRoot>/Figures/Figure_<Type>_Generalization.eps

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
fprintf('makeFigure5 loading data from: %s\n', resultsRoot);

% ------------------- Regime dirs and figure metadata -------------------
baseDirs = struct( ...
    'SSWM',       fullfile(resultsRoot,'Generalization','SSWM'), ...
    'CM_Asexual', fullfile(resultsRoot,'Generalization','CM_Asexual'), ...
    'CM_Sexual',  fullfile(resultsRoot,'Generalization','CM_Sexual'));

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
subplot_labels    = {'A','B','C','D','E','F'};
markers           = {'d','^','v','>','<','o'};
customcolor       = [0.5, 0.5, 0.5];

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
figDir = fullfile(resultsRoot, 'Figures');
if ~isfolder(figDir), mkdir(figDir); end
outPath = fullfile(figDir, [output_file '.eps']);
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
    commonDirs = {'SSWM','CM_Asexual','CM_Sexual'};
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

function addSubplotLabel(label, pos)
    annotation('textbox', [pos(1)-0.04, pos(2)+pos(4), 0.03, 0.03], ...
        'String', label, 'FontSize', 12, 'FontWeight', 'bold', 'EdgeColor', 'none');
end

function avgTraj = getAverageTrajectory(data)
    if isfield(data,'averageTrajectory')
        avgTraj = data.averageTrajectory;
    elseif isfield(data,'ave')
        avgTraj = data.ave;
    else
        error('Neither "averageTrajectory" nor "ave" found.');
    end
end

function plotFirstThreeSubplots(~, file_name, markers, color, customcolor, figureType)
    data      = load(file_name);
    simParams = data.simParams;
    avgTraj   = getAverageTrajectory(data);
    [~, inflection_pdf_level] = defineGaussianPDF(simParams);
    plotPDFContour(simParams, inflection_pdf_level, [0.7, 0.7, 0.7]); hold on;

    if strcmp(figureType,'NestedFGM')
        plotAverageTrajectories(avgTraj, markers, color);
    elseif isfield(data,'analyticalTrajectories')
        plotTrajectories(data.analyticalTrajectories, avgTraj, markers, color);
    else
        plotAverageTrajectories(avgTraj, markers, color);
    end

    plotReferenceLines(simParams, customcolor, figureType);
    customizeAxes(1.2.*[-2.8, 0.05], 1.2.*[-2.1875, 0.05]);
    text(-3, 0.2, 'Module 1 Performance', 'FontName','Helvetica', 'FontSize',10);
    text(0.2, -0.1, 'Module 2 Performance', 'FontName','Helvetica', 'FontSize',10, 'Rotation',270);
    text(0.08, 0.12, '0', 'FontName','Helvetica', 'FontSize',10);
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

    %col_equal = [0.8, 0.3, 0];       % s1 = s2
    col_msl = [0.8, 0.3, 0];  % s1/K1 = s2/K2
    %y_esbl = log((simParams.ellipseParams(2)^2)/(simParams.ellipseParams(1)^2));
    %yline(y_esbl, '-', 'LineWidth', 2, 'Color', col_equal);
    if isfield(simParams,'geneticTargetSize')
        y_msl = log((simParams.geneticTargetSize(2)*simParams.ellipseParams(2)^2) / ...
                    (simParams.geneticTargetSize(1)*simParams.ellipseParams(1)^2));
        yline(y_msl, '-', 'LineWidth', 2, 'Color', col_msl);
    end

    nt = length(avgTraj.logRatio{1});
    xticks(linspace(1, nt, 5)); xticklabels(linspace(0, 1, 5));
    xlim([1, nt]); ylim([-2, 2]);
    set(gca,'TickLabelInterpreter','latex','FontSize',8);
    xlabel('Normalized Time','FontName','Helvetica','FontSize',10);
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

function plotPDFContour(simParams, inflection_pdf_level, lineColor)
    fc = fcontour(@(x1,x2) exp(-sqrt((x1./simParams.ellipseParams(1)).^2 + ...
                                     (x2./simParams.ellipseParams(2)).^2).^2 / ...
                            (2 * simParams.landscapeStdDev^2)), ...
                  'LineColor', lineColor, 'LineStyle','-', 'LineWidth',1.4);
    fc.LevelList = [0.99, inflection_pdf_level];
end

function plotReferenceLines(simParams, ~, figureType)
    x1 = linspace(-3, 0.5, 1000);
    %col_equal = [0.8, 0.3, 0];
    col_msl = [0.8, 0.3, 0];
    %x2_esbl = (simParams.ellipseParams(2)^2)/(simParams.ellipseParams(1)^2) * x1;
    %plot(x1, x2_esbl, '-', 'Color', col_equal, 'LineWidth', 2);
    %text(-3.3, -1.8, '$\mathbf{s_1 = s_2}$', 'Interpreter','latex', ...
         %'FontSize',10, 'Color', col_equal);
    if isfield(simParams,'geneticTargetSize')
        R_bar = (simParams.ellipseParams(2)^2 * simParams.geneticTargetSize(2)) / ...
                (simParams.ellipseParams(1)^2 * simParams.geneticTargetSize(1));
        x2_msl = R_bar * x1;
        plot(x1, x2_msl, '-', 'Color', col_msl, 'LineWidth', 2);
        text(-3.2, -2.4, '$\mathbf{s_1/K_1 = s_2/K_2}$', 'Interpreter','latex', ...
         'FontSize',10, 'Color', col_msl);
    end
end

function plotAverageTrajectories(avgTraj, markers, color)
    for j = 1:length(avgTraj.averageTimeStamp)
        ts = avgTraj.averageTimeStamp{j};
        plot(ts(1,:), ts(2,:), '-', 'Color', color, 'LineWidth', 1);
        scatter(ts(1,1), ts(2,1), 30, 'Marker', markers{j}, ...
            'MarkerEdgeColor', color, 'MarkerFaceColor', color);
    end
end

function plotTrajectories(analyticalTrajectories, avgTraj, markers, color)
    for j = 1:length(analyticalTrajectories)
        Record = analyticalTrajectories{j,1}';
        plot(Record(1,:), Record(2,:), '-', 'Color', color, 'LineWidth', 2);
        ts = avgTraj.averageTimeStamp{j};
        plot(ts(1,:), ts(2,:), '-', 'Color', '#2776A8', 'LineWidth', 1);
        scatter(ts(1,1), ts(2,1), 30, 'Marker', markers{j}, ...
            'MarkerEdgeColor', '#2776A8', 'MarkerFaceColor', '#2776A8');
    end
end

function saveFigure(file_name, width, height)
    set(gcf,'PaperUnits','centimeters');
    set(gcf,'PaperPosition',[0,0,width,height]);
    set(gcf,'PaperSize',[width,height]);
    print(file_name,'-depsc','-r300');
end

function customizeAxes(x_limits,y_limits)
    ax = gca;
    ax.XAxisLocation='origin'; ax.YAxisLocation='origin';
    ax.Box='off'; ax.XColor='k'; ax.YColor='k';
    ax.LineWidth=1; xlim(x_limits); ylim(y_limits);
    ax.XTick=[]; ax.YTick=[];
end

function out = tern(cond,a,b)
    if cond, out = a; else, out = b; end
end
