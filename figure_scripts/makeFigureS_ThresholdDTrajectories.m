function makeFigureS_ThresholdDTrajectories(simDataFile, theoryDataFile, varargin)
% makeFigureS_ThresholdDTrajectories - Visualize threshold D sensitivity analysis.
%
% Description:
%   Creates supplementary figure (Figure S4) showing simulated mean trajectories
%   compared to theoretical predictions with different rate-ratio threshold D values.
%   The same simulation data (CM linked modules / clonal interference regime) is 
%   shown in all 4 panels, with different theoretical predictions (D = 10000, 1000, 
%   100, 10) overlaid. 
%
% Inputs:
%   simDataFile    - Path to CM Asexual (linked modules) simulation data from Run_modularFGM
%                    (e.g., 'ModularFGM_CM_Asexual_M2.0e-04.mat')
%   theoryDataFile - Path to theory data from Run_ThresholdDAnalysis
%                    (e.g., 'ThresholdD_Analysis.mat')
%   varargin       - Optional name-value pairs:
%       'outputDir' - Directory for output figure (default: current)
%
% Outputs:
%   Generates EPS file: FigureS_ThresholdDTrajectories.eps
%       2x2 panel figure showing:
%       - Same simulation trajectories (blue) in all panels
%       - Different theoretical predictions (yellow) for D = 10000, 1000, 100, 10
%
% Example:
%   makeFigureS_ThresholdDTrajectories('CM_Asexual/ModularFGM_CM_Asexual_M2.0e-04.mat', ...
%                                       'ThresholdD_Analysis.mat');
%
% Reference:
%   Kim, M., Ardell, S. M., & Kryazhimskiy, S. (2025).
%   "Module-Selection Balance in the Evolution of Modular Organisms."
%
% See also: Run_ThresholdDAnalysis, Run_modularFGM, predictModularCM
%
% Copyright (c) 2025 Minkyu Kim, Cornell University
% Licensed under MIT License

%% Parse inputs
p = inputParser;
addRequired(p, 'simDataFile', @ischar);
addRequired(p, 'theoryDataFile', @ischar);
addParameter(p, 'outputDir', '.', @ischar);
parse(p, simDataFile, theoryDataFile, varargin{:});

outputDir = p.Results.outputDir;

if ~isfolder(outputDir)
    mkdir(outputDir);
end

%% Load simulation data (CM Asexual, linked modules)
fprintf('Loading simulation data from %s...\n', simDataFile);
simData = load(simDataFile);

% Get average trajectory from simulation
if isfield(simData, 'averageTrajectory')
    avgTraj = simData.averageTrajectory;
elseif isfield(simData, 'ave')
    avgTraj = simData.ave;
else
    error('Could not find averageTrajectory in simulation data');
end

% Get simulation parameters for plotting contours
simParams = simData.simParams;

%% Load theory data
fprintf('Loading theory data from %s...\n', theoryDataFile);
theoryData = load(theoryDataFile);

% Extract analytical trajectories - handle both variable names
if isfield(theoryData, 'Theory_D_check')
    analyticalTrajectories = theoryData.Theory_D_check;
elseif isfield(theoryData, 'analyticalTrajectories')
    analyticalTrajectories = theoryData.analyticalTrajectories;
else
    error('Could not find trajectory data in theory file');
end

if isfield(theoryData, 'D_values')
    D_values = theoryData.D_values;
else
    D_values = [10000, 1000, 100, 10];
end

%% Figure setup
subplot_labels = {'A', 'B', 'C', 'D'};
markers = {'d', '^', 'v', '>', '<', 'o'};
figureWidth = 17.8;
figureHeight = 15;

% Colors matching Figure 3
theory_color = '#EDB120';  % Yellow for theory
sim_color = '#2776A8';     % Blue for simulation
customcolor = [0.4, 0.4, 0.4];

% Determine number of panels based on available D values
numPanels = min(4, length(D_values));
if numPanels < 4
    fprintf('Note: Only %d D values available, generating %d panels\n', length(D_values), numPanels);
end

figure('Units', 'centimeters', 'Position', [1, 1, figureWidth, figureHeight]);

%% Compute fitness contour levels
[~, inflection_pdf_level] = defineGaussianPDF(simParams);

%% Subplot positions (2x2 grid) — leave room at top for D titles
subplot_positions = {
    [0.08, 0.54, 0.38, 0.38];   % A: top-left  (D=10000)
    [0.56, 0.54, 0.38, 0.38];   % B: top-right (D=1000)
    [0.08, 0.08, 0.38, 0.38];   % C: bottom-left  (D=100)
    [0.56, 0.08, 0.38, 0.38];   % D: bottom-right (D=10)
};

%% Generate subplots
for i = 1:numPanels
    subplot('Position', subplot_positions{i});
    
    % Subplot label
    annotation('textbox', [subplot_positions{i}(1)-0.06, subplot_positions{i}(2)+subplot_positions{i}(4)-0.02, 0.03, 0.03], ...
               'String', subplot_labels{i}, 'FontSize', 14, 'FontWeight', 'bold', 'EdgeColor', 'none');

    % Plot fitness contours
    plotPDFContour(simParams, inflection_pdf_level, [0.7, 0.7, 0.7]);
    hold on;

    % Plot theoretical trajectories for this D value (yellow)
    for j = 1:length(avgTraj.averageTimeStamp)
        Record = analyticalTrajectories{j, i};
        if ~isempty(Record) && size(Record, 1) > 1
            plot(Record(:,1), Record(:,2), '-', 'Color', theory_color, 'LineWidth', 1.8);
        end
    end
    
    % Plot simulation trajectories (blue) - same in all panels
    for j = 1:length(avgTraj.averageTimeStamp)
        ts = avgTraj.averageTimeStamp{j};
        plot(ts(1,:), ts(2,:), '-', 'Color', sim_color, 'LineWidth', 1);
        scatter(ts(1,1), ts(2,1), 30, 'Marker', markers{j}, ...
            'MarkerEdgeColor', sim_color, 'MarkerFaceColor', sim_color);
    end

    % Reference lines
    plotReferenceLines(simParams, customcolor);

    % Axes setup — consistent with customizeAxes in main figures
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    ax.Box = 'off';
    ax.XColor = 'k';
    ax.YColor = 'k';
    ax.LineWidth = 1;
    xlim(1.2 * [-2.8, 0.05]);
    ylim(1.2 * [-2.1875, 0.05]);
    ax.XTick = [-2, -1];
    ax.YTick = [-2, -1];
    ax.TickLength = [0.015, 0.015];
    ax.TickDir = 'in';
    ax.FontSize = 8;
    ax.TickLabelInterpreter = 'latex';

    % Origin label only
    text(0.08, 0.12, '0', 'FontName', 'Helvetica', 'FontSize', 10);

    % Threshold D as panel title 
    titleX = subplot_positions{i}(1) + subplot_positions{i}(3)/2;
    titleY = subplot_positions{i}(2) + subplot_positions{i}(4) + 0.015;
    annotation('textbox', [titleX - 0.12, titleY, 0.24, 0.04], ...
        'String', sprintf('$D = %d$', D_values(i)), ...
        'Interpreter', 'latex', ...
        'FontSize', 12, ...
        'HorizontalAlignment', 'center', ...
        'EdgeColor', 'none', ...
        'FitBoxToText', 'off');
end

%% Save figure
set(gcf, 'Color', 'w');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 figureWidth figureHeight]);
set(gcf, 'PaperSize', [figureWidth figureHeight]);

outputFile = fullfile(outputDir, 'FigureS_ThresholdDTrajectories.eps');
print(outputFile, '-depsc', '-r300');
fprintf('Figure saved to %s\n', outputFile);

close(gcf);

end

%% Helper functions

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
                  'LineColor', lineColor, 'LineStyle', '-', 'LineWidth', 1.4);
    fc.LevelList = [0.99, inflection_pdf_level];
end

function plotReferenceLines(simParams, customcolor)
    x1 = linspace(-3.5, 0.5, 100);
    
    % x1 = x2 diagonal (gray)
    plot(x1, x1, '-', 'Color', customcolor, 'LineWidth', 1);
    text(-2.3, -2.4, '$\mathbf{x_1 = x_2}$', 'Interpreter', 'latex', 'FontSize', 10, 'Color', customcolor);
    
    % s1 = s2 line (orange) - module-selection balance
    if isfield(simParams, 'geneticTargetSize')
        R_bar = (simParams.ellipseParams(2)^2 * simParams.geneticTargetSize(2)) / ...
                (simParams.ellipseParams(1)^2 * simParams.geneticTargetSize(1));
        x2_MSL = R_bar * x1;
        plot(x1, x2_MSL, '-', 'Color', [0.8, 0.3, 0], 'LineWidth', 2);
        text(-3.3, -1.8, '$\mathbf{s_1 = s_2}$', 'Interpreter', 'latex', 'FontSize', 10, 'Color', [0.8, 0.3, 0]);
    end
end