function makeFigureS_ThresholdDTrajectories(simDataFile, theoryDataFile, varargin)
% makeFigureS_ThresholdDTrajectories - Visualize evolutionary trajectories for different threshold D values.
%
% Description:
%   Creates supplementary figure showing how evolutionary trajectories change
%   with different rate-ratio threshold D values. Compares theoretical predictions
%   (blue) against mean simulation trajectories (black).
%
% Inputs:
%   simDataFile    - Path to simulation data (e.g., 'DynSim_2.mat')
%   theoryDataFile - Path to theory data (e.g., 'Theory_D_check.mat')
%   varargin       - Optional name-value pairs:
%       'D_values'  - Threshold values to display (default: [10000, 1000, 100, 10])
%       'outputDir' - Directory for output figure (default: current)
%
% Outputs:
%   Generates EPS file: FigureS_ThresholdDTrajectories.eps
%       4-panel figure showing trajectories for D = 10000, 1000, 100, 10
%
% Example:
%   makeFigureS_ThresholdDTrajectories('DynSim_2.mat', 'Theory_D_check.mat');
%
% Reference:
%   Kim, M., Ardell, S. M., & Kryazhimskiy, S. (2025).
%   "Module-Selection Balance in the Evolution of Modular Organisms."
%
% See also: predictCM_ThresholdD, Run_ThresholdDAnalysis
%
% Copyright (c) 2025 Minkyu Kim, Cornell University
% Licensed under MIT License

%% Parse inputs
p = inputParser;
addRequired(p, 'simDataFile', @ischar);
addRequired(p, 'theoryDataFile', @ischar);
addParameter(p, 'D_values', [10000, 1000, 100, 10], @isnumeric);
addParameter(p, 'outputDir', '.', @ischar);
parse(p, simDataFile, theoryDataFile, varargin{:});

D_values = p.Results.D_values;
outputDir = p.Results.outputDir;

if ~isfolder(outputDir)
    mkdir(outputDir);
end

%% Load data
fprintf('Loading simulation data from %s...\n', simDataFile);
simData = load(simDataFile);
summary = simData.summary;
MeanTrajectory = simData.MeanTrajectory;

fprintf('Loading theory data from %s...\n', theoryDataFile);
theoryData = load(theoryDataFile);
Theory_D_check = theoryData.Theory_D_check;

%% Figure setup
subplot_labels = {'A', 'B', 'C', 'D'};
markers = {'d', '^', 'v', '>', '<', 'o'};
figureWidth = 17.8;
figureHeight = 15;

figure('Units', 'centimeters', 'Position', [1, 1, figureWidth, figureHeight]);

%% Define Gaussian PDF for contours
syms x1 x2
f = exp(-sqrt(summary.SelectionBias(1)*x1^2 + summary.SelectionBias(2)*x2^2)^2 / (2 * summary.sigW^2));
f_x1 = diff(f, x1); f_x2 = diff(f, x2);
f_x1x1 = diff(f_x1, x1); f_x2x2 = diff(f_x2, x2); f_x1x2 = diff(f_x1, x2);
H = [f_x1x1, f_x1x2; f_x1x2, f_x2x2];
det_H_func = matlabFunction(det(H), 'Vars', [x1, x2]);
inflection_level = fminsearch(@(x) abs(det_H_func(x(1), x(2))), [0, 0]);
inflection_pdf_level = double(subs(f, [x1, x2], inflection_level));

%% Subplot positions
subplot_positions = {
    [0.08, 0.55, 0.38, 0.35];
    [0.56, 0.55, 0.38, 0.35];
    [0.08, 0.1, 0.38, 0.35];
    [0.56, 0.1, 0.38, 0.35];
};

customcolor = [0.4, 0.4, 0.4];
CM_color = '#56B4E9';  % Blue for theory

%% Generate subplots
for i = 1:4
    subplot('Position', subplot_positions{i});
    
    % Subplot label
    annotation('textbox', [subplot_positions{i}(1)-0.06, subplot_positions{i}(2)+subplot_positions{i}(4)-0.0, 0.03, 0.03], ...
               'String', subplot_labels{i}, 'FontSize', 20, 'FontWeight', 'bold', 'EdgeColor', 'none');

    % Plot fitness contours
    fc = fcontour(exp(-sqrt(summary.SelectionBias(1)*x1^2 + summary.SelectionBias(2)*x2^2)^2 / ...
         (2 * summary.sigW^2)), 'LineColor', customcolor);
    fc.LevelList = [0.99, inflection_pdf_level];
    fc.LineStyle = '-';
    fc.LineWidth = 0.5;
    hold on;

    % Plot trajectories
    for j = 1:6
        % Theory trajectory (blue)
        Record = Theory_D_check{j, i};
        plot(Record(:,1), Record(:,2), '-', 'Color', CM_color, 'LineWidth', 1.8);
        
        % Mean simulation trajectory (black)
        Mean_TimeStamp = MeanTrajectory.Mean_TimeStamp{j};
        plot(Mean_TimeStamp(1,:), Mean_TimeStamp(2,:), '-', 'Color', 'k', 'LineWidth', 1);

        % Initial point marker
        scatter(Mean_TimeStamp(1,1), Mean_TimeStamp(2,1), 30, 'Marker', markers{j}, ...
            'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
    end

    % Reference lines
    ESL = [-2*summary.SelectionBias(2), -2*summary.SelectionBias(1)]';
    plot([ESL(1), 0], [ESL(2), 0], '-', 'Color', customcolor, 'LineWidth', 1);
    plot([-3, 0], [-3, 0], '--', 'Color', customcolor, 'LineWidth', 1);

    % Axes setup
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    ax.Box = 'off';
    ax.XColor = 'k'; 
    ax.YColor = 'k';
    ax.LineWidth = 1;
    xlim([-3.2, 0.0]);
    ylim([-2.5, 0.0]);
    ax.XTick = [];
    ax.YTick = [];

    % Labels and annotations
    text(-1.9, -2, '\boldmath{$x_1=x_2$}', 'Interpreter', 'latex', 'FontSize', 10, 'Color', 'k');
    text(-2.8, -1.5, '\boldmath{$s_1=s_2$}', 'Interpreter', 'latex', 'FontSize', 10, 'Color', 'k');
    text(-2.4, 0.1, '\boldmath{Module 1 performance $x_1$}', 'Interpreter', 'latex', 'FontSize', 10, 'Color', 'k');
    text(0.15, -0.25, '\boldmath{Module 2 performance $x_2$}', 'Interpreter', 'latex', 'FontSize', 10, 'Color', 'k', 'Rotation', 270);
    text(0.05, -0.05, '\boldmath{0}', 'Interpreter', 'latex', 'FontSize', 10);

    % Threshold label
    text(-1.75, -2.35, sprintf('Threshold $D = %d$', D_values(i)), 'Interpreter', 'latex', 'FontSize', 10, 'Color', 'k');
end

%% Save figure
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 figureWidth figureHeight]);
set(gcf, 'PaperSize', [figureWidth figureHeight]);

outputFile = fullfile(outputDir, 'FigureS_ThresholdDTrajectories.eps');
print(outputFile, '-depsc', '-r300');
fprintf('Figure saved to %s\n', outputFile);

end
