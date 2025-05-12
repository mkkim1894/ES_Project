clearvars
clc

%% Figure 2
subplot_labels = {'A', 'B', 'C', 'D'};
markers = {'d', '^', 'v', '>', '<', 'o'};
figureWidth = 17.8;
figureHeight = 15;

figure('Units', 'centimeters', 'Position', [1, 1, figureWidth, figureHeight]);

fname = 'DynSim_2.mat';
load(fname)

% Define the PDF
syms x1 x2
f = exp(-sqrt(summary.SelectionBias(1)*x1^2 + summary.SelectionBias(2)*x2^2)^2 / (2 * summary.sigW^2));
f_x1 = diff(f, x1); f_x2 = diff(f, x2);
f_x1x1 = diff(f_x1, x1); f_x2x2 = diff(f_x2, x2); f_x1x2 = diff(f_x1, x2);
H = [f_x1x1, f_x1x2; f_x1x2, f_x2x2];
det_H_func = matlabFunction(det(H), 'Vars', [x1, x2]);
inflection_level = fminsearch(@(x) abs(det_H_func(x(1), x(2))), [0, 0]);
inflection_pdf_level = double(subs(f, [x1, x2], inflection_level));

subplot_positions = {
    [0.08, 0.55, 0.38, 0.35];
    [0.56, 0.55, 0.38, 0.35];
    [0.08, 0.1, 0.38, 0.35];
    [0.56, 0.1, 0.38, 0.35];
};

customcolor = [0.4, 0.4, 0.4];
CM_color = '#56B4E9';  % Blue for CM

for i = 1:4
    subplot('Position', subplot_positions{i});
    
    annotation('textbox', [subplot_positions{i}(1)-0.06, subplot_positions{i}(2)+subplot_positions{i}(4)-0.0, 0.03, 0.03], ...
               'String', subplot_labels{i}, 'FontSize', 20, 'FontWeight', 'bold', 'EdgeColor', 'none');

    load(fname)
    CM_Theory = load('Theory_D_check.mat', 'Theory_D_check');

    fc = fcontour(exp(-sqrt(summary.SelectionBias(1)*x1^2 + summary.SelectionBias(2)*x2^2)^2 / ...
         (2 * summary.sigW^2)), 'LineColor', customcolor);
    fc.LevelList = [0.99, inflection_pdf_level];
    fc.LineStyle = '-';
    fc.LineWidth = 0.5;

    for j = 1:6
        Record = CM_Theory.Theory_D_check{j, i};
        hold on
        plot(Record(:,1), Record(:,2), '-', 'Color', CM_color, 'LineWidth', 1.8);  % Blue theory
        Mean_TimeStamp = MeanTrajectory.Mean_TimeStamp{j};
        plot(Mean_TimeStamp(1,:), Mean_TimeStamp(2,:), '-', 'Color', 'k', 'LineWidth', 1);  % Black mean

        scatter(Mean_TimeStamp(1,1), Mean_TimeStamp(2,1), 30, 'Marker', markers{j}, ...
            'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
    end

    % Reference lines (gray)
    ESL = [-2*summary.SelectionBias(2), -2*summary.SelectionBias(1)]';
    plot([ESL(1),0], [ESL(2),0], '-', 'Color', customcolor, 'LineWidth', 1);
    plot([-3, 0], [-3, 0], '--', 'Color', customcolor, 'LineWidth', 1);

    % Axes setup
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    ax.Box = 'off';
    ax.XColor = 'k'; ax.YColor = 'k';
    ax.LineWidth = 1;

    xlim([-3.2, 0.0])
    ylim([-2.5, 0.0])

    % Remove all tickmarks
    ax.XTick = [];
    ax.YTick = [];

    % Labels and annotations
    text(-1.9, -2, '\boldmath{$x_1=x_2$}', 'Interpreter', 'latex', 'FontSize', 10, 'Color', 'k');
    text(-2.8, -1.5, '\boldmath{$s_1=s_2$}', 'Interpreter', 'latex', 'FontSize', 10, 'Color', 'k');
    text(-2.4, 0.1, '\boldmath{Module 1 performance $x_1$}', 'Interpreter', 'latex', 'FontSize', 10, 'Color', 'k');
    text(0.15, -0.25, '\boldmath{Module 2 performance $x_2$}', 'Interpreter', 'latex', 'FontSize', 10, 'Color', 'k', 'Rotation', 270);
    text(0.05, -0.05, '\boldmath{0}', 'Interpreter', 'latex', 'FontSize', 10);

    thresholds = {'$D = 10000$', '$D = 1000$', '$D = 100$', '$D = 10$'};
    text(-1.75, -2.35, ['Threshold ', thresholds{i}], 'Interpreter', 'latex', 'FontSize', 10, 'Color', 'k');
end

% Save as EPS
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 figureWidth figureHeight]);
set(gcf, 'PaperSize', [figureWidth figureHeight]);

print('supp_fig_Threshold_D_check', '-depsc', '-r300');  % Publication-ready EPS

