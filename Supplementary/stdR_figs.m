clearvars;
clc;

% Load data
data1 = load('Trait1.mat');
data2 = load('Trait2.mat');
data3 = load('Two_Traits.mat');

% Set up figure size
figureWidth = 17.8;
figureHeight = 8;
figure('Units', 'centimeters', 'Position', [1, 1, figureWidth, figureHeight]);

%% Subplot A ? Initial Points in Phase Plane 
subplot('Position', [0.08, 0.15, 0.38, 0.75]);

% Parameter Setup
summary.step = 0.1;
summary.sigW = 2;
summary.SelectionBias = [1, 2];
summary.InitialAngle = [atan2(1, 6.4), atan2(1, 3.2), atan2(1, 1.6), atan2(1, 0.8), atan2(1, 0.4), atan2(1, 0.2)];
summary.d_ref_List = {-1*[1; 1], -2*[1; 1]};

syms x1 x2
f = exp(-sqrt(summary.SelectionBias(1)*x1^2 + summary.SelectionBias(2)*x2^2)^2 / (2 * summary.sigW^2));
H = hessian(f, [x1, x2]);
det_H_func = matlabFunction(det(H), 'Vars', [x1, x2]);
inflection_level = fminsearch(@(x) abs(det_H_func(x(1), x(2))), [0, 0]);
inflection_pdf_level = double(subs(f, [x1, x2], inflection_level));

fc = fcontour(f, 'LineColor', [0.4, 0.4, 0.4]);  % gray contour
fc.LevelList = [0.99, inflection_pdf_level];
fc.LineStyle = '-';
fc.LineWidth = 0.5;

for idx = 1:2
    summary.d_ref = summary.d_ref_List{idx};
    summary.W_ref = exp(-(vecnorm(summary.d_ref))^2 / (2 * summary.sigW^2));

    for i_pos = 1:6
        theta = summary.InitialAngle(i_pos);
        d = Find_InitialPoint(theta, summary.SelectionBias, summary.W_ref, summary.sigW);
        d = d(d > 0);
        WT = d .* [-cos(theta); -sin(theta)];
        hold on
        scatter(WT(1), WT(2), 30, 'k', 'filled');
        text(WT(1) - 0.1, WT(2) - 0.1, num2str((idx - 1) * 6 + i_pos), ...
            'FontSize', 8, 'Color', 'k');
    end
end

% Reference lines (gray)
ESL = [-2*summary.SelectionBias(2), -2*summary.SelectionBias(1)]';
plot([ESL(1),0], [ESL(2),0],'-', 'Color', [0.4, 0.4, 0.4], 'LineWidth', 1);
plot([-3, 0], [-3, 0], '--', 'Color', [0.4, 0.4, 0.4], 'LineWidth', 1);

xlim([-3.2, 0.0]);
ylim([-2.5, 0.0]);

% Axes setup
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
ax.Box = 'off';
ax.XColor = 'k';
ax.YColor = 'k';
ax.LineWidth = 1;
ax.XTick = [];
ax.YTick = [];
pbaspect([1 1 1]);

% Labels
text(-1.9, -2, '\boldmath{$x_1=x_2$}', 'Interpreter', 'latex', 'FontSize', 10, 'Color', 'k');
text(-2.8, -1.5, '\boldmath{$s_1=s_2$}', 'Interpreter', 'latex', 'FontSize', 10, 'Color', 'k');
text(-2.6, 0.1, '\boldmath{Module 1 performance $x_1$}', 'Interpreter', 'latex', 'FontSize', 10, 'Color', 'k');
text(0.15, -0.35, '\boldmath{Module 2 performance $x_2$}', 'Interpreter', 'latex', 'FontSize', 10, 'Color', 'k', 'Rotation', 270);
text(0.05, -0.05, '\boldmath{0}', 'Interpreter', 'latex', 'FontSize', 10);

% Subpanel label
text(-0.1, 0.95, 'A', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 14, ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');

%% Subplot B ? Trait Variance Comparison
subplot('Position', [0.58, 0.15, 0.38, 0.75]);
hold on;

for i = 1:12
    m1 = mean(data1.summary.var1(i, :));
    m2 = mean(data3.summary.var1(i, :));
    h1 = plot(m1, m2, 'o', 'MarkerFaceColor', '#0072BD', 'MarkerEdgeColor', '#0072BD', ...
              'MarkerSize', 5, 'LineWidth', 1, 'DisplayName', 'Trait 1');

    m3 = mean(data2.summary.var2(i, :));
    m4 = mean(data3.summary.var2(i, :));
    h2 = plot(m3, m4, 'o', 'MarkerFaceColor', '#D95319', 'MarkerEdgeColor', '#D95319', ...
              'MarkerSize', 5, 'LineWidth', 1, 'DisplayName', 'Trait 2');

    % Diagonal line
    x_limits = xlim;
    y_limits = ylim;
    min_limit = min([x_limits, y_limits]);
    max_limit = max([x_limits, y_limits]);
    plot([min_limit max_limit], [min_limit max_limit], 'k-', 'LineWidth', 0.5);
end

pbaspect([1 1 1]);

% Shrink tick font size
ax = gca;
ax.FontSize = 8;  % Reduce from default (typically 10?12)

xlabel('Single-trait Model Variance', 'Interpreter', 'latex', 'FontSize', 10);
ylabel('Two-trait Model Variance', 'Interpreter', 'latex', 'FontSize', 10);

legend({'Trait 1', 'Trait 2'}, 'Interpreter', 'latex', 'Location', 'southeast', 'Box', 'off');

% Subpanel label
text(-0.1, 0.95, 'B', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 14, ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');

%% Save as EPS
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 figureWidth figureHeight]);
set(gcf, 'PaperSize', [figureWidth figureHeight]);
print('std_state_R_combined', '-depsc', '-r300');

%% Helper function
function d = Find_InitialPoint(theta, ShapeParameter, Default_fitness, sigW)
if nargin < 3
    Default_fitness = 0.000123409804086679;
    sigW = 1;
end
syms r
f = exp(-sqrt(ShapeParameter(1)*(r*cos(theta))^2 + ShapeParameter(2)*(r*sin(theta))^2)^2 / (2*sigW^2));
eqn = Default_fitness == f;
sol = solve(eqn, r);
d = double(sol);
end
