% MATLAB script to visualize true phenotypic and fitness effects in nested FGM
% Supplementary figure for two modules with n_1 = 10 and n_2 = 40

clear; close all;

% Parameters
m = 0.1;           % Mutation magnitude (deltaTrait)
n_values = [10, 40]; % Dimensionalities for two modules
x_i_fixed = -1;     % Fixed x_i value for distribution comparison

% Range of x_i values for fitness effects
x_i_values = linspace(-5, -0.01, 100);
distance = abs(x_i_values);

% Range for phenotypic effect distribution
delta_x_range = linspace(-0.5, 0.5, 200); % Range for Delta x_i

% Preallocate arrays
pdf_delta_x = zeros(length(n_values), length(delta_x_range)); % True PDF of Delta x_i
fracBeneficial = zeros(length(n_values), length(x_i_values)); % True fraction beneficial

% Compute true distributions and fractions
for k = 1:length(n_values)
    n = n_values(k);
    
    % True distribution of phenotypic effects at fixed x_i
    variance_fixed = m * sqrt(2 * abs(x_i_fixed) / n);
    mu = -m^2 / 2;
    pdf_delta_x(k, :) = normpdf(delta_x_range, mu, variance_fixed);
    
    % True fraction of beneficial mutations across x_i range
    for i = 1:length(x_i_values)
        x_i = x_i_values(i);
        variance = m * sqrt(2 * abs(x_i) / n);
        % P(Delta x_i > 0) = 1 - CDF(0) for N(-m^2/2, variance)
        fracBeneficial(k, i) = 1 - normcdf(0, mu, variance);
    end
end

% Create supplementary figure with two subplots
figure('Position', [100, 100, 800, 300]);

% Subplot 1: True distribution of phenotypic effects (Delta x_i)
subplot(1, 2, 1);
hold on;
plot(delta_x_range, pdf_delta_x(1, :), 'k-', 'LineWidth', 2, ...
    'DisplayName', 'Module 1 (n = 10)');
plot(delta_x_range, pdf_delta_x(2, :), 'r-', 'LineWidth', 2, ...
    'DisplayName', 'Module 2 (n = 40)');
xlabel('Phenotypic Effect (\Delta x_i)', 'FontSize', 12);
ylabel('Probability Density', 'FontSize', 12);
legend('Location', 'best', 'FontSize', 10);
grid off;
set(gca, 'FontSize', 10);

% Add label "A" outside upper-left corner in axes coordinates
text(-0.125, 0.95, 'A', 'FontSize', 14, 'FontWeight', 'bold', ...
    'Units', 'normalized', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
hold off;

% Subplot 2: True fraction of beneficial mutations vs. distance
subplot(1, 2, 2);
hold on;
plot(distance, fracBeneficial(1, :), 'k-', 'LineWidth', 2, ...
    'DisplayName', 'Module 1 (n = 10)');
plot(distance, fracBeneficial(2, :), 'r-', 'LineWidth', 2, ...
    'DisplayName', 'Module 2 (n = 40)');
xlabel('Distance from Optimum (|x_i|)', 'FontSize', 12);
ylabel('Fraction of Beneficial Mutations', 'FontSize', 12);
legend('Location', 'best', 'FontSize', 10);
grid off;
set(gca, 'FontSize', 10);

% Add label "B" outside upper-left corner in axes coordinates
text(-0.125, 0.95, 'B', 'FontSize', 14, 'FontWeight', 'bold', ...
    'Units', 'normalized', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
hold off;

% Set figure background
set(gcf, 'Color', 'w');

% Save figure (optional)
print('-depsc', 'supp_fig_nestedFGM_modules_true.eps'); % For EPS format