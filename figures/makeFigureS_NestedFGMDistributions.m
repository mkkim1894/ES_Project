function makeFigureS_NestedFGMDistributions(varargin)
% makeFigureS_NestedFGMDistributions - Visualize phenotypic/fitness effects in nested FGM.
%
% Description:
%   Creates supplementary figure showing true distributions of phenotypic effects
%   and fraction of beneficial mutations for two modules with different dimensionalities
%   (n1 = 10, n2 = 40) in the nested Fisher's Geometric Model.
%
% Inputs:
%   varargin - Optional name-value pairs:
%       'n_values'   - Module dimensionalities [n1, n2] (default: [10, 40])
%       'm'          - Mutation magnitude (default: 0.1)
%       'x_i_fixed'  - Fixed x_i for distribution comparison (default: -1)
%       'outputDir'  - Directory for output figure (default: current)
%
% Outputs:
%   Generates EPS file: FigureS_NestedFGMDistributions.eps
%       Subplot A: Distribution of phenotypic effects (Delta x_i)
%       Subplot B: Fraction of beneficial mutations vs distance from optimum
%
% Example:
%   makeFigureS_NestedFGMDistributions();
%   makeFigureS_NestedFGMDistributions('n_values', [10, 20], 'outputDir', 'figures/');
%
% Reference:
%   Kim, M., Ardell, S. M., & Kryazhimskiy, S. (2025).
%   "Module-Selection Balance in the Evolution of Modular Organisms."
%
% See also: simulateNestedSSWM, simulateNestedCM
%
% Copyright (c) 2025 Minkyu Kim, Cornell University
% Licensed under MIT License

%% Parse inputs
p = inputParser;
addParameter(p, 'n_values', [10, 40], @isnumeric);
addParameter(p, 'm', 0.1, @isnumeric);
addParameter(p, 'x_i_fixed', -1, @isnumeric);
addParameter(p, 'outputDir', '.', @ischar);
parse(p, varargin{:});

n_values = p.Results.n_values;
m = p.Results.m;
x_i_fixed = p.Results.x_i_fixed;
outputDir = p.Results.outputDir;

if ~isfolder(outputDir)
    mkdir(outputDir);
end

%% Compute distributions

% Range of x_i values for fitness effects
x_i_values = linspace(-5, -0.01, 100);
distance = abs(x_i_values);

% Range for phenotypic effect distribution
delta_x_range = linspace(-0.5, 0.5, 200);

% Preallocate arrays
pdf_delta_x = zeros(length(n_values), length(delta_x_range));
fracBeneficial = zeros(length(n_values), length(x_i_values));

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

%% Create figure
figure('Position', [100, 100, 800, 300]);

%% Subplot A: Distribution of phenotypic effects
subplot(1, 2, 1);
hold on;

plot(delta_x_range, pdf_delta_x(1, :), 'k-', 'LineWidth', 2, ...
    'DisplayName', sprintf('Module 1 (n = %d)', n_values(1)));
plot(delta_x_range, pdf_delta_x(2, :), 'r-', 'LineWidth', 2, ...
    'DisplayName', sprintf('Module 2 (n = %d)', n_values(2)));

xlabel('Phenotypic Effect ($\Delta x_i$)', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Probability Density', 'FontSize', 12);
legend('Location', 'best', 'FontSize', 10);
grid off;
set(gca, 'FontSize', 10);

% Add label
text(-0.125, 0.95, 'A', 'FontSize', 14, 'FontWeight', 'bold', ...
    'Units', 'normalized', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
hold off;

%% Subplot B: Fraction of beneficial mutations
subplot(1, 2, 2);
hold on;

plot(distance, fracBeneficial(1, :), 'k-', 'LineWidth', 2, ...
    'DisplayName', sprintf('Module 1 (n = %d)', n_values(1)));
plot(distance, fracBeneficial(2, :), 'r-', 'LineWidth', 2, ...
    'DisplayName', sprintf('Module 2 (n = %d)', n_values(2)));

xlabel('Distance from Optimum ($|x_i|$)', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Fraction of Beneficial Mutations', 'FontSize', 12);
legend('Location', 'best', 'FontSize', 10);
grid off;
set(gca, 'FontSize', 10);

% Add label
text(-0.125, 0.95, 'B', 'FontSize', 14, 'FontWeight', 'bold', ...
    'Units', 'normalized', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
hold off;

%% Save figure
set(gcf, 'Color', 'w');
outputFile = fullfile(outputDir, 'FigureS_NestedFGMDistributions.eps');
print(outputFile, '-depsc', '-r300');
fprintf('Figure saved to %s\n', outputFile);

end
