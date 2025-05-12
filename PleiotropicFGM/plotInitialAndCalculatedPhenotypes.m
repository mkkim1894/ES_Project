%% Functions
function plotInitialAndCalculatedPhenotypes(simParams, currentPhenotypes)
% plotInitialAndCalculatedPhenotypes - Plots initial and calculated phenotypes
% with customized aesthetics according to publication guidelines.
%
% Inputs:
%   simParams - Structure containing the simulation parameters, including:
%       initialPhenotypes - Matrix of initial phenotype coordinates.
%   currentPhenotypes - Matrix representing the resulting phenotype after optimization.
%
% Example Usage:
%   plotInitialAndCalculatedPhenotypes(simParams, currentPhenotypes);

    %% Figure Setup
    % Set figure size according to publication guidelines (in centimeters)
    figureWidth = 17.8;  % Full page width in cm
    figureHeight = 12;

    % Create the figure
    figure('Units', 'centimeters', 'Position', [1, 1, figureWidth, figureHeight]);
    hold on;

    % Define markers and colors
    markerSize = 20;
    targetColor = [0.8500, 0.3250, 0.0980]; % Red color for target phenotypes
    calculatedColor = [0.0000, 0.4470, 0.7410]; % Blue color for calculated phenotypes
    labelFontSize = 10;
    labelFontName = 'Helvetica';

    % Plot the original initial phenotypes from simParams
    scatter(simParams.initialPhenotypes(:, 1), simParams.initialPhenotypes(:, 2), markerSize, targetColor, 'filled', 'DisplayName', 'Target Initial Phenotypes');

    % Plot the current phenotypes from the algorithm
    scatter(currentPhenotypes(:, 1), currentPhenotypes(:, 2), markerSize, calculatedColor, 'filled', 'DisplayName', 'Calculated Current Phenotypes');

    % Add labels to the points for clarity
    for i = 1:size(simParams.initialPhenotypes, 1)
        % Add label to target phenotype
        text(simParams.initialPhenotypes(i, 1), simParams.initialPhenotypes(i, 2), sprintf('Target %d', i), 'FontSize', labelFontSize, 'FontName', labelFontName, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
        % Add label to calculated phenotype
        text(currentPhenotypes(i, 1), currentPhenotypes(i, 2), sprintf('Calc %d', i), 'FontSize', labelFontSize, 'FontName', labelFontName, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
    end

    % Customize plot appearance
    xlabel('Trait 1 (x1)', 'FontSize', 12, 'FontName', labelFontName);
    ylabel('Trait 2 (x2)', 'FontSize', 12, 'FontName', labelFontName);
    title('Comparison of Target and Calculated Phenotypes', 'FontSize', 14, 'FontName', labelFontName);
    legend('Location', 'best', 'FontSize', 10);
    grid on;

    % Adjust axis properties
    ax = gca;
    ax.LineWidth = 1;
    ax.FontSize = 12;
    ax.FontName = labelFontName;
    ax.XColor = 'k';
    ax.YColor = 'k';
    ax.Box = 'off';

    % Finalize the figure
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 figureWidth figureHeight]);
    set(gcf, 'PaperSize', [figureWidth figureHeight]);

    hold off;
end
