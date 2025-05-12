function makeFigure6(figureType, SSWMRate, CMRate)
    % Main function to create and save simulation figures
    % figureType: 'ModularFGM', 'PleiotropicFGM', or 'NestedFGM'
    % SSWMRate: Mutation rate for the SSWM regime
    % CMRate: Mutation rate for the CM (asexual and sexual) regimes

    % Define global figure settings
    subplot_labels = {'A', 'B', 'C', 'D', 'E', 'F'};
    markers = {'d', '^', 'v', '>', '<', 'o'};
    customcolor = [0.4, 0.4, 0.4];

    % Directories for different regimes
    baseDirs = struct(...
        'SSWM', 'results/Generalization/nestedFGM/SSWM', ...
        'CM_Asexual', 'results/Generalization/nestedFGM/CM_Asexual', ...
        'CM_Sexual', 'results/Generalization/nestedFGM/CM_Sexual');

    % Titles and output based on figureType
    switch figureType
        case 'ModularFGM'
            dirs = {'SSWM', 'CM_Asexual', 'CM_Sexual'};
            column_titles = {'Successive mutations', ...
                             'Concurrent mutations, linked modules', ...
                             'Concurrent mutations, unlinked modules'};
            output_file = 'Figure_ModularFGM';
            regime_colors = {'#EDB120', '#56B4E9', '#D55E00'}; % Yellow, Blue, Orange
        case 'PleiotropicFGM'
            dirs = {'SSWM', 'CM_Asexual', 'CM_Sexual'};
            column_titles = {'Successive mutations', ...
                             'Concurrent mutations, linked genome', ...
                             'Concurrent mutations, unlinked genome'};
            output_file = 'Figure_PleiotropicFGM';
            regime_colors = {'#EDB120', '#EDB120', '#EDB120'}; % Yellow for all
        case 'NestedFGM'
            dirs = {'SSWM', 'CM_Asexual', 'CM_Sexual'};
            column_titles = {'Successive mutations', ...
                             'Concurrent mutations, linked genome', ...
                             'Concurrent mutations, unlinked genome'};
            output_file = 'Figure_NestedFGM';
            regime_colors = {'#EDB120', '#56B4E9', '#D55E00'}; % Adjust colors if necessary
        otherwise
            error('Invalid figureType. Choose either "ModularFGM", "PleiotropicFGM", or "NestedFGM".');
    end

    % Locate relevant files for each regime
    file_names = cell(1, 3); % One file per regime
    for i = 1:3
        dirPath = baseDirs.(dirs{i});
        if strcmp(dirs{i}, 'SSWM')
            % SSWM uses its own mutation rate
            filePattern = sprintf('%s_%s_M%.1e*.mat', figureType, dirs{i}, SSWMRate);
        else
            % CM regimes use the CM mutation rate
            filePattern = sprintf('%s_%s_M%.1e*.mat', figureType, dirs{i}, CMRate);
        end
        files = dir(fullfile(dirPath, filePattern));
        if isempty(files)
            error('No files found in %s for %s with mutation rate %.1e', dirPath, figureType, ...
                strcmp(dirs{i}, 'SSWM') * SSWMRate + ~strcmp(dirs{i}, 'SSWM') * CMRate); % Select appropriate rate
        end
        % Pick the first file or sort files if needed
        file_names{i} = fullfile(files(1).folder, files(1).name);
    end

    % Set up figure dimensions
    figureWidth = 17.8; % Full page width in cm
    figureHeight = 12;
    scaling_factor = 10 / 12; % Adjust subplot scaling

    % Subplot positions
    subplot_positions = defineSubplotPositions(scaling_factor);

    % Create figure
    figure('Units', 'centimeters', 'Position', [1, 1, figureWidth, figureHeight]);
    addColumnAnnotations(column_titles);

    % Loop through subplots
    for i = 1:6
        subplot('Position', subplot_positions{i});
        addSubplotLabel(subplot_labels{i}, subplot_positions{i});

        if i <= 3
            plotFirstThreeSubplots(i, file_names{i}, markers, regime_colors{i}, customcolor, figureType);
            hold on;
        else
            plotLastThreeSubplots(i, file_names{i - 3}, markers, customcolor, figureType);
            if strcmp(figureType, 'PleiotropicFGM') 
                ylim([-4, 2]);   
            end
            if strcmp(figureType, 'NestedFGM')               
                ylim([-2, 2]);   
            end
        end
    end

    % Save the figure
    saveFigure('figure6', figureWidth, figureHeight);
end



%% Define Subplot Positions
function positions = defineSubplotPositions(scaling_factor)
    positions = {
        [0.06, 0.56 * scaling_factor, 0.25, 0.38 * scaling_factor];
        [0.38, 0.56 * scaling_factor, 0.25, 0.38 * scaling_factor];
        [0.7,  0.56 * scaling_factor, 0.25, 0.38 * scaling_factor];
        [0.06, 0.08 * scaling_factor, 0.25, 0.38 * scaling_factor];
        [0.38, 0.08 * scaling_factor, 0.25, 0.38 * scaling_factor];
        [0.7,  0.08 * scaling_factor, 0.25, 0.38 * scaling_factor];
    };
end

%% Add Column Annotations
function addColumnAnnotations(column_titles)
    for i = 1:length(column_titles)
        x_pos = [0.06, 0.38, 0.7];
        annotation('textbox', [x_pos(i), 0.95, 0.25, 0.05], 'String', column_titles{i}, ...
                   'FontSize', 14, 'FontName', 'Helvetica', 'HorizontalAlignment', 'center', 'EdgeColor', 'none');
    end
end

%% Add Subplot Label
function addSubplotLabel(label, position)
    annotation('textbox', [position(1) - 0.04, position(2) + position(4) - 0.0, 0.03, 0.03], ...
        'String', label, 'FontSize', 12, 'FontWeight', 'bold', 'EdgeColor', 'none');
end

%% Plot First Three Subplots (A, B, C)
function plotFirstThreeSubplots(index, file_name, markers, color, customcolor, figureType)
    % Load data file
    data = load(file_name);

    % Extract variables
    simParams = data.simParams;
    averageTrajectory = data.averageTrajectory;

    % Define and process Gaussian PDF    
    [f, inflection_pdf_level] = defineGaussianPDF(simParams);     
    plotPDFContour(f, simParams, inflection_pdf_level, [0.5, 0.5, 0.5]);

    hold on;

    % For 'NestedFGM', only plot simulation data (averageTrajectory)
    if strcmp(figureType, 'NestedFGM')
        plotAverageTrajectories(averageTrajectory, markers, 'k');
    else
        % Plot analytical trajectories if available
        if isfield(data, 'analyticalTrajectories')
            plotTrajectories(data.analyticalTrajectories, averageTrajectory, markers, color);
        else
            plotAverageTrajectories(averageTrajectory, markers, 'k');
        end
    end

    % Plot reference lines
    plotReferenceLines(simParams, customcolor, figureType);

    % Customize axes
    customizeAxes(1.2.*[-2.8, 0.05], 1.2.*[-2.1875, 0.05]);

    % Non-math text with Helvetica
    text(-3, 0.2, 'Module 1 Performance', 'FontName', 'Helvetica', 'FontSize', 10, 'Color', 'k');

    % Non-math text with Helvetica
    text(0.2, -0.1, 'Module 2 Performance', 'FontName', 'Helvetica', 'FontSize', 10, 'Color', 'k', 'Rotation', 270);

    % Use Helvetica for the '0' text
    text(0.08, 0.12, '0', 'FontName', 'Helvetica', 'FontSize', 10, 'Color', 'k');
end

%% Function to plot the PDF contour
function plotPDFContour(f, simParams, inflection_pdf_level, lineColor)
    fc = fcontour(@(x1, x2) exp(-sqrt((x1./simParams.ellipseParams(1)).^2 + (x2./simParams.ellipseParams(2)).^2).^2 / ...
        (2 * simParams.landscapeStdDev^2)), 'LineColor', lineColor);
    fc.LevelList = [0.99, inflection_pdf_level];
    fc.LineStyle = '-';
    fc.LineWidth = 0.001;
    hold on;
end

%% Function to plot the Reference lines
function plotReferenceLines(simParams, color, figureType)
    % Generate the diagonal line (always plotted)
    x1 = linspace(1.2.*[-2.8], 0.05, 1000); % Adjust range as needed

    if strcmp(figureType, 'ModularFGM') || strcmp(figureType, 'NestedFGM')
        % For 'ModularFGM' and 'NestedFGM', plot reference lines

        % Calculate R_bar and plot MSL line if geneticTargetSize is available
        if isfield(simParams, 'geneticTargetSize')
            if simParams.geneticTargetSize(1) == simParams.geneticTargetSize(2)
                x2_diagonal = x1;
                plot(x1, x2_diagonal, '--', 'Color', color, 'LineWidth', 0.6); % Diagonal line
                
                if strcmp(figureType, 'NestedFGM')
                    x2_ESBL = (simParams.ellipseParams(2)^2) / (simParams.ellipseParams(1)^2) * x1;
                    plot(x1, x2_ESBL, '-', 'Color', color, 'LineWidth', 0.6); % ESBL line       
                else                   
                    R_bar = (simParams.ellipseParams(2)^2 * simParams.geneticTargetSize(2)) / ...
                        (simParams.ellipseParams(1)^2 * simParams.geneticTargetSize(1));                
                    x2_MSL = R_bar * x1;                
                    plot(x1, x2_MSL, '-', 'Color', color, 'LineWidth', 0.6); % MSL line
                end
                
                text(-2.2, -2.3, '$\boldmath{x_1 = x_2}$', 'Interpreter', 'latex', ...
                     'FontSize', 10, 'Color', 'k');  % LaTeX for math symbols
                text(-3.3, -1.8, '$\boldmath{s_1 = s_2}$', 'Interpreter', 'latex', ...
                     'FontSize', 10, 'Color', 'k');  % LaTeX for math symbols
            else
                R_bar = (simParams.ellipseParams(2)^2 * simParams.geneticTargetSize(2)) / ...
                        (simParams.ellipseParams(1)^2 * simParams.geneticTargetSize(1));
                x2_MSL = R_bar * x1;
                plot(x1, x2_MSL, '-', 'Color', color, 'LineWidth', 0.6); % MSL line

                x2_ESBL = (simParams.ellipseParams(2)^2) / (simParams.ellipseParams(1)^2) * x1;
                plot(x1, x2_ESBL, '--', 'Color', color, 'LineWidth', 0.6); % ESBL line

                text(-2.4, -2.5, '$\boldmath{s_1/K_1 = s_2/K_2}$', 'Interpreter', 'latex', ...
                     'FontSize', 10, 'Color', 'k');  % LaTeX for math symbols
                text(-3.3, -1.8, '$\boldmath{s_1 = s_2}$', 'Interpreter', 'latex', ...
                     'FontSize', 10, 'Color', 'k');  % LaTeX for math symbols
            end
        end
        
    elseif strcmp(figureType, 'PleiotropicFGM')
        % Only plot diagonal for PleiotropicFGM
        x2_diagonal = x1;
        plot(x1, x2_diagonal, '--', 'Color', color, 'LineWidth', 0.6); % Diagonal line

        text(-2.2, -2.3, '$\boldmath{x_1 = x_2}$', 'Interpreter', 'latex', ...
             'FontSize', 10, 'Color', 'k');  % LaTeX for math symbols
    end
end

%% Function to plot only the mean trajectories
function plotAverageTrajectories(averageTrajectory, markers, color)
    for j = 1:length(averageTrajectory.averageTimeStamp)
        averageTimeStamp = averageTrajectory.averageTimeStamp{j};
        plot(averageTimeStamp(1, :), averageTimeStamp(2, :), '-', 'Color', color, 'LineWidth', 1);

        % Plot starting point marker
        scatter(averageTimeStamp(1, 1), averageTimeStamp(2, 1), 30, 'Marker', markers{j}, ...
            'MarkerEdgeColor', color, 'MarkerFaceColor', color);
    end
end

%% Function to plot the theoretical and mean trajectories
function plotTrajectories(analyticalTrajectories, averageTrajectory, markers, color)
    for j = 1:length(analyticalTrajectories)
        % Plot analytical trajectory
        Record = analyticalTrajectories{j, 1}';
        plot(Record(1, :), Record(2, :), '-', 'Color', color, 'LineWidth', 1.8);

        % Plot mean trajectory
        averageTimeStamp = averageTrajectory.averageTimeStamp{j};
        plot(averageTimeStamp(1, :), averageTimeStamp(2, :), 'k-', 'LineWidth', 1);

        % Plot starting point marker
        scatter(averageTimeStamp(1, 1), averageTimeStamp(2, 1), 30, 'Marker', markers{j}, ...
            'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
    end
end

%% Plot Last Three Subplots (D, E, F)
function plotLastThreeSubplots(index, file_name, markers, customcolor, figureType)
    data = load(file_name);
    averageTrajectory = data.averageTrajectory;
    simParams = data.simParams;

    % Plot log ratios over time for each trajectory
    num_conditions = length(averageTrajectory.logRatio); % Number of conditions (6 in this case)
    for j = 1:num_conditions   
        log_ratios = averageTrajectory.logRatio{j}; % Extract log ratios for trajectory j    
        num_timepoints = length(log_ratios); % Number of timepoints in this trajectory

        % Plot log ratios over time
        plot(1:num_timepoints, log_ratios, '-', 'Color', 'k', ...
            'MarkerFaceColor', 'k', 'LineWidth', 1);
        hold on;

        % Plot starting point marker    
        plot(1, log_ratios(1), 'Marker', markers{j}, 'Color', 'k', ...
            'MarkerFaceColor', 'k', 'MarkerSize', 4);            
    end

    if strcmp(figureType, 'ModularFGM')
        if simParams.geneticTargetSize(1) == simParams.geneticTargetSize(2)
            yline(log(simParams.ellipseParams(2)^2 / simParams.ellipseParams(1)^2), '-', 'LineWidth', 1.5, 'Color', customcolor);                    
            yline(0, '--', 'LineWidth', 1.5, 'Color', customcolor); 
        else            
            yline(log(simParams.ellipseParams(2)^2 / simParams.ellipseParams(1)^2), '-', 'LineWidth', 1.5, 'Color', customcolor);      
            yline(log( ( simParams.geneticTargetSize(2)*simParams.ellipseParams(2)^2 ) /...
                ( simParams.geneticTargetSize(1)*simParams.ellipseParams(1)^2 ) ), '-', 'LineWidth', 1.5, 'Color', customcolor);        
        end
    elseif strcmp(figureType, 'PleiotropicFGM') || strcmp(figureType, 'NestedFGM')
        yline(0, '--', 'LineWidth', 1.5, 'Color', customcolor); 
    end

    % Scale x-axis from 0 to 1       
    xticks(linspace(1, num_timepoints, 5)); % Adjust the number of ticks as needed        
    xticklabels(linspace(0, 1, 5)); % Set the tick labels to normalized scale (0 to 1)

    % Set axis limits and other properties        
    xlim([1, num_timepoints]);        
    ylim([-2, 2]);        
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 8); 

    % Customize axes
    xlabel('Normalized Time', 'FontName', 'Helvetica', 'FontSize', 10);
    ylabel('log R', 'FontName', 'Helvetica', 'FontSize', 10);
end

%% Define Gaussian PDF
function [f, inflection_pdf_level] = defineGaussianPDF(simParams)
    syms x1 x2
    f = exp(-sqrt((x1./simParams.ellipseParams(1)).^2 + (x2./simParams.ellipseParams(2)).^2).^2 / ...
        (2 * simParams.landscapeStdDev^2));
    det_H = det(hessian(f, [x1, x2]));
    det_H_func = matlabFunction(det_H, 'Vars', [x1, x2]);
    inflection_level = fminsearch(@(x) abs(det_H_func(x(1), x(2))), [0, 0]);
    inflection_pdf_level = double(subs(f, [x1, x2], inflection_level));
end

%% Save the Figure
function saveFigure(file_name, width, height)
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0, 0, width, height]);
    set(gcf, 'PaperSize', [width, height]);
    print(file_name, '-depsc', '-r300');
end

%% Customize Axes
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
