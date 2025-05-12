function makeFigure4(mutationRate)
    % Function to create Figure 4 for a specified mutation rate
    % Input:
    %   mutationRate: Mutation rate (e.g., 5.0e-4) as a numeric value

    % Clear workspace and command window
    clearvars -except mutationRate
    clc

    % Labels for each subplot
    markers = {'d', 'o'};

    % Set the figure size according to PNAS guidelines (in centimeters)
    figureWidth = 8; % Full page width in cm
    figureHeight = 6;

    % Create a figure
    figure('Units', 'centimeters', 'Position', [1, 1, figureWidth, figureHeight]);

    customcolor = [0.4, 0.4, 0.4];
    R_list = zeros(1, 6);  % Pre-allocate R_list

    % Load files dynamically based on mutation rate
    files = dir(sprintf('results/CM_IntermediateRecomb/ModularFGM_CM_Recomb_R*_M%.1e.mat', mutationRate));
    if isempty(files)
        error('No files found for mutation rate %.1e', mutationRate);
    end
    
    % Iterate over the files
    for r_i = 1:length(files)    
        fname = fullfile(files(r_i).folder, files(r_i).name);   
        load(fname);  % Load the data file
    
        % Retrieve recombination rate from saved simParams    
        R_list(r_i) = simParams.recombinationRate;
   
        for ind = 1:length(simParams.initialAngles)                           
            logRatio_tmp(ind,:) = averageTrajectory.logRatio{ind};        
            logRatio(r_i,ind) = logRatio_tmp(ind,end);                      
        end
    end

    % Define ESL and plot y-lines
    yline(log(simParams.ellipseParams(2)^2 / simParams.ellipseParams(1)^2), '-', 'LineWidth', 1, 'Color', customcolor);

    % Apply pseudo-log transformation to handle zero
    epsilon = 5e-5;  % Small constant to avoid log(0)
    pseudo_log_R_list = log10(R_list + epsilon);  % Pseudo-log transformation

    % Sort R_list and corresponding pseudo_log_R_list in ascending order
    [pseudo_log_R_list_sorted, sortIdx] = sort(pseudo_log_R_list);
    R_list_sorted = R_list(sortIdx);

    % Compute Initial Ratios
    Initial_ln_R = log(tan(simParams.initialAngles));

    % Adjust y-axis limits to include both Initial_ln_R values
    y_min = min([Initial_ln_R, logRatio(:)']) - 0.5;  % Subtract a margin
    y_max = max([Initial_ln_R, logRatio(:)']) + 0.5;  % Add a margin
    ylim([y_min, y_max]);

    % Plotting means without error bars and connecting final ratios with lines
    for j_dep = 1:length(simParams.initialAngles)
        Dep_Variable = logRatio(:, j_dep);
        hold on
        % Plot final ratios (means) with markers and lines connecting them
        plot(pseudo_log_R_list_sorted, Dep_Variable(sortIdx), '-k', 'Marker', markers{j_dep}, ...
            'MarkerFaceColor', 'k', 'LineWidth', 0.5, 'MarkerSize', 3);
        % Plot initial ratios at the same x-positions, using the same markers
        plot(pseudo_log_R_list_sorted, repmat(Initial_ln_R(j_dep), size(pseudo_log_R_list_sorted)), ...
            'k', 'Marker', markers{j_dep}, 'MarkerFaceColor', 'w', 'LineStyle', 'none', ...
            'LineWidth', 0.5, 'MarkerSize', 3);

        % Now, for each data point, draw an arrow from initial ratio to final ratio
        for idx = 1:length(pseudo_log_R_list_sorted)
            x_data = pseudo_log_R_list_sorted(idx);
            y_start = Initial_ln_R(j_dep);
            y_end = Dep_Variable(sortIdx(idx));
            arrow_fraction = 0.8;  % Arrow length fraction
            delta_y = y_end - y_start;
            delta_y_adj = delta_y * arrow_fraction;
            y_offset = (delta_y - delta_y_adj) / 2;
            y_start_adj = y_start + y_offset;
            y_end_adj = y_start_adj + delta_y_adj;

            % Draw the arrow shaft as a line
            line([x_data, x_data], [y_start_adj, y_end_adj], 'Color', 'k', 'LineWidth', 0.5);

            % Draw the arrowhead with fixed size
            arrowhead_length = 0.08;
            arrowhead_width = 0.06;

            if delta_y_adj >= 0
                arrow_direction = 1;
            else
                arrow_direction = -1;
            end

            x_arrow = x_data;
            y_arrow = y_end_adj;
            x_vertices = [x_arrow, x_arrow - arrowhead_width / 2, x_arrow + arrowhead_width / 2];
            y_vertices = [y_arrow, y_arrow - arrow_direction * arrowhead_length, y_arrow - arrow_direction * arrowhead_length];
            fill(x_vertices, y_vertices, 'k');
        end
    end

    % Adjusting plot aesthetics
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 8);
    xlabel('$\rho$', 'Interpreter', 'latex', 'FontWeight', 'bold', 'FontSize', 10);
    ylabel('Ln R', 'FontName', 'Helvetica', 'FontSize', 10);

    % Set custom ticks to reflect original R_list values
    set(gca, 'XTick', pseudo_log_R_list_sorted);
    set(gca, 'XTickLabel', arrayfun(@(x) sprintf('%.2g', x), R_list_sorted, 'UniformOutput', false));
    xlim([min(pseudo_log_R_list_sorted) - 0.1, max(pseudo_log_R_list_sorted) + 0.1]);

    grid off

    % Adjust figure settings to match publication requirements
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 figureWidth figureHeight]);
    set(gcf, 'PaperSize', [figureWidth figureHeight]);

    saveFigure('figure4', figureWidth, figureHeight);

end

%% Save the Figure
function saveFigure(file_name, width, height)
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0, 0, width, height]);
    set(gcf, 'PaperSize', [width, height]);
    print(file_name, '-depsc', '-r300');
end
