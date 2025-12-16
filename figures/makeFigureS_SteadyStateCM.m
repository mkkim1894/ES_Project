function makeFigureS_SteadyStateCM(dataFile, varargin)
% makeFigureS_SteadyStateCM - Generate supplementary figures for CM steady-state validation.
%
% Description:
%   Creates publication-quality figures comparing Wright-Fisher simulation results
%   against analytical predictions for adaptation rates under concurrent mutations.
%
% Inputs:
%   dataFile - Path to .mat file from simulateSteadyStateCM (default: 'SteadyStateCM_N_10000.mat')
%   varargin - Optional name-value pairs:
%       'outputDir'  - Directory for output figures (default: current)
%       'figureSet'  - Which figures to generate: 'all', 'main', or figure number (default: 'all')
%
% Outputs:
%   Generates EPS files:
%       FigureS_CM_Main.eps         - Main validation (5 subplots)
%       FigureS_CM_Parameters.eps   - s and U parameter ranges
%       FigureS_CM_SteadyState.eps  - Fitness trajectory steady-state check
%       FigureS_CM_Ratios.eps       - v2/v1 ratio comparison
%       FigureS_CM_Weights.eps      - Weighted prediction comparisons
%       FigureS_CM_Summary.eps      - Mean squared distance summary
%
% Example:
%   makeFigureS_SteadyStateCM('results/SteadyStateCM_N_10000.mat');
%   makeFigureS_SteadyStateCM('data.mat', 'figureSet', 'main');
%
% Reference:
%   Kim, M., Ardell, S. M., & Kryazhimskiy, S. (2025).
%   "Module-Selection Balance in the Evolution of Modular Organisms."
%
% See also: simulateSteadyStateCM
%
% Copyright (c) 2025 Minkyu Kim, Cornell University
% Licensed under MIT License

%% Parse inputs
p = inputParser;
addRequired(p, 'dataFile', @ischar);
addParameter(p, 'outputDir', '.', @ischar);
addParameter(p, 'figureSet', 'all', @(x) ischar(x) || isnumeric(x));
parse(p, dataFile, varargin{:});

outputDir = p.Results.outputDir;
figureSet = p.Results.figureSet;

if ~isfolder(outputDir)
    mkdir(outputDir);
end

%% Load data
fprintf('Loading data from %s...\n', dataFile);
load(dataFile, 'WFsim', 'WFsim_vi', 'Parameter', 'Prediction', 'Prediction_U_weights', 'Prediction_s_weights');

%% Figure parameters
figureWidth = 17.8;
figureHeight = 13.8;

alpha_val = 0.4;
color = 0.7;
color_palette = [230, 97, 0; 139, 0, 139; 75, 83, 32; 26, 133, 255; 0, 0, 0] / 255;

map_index_to_color = @(i) (any([1,2,7]==i)*1 + any([3,8,12]==i)*2 + ...
    any([4,9,13,16]==i)*3 + any([5,10,14,17,19]==i)*4 + any([6,11,15,18,20,21]==i)*5);

%% Generate requested figures
if strcmp(figureSet, 'all') || strcmp(figureSet, 'main') || isequal(figureSet, 1)
    generateMainFigure();
end

if strcmp(figureSet, 'all') || isequal(figureSet, 2)
    generateParameterFigure();
end

if strcmp(figureSet, 'all') || isequal(figureSet, 3)
    generateSteadyStateFigure();
end

if strcmp(figureSet, 'all') || isequal(figureSet, 4)
    generateRatioFigure();
end

if strcmp(figureSet, 'all') || isequal(figureSet, 5)
    generateWeightsFigure();
end

fprintf('Figures saved to %s\n', outputDir);

%% Nested figure generation functions

    function generateMainFigure()
        figure('Units', 'centimeters', 'Position', [1, 1, figureWidth, figureHeight]);
        
        subplotPositions = [
            0.15 0.55 0.35 0.35;
            0.55 0.55 0.35 0.35;
            0.08 0.1 0.22 0.32;
            0.42 0.1 0.22 0.32;
            0.76 0.1 0.22 0.32;
        ];

        % Subplot A: v1' vs v2'
        subplot('Position', subplotPositions(1, :));
        hold on
        [ID_1, ID_2] = meshgrid(1:numel(Parameter.a), 1:numel(Parameter.a));
        
        v_sim = WFsim_vi.v1;
        v_sim = reshape(v_sim, 4, 6)';
        replicated_v1_sim = replicateSimData(v_sim);
        replicated_v2_sim = replicateSimDataFlip(v_sim);
        
        for i = 1:21  
            color_idx = map_index_to_color(i);
            v1 = zeros(numel(ID_1), 1);
            v1_sim = zeros(numel(ID_1), 1);
            v2 = zeros(numel(ID_1), 1);
            v2_sim = zeros(numel(ID_1), 1);
            
            for j = 1:numel(ID_1)
                v1(j,1) = Parameter.v1_check(i, ID_1(j))';
                v1_sim(j,1) = replicated_v1_sim{i, ID_1(j)}';
                v2(j,1) = Parameter.v2_check(i, ID_2(j))';
                v2_sim(j,1) = replicated_v2_sim{i, ID_2(j)}';
            end
            
            scatter(v1_sim, v2_sim, 'Marker', 'o', 'MarkerEdgeColor', color_palette(color_idx, :), ...
                'MarkerFaceColor', color_palette(color_idx, :), 'MarkerEdgeAlpha', 0, ...
                'MarkerFaceAlpha', alpha_val, 'SizeData', 10);
            scatter(v1, v2, 'Marker', 'o', 'MarkerEdgeColor', color_palette(color_idx, :), ...
                'MarkerFaceColor', color_palette(color_idx, :), 'MarkerEdgeAlpha', 1, ...
                'MarkerFaceAlpha', 0, 'SizeData', 15);
        end
        
        plot([10^-7, 10^0], [10^-7, 10^0], 'k--', 'LineWidth', 1);
        set(gca, 'XScale', 'log', 'YScale', 'log');
        xlim([10^-7, 10^-1]); ylim([10^-7, 10^-1]);
        xlabel('$v_1^\prime$', 'Interpreter', 'latex', 'FontSize', 14);
        ylabel('$v_2^\prime$', 'Interpreter', 'latex', 'FontSize', 14);
        box on; set(gca, 'LineWidth', 0.8);
        text(0.1, 0.9, 'A', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 16);
        pbaspect([1 1 1]);
        
        % Subplot B: v1 vs v2 (simulated)
        subplot('Position', subplotPositions(2, :));
        hold on
        
        for i = 1:21
            v1_2 = WFsim.v1_2{i};
            v2_1 = WFsim.v2_1{i};
            
            s1 = zeros(numel(ID_1), 1);
            s2 = zeros(numel(ID_1), 1);
            for j = 1:numel(ID_1)
                s1(j,1) = Parameter.s1U1{i, ID_1(j)}(1);
                s2(j,1) = Parameter.s2U2{i, ID_2(j)}(1);
            end
            
            v1_min = (2.*s1) ./ (WFsim.generations - WFsim.t_burnin + 1);
            v1_2(v1_min > WFsim.v1_2{i}) = 3*10^-7;
            
            v2_min = (2.*s2) ./ (WFsim.generations - WFsim.t_burnin + 1);
            v2_1(v2_min > WFsim.v2_1{i}) = 3*10^-7;
           
            color_idx = map_index_to_color(i);
            scatter(v1_2, v2_1, 'Marker', 'o', 'MarkerEdgeColor', color_palette(color_idx, :), ...
                'MarkerFaceColor', color_palette(color_idx, :), 'MarkerEdgeAlpha', 0, ...
                'MarkerFaceAlpha', alpha_val, 'SizeData', 10);
        end
        
        plot([10^-7, 10^0], [10^-7, 10^0], 'k--', 'LineWidth', 1);
        xline(3*10^-7, 'k--', 'LineWidth', 0.8);
        set(gca, 'XScale', 'log', 'YScale', 'log');
        xlim([10^-7, 10^-1]); ylim([10^-7, 10^-1]);
        xlabel('$v_{1}$', 'Interpreter', 'latex', 'FontSize', 14);
        ylabel('$v_{2}$', 'Interpreter', 'latex', 'FontSize', 14);
        box on; set(gca, 'LineWidth', 0.8);
        text(0.1, 0.9, 'B', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 16);
        pbaspect([1 1 1]);
        
        % Subplots C, D, E: Theory vs simulation for different regimes
        exclude_elements = [1, 2, 3, 7, 8, 12];
        shortened_vector = setdiff(1:21, exclude_elements);
        
        % C: Dominant regime (v2/v1 > 100)
        subplot('Position', subplotPositions(3, :));
        plotTheoryVsSimSubplot(exclude_elements, Parameter, WFsim, color, alpha_val, 'C', '$\frac{v_{2}}{v_{1}} > 100$');
        
        % D: Intermediate regime
        subplot('Position', subplotPositions(4, :));
        plotTheoryVsSimSubplot(shortened_vector, Parameter, WFsim, color, alpha_val, 'D', '$\frac{v_{2}}{v_{1}} \leq 100$');
        
        % E: Weighted predictions
        subplot('Position', subplotPositions(5, :));
        plotWeightedPredictions(shortened_vector, Parameter, WFsim, Prediction_s_weights, color, alpha_val, 'E');
        
        saveFigure(fullfile(outputDir, 'FigureS_CM_Main.eps'), figureWidth, figureHeight);
    end

    function generateParameterFigure()
        figure('Units', 'centimeters', 'Position', [1, 1, 14, 7]);
        
        [ID_1, ID_2] = meshgrid(1:numel(Parameter.a), 1:numel(Parameter.a));
        
        % s1 vs s2
        subplot(1, 2, 1);
        hold on
        for i = 1:21
            s1 = zeros(numel(ID_1), 1);
            s2 = zeros(numel(ID_1), 1);
            for j = 1:numel(ID_1)
                s1(j,1) = Parameter.s1U1{i, ID_1(j)}(1);
                s2(j,1) = Parameter.s2U2{i, ID_2(j)}(1);
            end
            color_idx = map_index_to_color(i);
            scatter(s1, s2, 'Marker', 'o', 'MarkerEdgeColor', color_palette(color_idx, :), ...
                'MarkerFaceColor', color_palette(color_idx, :), 'MarkerEdgeAlpha', alpha_val, ...
                'MarkerFaceAlpha', alpha_val, 'SizeData', 25);
        end
        plot([10^-3, 4*10^-1], [10^-3, 4*10^-1], 'k--', 'LineWidth', 1);
        set(gca, 'XScale', 'log', 'YScale', 'log');
        xlim([10^-3, 4*10^-1]); ylim([10^-3, 4*10^-1]);
        xlabel('$s_{1}$', 'Interpreter', 'latex', 'FontSize', 14);
        ylabel('$s_{2}$', 'Interpreter', 'latex', 'FontSize', 14);
        box on; text(0.1, 0.9, 'A', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 16);
        pbaspect([1 1 1]);
        
        % U1 vs U2
        subplot(1, 2, 2);
        hold on
        for i = 1:21
            U1 = zeros(numel(ID_1), 1);
            U2 = zeros(numel(ID_1), 1);
            for j = 1:numel(ID_1)
                U1(j,1) = Parameter.s1U1{i, ID_1(j)}(2);
                U2(j,1) = Parameter.s2U2{i, ID_2(j)}(2);
            end
            color_idx = map_index_to_color(i);
            scatter(U1, U2, 'Marker', 'o', 'MarkerEdgeColor', color_palette(color_idx, :), ...
                'MarkerFaceColor', color_palette(color_idx, :), 'MarkerEdgeAlpha', alpha_val, ...
                'MarkerFaceAlpha', alpha_val, 'SizeData', 25);
        end
        plot([10^-5, 2*10^-2], [10^-5, 2*10^-2], 'k--', 'LineWidth', 1);
        set(gca, 'XScale', 'log', 'YScale', 'log');
        xlim([2*10^-5, 2*10^-2]); ylim([2*10^-5, 2*10^-2]);
        xlabel('$U_{1}$', 'Interpreter', 'latex', 'FontSize', 14);
        ylabel('$U_{2}$', 'Interpreter', 'latex', 'FontSize', 14);
        box on; text(0.1, 0.9, 'B', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 16);
        pbaspect([1 1 1]);
        
        saveFigure(fullfile(outputDir, 'FigureS_CM_Parameters.eps'), 14, 7);
    end

    function generateSteadyStateFigure()
        figure('Units', 'centimeters', 'Position', [1, 1, 14, 7]);
        
        alph = 0.05;
        
        subplot(1, 2, 1);
        hold on
        for i = 1:numel(WFsim.F1_2)
            Mat = WFsim.F1_2{i};
            for j = 1:size(Mat, 1)
                plt1 = plot(Mat(j, :), 'k');
                plt1.Color(4) = alph;
            end
        end
        xline(WFsim.t_burnin, 'k--');
        set(gca, 'YScale', 'log');
        xlabel('$t$, generations', 'Interpreter', 'latex');
        ylabel('$F_{1}$, log fitness', 'Interpreter', 'latex');
        text(0.1, 0.9, 'A', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 16);
        ylim([10^-8, 10^4]);
        box on; pbaspect([1 1 1]);
        
        subplot(1, 2, 2);
        hold on
        for i = 1:numel(WFsim.F2_1)
            Mat = WFsim.F2_1{i};
            for j = 1:size(Mat, 1)
                plt2 = plot(Mat(j, :), 'k');
                plt2.Color(4) = alph;
            end
        end
        xline(WFsim.t_burnin, 'k--');
        set(gca, 'YScale', 'log');
        xlabel('$t$, generations', 'Interpreter', 'latex');
        ylabel('$F_{2}$, log fitness', 'Interpreter', 'latex');
        text(0.1, 0.9, 'B', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 16);
        ylim([10^-8, 10^4]);
        box on; pbaspect([1 1 1]);
        
        saveFigure(fullfile(outputDir, 'FigureS_CM_SteadyState.eps'), 14, 7);
    end

    function generateRatioFigure()
        figure('Units', 'centimeters', 'Position', [1, 1, 8, 8]);
        hold on
        
        [ID_1, ID_2] = meshgrid(1:numel(Parameter.a), 1:numel(Parameter.a));
        
        for i = 1:21
            v1_2 = WFsim.v1_2{i};
            v2_1 = WFsim.v2_1{i};
            v1 = zeros(numel(ID_1), 1);
            v2 = zeros(numel(ID_1), 1);
            s1 = zeros(numel(ID_1), 1);
            s2 = zeros(numel(ID_1), 1);
            
            for j = 1:numel(ID_1)
                s1(j,1) = Parameter.s1U1{i, ID_1(j)}(1);
                v1(j,1) = Parameter.v1_check(i, ID_1(j))';
                s2(j,1) = Parameter.s2U2{i, ID_2(j)}(1);
                v2(j,1) = Parameter.v2_check(i, ID_2(j))';
            end
            
            v1_min = (2.*s1) ./ (WFsim.generations - WFsim.t_burnin + 1);
            v1(v1_min > WFsim.v1_2{i}) = NaN;
            v1_2(v1_min > WFsim.v1_2{i}) = NaN;
            
            v2_min = (2.*s2) ./ (WFsim.generations - WFsim.t_burnin + 1);
            v2(v2_min > WFsim.v2_1{i}) = NaN;
            v2_1(v2_min > WFsim.v2_1{i}) = NaN;
           
            color_idx = map_index_to_color(i);
            scatter(v2./v1, v2_1./v1_2, 'Marker', 'o', 'MarkerEdgeColor', color_palette(color_idx, :), ...
                'MarkerFaceColor', color_palette(color_idx, :), 'MarkerEdgeAlpha', alpha_val*2, ...
                'MarkerFaceAlpha', alpha_val, 'SizeData', 12);
        end
        
        plot([10^-2, 10^6], [10^-2, 10^6], 'k--', 'LineWidth', 1);
        set(gca, 'XScale', 'log', 'YScale', 'log');
        xlim([10^-2, 10^5]); ylim([10^-2, 10^5]);
        xlabel('$v_2^\prime/v_1^\prime$', 'Interpreter', 'latex', 'FontSize', 14);
        ylabel('$v_2/v_1$', 'Interpreter', 'latex', 'FontSize', 14);
        box on; pbaspect([1 1 1]);
        
        saveFigure(fullfile(outputDir, 'FigureS_CM_Ratios.eps'), 8, 8);
    end

    function generateWeightsFigure()
        figure('Units', 'centimeters', 'Position', [1, 1, 14, 14]);
        
        subplotPositions = [
            0.1  0.55  0.4  0.35; 
            0.55  0.55  0.4  0.35;  
            0.1  0.1  0.4  0.35;  
            0.55  0.1  0.4  0.35; 
        ];
        
        shortened_vector = setdiff(1:21, [1, 2, 3, 7, 8, 12]);
        [ID_1, ID_2] = meshgrid(1:numel(Parameter.a), 1:numel(Parameter.a));
        
        % A: v' predictions
        subplot('Position', subplotPositions(1, :));
        plotWeightedComparison('vprime', shortened_vector, Parameter, WFsim, [], color, alpha_val, 'A');
        
        % B: Equal weights
        subplot('Position', subplotPositions(2, :));
        plotWeightedComparison('equal', shortened_vector, Parameter, WFsim, Prediction, color, alpha_val, 'B');
        
        % C: U weights
        subplot('Position', subplotPositions(3, :));
        plotWeightedComparison('U', shortened_vector, Parameter, WFsim, Prediction_U_weights, color, alpha_val, 'C');
        
        % D: s weights
        subplot('Position', subplotPositions(4, :));
        plotWeightedComparison('s', shortened_vector, Parameter, WFsim, Prediction_s_weights, color, alpha_val, 'D');
        
        saveFigure(fullfile(outputDir, 'FigureS_CM_Weights.eps'), 14, 14);
    end

end

%% Helper functions

function rep = replicateSimData(v_sim)
    rep = cell(21, 4);
    rep(1:6, :) = repmat(v_sim(1, :), 6, 1);
    rep(7:11, :) = repmat(v_sim(2, :), 5, 1);
    rep(12:15, :) = repmat(v_sim(3, :), 4, 1);
    rep(16:18, :) = repmat(v_sim(4, :), 3, 1);
    rep(19:20, :) = repmat(v_sim(5, :), 2, 1);
    rep(21, :) = v_sim(6, :);
end

function rep = replicateSimDataFlip(v_sim)
    rep = cell(21, 4);
    rep(1:6, :) = flip(v_sim(1:6, :));
    rep(7:11, :) = flip(v_sim(2:6, :));
    rep(12:15, :) = flip(v_sim(3:6, :));
    rep(16:18, :) = flip(v_sim(4:6, :));
    rep(19:20, :) = flip(v_sim(5:6, :));
    rep(21, :) = v_sim(6, :);
end

function plotTheoryVsSimSubplot(indices, Parameter, WFsim, color, alpha_val, label, regime_text)
    hold on
    [ID_1, ID_2] = meshgrid(1:numel(Parameter.a), 1:numel(Parameter.a));
    
    for i = indices
        v1_2 = WFsim.v1_2{i};
        v2_1 = WFsim.v2_1{i};
        v1 = zeros(numel(ID_1), 1);
        v2 = zeros(numel(ID_1), 1);
        s1 = zeros(numel(ID_1), 1);
        s2 = zeros(numel(ID_1), 1);
        
        for j = 1:numel(ID_1)
            s1(j,1) = Parameter.s1U1{i, ID_1(j)}(1);
            v1(j,1) = Parameter.v1_check(i, ID_1(j))';
            s2(j,1) = Parameter.s2U2{i, ID_2(j)}(1);
            v2(j,1) = Parameter.v2_check(i, ID_2(j))';
        end
        
        v1_min = (2.*s1) ./ (WFsim.generations - WFsim.t_burnin + 1);
        v1(v1_min > WFsim.v1_2{i}) = NaN;
        v1_2(v1_min > WFsim.v1_2{i}) = NaN;
        
        v2_min = (2.*s2) ./ (WFsim.generations - WFsim.t_burnin + 1);
        v2(v2_min > WFsim.v2_1{i}) = NaN;
        v2_1(v2_min > WFsim.v2_1{i}) = NaN;
               
        scatter1 = scatter(v1, v1_2, 'Marker', 'x', 'MarkerEdgeColor', [color 0 0], ...
            'MarkerEdgeAlpha', alpha_val, 'SizeData', 9);
        scatter2 = scatter(v2, v2_1, 'Marker', '+', 'MarkerEdgeColor', [0 0 color], ...
            'MarkerEdgeAlpha', alpha_val, 'SizeData', 9);
    end
    
    plot([10^-7, 1], [10^-7, 1], 'k--', 'LineWidth', 1);
    set(gca, 'XScale', 'log', 'YScale', 'log');
    xlim([2*10^-7, 10^-1]); ylim([2*10^-7, 10^-1]);
    xlabel('$v_i^\prime$', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('$v_i$', 'Interpreter', 'latex', 'FontSize', 14);
    legend([scatter1, scatter2], {'$v_1$', '$v_2$'}, 'Location', 'southeast', 'Interpreter', 'latex');
    box on;
    text(0.1, 0.9, label, 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 16);
    text(0.1, 0.7, regime_text, 'Interpreter', 'latex', 'Units', 'normalized', 'FontSize', 10);
    pbaspect([1 1 1]);
end

function plotWeightedPredictions(indices, Parameter, WFsim, Prediction, color, alpha_val, label)
    hold on
    [ID_1, ID_2] = meshgrid(1:numel(Parameter.a), 1:numel(Parameter.a));
    
    for i = indices
        y_axis = WFsim.v1_2{i};
        s1 = zeros(numel(ID_1), 1);
        for j = 1:numel(ID_1)
            s1(j,1) = Parameter.s1U1{i, ID_1(j)}(1);
        end
        v1_min = (2.*s1) ./ (WFsim.generations - WFsim.t_burnin + 1);
        y_axis(v1_min > WFsim.v1_2{i}) = NaN;
        
        predict = Prediction.v1_2(i, :)';
        predict(v1_min > WFsim.v1_2{i}) = NaN;
        scatter1 = scatter(predict, y_axis, 'Marker', 'x', 'MarkerEdgeColor', [color 0 0], ...
            'MarkerEdgeAlpha', alpha_val, 'SizeData', 9);

        y_axis = WFsim.v2_1{i};
        s2 = zeros(numel(ID_1), 1);
        for j = 1:numel(ID_1)
            s2(j,1) = Parameter.s2U2{i, ID_2(j)}(1);
        end
        v2_min = (2.*s2) ./ (WFsim.generations - WFsim.t_burnin + 1);
        y_axis(v2_min > WFsim.v2_1{i}) = NaN;
        
        predict = Prediction.v2_1(i, :)';
        predict(v2_min > WFsim.v2_1{i}) = NaN;
        scatter2 = scatter(predict, y_axis, 'Marker', '+', 'MarkerEdgeColor', [0 0 color], ...
            'MarkerEdgeAlpha', alpha_val, 'SizeData', 9);
    end
    
    plot([10^-8, 1], [10^-8, 1], 'k--', 'LineWidth', 1);
    set(gca, 'XScale', 'log', 'YScale', 'log');
    xlim([10^-7, 10^-1]); ylim([10^-7, 10^-1]);
    xlabel('$v_i^*$', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('$v_i$', 'Interpreter', 'latex', 'FontSize', 14);
    legend([scatter1, scatter2], {'$v_1$', '$v_2$'}, 'Location', 'southeast', 'Interpreter', 'latex');
    box on;
    text(0.1, 0.9, label, 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 16);
    text(0.1, 0.7, '$\frac{v_{2}}{v_{1}} \leq 100$', 'Interpreter', 'latex', 'Units', 'normalized', 'FontSize', 10);
    pbaspect([1 1 1]);
end

function plotWeightedComparison(weightType, indices, Parameter, WFsim, Prediction, color, alpha_val, label)
    hold on
    [ID_1, ID_2] = meshgrid(1:numel(Parameter.a), 1:numel(Parameter.a));
    
    for i = indices
        y_axis = WFsim.v1_2{i};
        s1 = zeros(numel(ID_1), 1);
        v1 = zeros(numel(ID_1), 1);
        for j = 1:numel(ID_1)
            s1(j,1) = Parameter.s1U1{i, ID_1(j)}(1);
            v1(j,1) = Parameter.v1_check(i, ID_1(j))';
        end
        v1_min = (2.*s1) ./ (WFsim.generations - WFsim.t_burnin + 1);
        v1(v1_min > WFsim.v1_2{i}) = NaN;
        y_axis(v1_min > WFsim.v1_2{i}) = NaN;
        
        if strcmp(weightType, 'vprime')
            x_data = v1;
        else
            predict = Prediction.v1_2(i, :)';
            predict(v1_min > WFsim.v1_2{i}) = NaN;
            x_data = predict;
        end
        scatter1 = scatter(x_data, y_axis, 'Marker', 'x', 'MarkerEdgeColor', [color 0 0], ...
            'MarkerEdgeAlpha', alpha_val, 'SizeData', 10);

        y_axis = WFsim.v2_1{i};
        s2 = zeros(numel(ID_1), 1);
        v2 = zeros(numel(ID_1), 1);
        for j = 1:numel(ID_1)
            s2(j,1) = Parameter.s2U2{i, ID_2(j)}(1);
            v2(j,1) = Parameter.v2_check(i, ID_2(j))';
        end
        v2_min = (2.*s2) ./ (WFsim.generations - WFsim.t_burnin + 1);
        v2(v2_min > WFsim.v2_1{i}) = NaN;
        y_axis(v2_min > WFsim.v2_1{i}) = NaN;
        
        if strcmp(weightType, 'vprime')
            x_data = v2;
        else
            predict = Prediction.v2_1(i, :)';
            predict(v2_min > WFsim.v2_1{i}) = NaN;
            x_data = predict;
        end
        scatter2 = scatter(x_data, y_axis, 'Marker', '+', 'MarkerEdgeColor', [0 0 color], ...
            'MarkerEdgeAlpha', alpha_val, 'SizeData', 10);
    end
    
    plot([10^-8, 1], [10^-8, 1], 'k--', 'LineWidth', 1);
    set(gca, 'XScale', 'log', 'YScale', 'log');
    xlim([10^-7, 10^-1]); ylim([10^-7, 10^-1]);
    
    switch weightType
        case 'vprime'
            xlabel('$v_i^\prime$', 'Interpreter', 'latex', 'FontSize', 14);
        case 'equal'
            xlabel('$v_i^*, w_i = \frac{1}{2}$', 'Interpreter', 'latex', 'FontSize', 14);
        case 'U'
            xlabel('$v_i^*, w_i = \frac{U_i}{U}$', 'Interpreter', 'latex', 'FontSize', 14);
        case 's'
            xlabel('$v_i^*, w_i = \frac{s_i}{s_1+s_2}$', 'Interpreter', 'latex', 'FontSize', 14);
    end
    ylabel('$v_i$', 'Interpreter', 'latex', 'FontSize', 14);
    legend([scatter1, scatter2], {'$v_1$', '$v_2$'}, 'Location', 'southeast', 'Interpreter', 'latex');
    box on;
    text(0.1, 0.9, label, 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 16);
    text(0.1, 0.7, '$\frac{v_{2}}{v_{1}} \leq 100$', 'Interpreter', 'latex', 'Units', 'normalized', 'FontSize', 10);
    pbaspect([1 1 1]);
end

function saveFigure(filename, width, height)
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 width height]);
    set(gcf, 'PaperSize', [width height]);
    print(filename, '-depsc', '-r300');
    close(gcf);
end
