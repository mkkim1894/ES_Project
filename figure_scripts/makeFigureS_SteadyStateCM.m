function makeFigureS_SteadyStateCM(dataFile, varargin)
% makeFigureS_SteadyStateCM - Generate supplementary figures for CM steady-state validation.
%
% Description:
%   Creates figures comparing Wright-Fisher simulation results against 
%   analytical predictions for adaptation rates under concurrent mutations.
%
% Inputs:
%   dataFile - Path to .mat file from simulateSteadyStateCM (default: 'SteadyStateCM_N_10000.mat')
%   varargin - Optional name-value pairs:
%       'outputDir'  - Directory for output figures (default: current)
%       'figureSet'  - Which figures to generate: 'all', 'main', or figure number (default: 'all')
%
% Outputs:
%   Generates EPS files:
%       FigureS_CM_Main.eps         - Main validation (6 subplots, 3x2 layout)
%       FigureS_CM_Weights.eps      - Weighted prediction comparisons (3 panels)
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
    generateWeightsFigure();
end

fprintf('Figures saved to %s\n', outputDir);

%% Nested figure generation functions

    function generateMainFigure()
        figure('Units', 'centimeters', 'Position', [1, 1, figureWidth, figureHeight*1.3]);

        % 3x2 layout: row1 = A,B; row2 = C,D; row3 = E,F (parameter panels)
        subplotPositions = [
            0.08 0.72 0.38 0.22;   % A
            0.51 0.72 0.38 0.22;   % B
            0.08 0.41 0.38 0.22;   % C
            0.51 0.41 0.38 0.22;   % D
            0.08 0.10 0.38 0.22;   % E (s1 vs s2)
            0.51 0.10 0.38 0.22;   % F (U1 vs U2)
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
                'MarkerFaceAlpha', alpha_val, 'SizeData', 15);
            scatter(v1, v2, 'Marker', 'o', 'MarkerEdgeColor', color_palette(color_idx, :), ...
                'MarkerFaceColor', color_palette(color_idx, :), 'MarkerEdgeAlpha', 1, ...
                'MarkerFaceAlpha', 0, 'SizeData', 15);
        end

        plot([10^-7, 10^0], [10^-7, 10^0], 'k--', 'LineWidth', 1);
        set(gca, 'XScale', 'log', 'YScale', 'log');
        xlim([10^-7, 10^-1]); ylim([10^-7, 10^-1]);
        hx = xlabel('$v_1^\prime$', 'Interpreter', 'latex', 'FontSize', 12);
        hx.Units = 'normalized'; hx.Position(2) = -0.08;
        ylabel('$v_2^\prime$', 'Interpreter', 'latex', 'FontSize', 12);
        box on; set(gca, 'LineWidth', 1);
        addSubplotLabel('A', subplotPositions(1,:));
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
                'MarkerFaceAlpha', alpha_val, 'SizeData', 12);
        end

        plot([10^-7, 10^0], [10^-7, 10^0], 'k--', 'LineWidth', 1);
        xline(3*10^-7, 'k--', 'LineWidth', 0.8);
        set(gca, 'XScale', 'log', 'YScale', 'log');
        xlim([10^-7, 10^-1]); ylim([10^-7, 10^-1]);
        hx = xlabel('$v_{1}$', 'Interpreter', 'latex', 'FontSize', 12);
        hx.Units = 'normalized'; hx.Position(2) = -0.08;
        ylabel('$v_{2}$', 'Interpreter', 'latex', 'FontSize', 12);
        box on; set(gca, 'LineWidth', 1);
        addSubplotLabel('B', subplotPositions(2,:));
        pbaspect([1 1 1]);

        % Subplots C, D
        exclude_elements = [1, 2, 3, 7, 8, 12];
        shortened_vector = setdiff(1:21, exclude_elements);

        % C: Stalling regime (v2/v1 > 100)
        subplot('Position', subplotPositions(3, :));
        plotTheoryVsSimSubplot(exclude_elements, Parameter, WFsim, color_palette, map_index_to_color, alpha_val, 'C', '$\frac{v_{2}}{v_{1}} > 100$', subplotPositions(3,:));

        % D: Simultaneous improvement regime (v2/v1 <= 100)
        subplot('Position', subplotPositions(4, :));
        plotTheoryVsSimSubplot(shortened_vector, Parameter, WFsim, color_palette, map_index_to_color, alpha_val, 'D', '$\frac{v_{2}}{v_{1}} \leq 100$', subplotPositions(4,:));

        % E: s1 vs s2 parameter ranges
        subplot('Position', subplotPositions(5, :));
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
        xlabel('$s_{1}$', 'Interpreter', 'latex', 'FontSize', 12);
        ylabel('$s_{2}$', 'Interpreter', 'latex', 'FontSize', 12);
        box on; set(gca, 'LineWidth', 1);
        addSubplotLabel('E', subplotPositions(5,:));
        pbaspect([1 1 1]);

        % F: U1 vs U2 parameter ranges
        subplot('Position', subplotPositions(6, :));
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
        xlabel('$U_{1}$', 'Interpreter', 'latex', 'FontSize', 12);
        ylabel('$U_{2}$', 'Interpreter', 'latex', 'FontSize', 12);
        box on; set(gca, 'LineWidth', 1);
        addSubplotLabel('F', subplotPositions(6,:));
        pbaspect([1 1 1]);

        saveFigure(fullfile(outputDir, 'FigureS_CM_Main.eps'), figureWidth, figureHeight*1.2);
    end


    function generateWeightsFigure()
        % 3-panel figure: removed redundant vprime panel (now S1D), relabeled B->A, C->B, D->C
        figure('Units', 'centimeters', 'Position', [1, 1, 14, 10]);

        subplotPositions = [
            0.10  0.1  0.22  0.75;   % A (equal weights)
            0.43  0.1  0.22  0.75;   % B (U weights)
            0.76  0.1  0.22  0.75;   % C (s weights)
        ];

        shortened_vector = setdiff(1:21, [1, 2, 3, 7, 8, 12]);

        % A: Equal weights
        subplot('Position', subplotPositions(1, :));
        plotWeightedComparison('equal', shortened_vector, Parameter, WFsim, Prediction, color_palette, map_index_to_color, alpha_val, 'A', subplotPositions(1,:));

        % B: U weights
        subplot('Position', subplotPositions(2, :));
        plotWeightedComparison('U', shortened_vector, Parameter, WFsim, Prediction_U_weights, color_palette, map_index_to_color, alpha_val, 'B', subplotPositions(2,:));

        % C: s weights
        subplot('Position', subplotPositions(3, :));
        plotWeightedComparison('s', shortened_vector, Parameter, WFsim, Prediction_s_weights, color_palette, map_index_to_color, alpha_val, 'C', subplotPositions(3,:));

        saveFigure(fullfile(outputDir, 'FigureS_CM_Weights.eps'), 14, 10);
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

function plotTheoryVsSimSubplot(indices, Parameter, WFsim, color_palette, map_index_to_color, alpha_val, label, regime_text, pos)
    hold on
    [ID_1, ID_2] = meshgrid(1:numel(Parameter.a), 1:numel(Parameter.a));

    scatter1 = []; scatter2 = [];
    for i = indices
        c = color_palette(map_index_to_color(i), :);
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

        scatter1 = scatter(v1, v1_2, 'Marker', 'x', 'MarkerEdgeColor', c, ...
            'MarkerEdgeAlpha', alpha_val, 'SizeData', 15);
        scatter2 = scatter(v2, v2_1, 'Marker', 'o', 'MarkerEdgeColor', c, ...
            'MarkerFaceColor', c, 'MarkerEdgeAlpha', alpha_val, ...
            'MarkerFaceAlpha', alpha_val, 'SizeData', 15);
    end

    plot([10^-7, 1], [10^-7, 1], 'k--', 'LineWidth', 1);
    set(gca, 'XScale', 'log', 'YScale', 'log');
    xlim([2*10^-7, 10^-1]); ylim([2*10^-7, 10^-1]);
    set(gca, 'XTick', [10^-5, 10^-3, 10^-1]);
    set(gca, 'YTick', [10^-5, 10^-3, 10^-1]);
    hx = xlabel('$v_i^\prime$', 'Interpreter', 'latex', 'FontSize', 12);
    hx.Units = 'normalized'; hx.Position(2) = -0.08;
    ylabel('$v_i$', 'Interpreter', 'latex', 'FontSize', 12);
    gray = [0.5 0.5 0.5];
    leg1 = scatter(nan, nan, 'Marker', 'x', 'MarkerEdgeColor', gray, 'SizeData', 15);
    leg2 = scatter(nan, nan, 'Marker', 'o', 'MarkerEdgeColor', gray, 'MarkerFaceColor', gray, 'SizeData', 15);
    legend([leg1, leg2], {'$v_1$', '$v_2$'}, 'Location', 'southeast', 'Interpreter', 'latex');
    box on; set(gca, 'LineWidth', 1);
    addSubplotLabel(label, pos);
    text(0.1, 0.7, regime_text, 'Interpreter', 'latex', 'Units', 'normalized', 'FontSize', 10);
    pbaspect([1 1 1]);
end

function plotWeightedComparison(weightType, indices, Parameter, WFsim, Prediction, color_palette, map_index_to_color, alpha_val, label, pos)
    hold on
    [ID_1, ID_2] = meshgrid(1:numel(Parameter.a), 1:numel(Parameter.a));

    scatter1 = []; scatter2 = [];
    for i = indices
        c = color_palette(map_index_to_color(i), :);
        y_axis = WFsim.v1_2{i};
        s1 = zeros(numel(ID_1), 1);
        for j = 1:numel(ID_1)
            s1(j,1) = Parameter.s1U1{i, ID_1(j)}(1);
        end
        v1_min = (2.*s1) ./ (WFsim.generations - WFsim.t_burnin + 1);
        y_axis(v1_min > WFsim.v1_2{i}) = NaN;

        if strcmp(weightType, 'vprime')
            v1 = zeros(numel(ID_1), 1);
            for j = 1:numel(ID_1)
                v1(j,1) = Parameter.v1_check(i, ID_1(j))';
            end
            v1(v1_min > WFsim.v1_2{i}) = NaN;
            x_data = v1;
        else
            predict = Prediction.v1_2(i, :)';
            predict(v1_min > WFsim.v1_2{i}) = NaN;
            x_data = predict;
        end
        scatter1 = scatter(x_data, y_axis, 'Marker', 'x', 'MarkerEdgeColor', c, ...
            'MarkerEdgeAlpha', alpha_val, 'SizeData', 10);

        y_axis = WFsim.v2_1{i};
        s2 = zeros(numel(ID_1), 1);
        for j = 1:numel(ID_1)
            s2(j,1) = Parameter.s2U2{i, ID_2(j)}(1);
        end
        v2_min = (2.*s2) ./ (WFsim.generations - WFsim.t_burnin + 1);
        y_axis(v2_min > WFsim.v2_1{i}) = NaN;

        if strcmp(weightType, 'vprime')
            v2 = zeros(numel(ID_1), 1);
            for j = 1:numel(ID_1)
                v2(j,1) = Parameter.v2_check(i, ID_2(j))';
            end
            v2(v2_min > WFsim.v2_1{i}) = NaN;
            x_data = v2;
        else
            predict = Prediction.v2_1(i, :)';
            predict(v2_min > WFsim.v2_1{i}) = NaN;
            x_data = predict;
        end
        scatter2 = scatter(x_data, y_axis, 'Marker', 'o', 'MarkerEdgeColor', c, ...
            'MarkerFaceColor', c, 'MarkerEdgeAlpha', alpha_val, ...
            'MarkerFaceAlpha', alpha_val, 'SizeData', 10);
    end

    plot([10^-8, 1], [10^-8, 1], 'k--', 'LineWidth', 1);
    set(gca, 'XScale', 'log', 'YScale', 'log');
    xlim([10^-7, 10^-1]); ylim([10^-7, 10^-1]);
    set(gca, 'XTick', [10^-5, 10^-3, 10^-1]);
    set(gca, 'YTick', [10^-5, 10^-3, 10^-1]);

    switch weightType
        case 'vprime'
            hx = xlabel('$v_i^\prime$', 'Interpreter', 'latex', 'FontSize', 10);
            hx.Units = 'normalized'; hx.Position(2) = -0.12;
        case 'equal'
            hx = xlabel('$v_i^*, w_i = \frac{1}{2}$', 'Interpreter', 'latex', 'FontSize', 10);
            hx.Units = 'normalized'; hx.Position(2) = -0.12;
        case 'U'
            hx = xlabel('$v_i^*, w_i = \frac{U_i}{U}$', 'Interpreter', 'latex', 'FontSize', 10);
            hx.Units = 'normalized'; hx.Position(2) = -0.12;
        case 's'
            hx = xlabel('$v_i^*, w_i = \frac{s_i}{s_1+s_2}$', 'Interpreter', 'latex', 'FontSize', 10);
            hx.Units = 'normalized'; hx.Position(2) = -0.12;
    end
    ylabel('$v_i$', 'Interpreter', 'latex', 'FontSize', 10);
    gray = [0.5 0.5 0.5];
    leg1 = scatter(nan, nan, 'Marker', 'x', 'MarkerEdgeColor', gray, 'SizeData', 10);
    leg2 = scatter(nan, nan, 'Marker', 'o', 'MarkerEdgeColor', gray, 'MarkerFaceColor', gray, 'SizeData', 10);
    legend([leg1, leg2], {'$v_1$', '$v_2$'}, 'Location', 'southeast', 'Interpreter', 'latex');
    box on;
    text(-0.25, 1.15, label, 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold');
    text(0.1, 0.7, '$\frac{v_{2}}{v_{1}} \leq 100$', 'Interpreter', 'latex', 'Units', 'normalized', 'FontSize', 10);
    pbaspect([1 1 1]);
end

function saveFigure(filename, ~, ~)
    print(filename, '-depsc', '-r300');
end

function addSubplotLabel(label, pos)
    annotation('textbox', [pos(1)-0.02, pos(2)+pos(4), 0.03, 0.03], ...
        'String', label, 'FontSize', 12, 'FontWeight', 'bold', 'EdgeColor', 'none');
end