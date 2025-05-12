%% --- Main plts --- %%
%% 1. Prediction vs. WF
clearvars;
clc;
load('Stdstate_N_10000.mat')

% Set the figure size according to PNAS guidelines (in centimeters)
figureWidth = 17.8; % 2 columns wide in cm
figureHeight = 13.8; % Maximum height in cm

% Create a figure
figure('Units', 'centimeters', 'Position', [1, 1, figureWidth, figureHeight]);

% Define positions and dimensions of subplots
subplotPositions = [
    0.15 0.55 0.35 0.35;  % Position of 'a'
    0.55 0.55 0.35 0.35;  % Position of 'b'
    0.08 0.1 0.22 0.32;  % Position of 'c'
    0.42 0.1 0.22 0.32;  % Position of 'd'
    0.76 0.1 0.22 0.32;  % Position of 'e'
];


alpha_val = 0.4;
color = 0.7;
color_palette = [230, 97, 0; ... % Dark orange
                 139, 0, 139;  ... % Dark magenta
                 75, 83, 32; ... % Dark green
                 26, 133, 255;    ... % Blue
                 0, 0, 0] / 255;   % Black

% Mapping function to assign color index based on i
map_index_to_color = @(i) (any([1,2,7]==i)*1 + any([3,8,12]==i)*2 + ...
    any([4,9,13,16]==i)*3 + any([5,10,14,17,19]==i)*4 + any([6,11,15,18,20,21]==i)*5);

%%% subplot a %%%
subplot('Position', subplotPositions(1, :));
hold on
[ID_1, ID_2] = meshgrid(1:numel(Parameter.a), 1:numel(Parameter.a));

v_sim = WFsim_vi.v1;
v_sim = reshape(v_sim, 4, 6)';
replicated_v1_sim(1:6, :) = repmat(v_sim(1, :), 6, 1); % Replicate the first row 6 times
replicated_v1_sim(7:11, :) = repmat(v_sim(2, :), 5, 1); % Replicate the second row 5 times
replicated_v1_sim(12:15, :) = repmat(v_sim(3, :), 4, 1); % Replicate the third row 4 times
replicated_v1_sim(16:18, :) = repmat(v_sim(4, :), 3, 1); % Replicate the fourth row 3 times
replicated_v1_sim(19:20, :) = repmat(v_sim(5, :), 2, 1); % Replicate the fifth row 2 times
replicated_v1_sim(21, :) = v_sim(6,:);

replicated_v2_sim(1:6, :) = flip(v_sim(1:6, :));
replicated_v2_sim(7:11, :) = flip(v_sim(2:6, :));
replicated_v2_sim(12:15, :) = flip(v_sim(3:6, :));
replicated_v2_sim(16:18, :) = flip(v_sim(4:6, :));
replicated_v2_sim(19:20, :) = flip(v_sim(5:6, :));
replicated_v2_sim(21, :) = v_sim(6,:);

for i = 1:21  
    color_idx = map_index_to_color(i);
    
    for j = 1:numel(ID_1)
        v1(j,1) = Parameter.v1_check(i,ID_1(j))';
        v1_sim(j,1) = replicated_v1_sim{i,ID_1(j)}';

        v2(j,1) = Parameter.v2_check(i,ID_2(j))';
        v2_sim(j,1) = replicated_v2_sim{i,ID_2(j)}';
    end
    
    scatter_handle = scatter(v1_sim, v2_sim, 'Marker', 'o', 'MarkerEdgeColor', color_palette(color_idx, :), 'MarkerFaceColor', color_palette(color_idx, :), 'MarkerEdgeAlpha', 0, 'MarkerFaceAlpha', alpha_val, 'SizeData', 10);
    Theo_vals = scatter(v1, v2, 'Marker', 'o', 'MarkerEdgeColor', color_palette(color_idx, :), 'MarkerFaceColor', color_palette(color_idx, :), 'MarkerEdgeAlpha', 1, 'MarkerFaceAlpha', 0, 'SizeData', 15);
end
plot([10^-7, 10^0], [10^-7, 10^0], 'k--', 'LineWidth', 1); 

set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlim([10^-7, 10^-1])
ylim([10^-7, 10^-1])

xlabel('$v_1^\prime$', 'Interpreter', 'latex', 'FontWeight', 'bold', 'FontSize', 14);
ylabel('$v_2^\prime$', 'Interpreter', 'latex', 'FontWeight', 'bold', 'FontSize', 14);

box on
set(gca, 'LineWidth', 0.8);
text(0.8, 0.2, 'A', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 20);
pbaspect([1 1 1]) 

%%% subplot b %%%
subplot('Position', subplotPositions(2, :));
hold on
[ID_1, ID_2] = meshgrid(1:numel(Parameter.a), 1:numel(Parameter.a));
for i = 1:21
    v1_2 = WFsim.v1_2{i};
    for j = 1:numel(ID_1)
        s1(j,1) = Parameter.s1U1{i,ID_1(j)}(1);
        v1(j,1) = Parameter.v1_check(i,ID_1(j))';
    end
    v1_min = (2.*s1)./(WFsim.generations-WFsim.t_burnin+1);
    v1_2(v1_min>WFsim.v1_2{i}) = 3*10^-7;
    
    v2_1 = WFsim.v2_1{i};
    for j = 1:numel(ID_1)
        s2(j,1) = Parameter.s2U2{i,ID_2(j)}(1);
        v2(j,1) = Parameter.v2_check(i,ID_2(j))';
    end
    v2_min = (2.*s2)./(WFsim.generations-WFsim.t_burnin+1);
    v2_1(v2_min>WFsim.v2_1{i}) = 3*10^-7;
   
    color_idx = map_index_to_color(i);
    
    scatter_handle = scatter(v1_2, v2_1, 'Marker', 'o', 'MarkerEdgeColor', color_palette(color_idx, :), 'MarkerFaceColor', color_palette(color_idx, :), 'MarkerEdgeAlpha', 0, 'MarkerFaceAlpha', alpha_val, 'SizeData', 10);
    %Theo_vals = scatter(v1, v2, 'Marker', 'o', 'MarkerEdgeColor', color_palette(color_idx, :), 'MarkerFaceColor', color_palette(color_idx, :), 'MarkerEdgeAlpha', 1, 'MarkerFaceAlpha', 0, 'SizeData', 15);
end
plot([10^-7, 10^0], [10^-7, 10^0], 'k--', 'LineWidth', 1); 
xline(3*10^-7, 'k--', 'LineWidth', 0.8);
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlim([10^-7, 10^-1])
ylim([10^-7, 10^-1])

xlabel('$v_{1}$', 'Interpreter', 'latex', 'FontWeight', 'bold', 'FontSize', 14);
ylabel('$v_{2}$', 'Interpreter', 'latex', 'FontWeight', 'bold', 'FontSize', 14);

box on
set(gca, 'LineWidth', 0.8);
text(0.8, 0.2, 'B', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 20);
pbaspect([1 1 1]) 

%%% subplot c %%%
% Create a vector from 1 to 21
full_vector = 1:21;

% Define the elements to exclude
exclude_elements = [1, 2, 3, 7, 8, 12];

% Use setdiff to create the vector excluding the specified elements
shortened_vector = setdiff(full_vector, exclude_elements);

subplot('Position', subplotPositions(3, :));
hold on
[ID_1, ID_2] = meshgrid(1:numel(Parameter.a), 1:numel(Parameter.a));
for i = exclude_elements 
    v1_2 = WFsim.v1_2{i};
    v2_1 = WFsim.v2_1{i};
    
    for j = 1:numel(ID_1)
        s1(j,1) = Parameter.s1U1{i,ID_1(j)}(1);
        v1(j,1) = Parameter.v1_check(i,ID_1(j))';
        
        s2(j,1) = Parameter.s2U2{i,ID_2(j)}(1);
        v2(j,1) = Parameter.v2_check(i,ID_2(j))';
    end
    
    v1_min = (2.*s1)./(WFsim.generations-WFsim.t_burnin+1);
    v1(v1_min>WFsim.v1_2{i}) = NaN;
    v1_2(v1_min>WFsim.v1_2{i}) = NaN;
    
    v2_min = (2.*s2)./(WFsim.generations-WFsim.t_burnin+1);
    v2(v2_min>WFsim.v2_1{i}) = NaN;
    v2_1(v2_min>WFsim.v2_1{i}) = NaN;
           
    scatter1 = scatter(v1, v1_2, 'Marker', 'x', 'MarkerEdgeColor', [color 0 0], 'MarkerFaceColor', [color 0 0], 'MarkerEdgeAlpha', alpha_val, 'MarkerFaceAlpha', 0, 'SizeData', 9);
    scatter2 = scatter(v2, v2_1, 'Marker', '+', 'MarkerEdgeColor', [0 0 color], 'MarkerFaceColor', [0 0 color], 'MarkerEdgeAlpha', alpha_val, 'MarkerFaceAlpha', 0, 'SizeData', 9);
end

plot([10^-7, 1], [10^-7, 1], 'k--', 'LineWidth', 1); 

set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlim([2*10^-7, 10^-1])
ylim([2*10^-7, 10^-1])

xlabel('$v_i^\prime$', 'Interpreter', 'latex', 'FontWeight', 'bold', 'FontSize', 16);
ylabel('$v_i$', 'Interpreter', 'latex', 'FontWeight', 'bold', 'FontSize', 16);

% Legend labels
legendLabels = {'$v_1$', '$v_2$'};
% Create the legend
legend([scatter1,scatter2], legendLabels, 'Location', 'southeast', 'Interpreter', 'latex','FontSize', 8);

box on
set(gca, 'LineWidth', 0.8);
text(0.1, 0.85, 'C', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 20);
text(0.1, 0.6, '$\frac{v_{2}}{v_{1}} > 100$', 'Interpreter', 'latex', 'Units', 'normalized', 'FontSize', 12);
pbaspect([1 1 1]) 

%%% subplot d %%%
subplot('Position', subplotPositions(4, :));
hold on
[ID_1, ID_2] = meshgrid(1:numel(Parameter.a), 1:numel(Parameter.a));
for i = shortened_vector  
    v1_2 = WFsim.v1_2{i};
    v2_1 = WFsim.v2_1{i};
    
    for j = 1:numel(ID_1)
        s1(j,1) = Parameter.s1U1{i,ID_1(j)}(1);
        v1(j,1) = Parameter.v1_check(i,ID_1(j))';
        
        s2(j,1) = Parameter.s2U2{i,ID_2(j)}(1);
        v2(j,1) = Parameter.v2_check(i,ID_2(j))';
    end
    
    v1_min = (2.*s1)./(WFsim.generations-WFsim.t_burnin+1);
    v1(v1_min>WFsim.v1_2{i}) = NaN;
    v1_2(v1_min>WFsim.v1_2{i}) = NaN;
    
    v2_min = (2.*s2)./(WFsim.generations-WFsim.t_burnin+1);
    v2(v2_min>WFsim.v2_1{i}) = NaN;
    v2_1(v2_min>WFsim.v2_1{i}) = NaN;
           
    scatter1 = scatter(v1, v1_2, 'Marker', 'x', 'MarkerEdgeColor', [color 0 0], 'MarkerFaceColor', [color 0 0], 'MarkerEdgeAlpha', alpha_val, 'MarkerFaceAlpha', 0, 'SizeData', 9);
    scatter2 = scatter(v2, v2_1, 'Marker', '+', 'MarkerEdgeColor', [0 0 color], 'MarkerFaceColor', [0 0 color], 'MarkerEdgeAlpha', alpha_val, 'MarkerFaceAlpha', 0, 'SizeData', 9);
end

plot([10^-7, 1], [10^-7, 1], 'k--', 'LineWidth', 1); 

set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlim([5*10^-7, 10^-1])
ylim([5*10^-7, 10^-1])

xlabel('$v_i^\prime$', 'Interpreter', 'latex', 'FontWeight', 'bold', 'FontSize', 16);
ylabel('$v_i$', 'Interpreter', 'latex', 'FontWeight', 'bold', 'FontSize', 16);

% Legend labels
legendLabels = {'$v_1$', '$v_2$'};
% Create the legend
legend([scatter1,scatter2], legendLabels, 'Location', 'southeast', 'Interpreter', 'latex','FontSize', 8);

box on
set(gca, 'LineWidth', 0.8);
text(0.1, 0.85, 'D', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 20);
text(0.1, 0.6, '$\frac{v_{2}}{v_{1}} \leq 100$', 'Interpreter', 'latex', 'Units', 'normalized', 'FontSize', 12);
pbaspect([1 1 1]) 

%%% subplot e %%
subplot('Position', subplotPositions(5, :));
hold on
[ID_1, ID_2] = meshgrid(1:numel(Parameter.a), 1:numel(Parameter.a));
for i = shortened_vector 
    y_axis = WFsim.v1_2{i};
    for j = 1:numel(ID_1)
        s1(j,1) = Parameter.s1U1{i,ID_1(j)}(1);
        v1(j,1) = Parameter.v1_check(i,ID_1(j))';
    end
    v1_min = (2.*s1)./(WFsim.generations-WFsim.t_burnin+1);
    
    v1(v1_min>WFsim.v1_2{i}) = NaN;
    y_axis(v1_min>WFsim.v1_2{i}) = NaN;
    
    predict = Prediction_s_weights.v1_2(i,1:end)';
    predict(v1_min>WFsim.v1_2{i}) = NaN;
    scatter_handle01 = scatter(predict, y_axis, 'Marker', 'x', 'MarkerEdgeColor', [color 0 0], 'MarkerFaceColor', [color 0 0], 'MarkerEdgeAlpha', alpha_val, 'MarkerFaceAlpha', 0, 'SizeData', 9);

    y_axis = WFsim.v2_1{i};
    for j = 1:numel(ID_1)
        s2(j,1) = Parameter.s2U2{i,ID_2(j)}(1);
    end
    v2_min = (2.*s2)./(WFsim.generations-WFsim.t_burnin+1);
    
    y_axis(v2_min>WFsim.v2_1{i}) = NaN;
    
    predict = Prediction_s_weights.v2_1(i,1:end)';
    predict(v2_min>WFsim.v2_1{i}) = NaN;
   
    scatter_handle02 = scatter(predict, y_axis, 'Marker', '+', 'MarkerEdgeColor', [0 0 color], 'MarkerFaceColor', [0 0 color], 'MarkerEdgeAlpha', alpha_val, 'MarkerFaceAlpha', 0, 'SizeData', 9);
end
plot([10^-8 1], [10^-8 1], 'k--', 'LineWidth', 1);  

set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlim([10^-7, 10^-1])
ylim([10^-7, 10^-1])

xlabel('$v_i^*$', 'Interpreter', 'latex', 'FontWeight', 'bold', 'FontSize', 14);
ylabel('$v_i$', 'Interpreter', 'latex', 'FontWeight', 'bold', 'FontSize', 14);

% Legend labels
legendLabels = {'$v_1$', '$v_2$'};
% Create the legend
legend([scatter_handle01,scatter_handle02], legendLabels, 'Location', 'southeast', 'Interpreter', 'latex','FontSize', 8);

box on
set(gca, 'LineWidth', 0.8);
text(0.1, 0.85, 'E', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 20);
text(0.1, 0.6, '$\frac{v_{2}}{v_{1}} \leq 100$', 'Interpreter', 'latex', 'Units', 'normalized', 'FontSize', 12);
pbaspect([1 1 1]) 

% Adjust spacing
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 figureWidth figureHeight]);
set(gcf, 'PaperSize', [figureWidth figureHeight]);

% Save the figure as a high-resolution PNG file
print('SI1', '-depsc', '-r300');

%% --- Suplementary plts --- %%
clearvars;
clc;
load('Stdstate_N_10000.mat')

%% Weighted predictions
clf;
color = 0.4;
alpha_val = 0.4;

subplotPositions = [
    0.1  0.55  0.4  0.35; 
    0.55  0.55  0.4  0.35;  
    0.1  0.1  0.4  0.35;  
    0.55  0.1  0.4  0.35; 
];

% v'
subplot('Position', subplotPositions(1, :));
hold on

% Initialize sums and counts for averaging
sumSquaredDiff_v1 = 0;
count_v1 = 0;
sumSquaredDiff_v2 = 0;
count_v2 = 0;

[ID_1, ID_2] = meshgrid(1:numel(Parameter.a), 1:numel(Parameter.a));

% Create a vector from 1 to 21
full_vector = 1:21;

% Define the elements to exclude
exclude_elements = [1, 2, 3, 7, 8, 12];

% Use setdiff to create the vector excluding the specified elements
shortened_vector = setdiff(full_vector, exclude_elements);

for i = shortened_vector
    y_axis = WFsim.v1_2{i};
    for j = 1:numel(ID_1)
        s1(j,1) = Parameter.s1U1{i,ID_1(j)}(1);
        v1(j,1) = Parameter.v1_check(i,ID_1(j))';
    end
    v1_min = (2.*s1)./(WFsim.generations-WFsim.t_burnin+1);
    v1(v1_min>WFsim.v1_2{i}) = NaN;
    y_axis(v1_min>WFsim.v1_2{i}) = NaN;
    scatter_handle01 = scatter(v1, y_axis, 'Marker', 'x', 'MarkerEdgeColor', [color 0 0], 'MarkerFaceColor', [color 0 0], 'MarkerEdgeAlpha', alpha_val, 'MarkerFaceAlpha', 0, 'SizeData', 10);

    % Calculate and accumulate log squared differences for v1
    validIndices_v1 = ~isnan(v1) & ~isnan(y_axis); % Indices of valid (non-NaN) points
    logSquaredDiff_v1 = (log(v1(validIndices_v1)) - log(y_axis(validIndices_v1))).^2;
    sumSquaredDiff_v1 = sumSquaredDiff_v1 + sum(logSquaredDiff_v1);
    count_v1 = count_v1 + numel(logSquaredDiff_v1);
    
    % Reset y_axis for v2 calculation
    y_axis = WFsim.v2_1{i};
    for j = 1:numel(ID_1)
        s2(j,1) = Parameter.s2U2{i,ID_2(j)}(1);
        v2(j,1) = Parameter.v2_check(i,ID_2(j))';
    end
    v2_min = (2.*s2)./(WFsim.generations-WFsim.t_burnin+1);
    v2(v2_min>WFsim.v2_1{i}) = NaN;
    y_axis(v2_min>WFsim.v2_1{i}) = NaN;
    scatter_handle02 = scatter(v2, y_axis, 'Marker', '+', 'MarkerEdgeColor', [0 0 color], 'MarkerFaceColor', [0 0 color], 'MarkerEdgeAlpha', alpha_val, 'MarkerFaceAlpha', alpha_val, 'SizeData', 10);
    
    % Calculate and accumulate log squared differences for v2
    validIndices_v2 = ~isnan(v2) & ~isnan(y_axis); % Indices of valid (non-NaN) points
    logSquaredDiff_v2 = (log(v2(validIndices_v2)) - log(y_axis(validIndices_v2))).^2;
    sumSquaredDiff_v2 = sumSquaredDiff_v2 + sum(logSquaredDiff_v2);
    count_v2 = count_v2 + numel(logSquaredDiff_v2);
end

% Compute the average log squared differences
avgLogSquaredDiff_v1 = sumSquaredDiff_v1 / count_v1;
avgLogSquaredDiff_v2 = sumSquaredDiff_v2 / count_v2;

v_star_Array(1,:) = [avgLogSquaredDiff_v1, avgLogSquaredDiff_v2]; 

plot([10^-8 1], [10^-8 1], 'k--', 'LineWidth', 1);  

set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlim([10^-7, 10^-1])
ylim([10^-7, 10^-1])

xlabel('$v_i^\prime$', 'Interpreter', 'latex', 'FontWeight', 'bold', 'FontSize', 14);
ylabel('$v_i$', 'Interpreter', 'latex', 'FontWeight', 'bold', 'FontSize', 14);

% Legend labels
legendLabels = {'$v_1$', '$v_2$'};
% Create the legend
legend([scatter_handle01,scatter_handle02], legendLabels, 'Location', 'southeast', 'Interpreter', 'latex','FontSize', 11);

box on
set(gca, 'LineWidth', 0.8);
text(0.1, 0.85, 'A', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 20);
text(0.1, 0.6, '$\frac{v_{2}}{v_{1}} \leq 100$', 'Interpreter', 'latex', 'Units', 'normalized', 'FontSize', 12);
pbaspect([1 1 1])

%%% Equal weights
subplot('Position', subplotPositions(2, :));
hold on

% Initialize sums and counts for averaging
sumSquaredDiff_v1 = 0;
count_v1 = 0;
sumSquaredDiff_v2 = 0;
count_v2 = 0;

NaNcount = 0;
[ID_1, ID_2] = meshgrid(1:numel(Parameter.a), 1:numel(Parameter.a));
for i = shortened_vector
    y_axis = WFsim.v1_2{i};
    for j = 1:numel(ID_1)
        s1(j,1) = Parameter.s1U1{i,ID_1(j)}(1);
        v1(j,1) = Parameter.v1_check(i,ID_1(j))';
    end
    v1_min = (2.*s1)./(WFsim.generations-WFsim.t_burnin+1);
    
    v1(v1_min>WFsim.v1_2{i}) = NaN;
    y_axis(v1_min>WFsim.v1_2{i}) = NaN;
    
    NaNcount = NaNcount + sum(isnan(Prediction.v1_2(i,1:end)'));
    predict = Prediction.v1_2(i,1:end)';
    predict(v1_min>WFsim.v1_2{i}) = NaN;
    scatter_handle01 = scatter(predict, y_axis, 'Marker', 'x', 'MarkerEdgeColor', [color 0 0], 'MarkerFaceColor', [color 0 0], 'MarkerEdgeAlpha', alpha_val, 'MarkerFaceAlpha', 0, 'SizeData', 10);

    % Calculate and accumulate log squared differences for v1
    validIndices_v1 = ~isnan(predict) & ~isnan(y_axis); % Indices of valid (non-NaN) points
    logSquaredDiff_v1 = (log(predict(validIndices_v1)) - log(y_axis(validIndices_v1))).^2;
    sumSquaredDiff_v1 = sumSquaredDiff_v1 + sum(logSquaredDiff_v1);
    count_v1 = count_v1 + numel(logSquaredDiff_v1);
    
    y_axis = WFsim.v2_1{i};
    for j = 1:numel(ID_1)
        s2(j,1) = Parameter.s2U2{i,ID_2(j)}(1);
    end
    v2_min = (2.*s2)./(WFsim.generations-WFsim.t_burnin+1);
    
    y_axis(v2_min>WFsim.v2_1{i}) = NaN;
    
    predict = Prediction.v2_1(i,1:end)';
    predict(v2_min>WFsim.v2_1{i}) = NaN;
    scatter_handle02 = scatter(predict, y_axis, 'Marker', '+', 'MarkerEdgeColor', [0 0 color], 'MarkerFaceColor', [0 0 color], 'MarkerEdgeAlpha', alpha_val, 'MarkerFaceAlpha', alpha_val, 'SizeData', 10);
    
    % Calculate and accumulate log squared differences for v2
    validIndices_v2 = ~isnan(predict) & ~isnan(y_axis); % Indices of valid (non-NaN) points
    logSquaredDiff_v2 = (log(predict(validIndices_v2)) - log(y_axis(validIndices_v2))).^2;
    sumSquaredDiff_v2 = sumSquaredDiff_v2 + sum(logSquaredDiff_v2);
    count_v2 = count_v2 + numel(logSquaredDiff_v2);
end

Prediction.NaNcount = NaNcount;

% Compute the average log squared differences
avgLogSquaredDiff_v1 = sumSquaredDiff_v1 / count_v1;
avgLogSquaredDiff_v2 = sumSquaredDiff_v2 / count_v2;

v_star_Array(2,:) = [avgLogSquaredDiff_v1, avgLogSquaredDiff_v2]; 

plot([10^-8 1], [10^-8 1], 'k--', 'LineWidth', 1);  

set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlim([10^-7, 10^-1])
ylim([10^-7, 10^-1])

xlabel('$v_i^*, w_i = \frac{1}{2}$', 'Interpreter', 'latex', 'FontWeight', 'bold', 'FontSize', 14);
ylabel('$v_i$', 'Interpreter', 'latex', 'FontWeight', 'bold', 'FontSize', 14);

% Legend labels
legendLabels = {'$v_1$', '$v_2$'};
% Create the legend
legend([scatter_handle01,scatter_handle02], legendLabels, 'Location', 'southeast', 'Interpreter', 'latex','FontSize', 11);

box on
set(gca, 'LineWidth', 0.8);
text(0.1, 0.85, 'B', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 20);
text(0.1, 0.6, '$\frac{v_{2}}{v_{1}} \leq 100$', 'Interpreter', 'latex', 'Units', 'normalized', 'FontSize', 12);
pbaspect([1 1 1])

%%% U weights
subplot('Position', subplotPositions(3, :));
hold on

% Initialize sums and counts for averaging
sumSquaredDiff_v1 = 0;
count_v1 = 0;
sumSquaredDiff_v2 = 0;
count_v2 = 0;

NaNcount = 0;
[ID_1, ID_2] = meshgrid(1:numel(Parameter.a), 1:numel(Parameter.a));
for i = shortened_vector
    y_axis = WFsim.v1_2{i};
    for j = 1:numel(ID_1)
        s1(j,1) = Parameter.s1U1{i,ID_1(j)}(1);
        v1(j,1) = Parameter.v1_check(i,ID_1(j))';
    end
    v1_min = (2.*s1)./(WFsim.generations-WFsim.t_burnin+1);
    
    v1(v1_min>WFsim.v1_2{i}) = NaN;
    y_axis(v1_min>WFsim.v1_2{i}) = NaN;
        
    NaNcount = NaNcount + sum(isnan(Prediction_U_weights.v1_2(i,1:end)'));
    predict = Prediction_U_weights.v1_2(i,1:end)';
    predict(v1_min>WFsim.v1_2{i}) = NaN;
    scatter_handle01 = scatter(predict, y_axis, 'Marker', 'x', 'MarkerEdgeColor', [color 0 0], 'MarkerFaceColor', [color 0 0], 'MarkerEdgeAlpha', alpha_val, 'MarkerFaceAlpha', 0, 'SizeData', 10);
    
    % Calculate and accumulate log squared differences for v1
    validIndices_v1 = ~isnan(predict) & ~isnan(y_axis); % Indices of valid (non-NaN) points
    logSquaredDiff_v1 = (log(predict(validIndices_v1)) - log(y_axis(validIndices_v1))).^2;
    sumSquaredDiff_v1 = sumSquaredDiff_v1 + sum(logSquaredDiff_v1);
    count_v1 = count_v1 + numel(logSquaredDiff_v1);
    
    y_axis = WFsim.v2_1{i};
    for j = 1:numel(ID_1)
        s2(j,1) = Parameter.s2U2{i,ID_2(j)}(1);
    end
    v2_min = (2.*s2)./(WFsim.generations-WFsim.t_burnin+1);
    
    y_axis(v2_min>WFsim.v2_1{i}) = NaN;
    
    predict = Prediction_U_weights.v2_1(i,1:end)';    
    predict(v2_min>WFsim.v2_1{i}) = NaN;
    scatter_handle02 = scatter(predict, y_axis, 'Marker', '+', 'MarkerEdgeColor', [0 0 color], 'MarkerFaceColor', [0 0 color], 'MarkerEdgeAlpha', alpha_val, 'MarkerFaceAlpha', alpha_val, 'SizeData', 10);
    
    % Calculate and accumulate log squared differences for v2
    validIndices_v2 = ~isnan(predict) & ~isnan(y_axis); % Indices of valid (non-NaN) points
    logSquaredDiff_v2 = (log(predict(validIndices_v2)) - log(y_axis(validIndices_v2))).^2;
    sumSquaredDiff_v2 = sumSquaredDiff_v2 + sum(logSquaredDiff_v2);
    count_v2 = count_v2 + numel(logSquaredDiff_v2);
end

Prediction_U_weights.NaNcount = NaNcount;

% Compute the average log squared differences
avgLogSquaredDiff_v1 = sumSquaredDiff_v1 / count_v1;
avgLogSquaredDiff_v2 = sumSquaredDiff_v2 / count_v2;

v_star_Array(3,:) = [avgLogSquaredDiff_v1, avgLogSquaredDiff_v2]; 

plot([10^-8 1], [10^-8 1], 'k--', 'LineWidth', 1);  

set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlim([10^-7, 10^-1])
ylim([10^-7, 10^-1])

xlabel('$v_i^*, w_i = \frac{U_i}{U}$', 'Interpreter', 'latex', 'FontWeight', 'bold', 'FontSize', 14);
ylabel('$v_i$', 'Interpreter', 'latex', 'FontWeight', 'bold', 'FontSize', 14);

% Legend labels
legendLabels = {'$v_1$', '$v_2$'};
% Create the legend
legend([scatter_handle01,scatter_handle02], legendLabels, 'Location', 'southeast', 'Interpreter', 'latex','FontSize', 11);

box on
set(gca, 'LineWidth', 0.8);
text(0.1, 0.85, 'C', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 20);
text(0.1, 0.6, '$\frac{v_{2}}{v_{1}} \leq 100$', 'Interpreter', 'latex', 'Units', 'normalized', 'FontSize', 12);
pbaspect([1 1 1]) 

%%% s weights
subplot('Position', subplotPositions(4, :));
hold on

% Initialize sums and counts for averaging
sumSquaredDiff_v1 = 0;
count_v1 = 0;
sumSquaredDiff_v2 = 0;
count_v2 = 0;

NaNcount = 0;
[ID_1, ID_2] = meshgrid(1:numel(Parameter.a), 1:numel(Parameter.a));
for i = shortened_vector
    y_axis = WFsim.v1_2{i};
    for j = 1:numel(ID_1)
        s1(j,1) = Parameter.s1U1{i,ID_1(j)}(1);
        v1(j,1) = Parameter.v1_check(i,ID_1(j))';
    end
    v1_min = (2.*s1)./(WFsim.generations-WFsim.t_burnin+1);
    
    v1(v1_min>WFsim.v1_2{i}) = NaN;
    y_axis(v1_min>WFsim.v1_2{i}) = NaN;
    
    NaNcount = NaNcount + sum(isnan(Prediction_s_weights.v1_2(i,1:end)'));
    predict = Prediction_s_weights.v1_2(i,1:end)';
    predict(v1_min>WFsim.v1_2{i}) = NaN;
    scatter_handle01 = scatter(predict, y_axis, 'Marker', 'x', 'MarkerEdgeColor', [color 0 0], 'MarkerFaceColor', [color 0 0], 'MarkerEdgeAlpha', alpha_val, 'MarkerFaceAlpha', 0, 'SizeData', 10);
    
    % Calculate and accumulate log squared differences for v1
    validIndices_v1 = ~isnan(predict) & ~isnan(y_axis); % Indices of valid (non-NaN) points
    logSquaredDiff_v1 = (log(predict(validIndices_v1)) - log(y_axis(validIndices_v1))).^2;
    sumSquaredDiff_v1 = sumSquaredDiff_v1 + sum(logSquaredDiff_v1);
    count_v1 = count_v1 + numel(logSquaredDiff_v1);
    
    y_axis = WFsim.v2_1{i};
    for j = 1:numel(ID_1)
        s2(j,1) = Parameter.s2U2{i,ID_2(j)}(1);
    end
    v2_min = (2.*s2)./(WFsim.generations-WFsim.t_burnin+1);
    
    y_axis(v2_min>WFsim.v2_1{i}) = NaN;
    
    predict = Prediction_s_weights.v2_1(i,1:end)';    
    predict(v2_min>WFsim.v2_1{i}) = NaN;
    scatter_handle02 = scatter(predict, y_axis, 'Marker', '+', 'MarkerEdgeColor', [0 0 color], 'MarkerFaceColor', [0 0 color], 'MarkerEdgeAlpha', alpha_val, 'MarkerFaceAlpha', alpha_val, 'SizeData', 10);
    
    % Calculate and accumulate log squared differences for v2
    validIndices_v2 = ~isnan(predict) & ~isnan(y_axis); % Indices of valid (non-NaN) points
    logSquaredDiff_v2 = (log(predict(validIndices_v2)) - log(y_axis(validIndices_v2))).^2;
    sumSquaredDiff_v2 = sumSquaredDiff_v2 + sum(logSquaredDiff_v2);
    count_v2 = count_v2 + numel(logSquaredDiff_v2);
end

Prediction_s_weights.NaNcount = NaNcount;

% Compute the average log squared differences
avgLogSquaredDiff_v1 = sumSquaredDiff_v1 / count_v1;
avgLogSquaredDiff_v2 = sumSquaredDiff_v2 / count_v2;

v_star_Array(4,:) = [avgLogSquaredDiff_v1, avgLogSquaredDiff_v2]; 

plot([10^-8 1], [10^-8 1], 'k--', 'LineWidth', 1);  

set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlim([10^-7, 10^-1])
ylim([10^-7, 10^-1])

xlabel('$v_i^*, w_i = \frac{s_i}{s_1+s_2}$', 'Interpreter', 'latex', 'FontWeight', 'bold', 'FontSize', 14);
ylabel('$v_i$', 'Interpreter', 'latex', 'FontWeight', 'bold', 'FontSize', 14);

% Legend labels
legendLabels = {'$v_1$', '$v_2$'};
% Create the legend
legend([scatter_handle01,scatter_handle02], legendLabels, 'Location', 'southeast', 'Interpreter', 'latex','FontSize', 11);

box on
set(gca, 'LineWidth', 0.8);
text(0.1, 0.85, 'D', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 20);
text(0.1, 0.6, '$\frac{v_{2}}{v_{1}} \leq 100$', 'Interpreter', 'latex', 'Units', 'normalized', 'FontSize', 12);
pbaspect([1 1 1]) 

print('SI5', '-depsc', '-r300');

%% Add metainfo 
clf;
color = 0.4;

subplotPositions = [
    0.15, 0.2, 0.3, 0.5;  
    0.55, 0.2, 0.3, 0.5; 
];

subplot('Position', subplotPositions(1, :));
hold on

box on
set(gca, 'LineWidth', 0.8);
plt1 = plot(v_star_Array(:,1),'x','color',[color 0 0]);
text(-0.15, 1, 'A', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 14, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
plt2 = plot(v_star_Array(:,2),'+','color',[0 0 color]);

xlim([0.5,4.5])
ylim([0.05,1.5])

legendLabels = {'Module $1$', 'Module $2$'};
legend([plt1,plt2], legendLabels, 'Location', 'northeast', 'Interpreter', 'latex','FontSize', 11);

xticks(1:4); 
xticklabels({'A', 'B', 'C', 'D'});
xlabel('Subplot ID', 'Interpreter', 'latex', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('Mean Squared Log Distance', 'Interpreter', 'latex', 'FontWeight', 'bold', 'FontSize', 12);

pbaspect([1 1 1]) 

% NaNcount
subplot('Position', subplotPositions(2, :));
hold on

box on
set(gca, 'LineWidth', 0.8);

NaNcount_Array = [Prediction.NaNcount,Prediction_U_weights.NaNcount,Prediction_s_weights.NaNcount];

plt1 = plot(2:4, 100 * NaNcount_Array(:)./numel(Prediction.v1_2),'ko-');
text(-0.15, 1, 'B', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 14, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
xlim([0.5,4.5])

xticks(1:4); 
xticklabels({'A', 'B', 'C', 'D'}); 
xlabel('Subplot ID', 'Interpreter', 'latex', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('\% Non-unique Solutions', 'Interpreter', 'latex', 'FontWeight', 'bold', 'FontSize', 12);

pbaspect([1 1 1]) 

text(0.1, 0.15, '$\textbf{tol} = 10$', 'Interpreter', 'latex', 'Units', 'normalized', 'FontSize', 14);

print('SI6', '-depsc', '-r300');

%% v' and v ratios %%
clf
alpha_val = 0.2;
color_palette = [230, 97, 0; ... % Dark orange
                 139, 0, 139;  ... % Dark magenta
                 75, 83, 32; ... % Dark green
                 26, 133, 255;    ... % Blue
                 0, 0, 0] / 255;   % Black
             
% Mapping function to assign color index based on i
map_index_to_color = @(i) (any([1,2,7]==i)*1 + any([3,8,12]==i)*2 + ...
    any([4,9,13,16]==i)*3 + any([5,10,14,17,19]==i)*4 + any([6,11,15,18,20,21]==i)*5);
          
hold on
[ID_1, ID_2] = meshgrid(1:numel(Parameter.a), 1:numel(Parameter.a));
for i = 1:21
    v1_2 = WFsim.v1_2{i};
    for j = 1:numel(ID_1)
        s1(j,1) = Parameter.s1U1{i,ID_1(j)}(1);
        v1(j,1) = Parameter.v1_check(i,ID_1(j))';
    end
    v1_min = (2.*s1)./(WFsim.generations-WFsim.t_burnin+1);
    v1(v1_min>WFsim.v1_2{i}) = NaN;
    v1_2(v1_min>WFsim.v1_2{i}) = NaN;
    
    v2_1 = WFsim.v2_1{i};
    for j = 1:numel(ID_1)
        s2(j,1) = Parameter.s2U2{i,ID_2(j)}(1);
        v2(j,1) = Parameter.v2_check(i,ID_2(j))';
    end
    v2_min = (2.*s2)./(WFsim.generations-WFsim.t_burnin+1);
    v2(v2_min>WFsim.v2_1{i}) = NaN;
    v2_1(v2_min>WFsim.v2_1{i}) = NaN;
   
    color_idx = map_index_to_color(i);
    
    scatter_handle = scatter(v2./v1, v2_1./v1_2, 'Marker', 'o', 'MarkerEdgeColor', color_palette(color_idx, :), 'MarkerFaceColor', color_palette(color_idx, :), 'MarkerEdgeAlpha', alpha_val*2, 'MarkerFaceAlpha', alpha_val, 'SizeData', 12);
end
plot([10^-2, 10^6], [10^-2, 10^6], 'k--', 'LineWidth', 1); 

set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlim([10^-2, 10^5])
ylim([10^-2, 10^5])

xlabel('$v_2^\prime/v_1^\prime$', 'Interpreter', 'latex', 'FontWeight', 'bold', 'FontSize', 14);
ylabel('$v_2/v_1$', 'Interpreter', 'latex', 'FontWeight', 'bold', 'FontSize', 14);

box on
set(gca, 'LineWidth', 0.8);
pbaspect([1 1 1]) 

print('SI4', '-depsc', '-r300');

%% s and U range
clearvars;
clc;
load('Stdstate_N_10000.mat')

clf
alpha_val = 0.4;

color_palette = [230, 97, 0; ... % Dark orange
                 139, 0, 139;  ... % Dark magenta
                 75, 83, 32; ... % Dark green
                 26, 133, 255;    ... % Blue
                 0, 0, 0] / 255;   % Black
             
% Mapping function to assign color index based on i
map_index_to_color = @(i) (any([1,2,7]==i)*1 + any([3,8,12]==i)*2 + ...
    any([4,9,13,16]==i)*3 + any([5,10,14,17,19]==i)*4 + any([6,11,15,18,20,21]==i)*5);

%%% 1st subplot %%%
subplot(1, 2, 1);
hold on

[ID_1, ID_2] = meshgrid(1:numel(Parameter.a), 1:numel(Parameter.a));
for i = 1:21
    v1_2 = WFsim.v1_2{i};
    for j = 1:numel(ID_1)
        s1(j,1) = Parameter.s1U1{i,ID_1(j)}(1);
        v1(j,1) = Parameter.v1_check(i,ID_1(j))';
    end
    v1_min = (2.*s1)./(WFsim.generations-WFsim.t_burnin+1);
    v1(v1_min>WFsim.v1_2{i}) = NaN;
    v1_2(v1_min>WFsim.v1_2{i}) = NaN;
    
    v2_1 = WFsim.v2_1{i};
    for j = 1:numel(ID_1)
        s2(j,1) = Parameter.s2U2{i,ID_2(j)}(1);
        v2(j,1) = Parameter.v2_check(i,ID_2(j))';
    end
    v2_min = (2.*s2)./(WFsim.generations-WFsim.t_burnin+1);
    v2(v2_min>WFsim.v2_1{i}) = NaN;
    v2_1(v2_min>WFsim.v2_1{i}) = NaN;
   
    % Determine color for this loop iteration
    color_idx = map_index_to_color(i);
    
    scatter_handle = scatter(s1, s2, 'Marker', 'o', 'MarkerEdgeColor', color_palette(color_idx, :), 'MarkerFaceColor', color_palette(color_idx, :), 'MarkerEdgeAlpha', alpha_val, 'MarkerFaceAlpha', alpha_val, 'SizeData', 25);
end
plot([10^-3, 4*10^-1], [10^-3, 4*10^-1], 'k--', 'LineWidth', 1); 

set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlim([10^-3, 4*10^-1])
ylim([10^-3, 4*10^-1])

xlabel('$s_{1}$', 'Interpreter', 'latex', 'FontWeight', 'bold', 'FontSize', 14);
ylabel('$s_{2}$', 'Interpreter', 'latex', 'FontWeight', 'bold', 'FontSize', 14);

box on
set(gca, 'LineWidth', 0.8);
text(0.8, 0.2, 'A', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 20);

pbaspect([1 1 1]) 

%%% 2nd subplot %%%
subplot(1, 2, 2);
hold on

[ID_1, ID_2] = meshgrid(1:numel(Parameter.a), 1:numel(Parameter.a));
for i = 1:21
    v1_2 = WFsim.v1_2{i};
    for j = 1:numel(ID_1)
        s1(j,1) = Parameter.s1U1{i,ID_1(j)}(1);
        U1(j,1) = Parameter.s1U1{i,ID_1(j)}(2);
        v1(j,1) = Parameter.v1_check(i,ID_1(j))';
    end
    v1_min = (2.*s1)./(WFsim.generations-WFsim.t_burnin+1);
    v1(v1_min>WFsim.v1_2{i}) = NaN;
    v1_2(v1_min>WFsim.v1_2{i}) = NaN;
    
    v2_1 = WFsim.v2_1{i};
    for j = 1:numel(ID_1)
        s2(j,1) = Parameter.s2U2{i,ID_2(j)}(1);
        U2(j,1) = Parameter.s2U2{i,ID_2(j)}(2);
        v2(j,1) = Parameter.v2_check(i,ID_2(j))';
    end
    v2_min = (2.*s2)./(WFsim.generations-WFsim.t_burnin+1);
    v2(v2_min>WFsim.v2_1{i}) = NaN;
    v2_1(v2_min>WFsim.v2_1{i}) = NaN;
   
    % Determine color for this loop iteration
    color_idx = map_index_to_color(i);
    
    scatter_handle = scatter(U1, U2, 'Marker', 'o', 'MarkerEdgeColor', color_palette(color_idx, :), 'MarkerFaceColor', color_palette(color_idx, :), 'MarkerEdgeAlpha', alpha_val, 'MarkerFaceAlpha', alpha_val, 'SizeData', 25);
end
plot([10^-5, 2*10^-2], [10^-5, 2*10^-2], 'k--', 'LineWidth', 1); 

set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlim([2*10^-5, 2*10^-2])
ylim([2*10^-5, 2*10^-2])

xlabel('$U_{1}$', 'Interpreter', 'latex', 'FontWeight', 'bold', 'FontSize', 14);
ylabel('$U_{2}$', 'Interpreter', 'latex', 'FontWeight', 'bold', 'FontSize', 14);

box on
set(gca, 'LineWidth', 0.8);
text(0.8, 0.2, 'B', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 20);

pbaspect([1 1 1]) 

print('SI2', '-depsc', '-r300');

%% F trajectories: Steady-state check

clf

%%% 1st subplot %%%
subplot(1, 2, 1);
hold on

alph = 0.05;
for i = 1:numel(WFsim.F1_2)
    Mat = WFsim.F1_2{i};
    for j = 1:numel(Mat(1:end,1))
        plt1 = plot(Mat(j, 1:end), 'k');%,'k','Color',[color color color]);
        plt1.Color(4) = alph;
    end
end
xline(10000,'k--');
set(gca, 'YScale', 'log')
hold off
xlabel('$t$, generations', 'Interpreter', 'latex')
ylabel('$F_{1}$, log fitness', 'Interpreter', 'latex')
text(0.1, 0.9, 'A', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 20);
ylim([10^-8, 10^4])

box on
pbaspect([1 1 1]) 

%%% 2nd subplot %%%
subplot(1, 2, 2);
hold on

for i = 1:numel(WFsim.F2_1)
    Mat = WFsim.F2_1{i};
    for j = 1:numel(Mat(1:end,1))
        plt2 = plot(Mat(j, 1:end), 'k');
        plt2.Color(4) = alph;
    end
end
xline(10000,'k--');
set(gca, 'YScale', 'log')
hold off
xlabel('$t$, generations', 'Interpreter', 'latex')
ylabel('$F_{2}$, log fitness', 'Interpreter', 'latex')
text(0.1, 0.9, 'B', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 20);
ylim([10^-8, 10^4])

box on
pbaspect([1 1 1]) 
print('SI3', '-depsc', '-r300');
