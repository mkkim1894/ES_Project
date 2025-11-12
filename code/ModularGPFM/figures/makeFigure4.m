function makeFigure4(mutationRate, varargin)
% makeFigure4 — Generate Figure 4: Recombination spectrum across regimes
%
% Usage:
%   makeFigure4(2e-4)             % auto-detect full/demo (prefers full)
%   makeFigure4(2e-4, 'demo')     % force demo mode
%   makeFigure4(2e-4, 'full')     % force full mode
%
% Inputs:
%   mutationRate : numeric (e.g., 2e-4)
%   mode         : optional ('demo', 'full', 'auto'); default = 'auto'
%
% Outputs:
%   Publication-quality EPS figure saved to:
%       results/Figures/Figure_ModularFGM_Recombination.eps
%       or results_demo/Figures/... if in demo mode
%
% Reference:
%   Kim, M., Ardell, S. M., & Kryazhimskiy, S. (2025).
%   "Module-Selection Balance in the Evolution of Modular Organisms."

% ---------------------------------------------------------------------
% Parse optional mode argument
% ---------------------------------------------------------------------
if nargin < 2
    mode = 'auto';
else
    mode = lower(varargin{1});
    if ~ismember(mode, {'demo', 'full', 'auto'})
        error('Invalid mode: use ''demo'', ''full'', or ''auto''.');
    end
end

% ---------------------------------------------------------------------
% Resolve results directory robustly
% ---------------------------------------------------------------------
figureRoot = fileparts(mfilename('fullpath'));

switch mode
    case 'demo'
        resultsRoot = fullfile(figureRoot, '../results_demo');
    case 'full'
        resultsRoot = fullfile(figureRoot, '../results');
    otherwise  % 'auto' — prefer full results for publication
        if isfolder(fullfile(figureRoot, '../results'))
            resultsRoot = fullfile(figureRoot, '../results');
        elseif isfolder(fullfile(figureRoot, '../results_demo'))
            resultsRoot = fullfile(figureRoot, '../results_demo');
        else
            error('Neither results/ nor results_demo/ found.');
        end
end

recombDir = fullfile(resultsRoot, 'CM_IntermediateRecomb');
if ~isfolder(recombDir)
    error('Could not find CM_IntermediateRecomb in %s', resultsRoot);
end
fprintf('makeFigure4 loading data from: %s\n', recombDir);

% ---------------------------------------------------------------------
% Initialize figure parameters
% ---------------------------------------------------------------------
markers = {'d', 'o'};
customcolor = [0.4, 0.4, 0.4];
figureWidth = 8;     % cm
figureHeight = 6;    % cm
figure('Units', 'centimeters', 'Position', [1, 1, figureWidth, figureHeight]);

% ---------------------------------------------------------------------
% Load recombination files dynamically
% ---------------------------------------------------------------------
files = dir(sprintf('%s/ModularFGM_CM_Recomb_R*_M%.1e.mat', recombDir, mutationRate));
if isempty(files)
    error('No recombination files found for mutation rate %.1e in %s', mutationRate, recombDir);
end

R_list = zeros(1, numel(files));
logRatio = [];
for r_i = 1:numel(files)
    data = load(fullfile(files(r_i).folder, files(r_i).name));
    simParams = data.simParams;
    if isfield(data, 'averageTrajectory')
        ave = data.averageTrajectory;
    elseif isfield(data, 'ave')
        ave = data.ave;
    else
        error('Missing averageTrajectory/ave in %s', files(r_i).name);
    end
    R_list(r_i) = simParams.recombinationRate;
    for ind = 1:length(simParams.initialAngles)
        logRatio_tmp(ind,:) = ave.logRatio{ind};
        logRatio(r_i,ind) = logRatio_tmp(ind,end);
    end
end

% ---------------------------------------------------------------------
% Plot recombination results
% ---------------------------------------------------------------------
% Baseline and parameter setup
epsilon = 5e-5;
pseudo_log_R_list = log10(R_list + epsilon);
[pseudo_log_R_list_sorted, sortIdx] = sort(pseudo_log_R_list);
R_list_sorted = R_list(sortIdx);

Initial_ln_R = log(tan(simParams.initialAngles));
y_min = min([Initial_ln_R, logRatio(:)']) - 0.5;
y_max = max([Initial_ln_R, logRatio(:)']) + 0.5;
ylim([y_min, y_max]);

% Reference line: ESBL
yline(log(simParams.ellipseParams(2)^2 / simParams.ellipseParams(1)^2), ...
    '-', 'LineWidth', 2, 'Color', [0.8 0.3 0]);

% Plot results per initial condition
for j_dep = 1:length(simParams.initialAngles)
    Dep_Variable = logRatio(:, j_dep);
    hold on
    plot(pseudo_log_R_list_sorted, Dep_Variable(sortIdx), '-', 'Color', '#2776A8', 'Marker', markers{j_dep}, ...
        'MarkerFaceColor', 'w', 'LineWidth', 1, 'MarkerSize', 3);
    plot(pseudo_log_R_list_sorted, repmat(Initial_ln_R(j_dep), size(pseudo_log_R_list_sorted)), ...
        'Color','#2776A8', 'Marker', markers{j_dep}, 'MarkerFaceColor', '#2776A8', 'LineStyle', 'none', ...
        'LineWidth', 1.5, 'MarkerSize', 3);

% ---------------------------------------------------------------------
% Draw vertical arrows from initial to final ratios (with clean arrowheads)
% ---------------------------------------------------------------------
for idx = 1:length(pseudo_log_R_list_sorted)
    x_data = pseudo_log_R_list_sorted(idx);
    y_start = Initial_ln_R(j_dep);
    y_end   = Dep_Variable(sortIdx(idx));

    % Scale and center arrow shaft
    arrow_fraction = 0.8;
    delta_y = y_end - y_start;
    delta_y_adj = delta_y * arrow_fraction;
    y_offset = (delta_y - delta_y_adj) / 2;
    y_start_adj = y_start + y_offset;
    y_end_adj   = y_start_adj + delta_y_adj;

    % Arrow parameters
    arrowhead_length = 0.12;
    arrowhead_width  = 0.085;
    arrow_direction  = sign(delta_y_adj);
    color_body = [0.15, 0.46, 0.66];  % blue tone (#2776A8)

    % --- Draw arrow line slightly shorter so arrowhead protrudes clearly ---
    tip_extension = 0.5 * arrowhead_length * arrow_direction;
    line([x_data, x_data], [y_start_adj, y_end_adj - tip_extension], ...
        'Color', color_body, 'LineWidth', 1);

    % --- Define arrowhead vertices (tip slightly beyond line end) ---
    x_vertices = [x_data, ...
                  x_data - arrowhead_width/2, ...
                  x_data + arrowhead_width/2];
    y_vertices = [y_end_adj + arrow_direction * (arrowhead_length/3), ...
                  y_end_adj - arrow_direction * arrowhead_length, ...
                  y_end_adj - arrow_direction * arrowhead_length];

    % --- Draw filled arrowhead ---
    fill(x_vertices, y_vertices, color_body, 'EdgeColor', 'none');
end

end

% ---------------------------------------------------------------------
% Axes and labels
% ---------------------------------------------------------------------
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 8);
xlabel('$\rho$', 'Interpreter', 'latex', 'FontWeight', 'bold', 'FontSize', 10);
ylabel('$\log(x_2/x_1)$', 'Interpreter', 'latex', 'FontSize', 10);
set(gca, 'XTick', pseudo_log_R_list_sorted);
set(gca, 'XTickLabel', arrayfun(@(x) sprintf('%.2g', x), R_list_sorted, 'UniformOutput', false));
xlim([min(pseudo_log_R_list_sorted) - 0.1, max(pseudo_log_R_list_sorted) + 0.1]);
grid off

% ---------------------------------------------------------------------
% Save figure (EPS, publication-quality)
% ---------------------------------------------------------------------
figDir = fullfile(resultsRoot, 'Figures');
if ~isfolder(figDir), mkdir(figDir); end
outPath = fullfile(figDir, sprintf('Figure_ModularFGM_Recombination_M%.1e.eps', mutationRate));
saveFigure(outPath, figureWidth, figureHeight);
fprintf('Saved %s\n', outPath);
end

% ======================================================================
% Helper function
% ======================================================================
function saveFigure(file_name, width, height)
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0, 0, width, height]);
set(gcf, 'PaperSize', [width, height]);
print(file_name, '-depsc', '-r300');
end
