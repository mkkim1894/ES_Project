function makeFigureS_NestedFGMDistributions(dataFile, varargin)
% makeFigureS_NestedFGMDistributions - Supplementary figure for nested FGM mutation distributions.
%
% Description:
%   Panel A: Mean trajectory in (x1,x2) space for a selected initial condition
%   from the SSWM nested FGM simulation, with three beads marking positions
%   at 5%, 50%, and 95% of arc length along the trajectory.
%   Panels B: Distributions of phenotypic effects Delta x_i ~
%   N(-m^2/2, (m*sqrt(2|x_i|/n_i))^2) for both modules at each bead position.
%
% Inputs:
%   dataFile - Path to NestedFGM_SSWM_*.mat file from Run_nestedFGM
%   varargin - Optional name-value pairs:
%       'outputDir' - Directory for output figure (default: current directory)
%       'angleIdx'  - Which initial angle to show (default: 1)
%
% Outputs:
%   Saves FigureS_NestedFGMDistributions.pdf to outputDir
%
% Reference:
%   Kim, M., Ardell, S. M., & Kryazhimskiy, S. (2025).
%   "Module-Selection Balance in the Evolution of Modular Organisms."
%
% See also: simulateNestedSSWM, Run_nestedFGM
%
% Copyright (c) 2025 Minkyu Kim, Cornell University
% Licensed under MIT License

%% Parse inputs
p = inputParser;
addRequired(p, 'dataFile', @ischar);
addParameter(p, 'outputDir', '.', @ischar);
addParameter(p, 'angleIdx', 1, @isnumeric);   % which initial angle to show
parse(p, dataFile, varargin{:});

outputDir  = p.Results.outputDir;
angleIdx   = p.Results.angleIdx;

if ~isfolder(outputDir)
    mkdir(outputDir);
end

%% Load data
fprintf('Loading %s ...\n', dataFile);
d  = load(dataFile);

% Validate file contents
if ~isfield(d, 'simParams')
    error('makeFigureS_NestedFGMDistributions: file does not contain simParams.\nExpected a NestedFGM_SSWM_*.mat file from Run_nestedFGM.');
end
if ~isfield(d.simParams, 'moduleDimension')
    error('makeFigureS_NestedFGMDistributions: simParams.moduleDimension not found.\nExpected a nested FGM file, not a modular/pleiotropic one.');
end
if ~isfield(d, 'averageTrajectory') && ~isfield(d, 'ave')
    error('makeFigureS_NestedFGMDistributions: no trajectory data found.\nExpected averageTrajectory or ave field.');
end

sp = d.simParams;
av = getAverageTrajectory(d);

m  = sp.deltaTrait;          % mutation-vector magnitude (= sqrt(2*delta))
n1 = sp.moduleDimension(1);
n2 = sp.moduleDimension(2);

%% Extract mean trajectory for selected angle
ts  = av.averageTimeStamp{angleIdx};   % [2 x T]  rows = x1, x2
x1t = ts(1, :);
x2t = ts(2, :);
T   = size(ts, 2);

% Three time stamps sampled at 5%, 50%, and 95% of arc length
stampFracs = [0.05, 0.5, 0.95];
stampLabels = {'Early', 'Middle', 'Late'};
beadColor   = '#D95319';

%% Layout
% Figure: 17.8 x 11 cm
% Panel A : left half  [0.06  0.12  0.38  0.78]
% Panels B: right half, stacked 3 rows
%   B-top    [0.54  0.68  0.40  0.24]
%   B-middle [0.54  0.38  0.40  0.24]
%   B-bottom [0.54  0.08  0.40  0.24]

figW = 17.8;  figH = 11.0;
figure('Units','centimeters','Position',[1 1 figW figH]);

posA  = [0.06  0.12  0.38  0.78];
posB  = {[0.54  0.68  0.40  0.22], ...
         [0.54  0.40  0.40  0.22], ...
         [0.54  0.12  0.40  0.22]};

%% ---- Panel A: trajectory --------------------------------------------------
subplot('Position', posA);
hold on;

% Fitness contours (same as Fig 5)
[~, lvl] = defineGaussianPDF(sp);
plotPDFContour(sp, lvl, [0.7 0.7 0.7]);

% Mean trajectory
plot(x1t, x2t, '-', 'Color', '#2776A8', 'LineWidth', 1.0);

% Starting marker
scatter(x1t(1), x2t(1), 30, ...
    'MarkerEdgeColor','#2776A8','MarkerFaceColor','#2776A8');

% Time-stamp beads - sampled at 5/50/95% of arc length along trajectory
dx = diff(x1t); dy = diff(x2t);
arcLen = cumsum([0, sqrt(dx.^2 + dy.^2)]);
totalArc = arcLen(end);
stampIdx = arrayfun(@(f) find(arcLen >= f*totalArc, 1, 'first'), stampFracs);
stampIdx = max(1, min(T, stampIdx));

for k = 1:3
    si = stampIdx(k);
    scatter(x1t(si), x2t(si), 60, 'o', ...
        'MarkerEdgeColor', beadColor, ...
        'MarkerFaceColor', beadColor);
    % label centered below the bead
    text(x1t(si), x2t(si) - 0.12, stampLabels{k}, ...
        'FontSize', 7, 'Color', beadColor, 'FontName', 'Helvetica', ...
        'HorizontalAlignment', 'center');
end

customizeAxes(1.2.*[-2.8, 0.05], 1.2.*[-2.1875, 0.05]);
pbaspect([1 1 1]);
text(-2.5,   0.22,  'Module 1 Performance', 'FontName','Helvetica','FontSize',10);
text( 0.25, -0.3,  'Module 2 Performance', 'FontName','Helvetica','FontSize',10,'Rotation',270);
text( 0.08, 0.12,'0',                    'FontName','Helvetica','FontSize',10);

annotation('textbox',[posA(1)-0.04, posB{1}(2)+posB{1}(4)+0.05, 0.03, 0.03], ...
    'String','A','FontSize',12,'FontWeight','bold','EdgeColor','none');

%% ---- Panels B: mutation-effect distributions ------------------------------
delta_x = linspace(-0.6, 0.4, 400);
mu_eff  = -(m^2) / 2;                % shared mean for both modules

col1 = [0.2  0.2  0.2];   % dark grey  - module 1
col2 = [0.6  0.6  0.6];   % light grey - module 2

% Compute unified y-axis max across all three time stamps
ymax = 0;
for k = 1:3
    si = stampIdx(k);
    s1 = m * sqrt(max(0, 2*abs(x1t(si))/n1));
    s2 = m * sqrt(max(0, 2*abs(x2t(si))/n2));
    if s1 > 0, ymax = max(ymax, max(normpdf(delta_x, mu_eff, s1))); end
    if s2 > 0, ymax = max(ymax, max(normpdf(delta_x, mu_eff, s2))); end
end
ymax = ymax * 1.15;

for k = 1:3
    si  = stampIdx(k);
    x1k = x1t(si);
    x2k = x2t(si);

    sigma1 = m * sqrt(max(0, 2*abs(x1k)/n1));
    sigma2 = m * sqrt(max(0, 2*abs(x2k)/n2));

    subplot('Position', posB{k});
    hold on;

    if sigma1 > 0
        plot(delta_x, normpdf(delta_x, mu_eff, sigma1), '-', ...
            'Color', col1, 'LineWidth', 1.5, 'DisplayName', 'Module 1');
    end
    if sigma2 > 0
        plot(delta_x, normpdf(delta_x, mu_eff, sigma2), '--', ...
            'Color', col2, 'LineWidth', 1.5, 'DisplayName', 'Module 2');
    end

    % Mark zero (beneficial boundary)
    xline(0, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8, 'HandleVisibility', 'off');

    if k == 3
        xlabel('$\Delta x_i$', 'Interpreter','latex','FontSize',10);
    else
        set(gca, 'XTickLabel', {});
    end
    if k == 2
        ylabel('Probability density', 'FontSize', 9, 'FontName', 'Helvetica');
    end

    if k == 1
        legend('Interpreter','latex','FontSize',7,'Location','northwest','Box','off');
    end

    th = title(stampLabels{k}, 'FontName','Helvetica','FontSize',9,'FontWeight','normal');
    th.Units = 'normalized';
    th.Position(2) = th.Position(2) - 0.08;

    box off;
    set(gca, 'FontSize', 8, 'LineWidth', 0.8);
    xlim([-0.6, 0.4]);
    ylim([0, ymax]);

    % Only label first panel as B
    if k == 1
        annotation('textbox', ...
            [posB{k}(1)-0.04, posB{k}(2)+posB{k}(4)+0.05, 0.03, 0.03], ...
            'String', 'B', ...
            'FontSize',12,'FontWeight','bold','EdgeColor','none');
    end
end

%% Save
set(gcf,'Color','w');
outFile = fullfile(outputDir, 'FigureS_NestedFGMDistributions.pdf');
print(outFile, '-dpdf', '-vector');
fprintf('Saved to %s\n', outFile);
close(gcf);

end

%% ---- Helper functions -----------------------------------------------------

function av = getAverageTrajectory(d)
    if isfield(d,'averageTrajectory'), av = d.averageTrajectory;
    elseif isfield(d,'ave'),           av = d.ave;
    else, error('No trajectory data found.');
    end
end

function [f, lvl] = defineGaussianPDF(sp)
    syms x1 x2
    f = exp(-sqrt((x1./sp.ellipseParams(1)).^2 + ...
                  (x2./sp.ellipseParams(2)).^2).^2 / (2*sp.landscapeStdDev^2));
    dH  = det(hessian(f,[x1,x2]));
    dHf = matlabFunction(dH,'Vars',[x1,x2]);
    xyStar = fminsearch(@(x) abs(dHf(x(1),x(2))), [0,0]);
    lvl = double(subs(f,[x1,x2],xyStar));
end

function plotPDFContour(sp, lvl, lineColor)
    fc = fcontour(@(x1,x2) exp(-sqrt((x1./sp.ellipseParams(1)).^2 + ...
        (x2./sp.ellipseParams(2)).^2).^2 / (2*sp.landscapeStdDev^2)), ...
        'LineColor', lineColor, 'LineWidth', 1.4);
    fc.LevelList = [0.99, lvl];
end

function customizeAxes(xl, yl)
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    ax.Box = 'off';
    ax.XColor = 'k';
    ax.YColor = 'k';
    ax.LineWidth = 1;
    xlim(xl); ylim(yl);
    ax.XTick = [-2, -1];
    ax.YTick = [-2, -1];
    ax.TickLength = [0.015, 0.015];
    ax.FontSize = 8;
    ax.TickLabelInterpreter = 'latex';
end