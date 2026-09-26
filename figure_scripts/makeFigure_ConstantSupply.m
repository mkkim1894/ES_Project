function makeFigure_ConstantSupply(simParamsRef, varargin)
% makeFigure_ConstantSupply  Figure for the constant-supply modular GPFM.
%
% Layout is identical to makeFigure5_Generations (nested FGM):
%   A, B, C - phenotypic trajectories in the three evolutionary regimes
%   D, E, F - log(x_2/x_1) against generations in the same three regimes
%
% The dashed orange line marks the module-selection balance
% log[(L_2*a_2^2)/(L_1*a_1^2)] of the DECLINING-supply modular GPFM (Figure 3),
% i.e. the equal fitness benefits line. It is drawn as a reference, not as a
% prediction for this model.
%
% What the three regimes do differ in is whether R approaches anything at all.
% In the successive mutations regime and under free reassortment the two modules
% are uncoupled, log(x_2/x_1) drifts without bound, and the six predicted curves
% never converge. Under complete linkage they are coupled by clonal interference
% even though their mutational supplies are constant, and that coupling alone
% produces a fixed point: it sits where the two modules would adapt at the same
% rate in isolation, f_DF(s_1,U_1,N) = f_DF(s_2,U_2,N), which with U_1 =/= U_2 is
% NOT the line s_1 = s_2. So a balance exists in that regime, at a ratio below
% the dashed line, and the theory curves in panel E approach it from both sides.
%
% Inputs
%   simParamsRef  Struct used to disambiguate result files when several exist in
%                 a results directory. Pass [] to take the most recent.
%
% Name-value pairs
%   'resultsRoot'      Directory containing ConstantSupply/{SSWM, CM_Asexual,
%                      CM_Sexual}. Default: results/ inside this directory,
%                      falling back to results_test/.
%   'proximityCutoff'  Exclude ratio data once |x_i| < cutoff. Scalar (applied to
%                      all three panels) or a 3-vector, one per panel D, E, F.
%                      Default [0, 0, 0], i.e. off.
%
%                      Why it is off. An earlier version applied a cutoff of 0.25
%                      to the two CM panels, on the reasoning that the terminal
%                      part of a CM run is uninformative. That cutoff introduces a
%                      bias in the wrong direction. It drops a replicate as soon as
%                      EITHER |x_i| falls below the threshold, and module 2 is the
%                      one nearest its optimum, so the replicates dropped first are
%                      systematically the low-R ones and the surviving mean is
%                      pulled upward. In the linked panel this was enough to make
%                      the displayed mean plateau ABOVE the equal fitness benefits
%                      line while the trait-space panel showed the same
%                      trajectories crossing below it. Without the cutoff the two
%                      panels agree, and both agree with the theory curve: the
%                      simulated endpoints are log R = -1.15 ... -0.96 across the
%                      six initial conditions, all below log(a_2^2/a_1^2) = -0.693.
%
%                      The 80 percent retention rule in plotLogRatioPanel is what
%                      handles the late-time bias instead, and it does so without
%                      conditioning on trait values.
%   'yLimits'          y-limits of panels D-F. Default [-3, 2], as in Figures 3
%                      and 5.
%   'outputFile'       Output filename without extension. Default
%                      Figure_ModularConstantSupply_Generations.
%
% Outputs
%   Writes <resultsRoot>/Figures/<outputFile>.pdf
%
% Reference
%   Kim, M., Ardell, S. M., & Kryazhimskiy, S. (2025). Module-Selection
%   Balance in the Evolution of Modular Organisms.
%
% See also: makeFigure5_Generations, makeFigure3_Generations

% ----------------------------- Parse inputs -----------------------------
p = inputParser;
addParameter(p, 'resultsRoot',     '', @ischar);
% Default [] means "use simParams.deltaTrait", which is the rule stated in
% Methods: exclude data points where either x_i exceeds -delta. Pass a number
% here only to override it.
addParameter(p, 'proximityCutoff', [], ...
    @(v) isnumeric(v) && (isempty(v) || isscalar(v) || numel(v) == 3));
addParameter(p, 'yLimits',   [-3, 2], @(v) isnumeric(v) && numel(v) == 2);
addParameter(p, 'outputFile', 'Figure_ModularConstantSupply_Generations', @ischar);
parse(p, varargin{:});

proximityCutoff = p.Results.proximityCutoff;
if isscalar(proximityCutoff)
    proximityCutoff = repmat(proximityCutoff, 1, 3);
elseif isempty(proximityCutoff)
    proximityCutoff = [NaN, NaN, NaN];   % resolved to deltaTrait per panel
end
yLimits         = p.Results.yLimits;
outputFile      = p.Results.outputFile;

if nargin < 1 || isempty(simParamsRef)
    simParamsRef = struct();
end

% ----------------------------- Settings --------------------------------
subplot_labels = {'A','B','C','D','E','F'};
markers        = {'d','^','v','>','<','o'};

thisDir = fileparts(mfilename('fullpath'));
if isempty(thisDir), thisDir = pwd; end

resultsRoot = p.Results.resultsRoot;
if isempty(resultsRoot)
    projRoot   = fileparts(thisDir);
    candidates = {fullfile(projRoot, 'results'), fullfile(projRoot, 'results_test')};
    resultsRoot = '';
    for i = 1:numel(candidates)
        if isfolder(fullfile(candidates{i}, 'ConstantSupply'))
            resultsRoot = candidates{i};
            break;
        end
    end
    if isempty(resultsRoot)
        error('makeFigure_ConstantSupply:NoResults', ...
              ['No ConstantSupply results found under %s or %s.\n' ...
               'Run Run_modularConstantSupply first.'], candidates{1}, candidates{2});
    end
end
fprintf('makeFigure_ConstantSupply loading data from: %s\n', resultsRoot);

dirs          = {'SSWM', 'CM_Asexual', 'CM_Sexual'};
column_titles = {'Successive mutations', ...
                 'Concurrent mutations, linked modules', ...
                 'Concurrent mutations, unlinked modules'};

% --------------------------- Locate files ------------------------------
file_names = cell(1,3);
for i = 1:3
    dirPath = fullfile(resultsRoot, 'ConstantSupply', dirs{i});
    files = dir(fullfile(dirPath, sprintf('ModularConstSupplyFGM_%s_*.mat', dirs{i})));

    if isempty(files)
        error('makeFigure_ConstantSupply:NoResults', ...
              ['No constant-supply results found in %s.\n' ...
               'Run Run_modularConstantSupply first.'], dirPath);
    end

    if isscalar(files)
        file_names{i} = fullfile(files(1).folder, files(1).name);
    else
        file_names{i} = pickBestMatchingFile(files, simParamsRef);
    end

    fprintf('  Loading: %s\n', file_names{i});
end

% --------------------------- Figure layout -----------------------------
figure('Units','centimeters','Position',[1,1,17.8,12]);
addColumnAnnotations(column_titles);
subplot_positions = defineSubplotPositions(10/12);

for i = 1:6
    subplot('Position', subplot_positions{i});
    addSubplotLabel(subplot_labels{i}, subplot_positions{i});

    if i <= 3
        plotTrajectoryPanel(file_names{i}, markers, dirs{i});
    else
        plotLogRatioPanel(file_names{i-3}, markers, proximityCutoff(i-3), yLimits, dirs{i-3});
    end
end

figDir = fullfile(resultsRoot, 'Figures');
if ~isfolder(figDir)
    mkdir(figDir);
end

print(fullfile(figDir, [outputFile '.pdf']), '-dpdf', '-vector');
fprintf('Saved %s\n', fullfile(figDir, [outputFile '.pdf']));
end

% ======================================================================
% Helper subfunctions
% ======================================================================

function picked = pickBestMatchingFile(files, simParamsRef)
    tags = {};

    if isfield(simParamsRef, 'popSize')
        tags{end+1} = sprintf('N%.0e', simParamsRef.popSize); %#ok<AGROW>
    end
    if isfield(simParamsRef, 'deltaTrait')
        tags{end+1} = sprintf('d%.2f', simParamsRef.deltaTrait); %#ok<AGROW>
    end
    if isfield(simParamsRef, 'ellipseParams')
        tags{end+1} = sprintf('eR%.2f', simParamsRef.ellipseParams(1) / simParamsRef.ellipseParams(2)); %#ok<AGROW>
    end
    if isfield(simParamsRef, 'landscapeStdDev')
        tags{end+1} = sprintf('s%.2f', simParamsRef.landscapeStdDev); %#ok<AGROW>
    end
    if isfield(simParamsRef, 'geneticTargetSize')
        tags{end+1} = sprintf('L%d-%d', simParamsRef.geneticTargetSize(1), simParamsRef.geneticTargetSize(2)); %#ok<AGROW>
    end
    if isfield(simParamsRef, 'beneficialFraction')
        tags{end+1} = sprintf('f%.3g-%.3g', simParamsRef.beneficialFraction(1), simParamsRef.beneficialFraction(2)); %#ok<AGROW>
    end
    if isfield(simParamsRef, 'mutationRate')
        tags{end+1} = sprintf('M%.1e', simParamsRef.mutationRate); %#ok<AGROW>
    end

    scores = zeros(1, numel(files));
    for f = 1:numel(files)
        for t = 1:numel(tags)
            if contains(files(f).name, tags{t})
                scores(f) = scores(f) + 1;
            end
        end
    end

    bestIdx = find(scores == max(scores));
    if isscalar(bestIdx)
        picked = fullfile(files(bestIdx).folder, files(bestIdx).name);
        return;
    end

    % A tie used to fall through to "most recent" with a warning, which is easy
    % to miss in a long run log and had already let a superseded result file be
    % picked up silently. Superseded runs belong in a _relegated/ subfolder,
    % which dir() does not match; anything still tied here is a genuine
    % ambiguity the caller has to resolve.
    tied = {files(bestIdx).name};
    error('makeFigure_ConstantSupply:AmbiguousResults', ...
          ['%d result files match the requested parameters equally well:\n  %s\n' ...
           'Move the superseded ones into a _relegated/ subfolder.'], ...
          numel(tied), strjoin(tied, sprintf('\n  ')));
end

function pos = defineSubplotPositions(s)
    pos = {
        [0.06, 0.56*s, 0.25, 0.38*s];
        [0.38, 0.56*s, 0.25, 0.38*s];
        [0.70, 0.56*s, 0.25, 0.38*s];
        [0.06, 0.08*s, 0.25, 0.38*s];
        [0.38, 0.08*s, 0.25, 0.38*s];
        [0.70, 0.08*s, 0.25, 0.38*s];
    };
end

function addColumnAnnotations(titles)
    x = [0.06, 0.38, 0.70];
    for i = 1:numel(titles)
        annotation('textbox', [x(i), 0.95, 0.25, 0.05], ...
            'String', titles{i}, ...
            'FontSize', 16, ...
            'FontName', 'Helvetica', ...
            'HorizontalAlignment', 'center', ...
            'EdgeColor', 'none');
    end
end

function addSubplotLabel(label, pos)
    annotation('textbox', [pos(1)-0.055, pos(2)+pos(4)+0.012, 0.045, 0.045], ...
        'String', label, ...
        'FontSize', 14, ...
        'FontWeight', 'bold', ...
        'EdgeColor', 'none');
end

function plotTrajectoryPanel(file, markers, regime)
    d  = load(file);
    sp = d.simParams;
    av = getAverageTrajectory(d);

    [~, lvl] = defineGaussianPDF(sp);
    plotPDFContour(sp, lvl, [0.7 0.7 0.7]);
    hold on;

    plotReferenceLines(sp);

    % Drawn in two passes rather than one loop, so that the layer order is
    % global rather than per-trajectory. In a single interleaved loop the
    % analytical curve of initial condition j+1 is drawn after the simulated
    % mean of initial condition j and overdraws it. The simulated means must
    % end up on top of every analytical curve.

    % Pass 1 - all analytical predictions for the constant-supply model.
    %
    % The concurrent-mutations curves are RECOMPUTED here rather than read from
    % the .mat. The stored ones were produced before the predictors learned to
    % stop at the edge of the Desai-Fisher domain, so they run on past it, and
    % using them would put a different curve in this panel from the one panel
    % D-F draws for the same regime. Recomputing costs a second or two and
    % guarantees the two rows agree. The successive-mutations prediction has no
    % such boundary (it is a closed-form exponential drawn only as far as the
    % simulated mean went), so the stored curve is used there.
    analytic = analyticCurvesForPanel(d, sp, regime, av);
    for j = 1:numel(analytic)
        R = analytic{j};
        if isempty(R) || size(R, 1) < 2, continue; end
        R = R(all(isfinite(R), 2), :);
        if size(R, 1) < 2, continue; end
        plot(R(:,1), R(:,2), '-', 'Color', '#EDB120', 'LineWidth', 2);
    end

    % Pass 2 - replicate spread around each mean path, drawn as in Figures 2-5:
    % each replicate's deviation from the mean path is projected onto the axis
    % perpendicular to the local direction of travel, so the band shows
    % variation in the SHAPE of the path rather than in how far along it a
    % replicate happens to be. Under the analytical curves, under the means.
    resultTable = getResultTable(d);
    for j = 1:numel(av.averageTimeStamp)
        ts = av.averageTimeStamp{j};
        bandData = computeTrajectoryPerpendicularBand(resultTable(j,:), size(ts, 2));
        fill([bandData.upperX1, fliplr(bandData.lowerX1)], ...
             [bandData.upperX2, fliplr(bandData.lowerX2)], ...
             [0.6118 0.7333 0.8431], 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    end

    % Pass 3 - all simulated averages and start markers, on top.
    for j = 1:numel(av.averageTimeStamp)
        ts = av.averageTimeStamp{j};
        plot(ts(1,:), ts(2,:), '-', 'Color', '#2776A8', 'LineWidth', 1.0);
        scatter(ts(1,1), ts(2,1), 70, markers{min(j, numel(markers))}, ...
            'MarkerEdgeColor', '#2776A8', 'MarkerFaceColor', '#2776A8');
    end

    customizeAxes(1.2.*[-2.8, 0.05], 1.2.*[-2.1875, 0.05]);
    text(-3, 0.4, 'Module 1 Performance', 'FontName','Helvetica', 'FontSize',12);
    text(0.4, -0.1, 'Module 2 Performance', 'FontName','Helvetica', 'FontSize',12, 'Rotation',270);
    text(0.08, 0.12, '0', 'FontName','Helvetica', 'FontSize',12);
end

function av = getAverageTrajectory(d)
    if isfield(d, 'averageTrajectory')
        av = d.averageTrajectory;
    elseif isfield(d, 'ave')
        av = d.ave;
    else
        error('No trajectory data found.');
    end
end

function plotPDFContour(sp, lvl, lineColor)
    fc = fcontour(@(x1,x2) exp(-sqrt((x1./sp.ellipseParams(1)).^2 + (x2./sp.ellipseParams(2)).^2).^2 / ...
        (2 * sp.landscapeStdDev^2)), ...
        'LineColor', lineColor, 'LineWidth', 1.4);
    fc.LevelList = [0.99, lvl];
end

function plotReferenceLines(sp)
    x1 = linspace(-3, 0.5, 1000);
    plot(x1, x1, '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 1);
    text(-2.3, -2.4, '$\mathbf{x_1 = x_2}$', 'Interpreter','latex', 'FontSize',12, 'Color',[0.4 0.4 0.4]);

    % Equal fitness benefits ray of the DECLINING-supply modular GPFM, drawn
    % dashed for reference. The constant-supply model is not predicted to
    % approach it in any regime; under complete linkage it does reach a balance,
    % but at a lower ratio (see the header of this file).
    R_bar = (sp.ellipseParams(2)^2 * sp.geneticTargetSize(2)) / ...
            (sp.ellipseParams(1)^2 * sp.geneticTargetSize(1));
    plot(x1, R_bar * x1, '--', 'Color', '#D95319', 'LineWidth', 1.5);
    text(-3.3, -1.8, '$\mathbf{s_1 = s_2}$', 'Interpreter','latex', 'FontSize',12, 'Color','#D95319');
end

function plotLogRatioPanel(file, markers, proximityCutoff, yLimits, regime)
    d  = load(file);
    sp = d.simParams;
    if ~isfinite(proximityCutoff)
        proximityCutoff = sp.deltaTrait;   % Methods rule: exclude |x_i| < delta
    end
    resultTable = getResultTable(d);

    numTimeStamp = 200;
    maxGenAll = 0;
    numSims = size(resultTable, 2);
    % Stop displaying log R once more than 20% of replicates have dropped out
    % (either finished, or excluded by proximityCutoff). Replicates do not drop
    % out at random: the ones that persist longest are those whose trait 1 is
    % still far from the optimum, which biases the surviving mean upward and
    % produces an apparent late-time increase in R. Requiring 80% retention
    % shortens the curves but removes that bias.
    %
    % NOTE: the previous rule, min(40, floor(0.5*numSims)), was capped at 40 and
    % so retained only 40/250 = 16% at paper scale, not the 50% it appears to say.
    minRetainedFraction = 0.80;
    minSimsCutoff = max(2, ceil(minRetainedFraction * numSims));

    hold on;

    % Theory first, simulation on top, matching the layer order of the
    % trait-space panels. The horizon is the longest run in the panel; the
    % prediction functions truncate themselves where their approximation stops
    % holding, so this is an upper bound rather than the drawn extent.
    tHorizon = 0;
    for j = 1:size(resultTable, 1)
        for k = 1:numSims
            traj = resultTable{j, k};
            if ~isempty(traj)
                tHorizon = max(tHorizon, max(traj(:, 2)));
            end
        end
    end
    % No theory curve on this row. Of the three regimes only the
    % successive-mutations prediction is available in closed form, and a row in
    % which one panel of three carries a curve reads as a missing result rather
    % than as a deliberate omission. The trait-space row above is unaffected.
    % predictLogRatioTrajectory_ConstantSupply remains available to callers.

    for j = 1:size(resultTable, 1)
        maxGen = 0;
        for k = 1:numSims
            traj = resultTable{j, k};
            if ~isempty(traj)
                maxGen = max(maxGen, max(traj(:, 2)));
            end
        end

        if maxGen <= 1
            continue;
        end

        genStamps = floor(linspace(1, maxGen, numTimeStamp));
        logRatioAll = NaN(numSims, numTimeStamp);

        for k = 1:numSims
            traj = resultTable{j, k};
            if isempty(traj)
                continue;
            end

            simTime = traj(:, 2);
            trait1  = traj(:, 3);
            trait2  = traj(:, 4);
            simMaxGen = max(simTime);

            for t = 1:numTimeStamp
                if genStamps(t) <= simMaxGen
                    idx = find(simTime <= genStamps(t), 1, 'last');
                    if isempty(idx)
                        idx = 1;
                    end

                    x1val = trait1(idx);
                    x2val = trait2(idx);

                    if proximityCutoff > 0 && (abs(x1val) < proximityCutoff || abs(x2val) < proximityCutoff)
                        break;
                    end

                    % Guard against the floating-point residue that repeated
                    % lattice additions can leave in place of an exact zero.
                    if abs(x1val) > 1e-9 && abs(x2val) > 1e-9
                        logRatioAll(k, t) = log(abs(x2val) / abs(x1val));
                    end
                end
            end
        end

        numContributing = sum(~isnan(logRatioAll), 1);
        validTimeIdx = numContributing >= minSimsCutoff;

        meanLogRatio = mean(logRatioAll, 1, 'omitnan');
        stdLogRatio  = std(logRatioAll, 0, 1, 'omitnan');

        meanLogRatio(~validTimeIdx) = NaN;
        stdLogRatio(~validTimeIdx)  = NaN;

        validGens = genStamps(validTimeIdx);
        if ~isempty(validGens)
            maxGenAll = max(maxGenAll, max(validGens));
        end

        validIdx = ~isnan(meanLogRatio);
        if any(validIdx)
            xFill = [genStamps(validIdx), fliplr(genStamps(validIdx))];
            yFill = [meanLogRatio(validIdx) + stdLogRatio(validIdx), ...
                     fliplr(meanLogRatio(validIdx) - stdLogRatio(validIdx))];
            fill(xFill, yFill, [0.15, 0.46, 0.66], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        end

        plot(genStamps(validIdx), meanLogRatio(validIdx), '-', ...
            'Color', '#2776A8', 'LineWidth', 0.9);

        if any(validIdx)
            firstValid = find(validIdx, 1, 'first');
            plot(genStamps(firstValid), meanLogRatio(firstValid), ...
                'Marker', markers{min(j, numel(markers))}, ...
                'Color', '#2776A8', ...
                'MarkerFaceColor', '#2776A8', ...
                'MarkerSize', 7);
        end
    end

    if maxGenAll <= 1
        maxGenAll = 100;
    end

    % Module-selection balance of the declining-supply model, for reference.
    yline(log((sp.geneticTargetSize(2)*sp.ellipseParams(2)^2) / ...
              (sp.geneticTargetSize(1)*sp.ellipseParams(1)^2)), ...
          '--', 'Color', '#D95319', 'LineWidth', 1.5);

    xlim([1, maxGenAll]);
    ylim(yLimits);
    set(gca, 'TickLabelInterpreter','latex', 'FontSize', 10, 'Box', 'on');
    xlabel('Generations', 'FontName','Helvetica', 'FontSize', 12);
    ylabel('$\log(x_2/x_1)$', 'Interpreter', 'latex', 'FontSize', 12);
end

function [f, lvl] = defineGaussianPDF(sp)
    syms x1 x2
    f = exp(-sqrt((x1./sp.ellipseParams(1)).^2 + (x2./sp.ellipseParams(2)).^2).^2 / ...
        (2 * sp.landscapeStdDev^2));
    dH  = det(hessian(f, [x1, x2]));
    dHf = matlabFunction(dH, 'Vars', [x1, x2]);

    xyStar = fminsearch(@(x) abs(dHf(x(1), x(2))), [0, 0]);
    lvl = double(subs(f, [x1, x2], xyStar));
end

function customizeAxes(xl, yl)
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    ax.Box = 'off';
    ax.XColor = 'k';
    ax.YColor = 'k';
    ax.LineWidth = 1;
    xlim(xl);
    ylim(yl);
    ax.XTick = [-2, -1];
    ax.YTick = [-2, -1];
    ax.TickLength = [0.015, 0.015];
    ax.FontSize = 10;
    ax.TickLabelInterpreter = 'latex';
end

function resultTable = getResultTable(d)
    fieldNames = fieldnames(d);
    resultFields = fieldNames(contains(fieldNames, 'result', 'IgnoreCase', true));

    for i = 1:length(resultFields)
        fn = resultFields{i};
        if isstruct(d.(fn)) && isfield(d.(fn), 'resultTable')
            resultTable = d.(fn).resultTable;
            return;
        end
    end

    error('Could not find resultTable in data file.');
end

function curves = analyticCurvesForPanel(d, sp, regime, av)
% Analytical trait-space curves for one regime, matching what the log-ratio
% panel of the same regime will draw.

    curves = {};

    switch upper(regime)
        case 'CM_ASEXUAL'
            spd = sp;
            spd.initialPhenotypes = latticeInitial(sp);
            curves = predictModularCM_ConstantSupply(spd);

        case 'CM_SEXUAL'
            spd = sp;
            spd.initialPhenotypes = latticeInitial(sp);
            curves = predictFullRecomb_ConstantSupply(spd);

        otherwise
            % SSWM: use the stored closed-form curve.
            if isfield(d, 'analyticalTrajectories')
                curves = d.analyticalTrajectories;
            end
    end

    if isempty(curves) && isfield(d, 'analyticalTrajectories')
        curves = d.analyticalTrajectories;
    end

    % Normalise to a plain cell column of [x1, x2] matrices.
    if ~isempty(curves)
        curves = curves(:);
    end

    if nargin > 3 && numel(curves) < numel(av.averageTimeStamp)
        curves(end+1:numel(av.averageTimeStamp)) = {[]};
    end
end

function X0 = latticeInitial(sp)
% Discretize the initial phenotypes onto the delta lattice, as the simulations
% do, so the theory starts where the populations actually started.
    X0 = min(-sp.deltaTrait * round(-sp.initialPhenotypes / sp.deltaTrait), 0);
end
