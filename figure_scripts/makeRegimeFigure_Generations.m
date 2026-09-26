function makeRegimeFigure_Generations(regime, pleiotropicSimParamsRef, modularSimParamsRef, varargin)
% makeRegimeFigure_Generations — Generate a regime-organized main figure.
%
% Organizes the main figures by REGIME rather than by MODEL, per
% request: each figure now fixes one evolutionary regime and places the two
% core GPFMs (Pleiotropic, Modular) side by side as columns, so the
% cross-model contrast (the main point of the paper) is visible within a
% single figure instead of requiring the reader to flip between figures.
%
% Layout (2x2):
%   Columns: Pleiotropic GPFM | Modular GPFM
%   Row 1:   phase-plane trajectory plot
%   Row 2:   log(x2/x1) vs. generations
%
% Usage:
%   makeRegimeFigure_Generations('SSWM',       pleioParams, modularParams)
%   makeRegimeFigure_Generations('CM_Asexual', pleioParams, modularParams, 'test')
%   makeRegimeFigure_Generations('CM_Sexual',  pleioParams, modularParams, 'full', ...
%                                 'proximityCutoff', 0.1, 'outputFile', 'Figure4')
%
% Inputs:
%   regime                  - 'SSWM', 'CM_Asexual', or 'CM_Sexual'
%   pleiotropicSimParamsRef - simParams struct used to disambiguate multiple
%                             matching PleiotropicFGM files in the regime folder
%   modularSimParamsRef     - simParams struct used to disambiguate multiple
%                             matching ModularFGM files in the regime folder
%
% Optional positional input:
%   mode - 'test', 'full', or 'auto' (default 'auto')
%
% Optional name-value input:
%   'proximityCutoff' - Exclude ratio data when |x_i| < cutoff. Default: 0.
%   'outputFile'      - Override default output filename (without extension).
%
% Outputs:
%   Saves Figure_<regime>_Generations.pdf to results/Figures/
%   (or outputFile.pdf if outputFile override is specified)

% ----------------------------- Validate regime --------------------------
validRegimes = {'SSWM', 'CM_Asexual', 'CM_Sexual'};
if ~ismember(regime, validRegimes)
    error('regime must be one of: %s', strjoin(validRegimes, ', '));
end

% ----------------------------- Parse inputs -----------------------------
mode = 'auto';
proximityCutoff = 0;
outputFileOverride = '';

idx = 1;
if nargin >= 4 && ischar(varargin{1}) && ismember(lower(varargin{1}), {'test','full','auto'})
    mode = lower(varargin{1});
    idx = 2;
end

while idx <= length(varargin) - 1
    key = varargin{idx};
    val = varargin{idx + 1};
    switch lower(key)
        case 'proximitycutoff'
            proximityCutoff = val;
        case 'outputfile'
            outputFileOverride = val;
        otherwise
            error('Unknown parameter: %s', key);
    end
    idx = idx + 2;
end

% ----------------------------- Settings --------------------------------
% A,B = phase-plane row. C,D = log-ratio row.
subplot_labels = {'A','B','C','D'};
markers        = {'d','^','v','>','<','o'};

regimeTitles = struct( ...
    'SSWM',       'Successive mutations', ...
    'CM_Asexual', 'Concurrent mutations, linked chromosomes', ...
    'CM_Sexual',  'Concurrent mutations, unlinked chromosomes');

column_titles = {'Pleiotropic GPFM', 'Modular GPFM'};

% ----------------------- Resolve results root --------------------------
figureRoot  = fileparts(mfilename('fullpath'));
resultsRoot = resolveResultsRoot(figureRoot, mode);
fprintf('makeRegimeFigure_Generations (%s) loading data from: %s\n', regime, resultsRoot);

regimeDir = fullfile(resultsRoot, regime);

% --------------------------- Locate result files ------------------------
% Column 1: PleiotropicFGM, Column 2: ModularFGM — both drawn from the SAME
% regime folder, since the whole point of this figure is fixing the regime
% and letting the model vary.
file_names = cell(1,2);
file_names{1} = locateFile(regimeDir, 'PleiotropicFGM', regime, pleiotropicSimParamsRef);
file_names{2} = locateFile(regimeDir, 'ModularFGM',     regime, modularSimParamsRef);

fprintf('  Loading (Pleiotropic): %s\n', file_names{1});
fprintf('  Loading (Modular):     %s\n', file_names{2});

% --------------------------- Figure layout -----------------------------
figure('Units','centimeters','Position',[1,1,17.8,17.1]);
addColumnAnnotations(column_titles);
%addRegimeTitle(regimeTitles.(regime));
subplot_positions = defineSubplotPositions();

for i = 1:4
    subplot('Position', subplot_positions{i});
    addSubplotLabel(subplot_labels{i}, subplot_positions{i});

    col = mod(i-1, 2) + 1;   % 1 = Pleiotropic, 2 = Modular
    row = ceil(i / 2);       % 1 = phase-plane, 2 = log-ratio

    if row == 1
        if col == 1
            plotPhasePlane_Pleiotropic(file_names{1}, markers);
        else
            plotPhasePlane_Modular(file_names{2}, markers, regime);
        end
    else
        if col == 1
            plotLogRatio_Pleiotropic(file_names{1}, markers, proximityCutoff, regime);
        else
            plotLogRatio_Modular(file_names{2}, markers, proximityCutoff, regime);
        end
    end
end

figDir = fullfile(resultsRoot, 'Figures');
if ~isfolder(figDir)
    mkdir(figDir);
end

if ~isempty(outputFileOverride)
    saveFile = outputFileOverride;
else
    saveFile = ['Figure_' regime '_Generations'];
end

print(fullfile(figDir, [saveFile '.pdf']), '-dpdf', '-vector');
fprintf('Saved %s\n', fullfile(figDir, [saveFile '.pdf']));
end

% ======================================================================
% File location (shared across both model columns)
% ======================================================================

function fname = locateFile(regimeDir, figureType, regime, simParamsRef)
% Resolve the single results file matching the requested model and regime.
%
% Previously this filtered only on population size and, when that failed to
% narrow the set, fell back to "most recent" with a warning. Every result file
% in a regime folder shares the same population size, so the filter never
% discriminated and the choice among duplicates was effectively alphabetical.
% Superseded runs sitting alongside current ones could therefore be picked up
% silently. Retired runs now live in a _relegated/ subfolder (not matched by
% dir), and what remains is filtered on every tag the filename carries. If more
% than one file still matches, that is an ambiguity the caller has to resolve,
% so it raises an error instead of guessing.

    files = dir(fullfile(regimeDir, sprintf('%s_%s_*.mat', figureType, regime)));
    files = files(~[files.isdir]);
    if isempty(files)
        error('makeRegimeFigure:NoResults', ...
              'No %s results for regime %s in %s.', figureType, regime, regimeDir);
    end

    names = {files.name};

    % Successive narrowing on the tags buildFilename writes into the name.
    % A tag is only applied if it actually leaves something behind, so a file
    % written under an older naming convention is not silently excluded.
    names = narrow(names, sprintf('_N%.0e_', simParamsRef.popSize));
    if isfield(simParamsRef, 'mutationRate')
        names = narrow(names, sprintf('_M%.1e_', simParamsRef.mutationRate));
    end
    if isfield(simParamsRef, 'deltaTrait')
        names = narrow(names, sprintf('_d%.2f_', simParamsRef.deltaTrait));
    end
    if isfield(simParamsRef, 'geneticTargetSize') && numel(simParamsRef.geneticTargetSize) == 2
        names = narrow(names, sprintf('_L%d-%d', simParamsRef.geneticTargetSize(1), ...
                                                simParamsRef.geneticTargetSize(2)));
    end

    if isscalar(names)
        fname = fullfile(regimeDir, names{1});
        return;
    end

    if isempty(names)
        error('makeRegimeFigure:NoMatchingResults', ...
              ['No %s file in %s matches the requested parameters ' ...
               '(N = %g, U = %g, delta = %g).\nFiles present: %s'], ...
              figureType, regimeDir, simParamsRef.popSize, ...
              getfielddef(simParamsRef, 'mutationRate', NaN), ...
              getfielddef(simParamsRef, 'deltaTrait', NaN), ...
              strjoin({files.name}, ', '));
    end

    error('makeRegimeFigure:AmbiguousResults', ...
          ['%d %s files in %s match the requested parameters equally well:\n  %s\n' ...
           'Move the superseded ones into a _relegated/ subfolder, or tighten ' ...
           'the reference parameters passed to this function.'], ...
          numel(names), figureType, regimeDir, strjoin(names, sprintf('\n  ')));
end

function kept = narrow(names, tag)
% Keep only the names containing tag, unless that would discard everything.
    hit = contains(names, tag);
    if any(hit)
        kept = names(hit);
    else
        kept = names;
    end
end

function v = getfielddef(s, f, default)
    if isfield(s, f); v = s.(f); else; v = default; end
end

% ======================================================================
% Layout helpers
% ======================================================================

function resultsRoot = resolveResultsRoot(figureRoot, mode)
    switch mode
        case 'test'
            candidates = {fullfile(figureRoot,'../results_test'), fullfile(figureRoot,'../../results_test')};
        case 'full'
            candidates = {fullfile(figureRoot,'../results'), fullfile(figureRoot,'../../results')};
        otherwise
            candidates = {fullfile(figureRoot,'../results'), fullfile(figureRoot,'../../results'), ...
                          fullfile(figureRoot,'../results_test'), fullfile(figureRoot,'../../results_test')};
    end
    for i = 1:numel(candidates)
        if isfolder(candidates{i})
            resultsRoot = candidates{i};
            return;
        end
    end
    error('No results directory found for mode ''%s''.', mode);
end

function pos = defineSubplotPositions()
    % 2x2 grid. The canvas was widened from 12.5 to 17.8 cm (journal full
    % width) and the height scaled by the SAME factor, 12.0 -> 17.1 cm, so
    % that these normalized fractions keep the designed panel proportions.
    % Panels are now 6.4 x 5.4 cm, i.e. the original 4.5 x 3.8 cm scaled by
    % 17.8/12.5 = 1.424. Changing the canvas width alone stretches every
    % panel horizontally and visibly distorts the phase-plane row, where the
    % x1 = x2 line then no longer renders at its true slope.
    % Gap between rows and title clearance scale with it. Title clearance is 2.6 cm
    % (y=0.95, height=0.05 on a 12 cm canvas), matching
    % makeFigure5_Generations.m and the original pre-edit Figures 2-4 --
    % an earlier revision of this file only gave 1.0 cm of clearance here,
    % which was too tight and let the "Module 1/2 Performance" text
    % (drawn outside the axis limits on purpose) crowd the title.
    colW  = 0.36;    % 6.4 cm / 17.8 cm
    rowH  = 0.3167;  % 5.4 cm / 17.1 cm
    xCol1 = 0.13;
    xCol2 = 0.57;
    yRow2 = 0.0667;  % bottom row (log-ratio): 1.1 cm / 17.1 cm
    yRow1 = 0.4667;  % top row (phase-plane): 8.0 cm / 17.1 cm

    pos = {
        [xCol1, yRow1, colW, rowH];   % A: Pleiotropic phase-plane
        [xCol2, yRow1, colW, rowH];   % B: Modular phase-plane
        [xCol1, yRow2, colW, rowH];   % C: Pleiotropic log-ratio
        [xCol2, yRow2, colW, rowH];   % D: Modular log-ratio
    };
end

function addColumnAnnotations(titles)
    x = [0.13, 0.57];
    for i = 1:numel(titles)
        annotation('textbox', [x(i), 0.95, 0.36, 0.05], ...
            'String', titles{i}, ...
            'FontSize', 16, 'FontName', 'Helvetica', ...
            'HorizontalAlignment', 'center', 'EdgeColor', 'none');
    end
end

function addRegimeTitle(titleStr)
    annotation('textbox', [0.05, 0.95, 0.90, 0.05], ...
        'String', titleStr, ...
        'FontSize', 17, 'FontWeight', 'bold', 'FontName', 'Helvetica', ...
        'HorizontalAlignment', 'center', 'EdgeColor', 'none');
end

function addSubplotLabel(label, pos)
    annotation('textbox', [pos(1)-0.065, pos(2)+pos(4)+0.012, 0.045, 0.045], ...
        'String', label, 'FontSize', 14, 'FontWeight', 'bold', 'EdgeColor', 'none');
end

% ======================================================================
% Pleiotropic column (row 1: phase-plane, row 2: log-ratio)
% Logic ported unchanged from makeFigure2_Generations.m
% ======================================================================

function plotPhasePlane_Pleiotropic(file, markers)
    color = '#EDB120';
    d = load(file);
    sp = d.simParams;
    av = getAverageTrajectory(d);
    resultTable = getResultTable(d);

    [~, lvl] = defineGaussianPDF(sp);
    plotPDFContour(sp, lvl, [0.7 0.7 0.7]); hold on;

    if isfield(d, 'analyticalTrajectories')
        plotTrajectories(d.analyticalTrajectories, av, markers, color, resultTable);
    else
        plotAverageTrajectories(av, markers, color, resultTable);
    end

    plotReferenceLines_Pleiotropic();
    customizeAxes(1.2.*[-2.8, 0.05], 1.2.*[-2.1875, 0.05]);
    text(-3, 0.4, 'Module 1 Performance', 'FontName','Helvetica', 'FontSize',12);
    text(0.4, -0.1, 'Module 2 Performance', 'FontName','Helvetica', 'FontSize',12, 'Rotation',270);
    text(0.08, 0.12, '0', 'FontName','Helvetica', 'FontSize',12);
end

function plotReferenceLines_Pleiotropic()
    x1 = linspace(-3, 0.5, 1000);
    plot(x1, x1, '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 1);
    text(-2.3, -2.4, '$\mathbf{x_1 = x_2}$', 'Interpreter','latex', 'FontSize',12, 'Color',[0.4 0.4 0.4]);
end

function plotLogRatioTheory(sp, regime, model, tMax)
    % Overlay the time-resolved theoretical prediction for log R = log(x_2/x_1).
    %
    % Drawn in the successive-mutations regime only. Under concurrent mutations
    % a time-resolved prediction exists for the modular model but not for the
    % pleiotropic one, so drawing it would fill one column of the row and leave
    % the other permanently empty. predictLogRatioTrajectory still supplies the
    % concurrent-mutations curves for any caller that wants them.
    if ~strcmpi(regime, 'SSWM')
        return;
    end
    if isempty(which('predictLogRatioTrajectory'))
        warning('predictLogRatioTrajectory not on the path; skipping log R theory.');
        return;
    end
    [tTh, logRTh] = predictLogRatioTrajectory(sp, regime, model, tMax);
    for j = 1:numel(tTh)
        if isempty(tTh{j}) || isempty(logRTh{j}), continue; end
        plot(tTh{j}, logRTh{j}, '-', 'Color', '#EDB120', 'LineWidth', 2);
    end
end

function plotLogRatio_Pleiotropic(file, markers, proximityCutoff, regime)
    d = load(file);
    sp = d.simParams;
    resultTable = getResultTable(d);

    [genStampsByAngle, meanByAngle, stdByAngle, validByAngle, maxGenAll] = ...
        computeLogRatioSeries(sp, resultTable, proximityCutoff);

    hold on;

    if maxGenAll <= 1, maxGenAll = 100; end

    % Theory first, simulation on top (same layer order as the trait-space
    % panels). Returns empty for the concurrent-mutations regimes: the only
    % pleiotropic theory is the successive-mutations result, which is a
    % time-free shape reference and would sit on the wrong clock here.
    plotLogRatioTheory(sp, regime, 'Pleiotropic', maxGenAll);

    for j = 1:numel(genStampsByAngle)
        plotOneLogRatioSeries(genStampsByAngle{j}, meanByAngle{j}, stdByAngle{j}, ...
            validByAngle{j}, markers{j});
    end

    % No theoretical balance-ratio line for the pleiotropic model (unlike
    % Modular) -- only the y = 0 reference.
    yline(0, '-', 'LineWidth', 1.4, 'Color', [0.4 0.4 0.4]);

    finalizeLogRatioAxes(maxGenAll);
end

% ======================================================================
% Modular column (row 1: phase-plane, row 2: log-ratio)
% Logic ported unchanged from makeFigure3_Generations.m
% ======================================================================

function plotPhasePlane_Modular(file, markers, regime)
    color = '#EDB120';
    d = load(file);
    sp = d.simParams;
    av = getAverageTrajectory(d);
    resultTable = getResultTable(d);

    [~, lvl] = defineGaussianPDF(sp);
    plotPDFContour(sp, lvl, [0.7 0.7 0.7]); hold on;

    if isfield(d, 'analyticalTrajectories')
        plotTrajectories(d.analyticalTrajectories, av, markers, color, resultTable);
    else
        plotAverageTrajectories(av, markers, color, resultTable);
    end

    plotReferenceLines_Modular(sp, regime);
    customizeAxes(1.2.*[-2.8, 0.05], 1.2.*[-2.1875, 0.05]);
    text(-3, 0.4, 'Module 1 Performance', 'FontName','Helvetica', 'FontSize',12);
    text(0.4, -0.1, 'Module 2 Performance', 'FontName','Helvetica', 'FontSize',12, 'Rotation',270);
    text(0.08, 0.12, '0', 'FontName','Helvetica', 'FontSize',12);
end

function plotReferenceLines_Modular(sp, regime)
% Trait-space counterpart of plotModularBalanceLines: the ray the population is
% predicted to settle onto, which differs by regime (see that function for the
% reasoning). Under free reassortment the equal fitness benefits ray is kept as
% a dashed reference, since populations are expected to approach it once they
% leave the concurrent mutations regime.

    if nargin < 2, regime = ''; end

    x1 = linspace(-3, 0.5, 1000);
    plot(x1, x1, '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 1);
    text(-2.3, -2.4, '$\mathbf{x_1 = x_2}$', 'Interpreter','latex', 'FontSize',12, 'Color',[0.4 0.4 0.4]);

    if ~isfield(sp, 'geneticTargetSize')
        return;
    end

    R_bar = (sp.ellipseParams(2)^2 * sp.geneticTargetSize(2)) / ...
            (sp.ellipseParams(1)^2 * sp.geneticTargetSize(1));

    if strcmpi(regime, 'CM_Sexual')
        R_recomb = sp.ellipseParams(2) / sp.ellipseParams(1);
        plot(x1, R_bar * x1,    '--', 'Color', [0.8, 0.3, 0], 'LineWidth', 1.2);
        plot(x1, R_recomb * x1, '-',  'Color', [0.8, 0.3, 0], 'LineWidth', 2);
        % Only the solid (regime's own) ray is labelled in the panel; the dashed
        % equal fitness benefits ray is identified in the caption, to keep the
        % panel from becoming crowded.
        text(-3.3, -1.8, '$\mathbf{x_2/x_1 = a_2/a_1}$', 'Interpreter','latex', 'FontSize',11, 'Color',[0.8,0.3,0]);
    else
        plot(x1, R_bar * x1, '-', 'Color', [0.8, 0.3, 0], 'LineWidth', 2);
        text(-3.3, -1.8, '$\mathbf{s_1 = s_2}$', 'Interpreter','latex', 'FontSize',12, 'Color',[0.8,0.3,0]);
    end
end

function plotLogRatio_Modular(file, markers, proximityCutoff, regime)
    d = load(file);
    sp = d.simParams;
    resultTable = getResultTable(d);

    [genStampsByAngle, meanByAngle, stdByAngle, validByAngle, maxGenAll] = ...
        computeLogRatioSeries(sp, resultTable, proximityCutoff);

    hold on;

    if maxGenAll <= 1, maxGenAll = 100; end

    % Theory first, simulation on top (same layer order as the trait-space panels)
    plotLogRatioTheory(sp, regime, 'Modular', maxGenAll);

    for j = 1:numel(genStampsByAngle)
        plotOneLogRatioSeries(genStampsByAngle{j}, meanByAngle{j}, stdByAngle{j}, ...
            validByAngle{j}, markers{j});
    end

    yline(0, '-', 'LineWidth', 1.4, 'Color', [0.4 0.4 0.4]);
    plotModularBalanceLines(sp, regime);

    finalizeLogRatioAxes(maxGenAll);
end

function plotModularBalanceLines(sp, regime)
% Draw the module-selection balance that the theory predicts FOR THIS REGIME.
%
% The balance value is not the same in all three regimes. In the successive
% mutations regime and in the concurrent mutations regime with complete linkage,
% R converges to a_2^2/a_1^2 (the equal fitness benefits line, equation
% "equal fitness benefits line through traits"). Under free reassortment the two
% modules evolve independently and equation "recomb R" gives a different limit,
%
%     R(t) -> a_2/a_1,
%
% because the rate constants gamma_i and the offsets A_i enter as a_i^2 inside a
% logarithm rather than as a_i^2 directly. Drawing a_2^2/a_1^2 in the
% recombination panel therefore puts the simulated curves above a line they were
% never heading for, which is what made their late-time rise look like an
% artifact. It is not: they are converging on a_2/a_1 from both sides.
%
% In the recombination panel both lines are drawn, the regime's own limit solid
% and the equal fitness benefits line dashed, because the manuscript argues that
% populations eventually leave the concurrent mutations regime and then approach
% the latter.

    balanceColor = [0.8, 0.3, 0];

    % Equal fitness benefits line, weighted by chromosome length when the two
    % chromosomes differ (the L's cancel when L_1 = L_2).
    if isfield(sp,'geneticTargetSize') && sp.geneticTargetSize(1) ~= sp.geneticTargetSize(2)
        logRbarLinked = log( (sp.geneticTargetSize(2)*sp.ellipseParams(2)^2) / ...
                             (sp.geneticTargetSize(1)*sp.ellipseParams(1)^2) );
    else
        logRbarLinked = log(sp.ellipseParams(2)^2 / sp.ellipseParams(1)^2);
    end

    % Free-recombination limit, a_2/a_1. The ratio is unchanged by the
    % manuscript/code unit convention a_i = sqrt(2)*sigma*ellipseParams(i).
    logRbarRecomb = log(sp.ellipseParams(2) / sp.ellipseParams(1));

    if strcmpi(regime, 'CM_Sexual')
        yline(logRbarLinked, '--', 'LineWidth', 1.2, 'Color', balanceColor);
        yline(logRbarRecomb, '-',  'LineWidth', 2,   'Color', balanceColor);
    else
        yline(logRbarLinked, '-',  'LineWidth', 2,   'Color', balanceColor);
    end
end

% ======================================================================
% Shared plotting/data helpers (identical across both original files)
% ======================================================================

function av = getAverageTrajectory(d)
    if isfield(d,'averageTrajectory'), av = d.averageTrajectory;
    elseif isfield(d,'ave'), av = d.ave;
    else, error('No trajectory data found.');
    end
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

function plotPDFContour(sp, lvl, lineColor)
    fc = fcontour(@(x1,x2) exp(-sqrt((x1./sp.ellipseParams(1)).^2 + (x2./sp.ellipseParams(2)).^2).^2 / ...
        (2 * sp.landscapeStdDev^2)), 'LineColor', lineColor, 'LineWidth', 1.4);
    fc.LevelList = [0.99, lvl];
end

function plotAverageTrajectories(av, markers, color, resultTable)
    bandColor = hex2rgbLocal(color);
    % Two passes so the layer order is global: bands underneath, all
    % simulated means and their start markers on top.
    for j = 1:numel(av.averageTimeStamp)
        ts = av.averageTimeStamp{j};
        numTimeStamp = size(ts, 2);
        bandData = computeTrajectoryPerpendicularBand(resultTable(j,:), numTimeStamp);
        fill([bandData.upperX1, fliplr(bandData.lowerX1)], ...
             [bandData.upperX2, fliplr(bandData.lowerX2)], ...
             bandColor, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    end

    for j = 1:numel(av.averageTimeStamp)
        ts = av.averageTimeStamp{j};
        plot(ts(1,:), ts(2,:), '-', 'Color', color, 'LineWidth', 1);
        scatter(ts(1,1), ts(2,1), 70, markers{j}, 'MarkerEdgeColor', color, 'MarkerFaceColor', color);
    end
end

function rgb = hex2rgbLocal(hexStr)
    % fill() requires an RGB triplet, not a hex string, when passed
    % positionally (same issue as in Run_nModuleTest.m). color here is a
    % variable, not always the same literal, so converting generically
    % rather than hardcoding a specific triplet.
    hexStr = strrep(hexStr, '#', '');
    rgb = [hex2dec(hexStr(1:2)), hex2dec(hexStr(3:4)), hex2dec(hexStr(5:6))] / 255;
end

function plotTrajectories(analytic, av, markers, color, resultTable)
    avgLineColor = '#2776A8';
    bandColor = hex2rgbLocal(avgLineColor);

    % Drawn in three passes rather than one loop, so that the layer order is
    % global rather than per-trajectory. In a single interleaved loop the
    % analytical curve of initial condition j+1 is drawn after the simulated
    % mean of initial condition j and overdraws it.
    %
    % Layer order (bottom to top): bands, analytical curves, simulated means.
    % The simulated means and their start markers must end up on top.

    % Pass 1 - uncertainty bands. The band uses avgLineColor, matching the
    % SIMULATED average line, not 'color' (the analytical curve's colour). The
    % analytical curve is a deterministic ODE solution with no replicate
    % variance, so it never gets a band.
    for j = 1:numel(analytic)
        ts = av.averageTimeStamp{j};
        bandData = computeTrajectoryPerpendicularBand(resultTable(j,:), size(ts, 2));
        fill([bandData.upperX1, fliplr(bandData.lowerX1)], ...
             [bandData.upperX2, fliplr(bandData.lowerX2)], ...
             bandColor, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    end

    % Pass 2 - all analytical/theoretical curves.
    for j = 1:numel(analytic)
        if ~isempty(analytic{j,1}) && size(analytic{j,1}, 1) > 0
            R = analytic{j,1}';
            if size(R, 1) >= 2 && size(R, 2) > 0
                plot(R(1,:), R(2,:), '-', 'Color', color, 'LineWidth', 2);
            end
        end
    end

    % Pass 3 - all simulated averages and start markers, on top.
    for j = 1:numel(analytic)
        ts = av.averageTimeStamp{j};
        plot(ts(1,:), ts(2,:), '-', 'Color', avgLineColor, 'LineWidth', 1);
        scatter(ts(1,1), ts(2,1), 70, 'Marker', markers{j}, ...
            'MarkerEdgeColor', avgLineColor, 'MarkerFaceColor', avgLineColor);
    end
end

function [genStampsByAngle, meanByAngle, stdByAngle, validByAngle, maxGenAll] = ...
        computeLogRatioSeries(sp, resultTable, proximityCutoff)

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
    nAngles = length(sp.initialAngles);

    genStampsByAngle = cell(1, nAngles);
    meanByAngle = cell(1, nAngles);
    stdByAngle = cell(1, nAngles);
    validByAngle = cell(1, nAngles);

    for j = 1:nAngles
        maxGen = 0;
        for k = 1:numSims
            traj = resultTable{j, k};
            if ~isempty(traj)
                maxGen = max(maxGen, max(traj(:, 2)));
            end
        end

        if maxGen <= 1
            genStampsByAngle{j} = [];
            meanByAngle{j} = [];
            stdByAngle{j} = [];
            validByAngle{j} = [];
            continue;
        end

        genStamps = floor(linspace(1, maxGen, numTimeStamp));
        logRatioAll = NaN(numSims, numTimeStamp);

        for k = 1:numSims
            traj = resultTable{j, k};
            if isempty(traj), continue; end

            simTime = traj(:, 2);
            trait1 = traj(:, 3);
            trait2 = traj(:, 4);
            simMaxGen = max(simTime);

            for t = 1:numTimeStamp
                if genStamps(t) <= simMaxGen
                    idx = find(simTime <= genStamps(t), 1, 'last');
                    if isempty(idx), idx = 1; end

                    x1val = trait1(idx);
                    x2val = trait2(idx);

                    if proximityCutoff > 0 && (abs(x1val) < proximityCutoff || abs(x2val) < proximityCutoff)
                        break;
                    end

                    if x1val ~= 0 && x2val ~= 0
                        logRatioAll(k, t) = log(abs(x2val) / abs(x1val));
                    end
                end
            end
        end

        numContributing = sum(~isnan(logRatioAll), 1);
        validTimeIdx = numContributing >= minSimsCutoff;

        meanLogRatio = mean(logRatioAll, 1, 'omitnan');
        stdLogRatio = std(logRatioAll, 0, 1, 'omitnan');
        meanLogRatio(~validTimeIdx) = NaN;
        stdLogRatio(~validTimeIdx) = NaN;

        validGens = genStamps(validTimeIdx);
        if ~isempty(validGens)
            maxGenAll = max(maxGenAll, max(validGens));
        end

        genStampsByAngle{j} = genStamps;
        meanByAngle{j} = meanLogRatio;
        stdByAngle{j} = stdLogRatio;
        validByAngle{j} = ~isnan(meanLogRatio);
    end
end

function plotOneLogRatioSeries(genStamps, meanLogRatio, stdLogRatio, validIdx, marker)
    if isempty(genStamps)
        return;
    end

    if any(validIdx)
        xFill = [genStamps(validIdx), fliplr(genStamps(validIdx))];
        yFill = [meanLogRatio(validIdx) + stdLogRatio(validIdx), ...
                 fliplr(meanLogRatio(validIdx) - stdLogRatio(validIdx))];
        fill(xFill, yFill, [0.15, 0.46, 0.66], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    end

    plot(genStamps(validIdx), meanLogRatio(validIdx), '-', 'Color', '#2776A8', 'LineWidth', 0.9);

    if any(validIdx)
        firstValid = find(validIdx, 1, 'first');
        plot(genStamps(firstValid), meanLogRatio(firstValid), 'Marker', marker, 'Color', '#2776A8', ...
             'MarkerFaceColor', '#2776A8', 'MarkerSize', 7);
    end
end

function finalizeLogRatioAxes(maxGenAll)
    xlim([1, maxGenAll]);
    ylim([-3, 2]);
    set(gca, 'TickLabelInterpreter','latex','FontSize',10);
    xlabel('Generations', 'FontName','Helvetica','FontSize',12);
    ylabel('$\log(x_2/x_1)$', 'Interpreter', 'latex', 'FontSize', 12);
end

function [f, lvl] = defineGaussianPDF(sp)
    syms x1 x2
    f = exp(-sqrt((x1./sp.ellipseParams(1)).^2 + (x2./sp.ellipseParams(2)).^2).^2 / ...
        (2 * sp.landscapeStdDev^2));
    dH = det(hessian(f, [x1, x2]));
    dHf = matlabFunction(dH, 'Vars', [x1, x2]);
    lvl = double(subs(f, [x1,x2], fminsearch(@(x) abs(dHf(x(1),x(2))), [0,0])));
end

function customizeAxes(xl, yl)
    ax = gca;
    ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin';
    ax.Box = 'off'; ax.XColor = 'k'; ax.YColor = 'k'; ax.LineWidth = 1;
    xlim(xl); ylim(yl);
    ax.XTick = [-2, -1];
    ax.YTick = [-2, -1];
    ax.TickLength = [0.015, 0.015];
    ax.FontSize = 10;
    ax.TickLabelInterpreter = 'latex';
end