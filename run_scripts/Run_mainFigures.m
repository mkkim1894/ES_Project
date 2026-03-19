%% Run_mainFigures.m
% Generates main figures (2-5) from existing simulation data.
% No simulations are run.
%
% Usage:
%   Run_mainFigures                        % all figures, reproduce mode
%   Run_mainFigures('test')                % all figures, test mode
%   Run_mainFigures('reproduce', [3 4])    % figures 3 and 4 only
%   Run_mainFigures('reproduce', 2)        % figure 2 only
%
% Figure assignments:
%   Figure 2 — PleiotropicFGM
%   Figure 3 — ModularFGM
%   Figure 4 — DiscordantFGM
%   Figure 5 — NestedFGM (moduleDimension [10,10])
%
% Outputs:
%   .eps files in ./results/Figures/ (or ./results_test/Figures/)
%
% Reference:
%   Kim, M., Ardell, S. M., & Kryazhimskiy, S. (2025).
%   "Module-Selection Balance in the Evolution of Modular Organisms."

function Run_mainFigures(mode, figures)
    if nargin < 1 || isempty(mode)
        mode = 'reproduce';
    end
    if nargin < 2 || isempty(figures)
        figures = [2, 3, 4, 5];
    end

    validFigs = [2, 3, 4, 5];
    if ~all(ismember(figures, validFigs))
        error('figures must be a subset of [2, 3, 4, 5].');
    end

    % ----------------------------------------------------------------
    % Paths
    % ----------------------------------------------------------------
    thisDir  = fileparts(mfilename('fullpath'));
    projRoot = fileparts(thisDir);
    addpath(fullfile(projRoot, 'simulation_scripts'));
    addpath(fullfile(projRoot, 'analysis_scripts'));
    addpath(fullfile(projRoot, 'utils'));
    addpath(fullfile(projRoot, 'figure_scripts'));

    fprintf('=== Main Figures %s (%s mode) ===\n', mat2str(figures), mode);

    % ----------------------------------------------------------------
    % Parameters (used only to construct simParams for file matching,
    % not for running simulations)
    % ----------------------------------------------------------------
    switch lower(mode)
        case 'test'
            isTest = true;
            C.numIteration      = 4;
            C.initialAngles     = [atan2(1,6.4), atan2(1,1.6), atan2(1,0.4)];
            C.mutationRateSlow  = 1e-7;
            C.mutationRateFast  = 2e-4;
            C.popSize           = 1e3;
            C.ellipseRatio      = sqrt(2);
            C.deltaTrait        = 0.1;
            C.landscapeStdDev   = 2;
            C.geneticTargetSize = [10, 10];
            C.discordantAngles  = [pi/16, 7*pi/16];
            C.moduleDimension   = [10, 10];
            C.L_CM              = 400;

        case {'reproduce', 'full'}
            isTest = false;
            C.numIteration      = 250;
            C.initialAngles     = [atan2(1,6.4), atan2(1,3.2), atan2(1,1.6), ...
                                   atan2(1,0.8), atan2(1,0.4), atan2(1,0.2)];
            C.mutationRateSlow  = 1e-7;
            C.mutationRateFast  = 2e-4;
            C.popSize           = 1e4;
            C.ellipseRatio      = sqrt(2);
            C.deltaTrait        = 0.1;
            C.landscapeStdDev   = 2;
            C.geneticTargetSize = [10, 10];
            C.discordantAngles  = [pi/16, 7*pi/16];
            C.moduleDimension   = [10, 10];
            C.L_CM              = 400;

        otherwise
            error('Unknown mode ''%s''. Use ''test'' or ''reproduce''.', mode);
    end

    figMode          = tern(isTest, 'test', 'full');
    proximityCutoff  = C.deltaTrait;
    nestedM          = sqrt(2 * C.deltaTrait);

    % ----------------------------------------------------------------
    % Discordant admissible angles (needed for figure 4 file matching)
    % ----------------------------------------------------------------
    masterAngles    = C.initialAngles;
    thetaDiscordant = C.discordantAngles;

    M_disc = [cos(thetaDiscordant(1)), cos(thetaDiscordant(2));
              sin(thetaDiscordant(1)), sin(thetaDiscordant(2))];

    ellipseParams_disc = tern(C.ellipseRatio >= 1, [1, 1/C.ellipseRatio], [C.ellipseRatio, 1]);

    allInitialPhenotypes = findInitialPhenotypes(masterAngles, ellipseParams_disc, 0.25, C.landscapeStdDev);

    discordantMask = false(1, numel(masterAngles));
    for ii = 1:numel(masterAngles)
        y0 = M_disc \ allInitialPhenotypes(ii, :)';
        discordantMask(ii) = all(y0 <= -C.deltaTrait);
    end
    discordantInitialAngles = masterAngles(discordantMask);
    discordantAngleIdx      = find(discordantMask);

    fprintf('DiscordantFGM: %d / %d initial angles retained.\n', ...
        numel(discordantInitialAngles), numel(masterAngles));

    % ----------------------------------------------------------------
    % Figure 2 — PleiotropicFGM
    % ----------------------------------------------------------------
    if ismember(2, figures)
        fprintf('\nFigure 2 (PleiotropicFGM)...\n'); tic;
        pleiotropicFigParams = initializeSimParams( ...
            'numIteration',      C.numIteration, ...
            'initialAngles',     C.initialAngles, ...
            'popSize',           C.popSize, ...
            'ellipseRatio',      C.ellipseRatio, ...
            'deltaTrait',        C.deltaTrait, ...
            'landscapeStdDev',   C.landscapeStdDev, ...
            'mutationRate',      C.mutationRateFast, ...
            'recombinationRate', 1, ...
            'omitParams',        {'geneticTargetSize'});
        makeFigure2_Generations('PleiotropicFGM', pleiotropicFigParams, figMode, ...
                                'proximityCutoff', proximityCutoff);
        fprintf('  Done in %.2f s\n', toc);
    end

    % ----------------------------------------------------------------
    % Figure 3 — ModularFGM
    % ----------------------------------------------------------------
    if ismember(3, figures)
        fprintf('\nFigure 3 (ModularFGM)...\n'); tic;
        modularFigParams = initializeSimParams( ...
            'numIteration',      C.numIteration, ...
            'initialAngles',     C.initialAngles, ...
            'popSize',           C.popSize, ...
            'ellipseRatio',      C.ellipseRatio, ...
            'deltaTrait',        C.deltaTrait, ...
            'landscapeStdDev',   C.landscapeStdDev, ...
            'geneticTargetSize', C.geneticTargetSize, ...
            'mutationRate',      C.mutationRateFast);
        makeFigure3_Generations('ModularFGM', modularFigParams, figMode, ...
                                'proximityCutoff', proximityCutoff);
        fprintf('  Done in %.2f s\n', toc);
    end

    % ----------------------------------------------------------------
    % Figure 4 — DiscordantFGM
    % ----------------------------------------------------------------
    if ismember(4, figures)
        fprintf('\nFigure 4 (DiscordantFGM)...\n'); tic;
        discordantFigParams = initializeSimParams( ...
            'numIteration',      C.numIteration, ...
            'initialAngles',     discordantInitialAngles, ...
            'popSize',           C.popSize, ...
            'ellipseRatio',      C.ellipseRatio, ...
            'deltaTrait',        C.deltaTrait, ...
            'landscapeStdDev',   C.landscapeStdDev, ...
            'geneticTargetSize', C.geneticTargetSize, ...
            'mutationRate',      C.mutationRateFast, ...
            'discordantAngles',  C.discordantAngles);
        discordantFigParams.masterInitialAngles = masterAngles;
        discordantFigParams.initialAngleIdx     = discordantAngleIdx;
        makeFigure4_Generations('DiscordantFGM', discordantFigParams, figMode, ...
                                'proximityCutoff', proximityCutoff);
        fprintf('  Done in %.2f s\n', toc);
    end

    % ----------------------------------------------------------------
    % Figure 5 — NestedFGM
    % ----------------------------------------------------------------
    if ismember(5, figures)
        fprintf('\nFigure 5 (NestedFGM)...\n'); tic;
        nestedFigParams = initializeSimParams( ...
            'numIteration',    C.numIteration, ...
            'initialAngles',   C.initialAngles, ...
            'popSize',         C.popSize, ...
            'ellipseRatio',    C.ellipseRatio, ...
            'deltaTrait',      nestedM, ...
            'landscapeStdDev', C.landscapeStdDev, ...
            'mutationRate',    C.mutationRateFast);
        nestedFigParams.moduleDimension = C.moduleDimension;
        makeFigure5_Generations('NestedFGM', nestedFigParams, figMode, ...
                                'proximityCutoff', proximityCutoff);
        fprintf('  Done in %.2f s\n', toc);
    end

    fprintf('\nAll requested figures done.\n');
end

function y = tern(cond, a, b)
    if cond; y = a; else; y = b; end
end
