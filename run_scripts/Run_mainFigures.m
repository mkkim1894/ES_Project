%% Run_mainFigures.m
% Generates main figures (2-5) from existing simulation data.
% No simulations are run.
%
% Figures 2-4 are now organized by
% REGIME rather than by MODEL. Each fixes one evolutionary regime and
% shows Pleiotropic GPFM vs. Modular GPFM side by side (2x2), so the
% cross-model contrast — the paper's main point — is visible within a
% single figure instead of requiring the reader to flip between figures.
% Figure 5 (Nested FGM) is unchanged (still model-organized, 2x3;
% generalization check, not part of the core regime contrast).
%
% Usage:
%   Run_mainFigures                        % all figures, reproduce mode
%   Run_mainFigures('test')                % all figures, test mode
%   Run_mainFigures('reproduce', [3 4])    % figures 3 and 4 only
%   Run_mainFigures('reproduce', 2)        % figure 2 only
%
% Figure assignments:
%   Figure 2 — SSWM                (regime-organized: Pleiotropic | Modular)
%   Figure 3 — CM, linked           (regime-organized: Pleiotropic | Modular)
%   Figure 4 — CM, unlinked         (regime-organized: Pleiotropic | Modular)
%   Figure 5 — NestedFGM (moduleDimension [10,10])  [unchanged, model-organized]
%
% Outputs:
%   .pdf files in ./results/Figures/ (or ./results_test/Figures/)
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
    K = 10;   % reference genetic target size; defines per-locus rate mu = U_ref / (2*K)
    L = 200;  % number of loci per module
    switch lower(mode)
        case 'test'
            isTest = true;
            C.numIteration      = 4;
            C.initialAngles     = [atan2(1,6.4), atan2(1,1.6), atan2(1,0.4)];
            C.mutationRateSlow  = 1e-7 * (L / K);
            C.mutationRateFast  = 2e-4 * (L / K);
            C.popSize           = 1e3;
            C.ellipseRatio      = sqrt(2);
            C.deltaTrait        = 0.1;
            C.landscapeStdDev   = 2;
            C.geneticTargetSize = [L, L];
            C.moduleDimension   = [10, 10];
            C.L_CM              = 2*L;

        case {'reproduce', 'full'}
            isTest = false;
            C.numIteration      = 250;
            C.initialAngles     = [atan2(1,6.4), atan2(1,3.2), atan2(1,1.6), ...
                                   atan2(1,0.8), atan2(1,0.4), atan2(1,0.2)];
            C.mutationRateSlow  = 1e-7 * (L / K);
            C.mutationRateFast  = 2e-4 * (L / K);
            C.popSize           = 1e4;
            C.ellipseRatio      = sqrt(2);
            C.deltaTrait        = 0.1;
            C.landscapeStdDev   = 2;
            C.geneticTargetSize = [L, L];
            C.moduleDimension   = [10, 10];
            C.L_CM              = 2*L;

        otherwise
            error('Unknown mode ''%s''. Use ''test'' or ''reproduce''.', mode);
    end

    figMode          = tern(isTest, 'test', 'full');
    % Methods: "we exclude data points where either x_i exceeds -delta".
    % The cutoff is therefore delta exactly, not a fraction of it, and the
    % same value is used in every figure that plots log R.
    proximityCutoff  = C.deltaTrait;
    nestedM          = sqrt(2 * C.deltaTrait);

    % ----------------------------------------------------------------
    % simParams used for file-matching (Pleiotropic / Modular).
    %
    % These are REFERENCE structs only: they never feed a simulation, they just
    % tell locateFile which .mat to load. The mutation rate therefore has to
    % match the rate the data were actually generated at, which differs by
    % regime (slow in SSWM, fast in the concurrent-mutations regimes). Building
    % one shared struct at the fast rate used to leave the SSWM lookup unable to
    % discriminate between files, so it fell through to an alphabetical guess.
    % ----------------------------------------------------------------
    makePleiotropicRef = @(U) initializeSimParams( ...
        'numIteration',      C.numIteration, ...
        'initialAngles',     C.initialAngles, ...
        'popSize',           C.popSize, ...
        'ellipseRatio',      C.ellipseRatio, ...
        'deltaTrait',        C.deltaTrait, ...
        'landscapeStdDev',   C.landscapeStdDev, ...
        'mutationRate',      U, ...
        'recombinationRate', 1, ...
        'omitParams',        {'geneticTargetSize'});

    makeModularRef = @(U) initializeSimParams( ...
        'numIteration',      C.numIteration, ...
        'initialAngles',     C.initialAngles, ...
        'popSize',           C.popSize, ...
        'ellipseRatio',      C.ellipseRatio, ...
        'deltaTrait',        C.deltaTrait, ...
        'landscapeStdDev',   C.landscapeStdDev, ...
        'geneticTargetSize', C.geneticTargetSize, ...
        'mutationRate',      U);

    pleiotropicSlowRef = makePleiotropicRef(C.mutationRateSlow);
    modularSlowRef     = makeModularRef(C.mutationRateSlow);
    pleiotropicFastRef = makePleiotropicRef(C.mutationRateFast);
    modularFastRef     = makeModularRef(C.mutationRateFast);

    % ----------------------------------------------------------------
    % Figure 2 — SSWM (Pleiotropic | Modular)
    % ----------------------------------------------------------------
    if ismember(2, figures)
        fprintf('\nFigure 2 (SSWM, Pleiotropic vs. Modular)...\n'); tic;
        makeRegimeFigure_Generations('SSWM', pleiotropicSlowRef, modularSlowRef, figMode, ...
                                     'proximityCutoff', proximityCutoff, 'outputFile', 'Figure2_SSWM');
        fprintf('  Done in %.2f s\n', toc);
    end

    % ----------------------------------------------------------------
    % Figure 3 — CM, linked chromosomes (Pleiotropic | Modular)
    % ----------------------------------------------------------------
    if ismember(3, figures)
        fprintf('\nFigure 3 (CM linked, Pleiotropic vs. Modular)...\n'); tic;
        makeRegimeFigure_Generations('CM_Asexual', pleiotropicFastRef, modularFastRef, figMode, ...
                                     'proximityCutoff', proximityCutoff, 'outputFile', 'Figure3_CMLinked');
        fprintf('  Done in %.2f s\n', toc);
    end

    % ----------------------------------------------------------------
    % Figure 4 — CM, unlinked chromosomes (Pleiotropic | Modular)
    % ----------------------------------------------------------------
    if ismember(4, figures)
        fprintf('\nFigure 4 (CM unlinked, Pleiotropic vs. Modular)...\n'); tic;
        makeRegimeFigure_Generations('CM_Sexual', pleiotropicFastRef, modularFastRef, figMode, ...
                                     'proximityCutoff', proximityCutoff, 'outputFile', 'Figure4_CMUnlinked');
        fprintf('  Done in %.2f s\n', toc);
    end

    % ----------------------------------------------------------------
    % Figure 5 — NestedFGM (unchanged: model-organized, 2x3)
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