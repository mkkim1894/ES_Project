%% Run_supplementaryFigures.m
% Generates all supplementary figures from existing data.
% No simulations are run. Requires outputs from Run_supplementary,
% Run_ThresholdDAnalysis, Run_modularFGM, and Run_nestedFGM.
%
% Usage:
%   Run_supplementaryFigures                   % all figures, reproduce mode
%   Run_supplementaryFigures('test')           % all figures, test mode
%   Run_supplementaryFigures('reproduce', 1)   % Figure S1 only
%   Run_supplementaryFigures('reproduce', 6)   % Figure S6 only
%
% Figure assignments:
%   S1 — Steady-state CM main validation
%   S2 — CM parameter ranges
%   S3 — CM weight comparisons
%   S4 — Threshold D sensitivity
%   S5 — Nested FGM mutation-effect distributions
%   S6 — Asymmetric nested FGM dynamics (n1=10, n2=20)
%
% Outputs:
%   .eps files in ./results_supplementary/Figures/
%
% Reference:
%   Kim, M., Ardell, S. M., & Kryazhimskiy, S. (2025).
%   "Module-Selection Balance in the Evolution of Modular Organisms."

function Run_supplementaryFigures(mode, figNum)
    if nargin < 1 || isempty(mode)
        mode = 'reproduce';
    end
    if nargin < 2 || isempty(figNum)
        figNum = 0;  % 0 = all
    end

    isTest = strcmpi(mode, 'test');

    %% N follows the same test/reproduce split as Run_supplementary
    N = tern(isTest, 1e3, 1e4);

    %% Setup paths
    thisDir  = fileparts(mfilename('fullpath'));
    projRoot = fileparts(thisDir);
    addpath(fullfile(projRoot, 'simulation_scripts'));
    addpath(fullfile(projRoot, 'analysis_scripts'));
    addpath(fullfile(projRoot, 'utils'));
    addpath(fullfile(projRoot, 'figure_scripts'));

    suppDir    = fullfile(projRoot, 'results_supplementary');
    figuresDir = fullfile(suppDir, 'Figures');
    resultsDir = fullfile(projRoot, tern(isTest, 'results_test', 'results'));

    if ~isfolder(figuresDir), mkdir(figuresDir); end

    fprintf('Mode: %s  (N = %d)\n', mode, N);
    fprintf('Output: %s\n\n', figuresDir);

    %% Locate data files
    steadyFile  = fullfile(suppDir, sprintf('SteadyStateCM_N_%d.mat', N));
    threshFile  = fullfile(resultsDir, 'ThresholdD_Analysis.mat');
    cmFiles     = dir(fullfile(resultsDir, 'CM_Asexual', 'ModularFGM_CM_Asexual_*.mat'));
    nestedFiles = dir(fullfile(resultsDir, 'Generalization', 'NestedFGM', 'SSWM', ...
                               'NestedFGM_SSWM_*.mat'));

    runAll = (figNum == 0);

    % --- S1: Steady-state CM main ---
    if runAll || figNum == 1
        fprintf('Generating Figure S1 (CM Main)...\n');
        requireFile(steadyFile, 'Run_supplementary');
        makeFigureS_SteadyStateCM(steadyFile, 'outputDir', figuresDir, 'figureSet', 1);
    end

    % --- S2: CM parameter ranges ---
    if runAll || figNum == 2
        fprintf('Generating Figure S2 (CM Parameters)...\n');
        requireFile(steadyFile, 'Run_supplementary');
        makeFigureS_SteadyStateCM(steadyFile, 'outputDir', figuresDir, 'figureSet', 2);
    end

    % --- S3: CM weight comparisons ---
    if runAll || figNum == 3
        fprintf('Generating Figure S3 (CM Weights)...\n');
        requireFile(steadyFile, 'Run_supplementary');
        makeFigureS_SteadyStateCM(steadyFile, 'outputDir', figuresDir, 'figureSet', 3);
    end

    % --- S4: Threshold D sensitivity ---
    if runAll || figNum == 4
        fprintf('Generating Figure S4 (Threshold D)...\n');
        requireFile(threshFile, 'Run_ThresholdDAnalysis');
        if isempty(cmFiles)
            error('No ModularFGM_CM_Asexual_*.mat found in %s.\nRun Run_modularFGM first.', ...
                fullfile(resultsDir, 'CM_Asexual'));
        end
        [~, newest] = max([cmFiles.datenum]);
        cmFile = fullfile(cmFiles(newest).folder, cmFiles(newest).name);
        makeFigureS_ThresholdDTrajectories(cmFile, threshFile, 'outputDir', figuresDir);
    end

    % --- S5: Nested FGM mutation-effect distributions ---
    if runAll || figNum == 5
        fprintf('Generating Figure S5 (Nested FGM distributions)...\n');
        if isempty(nestedFiles)
            error('No NestedFGM_SSWM_*.mat found in %s.\nRun Run_nestedFGM first.', ...
                fullfile(resultsDir, 'Generalization', 'NestedFGM', 'SSWM'));
        end
        [~, newest] = max([nestedFiles.datenum]);
        nestedFile = fullfile(nestedFiles(newest).folder, nestedFiles(newest).name);
        makeFigureS_NestedFGMDistributions(nestedFile, 'outputDir', figuresDir);
    end

    % --- S6: Asymmetric nested FGM dynamics (n1=10, n2=20) ---
    if runAll || figNum == 6
        fprintf('Generating Figure S6 (Asymmetric Nested FGM)...\n');
        asymDir = fullfile(resultsDir, 'Generalization', 'NestedFGM', 'SSWM');
        asymFiles = dir(fullfile(asymDir, 'NestedFGM_SSWM_*n10-20*.mat'));
        if isempty(asymFiles)
            error('No NestedFGM_SSWM_*n10-20*.mat found in %s.\nRun Run_nestedFGM(''%s'', {}, [10,20]) first.', ...
                asymDir, tern(isTest, 'test', 'reproduce'));
        end
        [~, newest] = max([asymFiles.datenum]);
        tmp = load(fullfile(asymFiles(newest).folder, asymFiles(newest).name), 'simParams');
        % Use outputFile override to save under a distinct name,
        % avoiding any collision with the symmetric Figure_NestedFGM_Generations.eps
        makeFigure5_Generations('NestedFGM', tmp.simParams, tern(isTest, 'test', 'full'), ...
            'outputFile', 'Figure_NestedFGM_Asymmetric');
        % Move from results/Figures/ to supplementary figures directory
        src = fullfile(resultsDir, 'Figures', 'Figure_NestedFGM_Asymmetric.eps');
        dst = fullfile(figuresDir, 'FigureS_NestedFGM_Asymmetric.eps');
        if isfile(src)
            movefile(src, dst);
            fprintf('Moved to %s\n', dst);
        else
            warning('Expected output not found at %s', src);
        end
    end

    fprintf('\nDone. Figures saved to %s\n', figuresDir);
    figFiles = dir(fullfile(figuresDir, '*.eps'));
    for i = 1:length(figFiles)
        fprintf('  - %s\n', figFiles(i).name);
    end
end

function requireFile(fpath, runner)
    if ~isfile(fpath)
        error('Required file not found: %s\nRun %s first.', fpath, runner);
    end
end

function y = tern(cond, a, b)
    if cond; y = a; else; y = b; end
end