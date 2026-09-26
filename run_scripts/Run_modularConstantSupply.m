function Run_modularConstantSupply(varargin)
% Run_modularConstantSupply  Modular GPFM with a constant mutational supply.
%
% Runs the modular GPFM in which the supply of module-improving mutations does
% NOT decline as a module approaches its optimum, in all three evolutionary
% regimes, and generates the corresponding figure. This is the control model
% used as a control.
%
% The model is identical to Run_modularFGM except that a constant fraction f_i of
% the L_i loci of module i is module-improving, everywhere in trait space and at
% all times, so the beneficial supply is
%
%     U_i = mu_i * f_i * L_i = U * f_i / 2        (mu_i = U / (2*L_i))
%
% rather than the declining U_i = U*|x_i| / (2*delta*L_i) of Run_modularFGM.
% Population size, landscape, anisotropy, step size, genetic target sizes,
% mutation rates, recombination and initial conditions are all unchanged, so any
% difference in outcome is attributable to the supply alone.
%
% Simulation output is written to results/ConstantSupply and figures to
% results/Figures, both inside this directory.
%
% Name-value pairs
%   'mode'                'reproduce' (default) or 'test'. Selects the full
%                         paper-scale parameter set or a fast sanity run.
%   'beneficialFraction'  [f1, f2], the constant fraction of module-improving
%                         loci in each module. Default [0.05, 0.10]. The two
%                         values differ between modules but do not vary across
%                         the trait space, across initial conditions, or over
%                         time.
%   'stages'              Cell array selecting stages to run, any of 'sswm',
%                         'cm_asexual', 'cm_sexual', 'figures'. Default all
%                         four; use this to regenerate the figure without
%                         repeating the simulations.
%
% Examples
%   Run_modularConstantSupply('mode', 'test')
%   Run_modularConstantSupply
%   Run_modularConstantSupply('beneficialFraction', [0.05 0.10])
%   Run_modularConstantSupply('stages', {'figures'})
%
% Dependencies
%   Local to this directory: simulateModularSSWM_ConstantSupply,
%   simulateModularCM_ConstantSupply, predictModularSSWM_ConstantSupply,
%   predictModularCM_ConstantSupply, predictFullRecomb_ConstantSupply,
%   makeFigure_ConstantSupply.
%
%   Reused unchanged from the parent project: initializeSimParams,
%   findInitialPhenotypes (utils/) and computeAverageTrajectory
%   (analysis_scripts/). These are deliberately NOT copied here, so that the
%   constant-supply model and the modular GPFM of Figure 3 are guaranteed to
%   share the same parameter initialization and averaging code.
%
%   Also required: mybinornd, used by the SSWM simulation for parity with
%   simulateModularSSWM. See the note in README.md.
%
%   Toolboxes: Statistics and Machine Learning, Symbolic Math, Parallel
%   Computing (optional, for parfor).
%
% Reference
%   Kim, M., Ardell, S. M., & Kryazhimskiy, S. (2025). Module-Selection
%   Balance in the Evolution of Modular Organisms.

% ----------------------------- Parse inputs -----------------------------
p = inputParser;
addParameter(p, 'mode',               'reproduce',  @ischar);
addParameter(p, 'beneficialFraction', [0.05, 0.10], @(v) isnumeric(v) && numel(v) == 2);
addParameter(p, 'stages', {'sswm', 'cm_asexual', 'cm_sexual', 'figures'}, @iscellstr);
parse(p, varargin{:});

opts   = p.Results;
mode   = lower(opts.mode);
stages = lower(opts.stages);
f      = opts.beneficialFraction(:)';

if ~ismember(mode, {'reproduce', 'test'})
    error('Run_modularConstantSupply:UnknownMode', ...
          'Unknown mode ''%s''. Valid modes are ''reproduce'' and ''test''.', mode);
end
isTest = strcmp(mode, 'test');

validStages = {'sswm', 'cm_asexual', 'cm_sexual', 'figures'};
unknown     = setdiff(stages, validStages);
if ~isempty(unknown)
    error('Run_modularConstantSupply:UnknownStage', ...
          'Unknown stage: %s. Valid stages are %s.', ...
          strjoin(unknown, ', '), strjoin(validStages, ', '));
end

if any(f <= 0) || any(f > 1)
    error('Run_modularConstantSupply:BadFraction', ...
          'beneficialFraction entries must lie in (0, 1].');
end

% --------------------------- Path and dependencies ----------------------
thisDir = fileparts(mfilename('fullpath'));
if isempty(thisDir), thisDir = pwd; end
addpath(thisDir);

% Reuse the parent project's shared utilities rather than duplicating them.
projRoot = fileparts(thisDir);
for d = {'utils', 'analysis_scripts'}
    candidate = fullfile(projRoot, d{1});
    if isfolder(candidate)
        addpath(candidate);
    end
end

required = {'simulateModularSSWM_ConstantSupply', 'simulateModularCM_ConstantSupply', ...
            'predictModularSSWM_ConstantSupply',  'predictModularCM_ConstantSupply', ...
            'predictFullRecomb_ConstantSupply',   'makeFigure_ConstantSupply', ...
            'initializeSimParams', 'findInitialPhenotypes', 'computeAverageTrajectory'};
missing = required(cellfun(@(fn) isempty(which(fn)), required));
if ~isempty(missing)
    error('Run_modularConstantSupply:MissingDependencies', ...
          'Not on the MATLAB path: %s', strjoin(missing, ', '));
end

if isempty(which('mybinornd'))
    error('Run_modularConstantSupply:MissingMybinornd', ...
          ['mybinornd is not on the MATLAB path. It is required here for parity ' ...
           'with simulateModularSSWM, which also calls it, but it is not part of ' ...
           'the ES_Project repository - add the folder that contains it to the ' ...
           'path. See the note in README.md.']);
end

% ------------------------------ Output roots ----------------------------
resultsRoot = fullfile(projRoot, tern(isTest, 'results_test', 'results'));
csRoot      = fullfile(resultsRoot, 'ConstantSupply');

% ------------------------------ Config knobs ----------------------------
% Identical to Run_modularFGM except for beneficialFraction.
%   SSWM: U = 1e-7*(L/K) = 2e-6, mu = U/(2*L) = 5e-9  (matches Run_modularFGM
%         after its SSWM mutationRate fix; the old 1e-7 default gave mu = 2.5e-10)
%   CM:   U = 2e-4*(L/K) = 4e-3 at K = 10, L = 200, so mu = 1e-5
K = 10;   % reference genetic target size
L = 200;  % actual number of loci per module

testCfg.numIteration      = 8;
testCfg.numTimeStamp      = 20;
testCfg.initialAngles     = [atan2(1,6.4), atan2(1,3.2), atan2(1,1.6), ...
                             atan2(1,0.8), atan2(1,0.4), atan2(1,0.2)];
testCfg.mutationRateSlow  = 1e-7 * (L / K);
testCfg.mutationRateFast  = 2e-4 * (L / K);
testCfg.popSize           = 10^4;
testCfg.ellipseRatio      = sqrt(2);
testCfg.deltaTrait        = 0.1;
testCfg.landscapeStdDev   = 2;
testCfg.geneticTargetSize = [L, L];

reproduceCfg.numIteration      = 250;
reproduceCfg.numTimeStamp      = 20;
reproduceCfg.initialAngles     = [atan2(1,6.4), atan2(1,3.2), atan2(1,1.6), ...
                                  atan2(1,0.8), atan2(1,0.4), atan2(1,0.2)];
reproduceCfg.mutationRateSlow  = 1e-7 * (L / K);
reproduceCfg.mutationRateFast  = 2e-4 * (L / K);
reproduceCfg.popSize           = 10^4;
reproduceCfg.ellipseRatio      = sqrt(2);
reproduceCfg.deltaTrait        = 0.1;
reproduceCfg.landscapeStdDev   = 2;
reproduceCfg.geneticTargetSize = [L, L];

C = tern(isTest, testCfg, reproduceCfg);
C.beneficialFraction = f;

a = ellipseAxes(C.ellipseRatio);

fprintf('\n=== Modular GPFM, constant mutational supply (%s mode) ===\n', mode);
fprintf('Output directory: %s\n', resultsRoot);
fprintf('Beneficial fractions: f1 = %.4g, f2 = %.4g\n', f(1), f(2));
fprintf('Predicted SSWM rate ratio beta2/beta1 = (f2/f1)*(a1^2/a2^2) = %.4g\n', ...
    (f(2)/f(1)) * (a(1)^2/a(2)^2));
if abs((f(2)/f(1)) - (a(2)^2/a(1)^2)) < 1e-12
    fprintf(['NOTE: f2/f1 equals a2^2/a1^2, so beta1 = beta2 and log(x2/x1) is\n' ...
             '      frozen at its initial value. That is a knife-edge, not a\n' ...
             '      balance: each initial condition keeps its own ratio.\n']);
end
fprintf('\n');

% ---------------------------------- SSWM --------------------------------
if ismember('sswm', stages)
    fprintf('Running ModularConstSupply SSWM...\n');
    tic;
    simParams = baseParams(C, C.mutationRateSlow, 0);
    resultModularSSWMConst = simulateModularSSWM_ConstantSupply(simParams); %#ok<NASGU>
    ave = computeAverageTrajectory(C.numTimeStamp, simParams, resultModularSSWMConst.resultTable);
    fprintf('  SSWM completed in %.2f s\n', toc);
    reportTermination('SSWM', resultModularSSWMConst.terminationStatus);

    outDir = ensureDir(csRoot, 'SSWM');
    fname  = fullfile(outDir, buildFilename('ModularConstSupplyFGM', 'SSWM', simParams));
    save(fname, 'simParams', 'resultModularSSWMConst', 'ave');
    [analyticalTrajectories, decayRates] = predictModularSSWM_ConstantSupply( ...
        discretizeInitialPhenotypes(simParams), ave); %#ok<ASGLU>
    save(fname, 'analyticalTrajectories', 'decayRates', '-append');
    reportDivergence(decayRates);
    fprintf('  Saved: %s\n', fname);
end

% ------------------------------ CM asexual ------------------------------
if ismember('cm_asexual', stages)
    fprintf('Running ModularConstSupply CM Asexual...\n');
    tic;
    simParams = baseParams(C, C.mutationRateFast, 0);
    resultModularCMConst = simulateModularCM_ConstantSupply(simParams); %#ok<NASGU>
    ave = computeAverageTrajectory(C.numTimeStamp, simParams, resultModularCMConst.resultTable);
    fprintf('  CM Asexual completed in %.2f s\n', toc);
    reportTermination('CM Asexual', resultModularCMConst.terminationStatus);

    outDir = ensureDir(csRoot, 'CM_Asexual');
    fname  = fullfile(outDir, buildFilename('ModularConstSupplyFGM', 'CM_Asexual', simParams));
    save(fname, 'simParams', 'resultModularCMConst', 'ave');
    analyticalTrajectories = predictModularCM_ConstantSupply(discretizeInitialPhenotypes(simParams)); %#ok<NASGU>
    save(fname, 'analyticalTrajectories', '-append');
    fprintf('  Saved: %s\n', fname);
end

% ------------------- CM sexual (full reassortment) ----------------------
if ismember('cm_sexual', stages)
    fprintf('Running ModularConstSupply CM Sexual (full recombination)...\n');
    tic;
    simParams = baseParams(C, C.mutationRateFast, 1);
    resultModularCMConst = simulateModularCM_ConstantSupply(simParams); %#ok<NASGU>
    ave = computeAverageTrajectory(C.numTimeStamp, simParams, resultModularCMConst.resultTable);
    fprintf('  CM Sexual completed in %.2f s\n', toc);
    reportTermination('CM Sexual', resultModularCMConst.terminationStatus);

    outDir = ensureDir(csRoot, 'CM_Sexual');
    fname  = fullfile(outDir, buildFilename('ModularConstSupplyFGM', 'CM_Sexual', simParams));
    save(fname, 'simParams', 'resultModularCMConst', 'ave');
    analyticalTrajectories = predictFullRecomb_ConstantSupply(discretizeInitialPhenotypes(simParams), ave); %#ok<NASGU>
    save(fname, 'analyticalTrajectories', '-append');
    fprintf('  Saved: %s\n', fname);
end

% --------------------------------- Figure -------------------------------
if ismember('figures', stages)
    fprintf('Generating figure...\n');
    simParamsRef.popSize            = C.popSize;
    simParamsRef.deltaTrait         = C.deltaTrait;
    simParamsRef.ellipseParams      = a;
    simParamsRef.landscapeStdDev    = C.landscapeStdDev;
    simParamsRef.geneticTargetSize  = C.geneticTargetSize;
    simParamsRef.beneficialFraction = f;

    makeFigure_ConstantSupply(simParamsRef, 'resultsRoot', resultsRoot);
end

fprintf('\nAll done. Mode: %s\n', mode);
end

% ============================== Helpers ==================================

function simParams = baseParams(C, U, rho)
    args = {'numIteration', C.numIteration, ...
            'initialAngles', C.initialAngles, ...
            'popSize', C.popSize, ...
            'ellipseRatio', C.ellipseRatio, ...
            'deltaTrait', C.deltaTrait, ...
            'landscapeStdDev', C.landscapeStdDev, ...
            'geneticTargetSize', C.geneticTargetSize, ...
            'mutationRate', U};
    if rho > 0
        args = [args, {'recombinationRate', rho}];
    end
    simParams = initializeSimParams(args{:});

    % Constant-supply settings are attached after initializeSimParams so that the
    % shared utility needs no modification.
    simParams.beneficialFraction = C.beneficialFraction;
end

function a = ellipseAxes(ellipseRatio)
    if ellipseRatio >= 1
        a = [1, 1/ellipseRatio];
    else
        a = [ellipseRatio, 1];
    end
end

function name = buildFilename(model, regime, sp)
    % Parameter-encoded filename, matching the convention in Run_modularFGM.
    name = sprintf('%s_%s_N%.0e_M%.1e_d%.2f_eR%.2f_s%.2f', ...
        model, regime, sp.popSize, sp.mutationRate, sp.deltaTrait, ...
        sp.ellipseParams(1)/sp.ellipseParams(2), sp.landscapeStdDev);

    if isfield(sp, 'geneticTargetSize')
        name = sprintf('%s_L%d-%d', name, sp.geneticTargetSize(1), sp.geneticTargetSize(2));
    end

    if isfield(sp, 'beneficialFraction')
        name = sprintf('%s_f%.3g-%.3g', name, sp.beneficialFraction(1), sp.beneficialFraction(2));
    end

    if isfield(sp, 'recombinationRate') && sp.recombinationRate > 0
        name = sprintf('%s_R%.4f', name, sp.recombinationRate);
    end

    name = [name, '.mat'];
end

function reportTermination(label, statusMatrix)
    ts = statusMatrix(:);
    nTotal = numel(ts);
    fprintf('  %s termination: %d/%d reached fitness, %d no beneficial muts\n', ...
        label, sum(ts == 1), nTotal, sum(ts == 0));
end

function reportDivergence(decayRates)
    fprintf('  Predicted SSWM drift of log(x2/x1), slope = -(beta2 - beta1):\n');
    fprintf('    beta = [%.4e  %.4e], slope = %+.4e per generation\n', ...
        decayRates.beta(1,1), decayRates.beta(1,2), decayRates.logRatioSlope(1));
    if decayRates.logRatioSlope(1) < 0
        fprintf('    -> log(x2/x1) decreases without bound (x2/x1 -> 0). No balance.\n');
    elseif decayRates.logRatioSlope(1) > 0
        fprintf('    -> log(x2/x1) increases without bound (x2/x1 -> Inf). No balance.\n');
    else
        fprintf('    -> log(x2/x1) frozen at its initial value (knife-edge, still no attractor).\n');
    end
end

function out = ensureDir(base, rel)
    out = fullfile(base, rel);
    if ~exist(out, 'dir')
        mkdir(out);
    end
end

function sp = discretizeInitialPhenotypes(sp)
    delta = sp.deltaTrait;
    x0    = sp.initialPhenotypes;
    x0    = -delta * round(-x0 / delta);
    x0    = min(x0, 0);
    sp.initialPhenotypes = x0;
end

function y = tern(cond, a, b)
    if cond
        y = a;
    else
        y = b;
    end
end
