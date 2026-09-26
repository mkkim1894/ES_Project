function [resultModularSSWMConst] = simulateModularSSWM_ConstantSupply(simParams)
% simulateModularSSWM_ConstantSupply - Modular FGM under SSWM with a constant
%                                      supply of module-improving mutations.
%
% Description:
%   Identical to simulateModularSSWM in every respect except one: the supply of
%   module-improving mutations does not depend on the module's performance.
%
%   In simulateModularSSWM the number of loci that can still improve module i is
%       b_i(x_i) = |x_i| / delta,
%   so the beneficial supply per genome per generation is
%       U_i(x_i) = mu_i * b_i(x_i) = U * |x_i| / (2 * delta * L_i)      (declining)
%   with per-locus rate mu_i = U / (2*L_i).
%
%   Here a fixed fraction f_i of the L_i loci of module i is module-improving,
%       b_i = f_i * L_i,
%   giving
%       U_i = mu_i * f_i * L_i = U * f_i / 2                            (constant)
%   independent of x_i, of the initial condition, and of L_i.
%
%   Note that this deliberately severs the locus bookkeeping of the declining
%   model. There, |x_i| = b_i * delta by definition, so b_i cannot be held fixed:
%   a fixed b_i would pin |x_i| to a single value. Severing that link is the
%   assumption under test, not an incidental modelling choice. f_1 and f_2 are
%   properties of the genotype-phenotype map: they differ between the two modules
%   but do not vary across the trait space, across initial conditions, or over the
%   course of evolution.
%
%   The supply is therefore constant EVERYWHERE in trait space, including at and
%   beyond the optimum, and is never switched off. What stops a module at its
%   optimum is selection, not the supply: a +delta step proposed at x_i = 0 has
%   s < 0 and is rejected by the Kimura factor below (Pr_fix ~ 3e-14 for N = 1e4,
%   delta = 0.1, sigma = 2). No special case is needed, and none is used. The only
%   cost is that once a module has reached its optimum, the proposals it continues
%   to receive never fix, so the loop performs roughly (U_1+U_2)/U_1 times as many
%   iterations over the remainder of the run.
%
%   This is the control model used as a control. It
%   corresponds to a genotype-phenotype map in which the optimal module state has
%   no greater genotypic redundancy than any other state, so module-improving
%   mutations remain equally available throughout.
%
%   Consequence (derived in predictModularSSWM_ConstantSupply): the per-module
%   fixation flux becomes linear rather than quadratic in |x_i|,
%
%       dx_i/dt = -beta_i * x_i,   beta_i = N * U * f_i * delta^2 / (sigma^2 * a_i^2),
%
%   so each trait decays EXPONENTIALLY and
%
%       log( x_2(t)/x_1(t) ) = log( x_2(0)/x_1(0) ) - (beta_2 - beta_1) * t
%
%   drifts linearly and without bound, with
%
%       beta_2 / beta_1 = (f_2 / f_1) * (a_1^2 / a_2^2).
%
%   There is no module-selection balance: the ratio of module performances goes
%   to 0 or to infinity and never forgets its initial condition. Compare with
%   simulateModularSSWM, where the declining supply makes both traits decay as a
%   power law and x_2/x_1 converges to the attractor (L_2*a_2^2)/(L_1*a_1^2), at
%   which s_1 = s_2.
%
%   As in simulateModularSSWM, the initial phenotype is lattice-initialized to the
%   nearest multiple of delta and clamped to <= 0, only +delta mutations are
%   proposed, and the Kimura formula handles both stochastic loss of beneficial
%   mutations and rejection of steps that would overshoot the optimum.
%
%   Termination conditions:
%     1. Fitness reaches 0.99 (terminationStatus = 1)
%     2. No mutations remain (terminationStatus = 0). With a constant supply this
%        cannot occur; the branch is retained as a guard.
%
% Inputs:
%   simParams - Structure containing simulation parameters:
%       .initialAngles      - Vector of initial angles in phenotype space
%       .initialPhenotypes  - Matrix of initial phenotype coordinates [nAngles x 2]
%       .popSize            - Population size (N)
%       .mutationRate       - Per-genome mutation rate (U)
%       .deltaTrait         - Mutational step size (delta)
%       .landscapeStdDev    - Fitness landscape width parameter (sigma)
%       .ellipseParams      - Ellipse axes [a1, a2] for anisotropic selection
%       .geneticTargetSize  - Number of loci per module [L1, L2]
%       .numIteration       - Number of replicate simulations
%       .beneficialFraction - [f1, f2], the fixed fraction of module-improving
%                             loci in each module. Required.
%
% Outputs:
%   resultModularSSWMConst - Structure containing:
%       .resultTable        - Cell array {nAngles x nIterations}, each cell contains
%                             trajectory matrix [fitness, time, trait1, trait2]
%       .terminationStatus  - 1: fitness threshold reached
%                             0: no mutations remain (should not occur)
%       .beneficialFraction - Echo of [f1, f2]
%       .beneficialSupply   - [U_1, U_2], the constant per-genome supplies
%
% Reference:
%   Kim, M., Ardell, S. M., & Kryazhimskiy, S. (2025).
%   "Module-Selection Balance in the Evolution of Modular Organisms."
%
% See also: simulateModularSSWM, simulateModularCM_ConstantSupply,
%           predictModularSSWM_ConstantSupply
%
% Copyright (c) 2025 Minkyu Kim, Cornell University
% Licensed under MIT License

    %% Fixed beneficial fractions
    if ~isfield(simParams, 'beneficialFraction') || isempty(simParams.beneficialFraction)
        error('simulateModularSSWM_ConstantSupply:missingFraction', ...
              'simParams.beneficialFraction = [f1, f2] must be supplied.');
    end
    f = simParams.beneficialFraction(:)';
    if numel(f) ~= 2 || any(f <= 0) || any(f > 1)
        error('simulateModularSSWM_ConstantSupply:badFraction', ...
              'beneficialFraction must be a 2-vector with entries in (0, 1].');
    end

    nAngles = length(simParams.initialAngles);

    resultModularSSWMConst.resultTable        = cell(nAngles, simParams.numIteration);
    resultModularSSWMConst.terminationStatus  = zeros(nAngles, simParams.numIteration);
    resultModularSSWMConst.beneficialFraction = f;

    allInitialPhenotypes = simParams.initialPhenotypes;
    popSize              = simParams.popSize;
    mutationRate         = simParams.mutationRate;  % U (genome-wide)
    deltaTrait           = simParams.deltaTrait;
    landscapeStdDev      = simParams.landscapeStdDev;
    ellipseParams        = simParams.ellipseParams;

    % Constant supply of module-improving mutations, per genome per generation:
    %   U_i = mu_i * f_i * L_i = U * f_i / 2,   mu_i = U / (2*L_i)
    % L_i cancels: the supply depends on the beneficial FRACTION alone.
    beneficialSupply1 = mutationRate * f(1) / 2;
    beneficialSupply2 = mutationRate * f(2) / 2;
    resultModularSSWMConst.beneficialSupply = [beneficialSupply1, beneficialSupply2];

    % Both are fixed for the whole run and for every initial condition, so the
    % total supply and the module-choice probability are computed once, here,
    % rather than inside the loop. That the loop cannot modify them is the
    % structural guarantee that the supply does not vary across the trait space.
    totalBeneficialRate = beneficialSupply1 + beneficialSupply2;
    p_trait1            = beneficialSupply1 / totalBeneficialRate;

    if totalBeneficialRate <= 0
        error('simulateModularSSWM_ConstantSupply:zeroSupply', ...
              'Total mutational supply is zero; check mutationRate and beneficialFraction.');
    end

    for i_pos = 1:nAngles
        initialPhenotypes = allInitialPhenotypes(i_pos, 1:end);

        temporaryTable    = cell(1, simParams.numIteration);
        terminationStatus = zeros(1, simParams.numIteration);

        parfor i_repeat = 1:simParams.numIteration
            % Discretize initial phenotype to delta-lattice and clamp to <= 0
            initialPhenotypeSlice = min(-deltaTrait * round(-initialPhenotypes(1:2) / deltaTrait), 0);
            ellipseParamsSlice    = ellipseParams;

            Fitness    = -((initialPhenotypeSlice(1)/ellipseParamsSlice(1))^2 + ...
                           (initialPhenotypeSlice(2)/ellipseParamsSlice(2))^2);
            expFitness = exp(Fitness / (2*landscapeStdDev^2));

            EvolutionLog = struct('MetaInfo', zeros(1,4), 'FixationRecords', []);
            EvolutionLog.FixationRecords    = [expFitness, 1, initialPhenotypeSlice(1), initialPhenotypeSlice(2)];
            EvolutionLog.MetaInfo(1, 1:4)   = [initialPhenotypeSlice(1), initialPhenotypeSlice(2), expFitness, 1];

            simulationEnd         = false;
            reachedFitness        = 0;
            t                     = 1;
            changeTime            = 1;
            finalFitnessThreshold = 0.99;

            while ~simulationEnd
                x1cur = EvolutionLog.MetaInfo(changeTime, 1);
                x2cur = EvolutionLog.MetaInfo(changeTime, 2);

                %% Waiting time drawn from the constant supply
                waitingTime = floor(exprnd(1 / (popSize * totalBeneficialRate)));
                t = t + waitingTime;

                %% Choose which module mutates, proportional to its constant supply
                if mybinornd(1, p_trait1) == 1
                    mutDim = 1;
                else
                    mutDim = 2;
                end

                %% Propose +delta in that module
                % There is deliberately no test of whether the module is still
                % below its optimum. If it is not, the step has s < 0 and Pr_fix
                % below rejects it. Selection, not the supply, is what stops a
                % module from overshooting.
                direction = 1;
                newPhenotype         = [x1cur, x2cur];
                newPhenotype(mutDim) = newPhenotype(mutDim) + direction * deltaTrait;

                newFitness    = -((newPhenotype(1)/ellipseParamsSlice(1))^2 + ...
                                  (newPhenotype(2)/ellipseParamsSlice(2))^2);
                expFitnessNew = exp(newFitness / (2*landscapeStdDev^2));
                s = log(expFitnessNew) - log(EvolutionLog.MetaInfo(changeTime, 3));

                Pr_fix   = (1 - exp(-2*s)) / (1 - exp(-2*popSize*s));
                Fixation = mybinornd(1, Pr_fix);

                if Fixation == 1
                    changeTime = changeTime + 1;
                    EvolutionLog.MetaInfo(changeTime, 1:2) = newPhenotype;
                    EvolutionLog.MetaInfo(changeTime, 3)   = expFitnessNew;
                    EvolutionLog.MetaInfo(changeTime, 4)   = t;
                    EvolutionLog.FixationRecords = [EvolutionLog.FixationRecords; ...
                        expFitnessNew, t, newPhenotype(1), newPhenotype(2)];

                    if expFitnessNew >= finalFitnessThreshold
                        simulationEnd  = true;
                        reachedFitness = 1;
                    end
                end
            end

            temporaryTable{i_repeat}    = EvolutionLog.FixationRecords;
            terminationStatus(i_repeat) = reachedFitness;
        end

        resultModularSSWMConst.resultTable(i_pos, :)       = temporaryTable;
        resultModularSSWMConst.terminationStatus(i_pos, :) = terminationStatus;
    end
end
