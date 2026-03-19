function [resultModularSSWM] = simulateModularSSWM(simParams)
% simulateModularSSWM - Simulate modular FGM under Strong Selection Weak Mutation.
%
% Description:
%   Runs stochastic simulations of adaptive evolution in a two-module Fisher's
%   Geometric Model (FGM) under the SSWM regime, where mutations fix sequentially
%   and the population remains monomorphic between fixation events.
%
%   The clock runs at the beneficial mutation supply only. Both +delta and -delta
%   displacements are proposed with equal probability; the Kimura fixation formula
%   handles deleterious outcomes. After each fixation, trait values are rounded to
%   the nearest multiple of delta and clamped to <= 0, enforcing the discrete locus
%   structure of the underlying binary genotype model.
%
%   Termination conditions:
%     1. Fitness reaches 0.99 (terminationStatus = 1)
%     2. No beneficial mutations remain (terminationStatus = 0)
%
% Inputs:
%   simParams - Structure containing simulation parameters:
%       .initialAngles     - Vector of initial angles in phenotype space
%       .initialPhenotypes - Matrix of initial phenotype coordinates [nAngles x 2]
%       .popSize           - Population size (N)
%       .mutationRate      - Per-genome mutation rate (U)
%       .deltaTrait        - Mutational step size (delta)
%       .landscapeStdDev   - Fitness landscape width parameter (sigma)
%       .ellipseParams     - Ellipse axes [a1, a2] for anisotropic selection
%       .geneticTargetSize - Genetic target size per module [K1, K2]
%       .numIteration      - Number of replicate simulations
%
% Outputs:
%   resultModularSSWM - Structure containing:
%       .resultTable       - Cell array {nAngles x nIterations}, each cell contains
%                            trajectory matrix [fitness, time, trait1, trait2]
%       .terminationStatus - 1: fitness threshold reached
%                            0: no beneficial mutations remain
%
% Reference:
%   Kim, M., Ardell, S. M., & Kryazhimskiy, S. (2025).
%   "Module-Selection Balance in the Evolution of Modular Organisms."
%
% See also: simulateModularCM, predictModularSSWM, initializeSimParams
%
% Copyright (c) 2025 Minkyu Kim, Cornell University
% Licensed under MIT License

    resultModularSSWM.resultTable = cell(length(simParams.initialAngles), simParams.numIteration);
    resultModularSSWM.terminationStatus = zeros(length(simParams.initialAngles), simParams.numIteration);

    allInitialPhenotypes = simParams.initialPhenotypes;
    popSize              = simParams.popSize;
    mutationRate         = simParams.mutationRate;  % U (genome-wide)
    deltaTrait           = simParams.deltaTrait;
    landscapeStdDev      = simParams.landscapeStdDev;
    ellipseParams        = simParams.ellipseParams;
    geneticTargetSize    = simParams.geneticTargetSize;  % [K1, K2]

    % Per-locus rate: mu = U / (2*K_i) for each direction
    % Beneficial supply per module i: mu * b_i = U/(2*K_i) * max(-x_i,0)/delta
    % = mutationRate * max(-x_i, 0) / (2 * deltaTrait * K_i)
    invTwoK1 = 1 / (2 * deltaTrait * geneticTargetSize(1));
    invTwoK2 = 1 / (2 * deltaTrait * geneticTargetSize(2));

    for i_pos = 1:length(simParams.initialAngles)
        initialPhenotypes = allInitialPhenotypes(i_pos, 1:end);

        temporaryTable    = cell(1, simParams.numIteration);
        terminationStatus = zeros(1, simParams.numIteration);

        parfor i_repeat = 1:simParams.numIteration
            initialPhenotypeSlice = initialPhenotypes;
            ellipseParamsSlice    = ellipseParams;
            geneticTargetSizeSlice = geneticTargetSize;

            Fitness    = -((initialPhenotypeSlice(1)/ellipseParamsSlice(1))^2 + ...
                           (initialPhenotypeSlice(2)/ellipseParamsSlice(2))^2);
            expFitness = exp(Fitness / (2*landscapeStdDev^2));

            EvolutionLog = struct('MetaInfo', zeros(1,4), 'FixationRecords', []);
            EvolutionLog.FixationRecords    = [expFitness, 1, initialPhenotypeSlice(1), initialPhenotypeSlice(2)];
            EvolutionLog.MetaInfo(1, 1:4)   = [initialPhenotypeSlice(1), initialPhenotypeSlice(2), expFitness, 1];

            simulationEnd        = false;
            reachedFitness       = 0;
            t                    = 1;
            changeTime           = 1;
            finalFitnessThreshold = 0.99;

            while ~simulationEnd
                x1cur = EvolutionLog.MetaInfo(changeTime, 1);
                x2cur = EvolutionLog.MetaInfo(changeTime, 2);

                %% Beneficial supply: mu * b_i = U/(2*K_i) * max(-x_i,0)/delta
                beneficialRate1 = mutationRate * max(-x1cur, 0) * invTwoK1;
                beneficialRate2 = mutationRate * max(-x2cur, 0) * invTwoK2;
                totalBeneficialRate = beneficialRate1 + beneficialRate2;

                if totalBeneficialRate <= 0
                    reachedFitness = 0;
                    break;
                end

                %% Waiting time drawn from beneficial supply
                waitingTime = floor(exprnd(1 / (popSize * totalBeneficialRate)));
                t = t + waitingTime;

                %% Choose which module mutates proportional to beneficial supply
                p_trait1 = beneficialRate1 / totalBeneficialRate;
                if mybinornd(1, p_trait1) == 1
                    mutDim = 1;
                else
                    mutDim = 2;
                end

                %% Propose ±delta with equal probability
                % Clock runs at beneficial rate; occasional deleterious proposals handled by Kimura. 
                % Deleterious fixations are rare under SSWM and negligible.
                direction = 2 * (rand < 0.5) - 1;
                newPhenotype    = [x1cur, x2cur];
                newPhenotype(mutDim) = newPhenotype(mutDim) + direction * deltaTrait;

                %% Round to nearest delta-multiple and clamp to <= 0
                % Enforces discrete locus structure after fixation
                newPhenotype(mutDim) = -deltaTrait * round(-newPhenotype(mutDim) / deltaTrait);
                newPhenotype(mutDim) = min(newPhenotype(mutDim), 0);

                newFitness    = -((newPhenotype(1)/ellipseParamsSlice(1))^2 + ...
                                  (newPhenotype(2)/ellipseParamsSlice(2))^2);
                expFitnessNew = exp(newFitness / (2*landscapeStdDev^2));
                s = log(expFitnessNew) - log(EvolutionLog.MetaInfo(changeTime, 3));

                Pr_fix  = (1 - exp(-2*s)) / (1 - exp(-2*popSize*s));
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

        resultModularSSWM.resultTable(i_pos, :)        = temporaryTable;
        resultModularSSWM.terminationStatus(i_pos, :)  = terminationStatus;
    end
end
