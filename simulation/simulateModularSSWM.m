function [resultModularSSWM] = simulateModularSSWM( simParams )
% simulateModularSSWM - Simulate modular FGM under Strong Selection Weak Mutation.
%
% Description:
%   Runs stochastic simulations of adaptive evolution in a two-module Fisher's
%   Geometric Model (FGM) under the SSWM regime, where mutations fix sequentially
%   and the population remains monomorphic between fixation events.
%
% Inputs:
%   simParams - Structure containing simulation parameters:
%       .initialAngles     - Vector of initial angles in phenotype space
%       .initialPhenotypes - Matrix of initial phenotype coordinates [nAngles x 2]
%       .popSize           - Population size (N)
%       .mutationRate      - Per-locus mutation rate (mu)
%       .deltaTrait        - Mutational step size (delta)
%       .landscapeStdDev   - Fitness landscape width parameter (sigma)
%       .ellipseParams     - Ellipse axes [a1, a2] for anisotropic selection
%       .geneticTargetSize - Number of loci per module [K1, K2]
%       .numIteration      - Number of replicate simulations
%
% Outputs:
%   resultModularSSWM - Structure containing:
%       .resultTable - Cell array {nAngles x nIterations}, each cell contains
%                      trajectory matrix [fitness, time, trait1, trait2]
%
% Algorithm:
%   For each replicate:
%   1. Draw waiting time to next mutation from exponential distribution
%   2. Select module to mutate proportional to beneficial mutation supply
%   3. Apply mutation with random direction (+/- delta)
%   4. Calculate fixation probability using Kimura's formula
%   5. If fixed, update population state; repeat until fitness threshold
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

    allInitialPhenotypes = simParams.initialPhenotypes;
    popSize = simParams.popSize;
    mutationRate = simParams.mutationRate;
    deltaTrait = simParams.deltaTrait;
    landscapeStdDev = simParams.landscapeStdDev;
    ellipseParams = simParams.ellipseParams;
    geneticTargetSize = simParams.geneticTargetSize;

    for i_pos = 1:length(simParams.initialAngles)
        initialPhenotypes = allInitialPhenotypes(i_pos, 1:end);

        temporaryTable = cell(1, simParams.numIteration);
        parfor i_repeat = 1:simParams.numIteration
            initialPhenotypeSlice = initialPhenotypes;
            ellipseParamsSlice = ellipseParams;
            geneticTargetSizeSlice = geneticTargetSize;
            
            Fitness = -( (initialPhenotypeSlice(1)/ellipseParamsSlice(1))^2 + (initialPhenotypeSlice(2)/ellipseParamsSlice(2))^2 );
            expFitness = exp( Fitness./(2*landscapeStdDev^2) ); 

            EvolutionLog = struct('MetaInfo', zeros(1, 4), 'FixationRecords', []); 
            EvolutionLog.FixationRecords = [expFitness, 1, initialPhenotypeSlice(1), initialPhenotypeSlice(2)]; % Initial change point (time, trait1, trait2)
                            
            EvolutionLog.MetaInfo(1, 1:4) = [initialPhenotypeSlice(1), initialPhenotypeSlice(2), expFitness, 1]; % Initial values

            simulationEnd = false;
            t = 1;
            changeTime = 1; 
            finalFitnessThreshold = 0.99; % Ends the simulation at this fitness level.

            while ~simulationEnd
                %% Draw a waiting time
                fractionBeneficial1 = abs(EvolutionLog.MetaInfo(changeTime, 1)) * 1/(deltaTrait*geneticTargetSizeSlice(1));
                fractionBeneficial2 = abs(EvolutionLog.MetaInfo(changeTime, 2)) * 1/(deltaTrait*geneticTargetSizeSlice(2));
                waitingTime = floor(exprnd(1 / (popSize * mutationRate * (fractionBeneficial1 + fractionBeneficial2))));
                t = t + waitingTime;

                %% Mutation
                traitDisplacement = zeros(2, 1);

                % Decide which trait mutates
                % p_trait1 = probability of mutating trait 1
                p_trait1 = fractionBeneficial1 / (fractionBeneficial1 + fractionBeneficial2);
                
                % mybinornd(1, p) returns 1 with probability p, 0 with probability (1-p)
                % So if mybinornd returns 1, we mutate trait 1; if 0, we mutate trait 2
                if mybinornd(1, p_trait1) == 1
                    mutDim = 1;  % Mutate trait 1
                else
                    mutDim = 2;  % Mutate trait 2
                end

                % Randomly choose sign (+ or −) for the mutation step
                direction = 2 * (rand < 0.5) - 1;     % +1 or -1

                % Apply the mutation step
                traitDisplacement(mutDim) = deltaTrait * direction;
                
                newPhenotype = EvolutionLog.MetaInfo(changeTime, 1:2) + traitDisplacement';
                newFitness = -( (newPhenotype(1)/ellipseParamsSlice(1))^2 + (newPhenotype(2)/ellipseParamsSlice(2))^2 );            
                expFitness = exp( newFitness./(2*landscapeStdDev^2) ); 
                s = log(expFitness) - log(EvolutionLog.MetaInfo(changeTime, 3));

                Pr_fix = (1 - exp(-2 * s)) / (1 - exp(-2 * popSize * s));
                Fixation = mybinornd(1, Pr_fix);

                if Fixation == 1
                    changeTime = changeTime + 1;
                    EvolutionLog.MetaInfo(changeTime, 1:2) = newPhenotype;
                    EvolutionLog.MetaInfo(changeTime, 3) = expFitness;
                    EvolutionLog.MetaInfo(changeTime, 4) = t;
                    EvolutionLog.FixationRecords = [EvolutionLog.FixationRecords; expFitness, t, newPhenotype]; % Add change point

                    if expFitness >= finalFitnessThreshold
                        simulationEnd = true;
                    end
                end
            end
            temporaryTable{i_repeat} = EvolutionLog.FixationRecords;
        end
        resultModularSSWM.resultTable(i_pos, :) = temporaryTable;
    end
end
