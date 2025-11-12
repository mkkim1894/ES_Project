function [resultModularSSWM] = simulateModularSSWM( simParams )
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
                p_trait = fractionBeneficial1 / (fractionBeneficial1 + fractionBeneficial2);
                mutDim = mybinornd(1, p_trait) + 1;   % 1 → trait1, 2 → trait2

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
