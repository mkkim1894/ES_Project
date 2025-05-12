function [resultModularCM] = simulateModularCM( simParams )
    resultModularCM.resultTable = cell(length(simParams.initialAngles), simParams.numIteration);
        
    % Check if pre-run populationMatrices exist
    usePreRun = isfield(simParams, 'populationMatrices') && ~isempty(simParams.populationMatrices);
    if ~usePreRun
        warning('Pre-run population data (simParams.populationMatrices) not found. Using monomorphic initialization.');
    end
    
    allInitialPhenotypes = simParams.initialPhenotypes;
    popSize = simParams.popSize;
    mutationRate = simParams.mutationRate;
    deltaTrait = simParams.deltaTrait;
    landscapeStdDev = simParams.landscapeStdDev;
    ellipseParams = simParams.ellipseParams;
    geneticTargetSize = simParams.geneticTargetSize;
    recombinationRate = simParams.recombinationRate;
    finalFitnessThreshold = 0.99; % Threshold to end the simulation
    
    for i_pos = 1:length(simParams.initialAngles)
        if usePreRun
            % Use pre-run population matrix
            basePopulationMatrix = simParams.populationMatrices{i_pos};
        else
            % Monomorphic initialization
            initialPhenotype = allInitialPhenotypes(i_pos, 1:end);
            Fitness = -( (initialPhenotype(1) / ellipseParams(1))^2 + (initialPhenotype(2) / ellipseParams(2))^2 );
            expFitness = exp(Fitness / (2 * landscapeStdDev^2));
            
            basePopulationMatrix = zeros(popSize, 5);            
            basePopulationMatrix(:, 1) = initialPhenotype(1);
            basePopulationMatrix(:, 2) = initialPhenotype(2);
            basePopulationMatrix(:, 3) = expFitness;
            basePopulationMatrix(:, 4) = mutationRate * abs( basePopulationMatrix(:, 1) ) .* ( 1/(deltaTrait*geneticTargetSize(1)) );
            basePopulationMatrix(:, 5) = mutationRate * abs( basePopulationMatrix(:, 2) ) .* ( 1/(deltaTrait*geneticTargetSize(2)) );
        end

        temporaryTable = cell(1, simParams.numIteration);
        
        parfor i_repeat = 1:simParams.numIteration
            ellipseParamsSlice = ellipseParams;
            geneticTargetSizeSlice = geneticTargetSize;
            
            % Copy initial population matrix for this repeat
            populationMatrix = basePopulationMatrix;
                
            % Initialize EvolutionLog for recording
            EvolutionLog = struct('FixationRecords', []); 
            meanFitness = mean(populationMatrix(:, 3));
            meanPhenotype = [mean(populationMatrix(:, 1)), mean(populationMatrix(:, 2))];
            EvolutionLog.FixationRecords = [meanFitness, 1, meanPhenotype];

            simulationEnd = false;
            t = 1;
            while ~simulationEnd
                %% Mutation        
                numMutantTrait1 = poissrnd( sum(populationMatrix(1:end, 4)) ); 
                if numMutantTrait1 > 0
                    mutID_trait1 = datasample(1:popSize, numMutantTrait1, 'Weights', populationMatrix(1:end, 4), 'Replace', false);
                    mutMatrix1 = populationMatrix(mutID_trait1, 1:end);
                    mutMatrix1(1:end, 1) = mutMatrix1(1:end, 1) + deltaTrait;
                    
                    newFitness1 = -( (mutMatrix1(1:end, 1)./ellipseParamsSlice(1)).^2 + (mutMatrix1(1:end, 2)./ellipseParamsSlice(2)).^2 );
                    mutMatrix1(1:end, 3) = exp( newFitness1./(2*landscapeStdDev^2) ); 
                    mutMatrix1(1:end, 4) = mutationRate * abs( mutMatrix1(1:end, 1) ) .* ( 1/(deltaTrait*geneticTargetSizeSlice(1)) );

                    populationMatrix(mutID_trait1, 1:end) = mutMatrix1;
                end
                
                numMutantTrait2 = poissrnd( sum(populationMatrix(1:end, 5)) );
                if numMutantTrait2 > 0
                    mutID_trait2 = datasample(1:popSize, numMutantTrait2, 'Weights', populationMatrix(1:end, 5), 'Replace', false);
                    mutMatrix2 = populationMatrix(mutID_trait2, 1:end);
                    mutMatrix2(1:end, 2) = mutMatrix2(1:end, 2) + deltaTrait;
                    
                    newFitness2 = -( (mutMatrix2(1:end, 1)./ellipseParamsSlice(1)).^2 + (mutMatrix2(1:end, 2)./ellipseParamsSlice(2)).^2 );
                    mutMatrix2(1:end, 3) = exp( newFitness2./(2*landscapeStdDev^2) ); 
                    mutMatrix2(1:end, 5) = mutationRate * abs( mutMatrix2(1:end, 2) ) .* ( 1/(deltaTrait*geneticTargetSizeSlice(2)) );

                    populationMatrix(mutID_trait2, 1:end) = mutMatrix2;
                end
                
                %% Recombination
                if recombinationRate ~= 0
                    if recombinationRate == 1  % Full recombination
                        numPairs = popSize / 2;
                    else
                        numPairs = binornd(popSize, recombinationRate / 2);
                        numPairs = min(numPairs, floor(popSize / 2));
                    end

                    recombIndices = randperm(popSize, 2 * numPairs);
                    recMatrix = populationMatrix(recombIndices, 1:2);
                    % Swap trait 2 between pairs
                    recMatrix(1:numPairs, 2) = populationMatrix(recombIndices(numPairs+1:end), 2);
                    recMatrix(numPairs+1:end, 2) = populationMatrix(recombIndices(1:numPairs), 2);

                    recombFitness = -( (recMatrix(1:end, 1)./ellipseParamsSlice(1)).^2 + (recMatrix(1:end, 2)./ellipseParamsSlice(2)).^2 );
                    recMatrix(1:end, 3) = exp( recombFitness./(2*landscapeStdDev^2) ); 
                    recMatrix(1:end, 4) = mutationRate * abs( recMatrix(1:end, 1) ) .* 1/(deltaTrait*geneticTargetSizeSlice(1));
                    recMatrix(1:end, 5) = mutationRate * abs( recMatrix(1:end, 2) ) .* 1/(deltaTrait*geneticTargetSizeSlice(2));

                    populationMatrix(recombIndices, 1:end) = recMatrix;
                end
                
                %% Selection
                fitnessProbs = populationMatrix(1:end, 3) / sum(populationMatrix(1:end, 3));
                numOffsprings = mnrnd(popSize, fitnessProbs);
                
                % Reproduce individuals according to their fitness
                newPopulationMatrix = zeros(size(populationMatrix));
                idx = 1;
                for i = 1:popSize
                    numOffspring = numOffsprings(i);
                    if numOffspring > 0
                        newPopulationMatrix(idx:idx + numOffspring - 1, 1:end) = repmat(populationMatrix(i, 1:end), numOffspring, 1);
                        idx = idx + numOffspring;
                    end
                end
                populationMatrix = newPopulationMatrix(1:popSize, 1:end);     
                
                t = t + 1;
                newPhenotype = [mean(populationMatrix(1:popSize, 1)), mean(populationMatrix(1:popSize, 2))];
                newFitness = mean(populationMatrix(1:popSize, 3));
                EvolutionLog.FixationRecords = [EvolutionLog.FixationRecords; newFitness, t, newPhenotype]; 
                
                if newFitness >= finalFitnessThreshold
                    simulationEnd = true;    
                end    
            end
            temporaryTable{i_repeat} = EvolutionLog.FixationRecords;
        end
        resultModularCM.resultTable(i_pos, :) = temporaryTable;
    end
end
