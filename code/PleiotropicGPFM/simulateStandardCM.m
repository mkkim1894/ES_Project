function [resultModularCM] = simulateStandardCM( simParams )
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
            
            basePopulationMatrix = zeros(popSize, 3);            
            basePopulationMatrix(:, 1) = initialPhenotype(1);
            basePopulationMatrix(:, 2) = initialPhenotype(2);
            basePopulationMatrix(:, 3) = expFitness;
        end

        temporaryTable = cell(1, simParams.numIteration);
        
        parfor i_repeat = 1:simParams.numIteration
            ellipseParamsSlice = ellipseParams;
            
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
                numMutant = poissrnd( popSize*mutationRate );                
                if numMutant > 0
                    mutTheta = 2 * pi * rand(numMutant,1);  
                    
                    traitDisplacement = zeros(numMutant,2);
                    traitDisplacement(:, 1) = deltaTrait .* cos(mutTheta);                
                    traitDisplacement(:, 2) = deltaTrait .* sin(mutTheta); 
                    
                    mutID_trait = datasample(1:popSize, numMutant, 'Replace', false);
                    mutMatrix = populationMatrix(mutID_trait, 1:end);
                    mutMatrix(1:end, 1) = mutMatrix(1:end, 1) + traitDisplacement(:, 1);
                    mutMatrix(1:end, 2) = mutMatrix(1:end, 2) + traitDisplacement(:, 2);
                    
                    newFitness1 = -( (mutMatrix(1:end, 1)./ellipseParamsSlice(1)).^2 + (mutMatrix(1:end, 2)./ellipseParamsSlice(2)).^2 );
                    mutMatrix(1:end, 3) = exp( newFitness1./(2*landscapeStdDev^2) ); 

                    populationMatrix(mutID_trait, 1:end) = mutMatrix;
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
