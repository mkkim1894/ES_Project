function [resultPleiotropicCM] = simulatePleiotropicCM(simParams, genomeParams)
% simulatePleiotropicCM - Simulate pleiotropic FGM under Concurrent Mutations regime.
%
% Description:
%   Runs Wright-Fisher simulations where each mutation affects both traits with
%   a pleiotropic angle. Multiple mutations can segregate simultaneously.
%
% Inputs:
%   simParams    - Structure containing simulation parameters (see initializeSimParams)
%   genomeParams - Structure containing:
%       .genomeTheta    - Vector of pleiotropic angles for each locus
%       .initialGenomes - Matrix of initial allele states [nAngles x nLoci]
%
% Outputs:
%   resultPleiotropicCM - Structure containing:
%       .resultTable - Cell array of trajectory matrices [fitness, time, x1, x2]
%
% Reference:
%   Kim, M., Ardell, S. M., & Kryazhimskiy, S. (2025).
%   "Module-Selection Balance in the Evolution of Modular Organisms."
%
% See also: simulatePleiotropicSSWM, initializeGenomeTheta
%
% Copyright (c) 2025 Minkyu Kim, Cornell University
% Licensed under MIT License

    resultPleiotropicCM.resultTable = cell(length(simParams.initialAngles), simParams.numIteration);

    % Extract simulation parameters
    popSize = simParams.popSize;
    mutationRate = simParams.mutationRate; % Mutation rate per genome per generation
    deltaTrait = simParams.deltaTrait;
    landscapeStdDev = simParams.landscapeStdDev;
    ellipseParams = simParams.ellipseParams;
    recombinationRate = simParams.recombinationRate; % Recombination rate between genomes
    finalFitnessThreshold = 0.99; % Threshold to end the simulation

    % Extract genome parameters
    numLoci = length(genomeParams.genomeTheta);
    genomeTheta = genomeParams.genomeTheta; % Direction of effect for each locus
    initialGenomes = genomeParams.initialGenomes; % Initial genomes for each simulation

    for i_pos = 1:length(simParams.initialAngles)

        temporaryTable = cell(1, simParams.numIteration);

        parfor i_repeat = 1:simParams.numIteration
            initialGenomeSlice = initialGenomes;
            
            ellipseParamsSlice = ellipseParams;
            
            % Initialize genomes
            genomeMatrix = repmat(initialGenomeSlice(i_pos, :), popSize, 1); % popSize x numLoci

            % Initialize phenotypes based on genomes
            % Calculate initial phenotypes for each individual
            % Phenotype = sum over loci of allele * effect size
            % -deltaTrait is applied due to initialization process is done
            % by iteratively flipping 0 to 1 from all 0 genome.
            alleleEffects = - deltaTrait * [cos(genomeTheta); sin(genomeTheta)]'; % numLoci x 2
            currentPhenotypes = (genomeMatrix * alleleEffects); % popSize x 2

            % Calculate initial fitness
            Fitness = -((currentPhenotypes(:, 1) / ellipseParamsSlice(1)).^2 + (currentPhenotypes(:, 2) / ellipseParamsSlice(2)).^2);
            currentFitness = exp(Fitness / (2 * landscapeStdDev^2)); % popSize x 1

            % Initialize EvolutionLog for recording
            EvolutionLog = struct('Time', [], 'MeanPhenotype', [], 'MeanFitness', []);
            EvolutionLog.Time = 1;
            EvolutionLog.MeanPhenotype = mean(currentPhenotypes, 1);
            EvolutionLog.MeanFitness = mean(currentFitness);
                        
            EvolutionLog.FixationRecords = [EvolutionLog.MeanFitness, EvolutionLog.Time, EvolutionLog.MeanPhenotype(1), EvolutionLog.MeanPhenotype(2)];


            % Simulation parameters
            simulationEnd = false;
            t = 1; % Generation counter

            while ~simulationEnd
                t = t + 1;

                %% Mutation Phase
                % Total number of mutations in the population
                totalMutations = poissrnd(popSize * mutationRate);

                if totalMutations > 0
                    % Randomly select individuals to mutate
                    mutatedIndividuals = randi(popSize, totalMutations, 1);

                    % Randomly select loci to mutate
                    mutationLoci = randi(numLoci, totalMutations, 1);

                    % Initialize mutationalEffects to zeros
                    mutationalEffects = zeros(popSize, 2);

                    % Apply mutations and compute mutational effects
                    for idx = 1:totalMutations
                        individual = mutatedIndividuals(idx);
                        locus = mutationLoci(idx);
                        currentAllele = genomeMatrix(individual, locus);

                        % Calculate mutational effect based on current allele
                        if currentAllele == 0
                            % Allele changes from 0 to 1
                            mutationalEffect = [-deltaTrait * cos(genomeTheta(locus)), -deltaTrait * sin(genomeTheta(locus))];
                        else
                            % Allele changes from 1 to 0
                            mutationalEffect = [deltaTrait * cos(genomeTheta(locus)), deltaTrait * sin(genomeTheta(locus))];
                        end
                        
                        % Add mutational effect to individual's total mutational effect
                        mutationalEffects(individual, :) = mutationalEffects(individual, :) + mutationalEffect;

                        % Flip allele
                        genomeMatrix(individual, locus) = 1 - currentAllele;
                    end
                    
                    % Update phenotypes
                    currentPhenotypes = currentPhenotypes + mutationalEffects;
                end
                
                % Calculate fitness
                Fitness = -((currentPhenotypes(:, 1) / ellipseParamsSlice(1)).^2 + ...
                    (currentPhenotypes(:, 2) / ellipseParamsSlice(2)).^2);
                currentFitness = exp(Fitness / (2 * landscapeStdDev^2)); % popSize x 1
                
                %% Recombination Phase
                if recombinationRate > 0
                    % Determine the number of pairs participating in recombination
                    if recombinationRate == 1  % Full recombination
                        numPairs = floor(popSize / 2);
                    else
                        numPairs = binornd(floor(popSize / 2), recombinationRate); % Binomial sampling for pairs
                    end
                    
                    if numPairs > 0
                        % Randomly select individuals to form pairs
                        recombIndices = randperm(popSize, 2 * numPairs);
                        parent1 = recombIndices(1:numPairs);
                        parent2 = recombIndices(numPairs+1:end);
                        % Perform recombination for each pair
                                for i = 1:numPairs
                                    % Fixed midpoint crossover
                                    crossoverPoint = floor(numLoci / 2);
                                                                        
                                    % Randomly choose a crossover point if
                                    % needed
                                    %crossoverPoint = randi([1, numLoci - 1]);
                                    
                                    % Swap loci from the crossover point onward
                                    temp = genomeMatrix(parent1(i), crossoverPoint+1:end);
                                    genomeMatrix(parent1(i), crossoverPoint+1:end) = genomeMatrix(parent2(i), crossoverPoint+1:end);
                                    genomeMatrix(parent2(i), crossoverPoint+1:end) = temp;
                                    
                                    % Update phenotypes for the recombined individuals
                                    newPhenotype1 = genomeMatrix(parent1(i), :) * alleleEffects;
                                    newPhenotype2 = genomeMatrix(parent2(i), :) * alleleEffects;

                                    % Update fitness for the recombined individuals
                                    fitness1 = -((newPhenotype1(1) / ellipseParams(1))^2 + ...
                                        (newPhenotype1(2) / ellipseParams(2))^2);
                                    fitness2 = -((newPhenotype2(1) / ellipseParams(1))^2 + ...
                                        (newPhenotype2(2) / ellipseParams(2))^2);

                                    % Store updated phenotypes and fitness in their respective locations
                                    currentPhenotypes(parent1(i), :) = newPhenotype1;
                                    currentPhenotypes(parent2(i), :) = newPhenotype2;
                                    currentFitness(parent1(i)) = exp(fitness1 / (2 * landscapeStdDev^2));
                                    currentFitness(parent2(i)) = exp(fitness2 / (2 * landscapeStdDev^2));
                                end
                    end
                end

                %% Selection Phase
                % Normalize fitness to get selection probabilities
                selectionProbabilities = currentFitness / sum(currentFitness);
                numOffsprings = mnrnd(popSize, selectionProbabilities);

                % Build parent index list once
                parentIndices = repelem(1:popSize, numOffsprings);

                % Create next generation genomes directly
                newGenomeMatrix = genomeMatrix(parentIndices, :);

                % Update genomeMatrix for the next generation
                genomeMatrix = newGenomeMatrix(1:popSize, :);

                % Recalculate phenotypes for the new generation
                currentPhenotypes = (genomeMatrix * alleleEffects); % popSize x 2

                % Recalculate fitness for the new generation
                Fitness = -((currentPhenotypes(:, 1) / ellipseParamsSlice(1)).^2 + (currentPhenotypes(:, 2) / ellipseParamsSlice(2)).^2);
                currentFitness = exp(Fitness / (2 * landscapeStdDev^2)); % popSize x 1

                % Update EvolutionLog
                EvolutionLog.Time(end+1, 1) = t;
                EvolutionLog.MeanPhenotype(end+1, :) = mean(currentPhenotypes, 1);
                EvolutionLog.MeanFitness(end+1, 1) = mean(currentFitness);                    
                EvolutionLog.FixationRecords = [EvolutionLog.FixationRecords;...
                    EvolutionLog.MeanFitness(end, 1), EvolutionLog.Time(end, 1),...
                    EvolutionLog.MeanPhenotype(end, :)];
                
                % Check termination condition
                if EvolutionLog.MeanFitness(end) >= finalFitnessThreshold
                    simulationEnd = true;
                end
            end

            % Data recording
            temporaryTable{i_repeat} = EvolutionLog.FixationRecords;
        end

        resultPleiotropicCM.resultTable(i_pos, :) = temporaryTable;
    end
end
