function [resultPleiotropicCM] = simulatePleiotropicCM(simParams, genomeParams)
% simulatePleiotropicCM - Simulate pleiotropic FGM under Concurrent Mutations regime.
%
% Description:
%   Runs Wright-Fisher simulations where each mutation affects both traits with
%   a pleiotropic angle. Multiple mutations can segregate simultaneously.
%   Recombination is implemented as fixed midpoint crossover (single-point),
%   swapping the second half of the genome between randomly paired individuals.
%
%   Termination conditions:
%     1. Mean fitness reaches 0.99 (terminationStatus = 1)
%
% Inputs:
%   simParams    - Structure containing simulation parameters (see initializeSimParams)
%   genomeParams - Structure containing:
%       .genomeTheta    - Vector of pleiotropic angles for each locus
%       .initialGenomes - Matrix of initial allele states [nAngles x nLoci]
%
% Outputs:
%   resultPleiotropicCM - Structure containing:
%       .resultTable       - Cell array of trajectory matrices [fitness, time, x1, x2]
%       .terminationStatus - 1: fitness threshold reached
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
    resultPleiotropicCM.terminationStatus = zeros(length(simParams.initialAngles), simParams.numIteration);

    % Extract simulation parameters
    popSize = simParams.popSize;
    mutationRate = simParams.mutationRate;
    deltaTrait = simParams.deltaTrait;
    landscapeStdDev = simParams.landscapeStdDev;
    ellipseParams = simParams.ellipseParams;
    recombinationRate = simParams.recombinationRate;
    finalFitnessThreshold = 0.99;

    % Extract genome parameters
    numLoci = length(genomeParams.genomeTheta);
    genomeTheta = genomeParams.genomeTheta;
    initialGenomes = genomeParams.initialGenomes;

    % Precompute allele effects matrix (numLoci x 2)
    alleleEffects = -deltaTrait * [cos(genomeTheta); sin(genomeTheta)]';

    % Precompute crossover point for recombination
    crossoverPoint = floor(numLoci / 2);

    for i_pos = 1:length(simParams.initialAngles)

        temporaryTable = cell(1, simParams.numIteration);
        terminationStatus = zeros(1, simParams.numIteration);

        parfor i_repeat = 1:simParams.numIteration
            initialGenomeSlice = initialGenomes;
            ellipseParamsSlice = ellipseParams;
            alleleEffectsSlice = alleleEffects;

            % Initialize population: all individuals share the same starting genome
            genomeMatrix = repmat(initialGenomeSlice(i_pos, :), popSize, 1);

            % Compute initial phenotypes and fitness
            currentPhenotypes = genomeMatrix * alleleEffectsSlice;
            currentFitness = computeFitness(currentPhenotypes, ellipseParamsSlice, landscapeStdDev);

            % Initialize evolution log
            meanPhenotype = mean(currentPhenotypes, 1);
            meanFitness = mean(currentFitness);
            EvolutionLog = [meanFitness, 1, meanPhenotype(1), meanPhenotype(2)];

            simulationEnd = false;
            reachedFitness = 0;
            t = 1;

            while ~simulationEnd
                t = t + 1;

                %% Mutation Phase
                totalMutations = poissrnd(popSize * mutationRate);

                if totalMutations > 0
                    mutIndividuals = randi(popSize, totalMutations, 1);
                    mutLoci = randi(numLoci, totalMutations, 1);

                    % Accumulate phenotypic effects per individual
                    mutEffects = zeros(popSize, 2);

                    for idx = 1:totalMutations
                        ind = mutIndividuals(idx);
                        loc = mutLoci(idx);
                        oldAllele = genomeMatrix(ind, loc);

                        % Effect depends on current allele state
                        if oldAllele == 0
                            effect = alleleEffectsSlice(loc, :);   % [-delta*cos, -delta*sin]
                        else
                            effect = -alleleEffectsSlice(loc, :);  % [+delta*cos, +delta*sin]
                        end

                        mutEffects(ind, :) = mutEffects(ind, :) + effect;
                        genomeMatrix(ind, loc) = 1 - oldAllele;
                    end

                    currentPhenotypes = currentPhenotypes + mutEffects;
                end

                % Update fitness after mutations
                currentFitness = computeFitness(currentPhenotypes, ellipseParamsSlice, landscapeStdDev);

                %% Recombination Phase (fixed midpoint crossover)
                if recombinationRate > 0
                    if recombinationRate == 1
                        numPairs = floor(popSize / 2);
                    else
                        numPairs = binornd(floor(popSize / 2), recombinationRate);
                    end

                    if numPairs > 0
                        recombIndices = randperm(popSize, 2 * numPairs);
                        parent1 = recombIndices(1:numPairs);
                        parent2 = recombIndices(numPairs+1:end);

                        % Vectorized swap of loci after the crossover point
                        temp = genomeMatrix(parent1, crossoverPoint+1:end);
                        genomeMatrix(parent1, crossoverPoint+1:end) = genomeMatrix(parent2, crossoverPoint+1:end);
                        genomeMatrix(parent2, crossoverPoint+1:end) = temp;

                        % Recompute phenotypes and fitness for recombined individuals
                        recombined = [parent1, parent2];
                        currentPhenotypes(recombined, :) = genomeMatrix(recombined, :) * alleleEffectsSlice;
                        currentFitness(recombined) = computeFitness( ...
                            currentPhenotypes(recombined, :), ellipseParamsSlice, landscapeStdDev);
                    end
                end

                %% Selection Phase (Wright-Fisher sampling)
                selectionProbs = currentFitness / sum(currentFitness);
                numOffsprings = mnrnd(popSize, selectionProbs);

                parentIndices = repelem(1:popSize, numOffsprings);
                genomeMatrix = genomeMatrix(parentIndices(1:popSize), :);

                % Recompute phenotypes and fitness for new generation
                currentPhenotypes = genomeMatrix * alleleEffectsSlice;
                currentFitness = computeFitness(currentPhenotypes, ellipseParamsSlice, landscapeStdDev);

                % Record population mean state
                meanPhenotype = mean(currentPhenotypes, 1);
                meanFitness = mean(currentFitness);
                EvolutionLog = [EvolutionLog; meanFitness, t, meanPhenotype(1), meanPhenotype(2)];

                % Check termination: fitness threshold
                if meanFitness >= finalFitnessThreshold
                    simulationEnd = true;
                    reachedFitness = 1;
                end
            end

            temporaryTable{i_repeat} = EvolutionLog;
            terminationStatus(i_repeat) = reachedFitness;
        end

        resultPleiotropicCM.resultTable(i_pos, :) = temporaryTable;
        resultPleiotropicCM.terminationStatus(i_pos, :) = terminationStatus;
    end
end

%% Helper: compute fitness for a set of phenotypes
function w = computeFitness(phenotypes, ellipseParams, landscapeStdDev)
    logF = -((phenotypes(:,1) / ellipseParams(1)).^2 + ...
             (phenotypes(:,2) / ellipseParams(2)).^2);
    w = exp(logF / (2 * landscapeStdDev^2));
end