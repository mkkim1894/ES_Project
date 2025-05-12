function populationMatrices = preRunSimulation(initialParameters, simParams, numSteps, varargin)
% preRunSimulation - Simulates a steady-state additive model with fitness and phenotype normalization.
%
% Inputs:
%   initialParameters - A matrix where each row contains [U1, U2, s1, s2] for each phenotype.
%   simParams - Structure containing simulation parameters (e.g., population size, mutation rate).
%   numSteps - Number of discrete time steps to simulate.
%   varargin - Optional argument to specify a random seed.
%              Example: preRunSimulation(initialParameters, simParams, numSteps, 'rngSeed', 42)
%
% Outputs:
%   populationMatrices - A cell array where each cell contains the population state for a phenotype.
%                        Columns: [Trait1, Trait2, Fitness, U1, U2]

    % Parse optional arguments
    p = inputParser;
    addParameter(p, 'rngSeed', [], @(x) isnumeric(x) && isscalar(x)); % Optional RNG seed
    parse(p, varargin{:});
    rngSeed = p.Results.rngSeed;

    % Set RNG seed if specified
    if ~isempty(rngSeed)
        rng(rngSeed);
    end

    % Initialize cell array to store population matrices for each phenotype
    numPhenotypes = size(initialParameters, 1);
    populationMatrices = cell(numPhenotypes, 1);

    % Loop through each initial phenotype
    for i_pos = 1:numPhenotypes
        % Extract precomputed parameters
        U1 = initialParameters(i_pos, 1); % Mutation rate for Trait 1
        U2 = initialParameters(i_pos, 2); % Mutation rate for Trait 2
        s1 = initialParameters(i_pos, 3); % Selection coefficient for Trait 1
        s2 = initialParameters(i_pos, 4); % Selection coefficient for Trait 2
        
        % Initialize population matrix
        popSize = simParams.popSize;
        populationMatrix = zeros(popSize, 7); % Columns: [Trait1, Trait2, Fitness, U1, U2, MutCount1, MutCount2]

        % Initialize population state
        WT = simParams.initialPhenotypes(i_pos, :);
        populationMatrix(:, 1) = WT(1); % Trait 1 values
        populationMatrix(:, 2) = WT(2); % Trait 2 values
        populationMatrix(:, 3) = simParams.initialFitness; % Initial fitness
        populationMatrix(:, 4) = U1;   % Mutation rate for Trait 1
        populationMatrix(:, 5) = U2;   % Mutation rate for Trait 2
        populationMatrix(:, 6) = 0;    % Mutation count for Trait 1
        populationMatrix(:, 7) = 0;    % Mutation count for Trait 2

        % Simulate for the specified number of steps
        for t = 1:numSteps
            %% Mutation on Trait 1
            numMutantTrait1 = poissrnd(sum(populationMatrix(:, 4))); % Total mutations on Trait 1
            if numMutantTrait1 > 0
                mutID_trait1 = datasample(1:popSize, numMutantTrait1, 'Weights', populationMatrix(:, 4), 'Replace', false);
                populationMatrix(mutID_trait1, 1) = populationMatrix(mutID_trait1, 1) + simParams.deltaTrait; % Update Trait 1
                populationMatrix(mutID_trait1, 6) = populationMatrix(mutID_trait1, 6) + 1; % Track mutation count for Trait 1
                populationMatrix(mutID_trait1, 3) = populationMatrix(mutID_trait1, 3) * exp(s1); % Update fitness by s1
            end
            
            %% Mutation on Trait 2
            numMutantTrait2 = poissrnd(sum(populationMatrix(:, 5))); % Total mutations on Trait 2
            if numMutantTrait2 > 0
                mutID_trait2 = datasample(1:popSize, numMutantTrait2, 'Weights', populationMatrix(:, 5), 'Replace', false);
                populationMatrix(mutID_trait2, 2) = populationMatrix(mutID_trait2, 2) + simParams.deltaTrait; % Update Trait 2
                populationMatrix(mutID_trait2, 7) = populationMatrix(mutID_trait2, 7) + 1; % Track mutation count for Trait 2
                populationMatrix(mutID_trait2, 3) = populationMatrix(mutID_trait2, 3) * exp(s2); % Update fitness by s2
            end
            
            %% Recombination
            if simParams.recombinationRate > 0
                % Number of recombination pairs
                numPairs = floor(simParams.recombinationRate * popSize / 2);
                recombIndices = randperm(popSize, 2 * numPairs); % Randomly pick individuals for recombination
                
                % Swap Trait 2 values and mutation counts between pairs
                for pairIdx = 1:numPairs
                    ind1 = recombIndices(pairIdx);
                    ind2 = recombIndices(pairIdx + numPairs);
                    
                    % Swap Trait 2 values
                    tempTrait2 = populationMatrix(ind1, 2);
                    populationMatrix(ind1, 2) = populationMatrix(ind2, 2);
                    populationMatrix(ind2, 2) = tempTrait2;
                    
                    % Swap Mutation counts for Trait 2
                    tempMutCount2 = populationMatrix(ind1, 7);
                    populationMatrix(ind1, 7) = populationMatrix(ind2, 7);
                    populationMatrix(ind2, 7) = tempMutCount2;

                    % Recompute fitness after recombination
                    populationMatrix(ind1, 3) = simParams.initialFitness * exp(s1 * populationMatrix(ind1, 6) + s2 * populationMatrix(ind1, 7));
                    populationMatrix(ind2, 3) = simParams.initialFitness * exp(s1 * populationMatrix(ind2, 6) + s2 * populationMatrix(ind2, 7));
                end
            end
            
            %% Selection
            fitnessProbs = populationMatrix(:, 3) / sum(populationMatrix(:, 3)); % Fitness-based reproduction
            numOffsprings = mnrnd(popSize, fitnessProbs);
            newPopulationMatrix = zeros(size(populationMatrix));
            idx = 1;
            for i = 1:popSize
                numOffspring = numOffsprings(i);
                if numOffspring > 0
                    newPopulationMatrix(idx:idx + numOffspring - 1, :) = repmat(populationMatrix(i, :), numOffspring, 1);
                    idx = idx + numOffspring;
                end
            end
            populationMatrix = newPopulationMatrix(1:popSize, :); % Update population
        end

        %% Renormalize fitness to initial fitness
        meanFitness = mean(populationMatrix(:, 3));
        populationMatrix(:, 3) = populationMatrix(:, 3) / meanFitness * simParams.initialFitness;

        %% Center phenotypes around initial phenotype
        meanTrait1 = mean(populationMatrix(:, 1));
        meanTrait2 = mean(populationMatrix(:, 2));
        populationMatrix(:, 1) = populationMatrix(:, 1) - (meanTrait1 - WT(1));
        populationMatrix(:, 2) = populationMatrix(:, 2) - (meanTrait2 - WT(2));

        %% Trim unnecessary columns before storing
        populationMatrices{i_pos} = populationMatrix(:, 1:5);
    end
end
