function [resultNestedCM] = simulateNestedCM(simParams)
% simulateNestedCM - Simulate nested FGM under Concurrent Mutations regime.
%
% Description:
%   Runs Wright-Fisher simulations of the nested FGM where each module is itself
%   a multi-dimensional FGM. Supports recombination between modules.
%
%   Mutational effects on x_i are drawn from the full (unconditional) distribution:
%
%       Delta x_i ~ N( -m^2/2 , (m * sqrt(2*|x_i|/n_i))^2 )
%
%   where m = simParams.deltaTrait is the magnitude of the mutation vector
%   in the underlying y-space and n_i is the dimensionality of module i.
%   Both beneficial (Delta x_i > 0) and deleterious (Delta x_i < 0) mutations
%   are generated; selection acts on all of them via Wright-Fisher sampling.
%
%   The total mutation rate per genome per generation is U = simParams.mutationRate,
%   giving an expected number of mutations per module across the population of N*U/2.
%
% Inputs:
%   simParams - Structure containing simulation parameters:
%       .initialAngles     - Vector of initial angles in phenotype space
%       .initialPhenotypes - Matrix of initial phenotype coordinates [nAngles x 2]
%       .popSize           - Population size (N)
%       .mutationRate      - Per-genome mutation rate (U)
%       .deltaTrait        - Mutation vector magnitude m (= sqrt(2*delta))
%       .landscapeStdDev   - Fitness landscape width (sigma)
%       .ellipseParams     - Ellipse axes [a1, a2]
%       .recombinationRate - Recombination rate (rho)
%       .moduleDimension   - Module dimensionalities [n1, n2]
%       .numIteration      - Number of replicate simulations
%
% Outputs:
%   resultNestedCM - Structure containing:
%       .resultTable       - Cell array {nAngles x nIterations}, each cell contains
%                            trajectory matrix [fitness, time, meanTrait1, meanTrait2]
%       .terminationStatus - 1: fitness threshold reached
%
% Reference:
%   Kim, M., Ardell, S. M., & Kryazhimskiy, S. (2025).
%   "Module-Selection Balance in the Evolution of Modular Organisms."
%
% See also: simulateNestedSSWM, simulateModularCM
%
% Copyright (c) 2025 Minkyu Kim, Cornell University
% Licensed under MIT License

    resultNestedCM.resultTable = cell(length(simParams.initialAngles), simParams.numIteration);
    resultNestedCM.terminationStatus = zeros(length(simParams.initialAngles), simParams.numIteration);

    usePreRun = isfield(simParams, 'populationMatrices') && ~isempty(simParams.populationMatrices);
    if usePreRun
        warning(['simulateNestedCM: supplied populationMatrices will be ignored because ' ...
                 'they are not guaranteed to be consistent with the nested-FGM beneficial-supply approximation.']);
    end

    allInitialPhenotypes = simParams.initialPhenotypes;
    popSize = simParams.popSize;
    mutationRate = simParams.mutationRate;     % Per-module mutation opportunity rate U
    m = simParams.deltaTrait;                  % Nested FGM mutation-vector magnitude
    landscapeStdDev = simParams.landscapeStdDev;
    ellipseParams = simParams.ellipseParams;
    recombinationRate = simParams.recombinationRate;
    moduleDimension = simParams.moduleDimension;
    finalFitnessThreshold = 0.99;

    for i_pos = 1:length(simParams.initialAngles)
        % Monomorphic initialization only
        initialPhenotype = allInitialPhenotypes(i_pos, :);
        Fitness = -((initialPhenotype(1) / ellipseParams(1))^2 + ...
                    (initialPhenotype(2) / ellipseParams(2))^2);
        expFitness = exp(Fitness / (2 * landscapeStdDev^2));

        basePopulationMatrix = zeros(popSize, 5);
        basePopulationMatrix(:, 1) = initialPhenotype(1);
        basePopulationMatrix(:, 2) = initialPhenotype(2);
        basePopulationMatrix(:, 3) = expFitness;

        % Full mutation rate U stored in columns 4 and 5 (same for all individuals)

        basePopulationMatrix(:, 4) = mutationRate;   % Full mutation rate for module 1
        basePopulationMatrix(:, 5) = mutationRate;   % Full mutation rate for module 2

        temporaryTable = cell(1, simParams.numIteration);
        terminationStatus = zeros(1, simParams.numIteration);

        parfor i_repeat = 1:simParams.numIteration
            ellipseParamsSlice = ellipseParams;
            moduleDimensionSlice = moduleDimension;

            populationMatrix = basePopulationMatrix;

            meanFitness = mean(populationMatrix(:, 3));
            meanPhenotype = [mean(populationMatrix(:, 1)), mean(populationMatrix(:, 2))];
            EvolutionLog = [meanFitness, 1, meanPhenotype];

            simulationEnd = false;
            reachedFitness = 0;
            t = 1;

            while ~simulationEnd
                %% Mutation on module 1 (full distribution: beneficial and deleterious)
                % U/2 per module so that genome-wide total mutation rate remains NU
                numMutantTrait1 = poissrnd(popSize * mutationRate / 2);

                if numMutantTrait1 > 0
                    numMutantTrait1 = min(numMutantTrait1, popSize);
                    mutID_trait1 = datasample(1:popSize, numMutantTrait1, 'Replace', false);
                    mutMatrix1 = populationMatrix(mutID_trait1, :);

                    dx1 = sampleFullNestedEffectsVector(mutMatrix1(:,1), m, moduleDimensionSlice(1));
                    mutMatrix1(:,1) = mutMatrix1(:,1) + dx1;
                    mutMatrix1(:,1) = min(mutMatrix1(:,1), 0);  % enforce x_i <= 0 (FGM boundary)

                    newFitness1 = -((mutMatrix1(:,1) ./ ellipseParamsSlice(1)).^2 + ...
                                    (mutMatrix1(:,2) ./ ellipseParamsSlice(2)).^2);
                    mutMatrix1(:,3) = exp(newFitness1 ./ (2 * landscapeStdDev^2));

                    % Columns 4 and 5 remain U (unchanged); no update needed
                    populationMatrix(mutID_trait1, :) = mutMatrix1;
                end

                %% Mutation on module 2 (full distribution: beneficial and deleterious)
                % U/2 per module so that genome-wide total mutation rate remains NU
                numMutantTrait2 = poissrnd(popSize * mutationRate / 2);

                if numMutantTrait2 > 0
                    numMutantTrait2 = min(numMutantTrait2, popSize);
                    mutID_trait2 = datasample(1:popSize, numMutantTrait2, 'Replace', false);
                    mutMatrix2 = populationMatrix(mutID_trait2, :);

                    dx2 = sampleFullNestedEffectsVector(mutMatrix2(:,2), m, moduleDimensionSlice(2));
                    mutMatrix2(:,2) = mutMatrix2(:,2) + dx2;
                    mutMatrix2(:,2) = min(mutMatrix2(:,2), 0);  % enforce x_i <= 0 (FGM boundary)

                    newFitness2 = -((mutMatrix2(:,1) ./ ellipseParamsSlice(1)).^2 + ...
                                    (mutMatrix2(:,2) ./ ellipseParamsSlice(2)).^2);
                    mutMatrix2(:,3) = exp(newFitness2 ./ (2 * landscapeStdDev^2));

                    % Columns 4 and 5 remain U (unchanged); no update needed
                    populationMatrix(mutID_trait2, :) = mutMatrix2;
                end

                %% Recombination
                if recombinationRate ~= 0
                    if recombinationRate == 1
                        numPairs = floor(popSize / 2);
                    else
                        numPairs = binornd(popSize, recombinationRate / 2);
                        numPairs = min(numPairs, floor(popSize / 2));
                    end

                    if numPairs > 0
                        recombIndices = randperm(popSize, 2 * numPairs);
                        recMatrix = populationMatrix(recombIndices, :);

                        tempTrait2 = recMatrix(1:numPairs, 2);
                        recMatrix(1:numPairs, 2) = recMatrix(numPairs+1:end, 2);
                        recMatrix(numPairs+1:end, 2) = tempTrait2;

                        recombFitness = -((recMatrix(:,1) ./ ellipseParamsSlice(1)).^2 + ...
                                          (recMatrix(:,2) ./ ellipseParamsSlice(2)).^2);
                        recMatrix(:,3) = exp(recombFitness ./ (2 * landscapeStdDev^2));
                        % Columns 4 and 5 remain U (unchanged); no update needed

                        populationMatrix(recombIndices, :) = recMatrix;
                    end
                end

                %% Selection
                fitnessProbs = populationMatrix(:,3) / sum(populationMatrix(:,3));
                numOffsprings = mnrnd(popSize, fitnessProbs);

                parentIndices = repelem(1:popSize, numOffsprings);
                populationMatrix = populationMatrix(parentIndices(1:popSize), :);

                t = t + 1;
                newPhenotype = [mean(populationMatrix(:,1)), mean(populationMatrix(:,2))];
                newFitness = mean(populationMatrix(:,3));
                EvolutionLog = [EvolutionLog; newFitness, t, newPhenotype];

                if newFitness >= finalFitnessThreshold
                    simulationEnd = true;
                    reachedFitness = 1;
                end
            end

            temporaryTable{i_repeat} = EvolutionLog;
            terminationStatus(i_repeat) = reachedFitness;
        end

        resultNestedCM.resultTable(i_pos, :) = temporaryTable;
        resultNestedCM.terminationStatus(i_pos, :) = terminationStatus;
    end
end

% -------------------------------------------------------------------------
% Utility function — not called during simulation but retained for analysis.
% -------------------------------------------------------------------------
function [Ub1, Ub2] = computeBeneficialRatesFromX(x1, x2, U, m, moduleDimension)
    n1 = moduleDimension(1);
    n2 = moduleDimension(2);

    mu = -(m^2) / 2;

    sigma1 = m .* sqrt(max(0, 2 .* abs(x1) ./ n1));
    sigma2 = m .* sqrt(max(0, 2 .* abs(x2) ./ n2));

    pben1 = normalTailAboveZero(mu, sigma1);
    pben2 = normalTailAboveZero(mu, sigma2);

    Ub1 = U .* pben1;
    Ub2 = U .* pben2;
end

% -------------------------------------------------------------------------
% P(X > 0) for X ~ N(mu, sigma^2). Works for scalar or vector sigma.
% -------------------------------------------------------------------------
function p = normalTailAboveZero(mu, sigma)
    p = zeros(size(sigma));
    positiveMask = sigma > 0;
    z = zeros(size(sigma));
    z(positiveMask) = (0 - mu) ./ sigma(positiveMask);
    p(positiveMask) = 0.5 .* erfc(z(positiveMask) ./ sqrt(2));
    p(~positiveMask) = double(mu > 0);
end

% -------------------------------------------------------------------------
% Sample Delta x_i from the FULL (unconditional) normal distribution:
%     Delta x_i ~ N( -m^2/2,  (m * sqrt(2*|x_i|/n))^2 )
% xCurrent is a vector of current x_i values for one module.
% Both positive (beneficial) and negative (deleterious) effects are returned.
% -------------------------------------------------------------------------
function dx = sampleFullNestedEffectsVector(xCurrent, m, n)
    mu    = -(m^2) / 2;
    sigma = m .* sqrt(max(0, 2 .* abs(xCurrent) ./ n));

    dx = zeros(size(xCurrent));

    positiveMask = sigma > 0;
    if any(positiveMask)
        dx(positiveMask) = mu + sigma(positiveMask) .* randn(sum(positiveMask), 1);
    end
    % When sigma == 0 (x_i = 0, already at optimum), effect is deterministically mu < 0,
    % i.e., always deleterious; set dx = mu.
    dx(~positiveMask) = mu;
end