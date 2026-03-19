function [resultModularCM] = simulateModularCM(simParams)
% simulateModularCM - Simulate modular FGM under Concurrent Mutations regime.
%
% Description:
%   Runs Wright-Fisher simulations of adaptive evolution in a two-module
%   Fisher's Geometric Model under the CM regime.
%
%   Beneficial (1->0, +delta) and deleterious (0->1, -delta) mutations are
%   drawn separately each generation with correct rates from the binary locus
%   model:
%       Beneficial rate for individual j, module i:
%           mu * b_ji = U/(2*K_i) * max(-x_ji, 0) / delta
%       Deleterious rate for individual j, module i:
%           mu * (K_i - b_ji) = U/(2*K_i) * max(K_i*delta + x_ji, 0) / delta
%
%   After mutation and recombination, trait values are rounded to the nearest
%   multiple of delta and clamped to <= 0, enforcing the discrete locus
%   structure of the underlying binary genotype model.
%
%   Termination conditions:
%     1. Mean fitness reaches 0.99 (terminationStatus = 1)
%
% Inputs:
%   simParams - Structure containing simulation parameters:
%       .initialAngles      - Vector of initial angles in phenotype space
%       .initialPhenotypes  - Matrix of initial phenotype coordinates [nAngles x 2]
%       .popSize            - Population size (N)
%       .mutationRate       - Per-genome mutation rate (U)
%       .deltaTrait         - Mutational step size (delta)
%       .landscapeStdDev    - Fitness landscape width (sigma)
%       .ellipseParams      - Ellipse axes [a1, a2]
%       .geneticTargetSize  - Loci per module [K1, K2]
%       .recombinationRate  - Recombination rate (rho)
%       .numIteration       - Number of replicate simulations
%       .populationMatrices - (Optional) pre-initialized populations
%
% Outputs:
%   resultModularCM - Structure with fields:
%       .resultTable       - Cell {nAngles x nIter}, trajectory [fit,t,x1,x2]
%       .terminationStatus - 1: fitness threshold reached
%
% Reference:
%   Kim, M., Ardell, S. M., & Kryazhimskiy, S. (2025).
%   "Module-Selection Balance in the Evolution of Modular Organisms."
%
% See also: simulateModularSSWM, predictModularCM, initializeSimParams
%
% Copyright (c) 2025 Minkyu Kim, Cornell University
% Licensed under MIT License

    resultModularCM.resultTable       = cell(length(simParams.initialAngles), simParams.numIteration);
    resultModularCM.terminationStatus = zeros(length(simParams.initialAngles), simParams.numIteration);

    usePreRun = isfield(simParams, 'populationMatrices') && ~isempty(simParams.populationMatrices);

    allInitialPhenotypes = simParams.initialPhenotypes;
    popSize           = simParams.popSize;
    mutationRate      = simParams.mutationRate;   % U (genome-wide)
    deltaTrait        = simParams.deltaTrait;
    landscapeStdDev   = simParams.landscapeStdDev;
    ellipseParams     = simParams.ellipseParams;
    geneticTargetSize = simParams.geneticTargetSize;  % [K1, K2]
    recombinationRate = simParams.recombinationRate;
    finalFitnessThreshold = 0.99;

    if ~usePreRun && recombinationRate > 0
        warning('Pre-run population data not found for sexual CM. Using monomorphic initialization.');
    end

    % Per-locus rate: mu = U / (2*K_i)
    % Beneficial rate_ji  = mu * max(-x_ji, 0) / delta
    % Deleterious rate_ji = mu * max(K_i*delta + x_ji, 0) / delta
    mu1    = mutationRate / (2 * geneticTargetSize(1));  % per locus, module 1
    mu2    = mutationRate / (2 * geneticTargetSize(2));  % per locus, module 2
    K1del  = geneticTargetSize(1) * deltaTrait;          % K1 * delta (upper boundary module 1)
    K2del  = geneticTargetSize(2) * deltaTrait;          % K2 * delta (upper boundary module 2)
    invD   = 1 / deltaTrait;

    for i_pos = 1:length(simParams.initialAngles)
        if usePreRun
            basePopulationMatrix = simParams.populationMatrices{i_pos};
        else
            initialPhenotype = allInitialPhenotypes(i_pos, 1:end);
            Fitness    = -((initialPhenotype(1)/ellipseParams(1))^2 + (initialPhenotype(2)/ellipseParams(2))^2);
            expFitness = exp(Fitness / (2*landscapeStdDev^2));

            % Columns: [x1, x2, fitness]
            % Rates computed on the fly each generation from x1, x2
            basePopulationMatrix = zeros(popSize, 3);
            basePopulationMatrix(:, 1) = initialPhenotype(1);
            basePopulationMatrix(:, 2) = initialPhenotype(2);
            basePopulationMatrix(:, 3) = expFitness;
        end

        temporaryTable    = cell(1, simParams.numIteration);
        terminationStatus = zeros(1, simParams.numIteration);

        parfor i_repeat = 1:simParams.numIteration
            ellipseParamsSlice = ellipseParams;
            populationMatrix   = basePopulationMatrix;

            meanFitness  = mean(populationMatrix(:, 3));
            meanPhenotype = [mean(populationMatrix(:, 1)), mean(populationMatrix(:, 2))];
            EvolutionLog  = [meanFitness, 1, meanPhenotype];

            simulationEnd  = false;
            reachedFitness = 0;
            t = 1;

            while ~simulationEnd
                x1col = populationMatrix(:, 1);
                x2col = populationMatrix(:, 2);

                %% --- Module 1 mutations ---
                % Beneficial (1->0): rate = mu1 * max(-x1,0) / delta
                ben1_rates  = mu1 * max(-x1col, 0) * invD;   % [N x 1]
                totalBen1   = sum(ben1_rates);

                % Deleterious (0->1): rate = mu1 * max(K1*delta + x1, 0) / delta
                del1_rates  = mu1 * max(K1del + x1col, 0) * invD;
                totalDel1   = sum(del1_rates);

                % Beneficial mutations on module 1
                numBen1 = poissrnd(totalBen1);
                if numBen1 > 0 && totalBen1 > 0
                    mutIDs = datasample(1:popSize, min(numBen1, nnz(ben1_rates>0)), ...
                        'Weights', ben1_rates, 'Replace', false);
                    populationMatrix(mutIDs, 1) = populationMatrix(mutIDs, 1) + deltaTrait;
                    % Round and clamp
                    populationMatrix(mutIDs, 1) = min(-deltaTrait * round(-populationMatrix(mutIDs,1)/deltaTrait), 0);
                end

                % Deleterious mutations on module 1
                numDel1 = poissrnd(totalDel1);
                if numDel1 > 0 && totalDel1 > 0
                    mutIDs = datasample(1:popSize, min(numDel1, nnz(del1_rates>0)), ...
                        'Weights', del1_rates, 'Replace', false);
                    populationMatrix(mutIDs, 1) = populationMatrix(mutIDs, 1) - deltaTrait;
                    % Round and clamp
                    populationMatrix(mutIDs, 1) = min(-deltaTrait * round(-populationMatrix(mutIDs,1)/deltaTrait), 0);
                end

                %% --- Module 2 mutations ---
                x2col = populationMatrix(:, 2);

                ben2_rates  = mu2 * max(-x2col, 0) * invD;
                totalBen2   = sum(ben2_rates);

                del2_rates  = mu2 * max(K2del + x2col, 0) * invD;
                totalDel2   = sum(del2_rates);

                % Beneficial mutations on module 2
                numBen2 = poissrnd(totalBen2);
                if numBen2 > 0 && totalBen2 > 0
                    mutIDs = datasample(1:popSize, min(numBen2, nnz(ben2_rates>0)), ...
                        'Weights', ben2_rates, 'Replace', false);
                    populationMatrix(mutIDs, 2) = populationMatrix(mutIDs, 2) + deltaTrait;
                    populationMatrix(mutIDs, 2) = min(-deltaTrait * round(-populationMatrix(mutIDs,2)/deltaTrait), 0);
                end

                % Deleterious mutations on module 2
                numDel2 = poissrnd(totalDel2);
                if numDel2 > 0 && totalDel2 > 0
                    mutIDs = datasample(1:popSize, min(numDel2, nnz(del2_rates>0)), ...
                        'Weights', del2_rates, 'Replace', false);
                    populationMatrix(mutIDs, 2) = populationMatrix(mutIDs, 2) - deltaTrait;
                    populationMatrix(mutIDs, 2) = min(-deltaTrait * round(-populationMatrix(mutIDs,2)/deltaTrait), 0);
                end

                %% Recompute fitness after all mutations
                newFitness = -((populationMatrix(:,1)./ellipseParamsSlice(1)).^2 + ...
                               (populationMatrix(:,2)./ellipseParamsSlice(2)).^2);
                populationMatrix(:, 3) = exp(newFitness / (2*landscapeStdDev^2));

                %% Recombination (swap x2 between pairs)
                if recombinationRate ~= 0
                    if recombinationRate == 1
                        numPairs = floor(popSize / 2);
                    else
                        numPairs = min(binornd(popSize, recombinationRate/2), floor(popSize/2));
                    end

                    if numPairs > 0
                        recombIndices = randperm(popSize, 2*numPairs);
                        recMatrix = populationMatrix(recombIndices, :);

                        temp = recMatrix(1:numPairs, 2);
                        recMatrix(1:numPairs, 2)       = recMatrix(numPairs+1:end, 2);
                        recMatrix(numPairs+1:end, 2)   = temp;

                        % Round and clamp x2 after recombination
                        recMatrix(:, 2) = min(-deltaTrait * round(-recMatrix(:,2)/deltaTrait), 0);

                        % Recompute fitness
                        recFitness = -((recMatrix(:,1)./ellipseParamsSlice(1)).^2 + ...
                                       (recMatrix(:,2)./ellipseParamsSlice(2)).^2);
                        recMatrix(:, 3) = exp(recFitness / (2*landscapeStdDev^2));

                        populationMatrix(recombIndices, :) = recMatrix;
                    end
                end

                %% Selection (Wright-Fisher sampling)
                fitnessProbs  = populationMatrix(:, 3) / sum(populationMatrix(:, 3));
                numOffsprings = mnrnd(popSize, fitnessProbs);
                parentIndices = repelem(1:popSize, numOffsprings);
                populationMatrix = populationMatrix(parentIndices(1:popSize), :);

                t = t + 1;
                meanPhenotype = [mean(populationMatrix(:,1)), mean(populationMatrix(:,2))];
                meanFitness   = mean(populationMatrix(:, 3));
                EvolutionLog  = [EvolutionLog; meanFitness, t, meanPhenotype];

                if meanFitness >= finalFitnessThreshold
                    simulationEnd  = true;
                    reachedFitness = 1;
                end
            end

            temporaryTable{i_repeat}    = EvolutionLog;
            terminationStatus(i_repeat) = reachedFitness;
        end

        resultModularCM.resultTable(i_pos, :)       = temporaryTable;
        resultModularCM.terminationStatus(i_pos, :) = terminationStatus;
    end
end
