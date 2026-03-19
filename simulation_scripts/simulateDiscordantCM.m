function [resultDiscordantCM] = simulateDiscordantCM(simParams)
% simulateDiscordantCM - Simulate discordant-module GPM under Concurrent Mutations regime.
%
% Description:
%   Runs Wright-Fisher simulations of adaptive evolution in a two-chromosome
%   discordant-module genotype-phenotype map. Structurally identical to
%   simulateModularCM with the addition of the M-transformation from latent
%   chromosome states y to functional traits x.
%
%   Beneficial (1->0, +delta) and deleterious (0->1, -delta) mutations on
%   each chromosome are drawn separately each generation with correct rates:
%       Beneficial rate for individual j, chromosome k:
%           mu * b_jk = U/(2*K_k) * max(-y_jk, 0) / delta
%       Deleterious rate for individual j, chromosome k:
%           mu * (K_k - b_jk) = U/(2*K_k) * max(K_k*delta + y_jk, 0) / delta
%
%   After mutation and recombination, y values are rounded to the nearest
%   multiple of delta and clamped to <= 0, enforcing the discrete locus
%   structure. Fitness is always evaluated in x-space via x = M*y.
%
%   Trait mapping:
%       x1 = y1*cos(theta1) + y2*cos(theta2)
%       x2 = y1*sin(theta1) + y2*sin(theta2)
%
%   Termination conditions:
%     1. Mean fitness reaches 0.99 (terminationStatus = 1)
%     2. Max generations reached (terminationStatus = -1)
%
% Inputs:
%   simParams - Structure containing simulation parameters:
%       .initialPhenotypes  - Matrix of initial phenotype coordinates [nPos x 2]
%       .discordantAngles   - [theta1, theta2] chromosome effect angles
%       .popSize            - Population size (N)
%       .mutationRate       - Per-genome mutation rate (U)
%       .deltaTrait         - Mutational step size (delta)
%       .landscapeStdDev    - Fitness landscape width (sigma)
%       .ellipseParams      - Ellipse axes [a1, a2]
%       .geneticTargetSize  - Loci per chromosome [K1, K2]
%       .recombinationRate  - Recombination rate (rho)
%       .numIteration       - Number of replicate simulations
%       .maxGenerations     - (Optional) hard cap, default 1e6
%
% Outputs:
%   resultDiscordantCM - Structure with fields:
%       .resultTable       - Cell {nPos x nIter}, trajectory [fit,t,x1,x2]
%       .terminationStatus - 1: fitness threshold; -1: max generations
%
% Reference:
%   Kim, M., Ardell, S. M., & Kryazhimskiy, S. (2025).
%   "Module-Selection Balance in the Evolution of Modular Organisms."
%
% Copyright (c) 2025 Minkyu Kim, Cornell University
% Licensed under MIT License

    nPos = size(simParams.initialPhenotypes, 1);
    resultDiscordantCM.resultTable       = cell(nPos, simParams.numIteration);
    resultDiscordantCM.terminationStatus = zeros(nPos, simParams.numIteration);

    if ~isfield(simParams, 'discordantAngles') || numel(simParams.discordantAngles) ~= 2
        error('simulateDiscordantCM:MissingAngles', ...
            'simParams.discordantAngles = [theta1, theta2] must be provided.');
    end

    theta1 = simParams.discordantAngles(1);
    theta2 = simParams.discordantAngles(2);

    M = [cos(theta1), cos(theta2);
         sin(theta1), sin(theta2)];

    if abs(det(M)) < 1e-12
        error('simulateDiscordantCM:SingularMapping', ...
            'discordantAngles produce a singular mapping matrix.');
    end

    if isfield(simParams, 'populationMatrices') && ~isempty(simParams.populationMatrices)
        warning(['simulateDiscordantCM ignores simParams.populationMatrices. ' ...
                 'Discordant pre-run initialization has not been implemented.']);
    end

    allInitialPhenotypes = simParams.initialPhenotypes;
    popSize           = simParams.popSize;
    mutationRate      = simParams.mutationRate;   % U (genome-wide)
    deltaTrait        = simParams.deltaTrait;
    landscapeStdDev   = simParams.landscapeStdDev;
    ellipseParams     = simParams.ellipseParams;
    geneticTargetSize = simParams.geneticTargetSize;  % [K1, K2]
    recombinationRate = simParams.recombinationRate;
    finalFitnessThreshold = 0.99;

    if isfield(simParams, 'maxGenerations') && ~isempty(simParams.maxGenerations)
        maxGenerations = simParams.maxGenerations;
    else
        maxGenerations = 1e6;
    end

    % Per-locus rate: mu_k = U / (2*K_k)
    mu1   = mutationRate / (2 * geneticTargetSize(1));
    mu2   = mutationRate / (2 * geneticTargetSize(2));
    K1del = geneticTargetSize(1) * deltaTrait;   % upper boundary chr 1
    K2del = geneticTargetSize(2) * deltaTrait;   % upper boundary chr 2
    invD  = 1 / deltaTrait;

    for i_pos = 1:nPos
        initialPhenotype = allInitialPhenotypes(i_pos, :)';
        initialY         = M \ initialPhenotype;
        % initialY may not lie on the delta-lattice; snaps on first operation.

        x0         = M * initialY;
        Fitness    = -((x0(1)/ellipseParams(1))^2 + (x0(2)/ellipseParams(2))^2);
        expFitness = exp(Fitness / (2*landscapeStdDev^2));

        % Columns: [y1, y2, fitness]  (rates computed on the fly)
        basePopulationMatrix = zeros(popSize, 3);
        basePopulationMatrix(:, 1) = initialY(1);
        basePopulationMatrix(:, 2) = initialY(2);
        basePopulationMatrix(:, 3) = expFitness;

        temporaryTable    = cell(1, simParams.numIteration);
        terminationStatus = zeros(1, simParams.numIteration);

        parfor i_repeat = 1:simParams.numIteration
            ellipseParamsSlice = ellipseParams;
            populationMatrix   = basePopulationMatrix;

            meanY         = [mean(populationMatrix(:,1)); mean(populationMatrix(:,2))];
            meanPhenotype = M * meanY;
            meanFitness   = mean(populationMatrix(:, 3));
            EvolutionLog  = [meanFitness, 1, meanPhenotype'];

            simulationEnd  = false;
            reachedFitness = 0;
            t = 1;

            while ~simulationEnd
                if t >= maxGenerations
                    simulationEnd  = true;
                    reachedFitness = -1;
                    break;
                end

                y1col = populationMatrix(:, 1);
                y2col = populationMatrix(:, 2);

                %% --- Chromosome 1 mutations ---
                % Beneficial (1->0, +delta on y1): mu1 * max(-y1,0)/delta
                ben1_rates = mu1 * max(-y1col, 0) * invD;
                totalBen1  = sum(ben1_rates);

                % Deleterious (0->1, -delta on y1): mu1 * max(K1*delta+y1,0)/delta
                del1_rates = mu1 * max(K1del + y1col, 0) * invD;
                totalDel1  = sum(del1_rates);

                numBen1 = poissrnd(totalBen1);
                if numBen1 > 0 && totalBen1 > 0
                    mutIDs = datasample(1:popSize, min(numBen1, nnz(ben1_rates>0)), ...
                        'Weights', ben1_rates, 'Replace', false);
                    populationMatrix(mutIDs, 1) = populationMatrix(mutIDs, 1) + deltaTrait;
                    populationMatrix(mutIDs, 1) = min(-deltaTrait * round(-populationMatrix(mutIDs,1)/deltaTrait), 0);
                end

                numDel1 = poissrnd(totalDel1);
                if numDel1 > 0 && totalDel1 > 0
                    mutIDs = datasample(1:popSize, min(numDel1, nnz(del1_rates>0)), ...
                        'Weights', del1_rates, 'Replace', false);
                    populationMatrix(mutIDs, 1) = populationMatrix(mutIDs, 1) - deltaTrait;
                    populationMatrix(mutIDs, 1) = min(-deltaTrait * round(-populationMatrix(mutIDs,1)/deltaTrait), 0);
                end

                %% --- Chromosome 2 mutations ---
                y2col = populationMatrix(:, 2);

                ben2_rates = mu2 * max(-y2col, 0) * invD;
                totalBen2  = sum(ben2_rates);

                del2_rates = mu2 * max(K2del + y2col, 0) * invD;
                totalDel2  = sum(del2_rates);

                numBen2 = poissrnd(totalBen2);
                if numBen2 > 0 && totalBen2 > 0
                    mutIDs = datasample(1:popSize, min(numBen2, nnz(ben2_rates>0)), ...
                        'Weights', ben2_rates, 'Replace', false);
                    populationMatrix(mutIDs, 2) = populationMatrix(mutIDs, 2) + deltaTrait;
                    populationMatrix(mutIDs, 2) = min(-deltaTrait * round(-populationMatrix(mutIDs,2)/deltaTrait), 0);
                end

                numDel2 = poissrnd(totalDel2);
                if numDel2 > 0 && totalDel2 > 0
                    mutIDs = datasample(1:popSize, min(numDel2, nnz(del2_rates>0)), ...
                        'Weights', del2_rates, 'Replace', false);
                    populationMatrix(mutIDs, 2) = populationMatrix(mutIDs, 2) - deltaTrait;
                    populationMatrix(mutIDs, 2) = min(-deltaTrait * round(-populationMatrix(mutIDs,2)/deltaTrait), 0);
                end

                %% Recompute fitness in x-space after all mutations
                y1col = populationMatrix(:, 1);
                y2col = populationMatrix(:, 2);
                x1_all = y1col * cos(theta1) + y2col * cos(theta2);
                x2_all = y1col * sin(theta1) + y2col * sin(theta2);
                newFitness = -((x1_all./ellipseParamsSlice(1)).^2 + (x2_all./ellipseParamsSlice(2)).^2);
                populationMatrix(:, 3) = exp(newFitness / (2*landscapeStdDev^2));

                %% Recombination (swap y2 between pairs)
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
                        recMatrix(1:numPairs, 2)     = recMatrix(numPairs+1:end, 2);
                        recMatrix(numPairs+1:end, 2) = temp;

                        % Round and clamp y2 after recombination
                        recMatrix(:, 2) = min(-deltaTrait * round(-recMatrix(:,2)/deltaTrait), 0);

                        % Recompute fitness in x-space
                        x1_rec = recMatrix(:,1)*cos(theta1) + recMatrix(:,2)*cos(theta2);
                        x2_rec = recMatrix(:,1)*sin(theta1) + recMatrix(:,2)*sin(theta2);
                        recFitness = -((x1_rec./ellipseParamsSlice(1)).^2 + (x2_rec./ellipseParamsSlice(2)).^2);
                        recMatrix(:, 3) = exp(recFitness / (2*landscapeStdDev^2));

                        populationMatrix(recombIndices, :) = recMatrix;
                    end
                end

                %% Selection (Wright-Fisher sampling)
                totalFitness = sum(populationMatrix(:, 3));
                if totalFitness <= 0 || ~isfinite(totalFitness)
                    error('simulateDiscordantCM:InvalidFitnessWeights', ...
                          'Population fitness weights became non-finite or non-positive.');
                end

                fitnessProbs  = populationMatrix(:, 3) / totalFitness;
                numOffsprings = mnrnd(popSize, fitnessProbs);
                parentIndices = repelem((1:popSize)', numOffsprings);
                populationMatrix = populationMatrix(parentIndices, :);

                t = t + 1;
                y1col = populationMatrix(:, 1);
                y2col = populationMatrix(:, 2);
                meanY         = [mean(y1col); mean(y2col)];
                meanPhenotype = M * meanY;
                meanFitness   = mean(populationMatrix(:, 3));
                EvolutionLog  = [EvolutionLog; meanFitness, t, meanPhenotype'];

                if meanFitness >= finalFitnessThreshold
                    simulationEnd  = true;
                    reachedFitness = 1;
                end
            end

            temporaryTable{i_repeat}    = EvolutionLog;
            terminationStatus(i_repeat) = reachedFitness;
        end

        resultDiscordantCM.resultTable(i_pos, :)       = temporaryTable;
        resultDiscordantCM.terminationStatus(i_pos, :) = terminationStatus;
    end
end
