function [resultModularCMConst] = simulateModularCM_ConstantSupply(simParams)
% simulateModularCM_ConstantSupply - Modular FGM under Concurrent Mutations with a
%                                    constant supply of module-improving mutations.
%
% Description:
%   Wright-Fisher simulation identical to simulateModularCM except that the
%   fraction of loci in each module carrying a module-improving allele is held
%   FIXED, so that the beneficial mutation supply of a module does not decline as
%   the module approaches its optimum.
%
%   simulateModularCM (declining supply), per individual j and module i:
%       beneficial  rate = mu_i * max(-x_ji, 0) / delta
%       deleterious rate = mu_i * max(L_i*delta + x_ji, 0) / delta
%   with mu_i = U / (2*L_i), i.e. the b_ji = |x_ji|/delta loci that still improve
%   the module and the L_i - b_ji loci that do not.
%
%   Here, a fixed fraction f_i of the L_i loci of module i is module-improving and
%   the remaining fraction 1 - f_i is not, everywhere in trait space:
%       +delta rate = mu_i * f_i * L_i       = U * f_i / 2
%       -delta rate = mu_i * (1 - f_i) * L_i = U * (1 - f_i) / 2
%   Both are constants, independent of x_ji, of the initial condition, and of L_i,
%   and neither is ever switched off - not at the optimum, not beyond it. What
%   stops a module at its optimum is selection: a +delta step there has s < 0 and
%   is purged by Wright-Fisher sampling. No test on x_ji appears anywhere in the
%   mutation rates below.
%
%   Note that this deliberately severs the locus bookkeeping of the declining
%   model, in which |x_ji| = b_ji * delta by definition so that b_ji cannot be held
%   fixed. Severing that link is the assumption under test. f_1 and f_2 are
%   properties of the genotype-phenotype map: they differ between the two modules
%   but do not vary across the trait space or over the course of evolution.
%
%   Everything else - Poisson mutation sampling, weighted assignment of mutations
%   to individuals, recombination as a swap of the module-2 chromosome, and
%   Wright-Fisher multinomial resampling - is unchanged from simulateModularCM,
%   so any difference in outcome is attributable to the mutational supply alone.
%
%   The initial phenotype is lattice-initialized to the nearest multiple of delta
%   and clamped to <= 0, as in simulateModularCM.
%
%   Termination condition:
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
%       .geneticTargetSize  - Loci per module [L1, L2]
%       .recombinationRate  - Recombination rate (rho)
%       .numIteration       - Number of replicate simulations
%       .beneficialFraction - [f1, f2], constant fraction of module-improving loci
%       .populationMatrices - (Optional) pre-initialized populations
%
% Outputs:
%   resultModularCMConst - Structure with fields:
%       .resultTable        - Cell {nAngles x nIter}, trajectory [fit, t, x1, x2]
%       .terminationStatus  - 1: fitness threshold reached
%       .beneficialFraction - Echo of [f1, f2]
%       .beneficialSupply   - [U_1, U_2], the constant per-genome beneficial rates
%
% Reference:
%   Kim, M., Ardell, S. M., & Kryazhimskiy, S. (2025).
%   "Module-Selection Balance in the Evolution of Modular Organisms."
%
% See also: simulateModularCM, simulateModularSSWM_ConstantSupply,
%           predictModularCM_ConstantSupply
%
% Copyright (c) 2025 Minkyu Kim, Cornell University
% Licensed under MIT License

    %% Constant beneficial fractions
    if ~isfield(simParams, 'beneficialFraction') || isempty(simParams.beneficialFraction)
        error('simulateModularCM_ConstantSupply:missingFraction', ...
              'simParams.beneficialFraction = [f1, f2] must be supplied.');
    end
    f = simParams.beneficialFraction(:)';
    if numel(f) ~= 2 || any(f <= 0) || any(f > 1)
        error('simulateModularCM_ConstantSupply:badFraction', ...
              'beneficialFraction must be a 2-vector with entries in (0, 1].');
    end

    resultModularCMConst.resultTable        = cell(length(simParams.initialAngles), simParams.numIteration);
    resultModularCMConst.terminationStatus  = zeros(length(simParams.initialAngles), simParams.numIteration);
    resultModularCMConst.beneficialFraction = f;

    usePreRun = isfield(simParams, 'populationMatrices') && ~isempty(simParams.populationMatrices);

    allInitialPhenotypes = simParams.initialPhenotypes;
    popSize           = simParams.popSize;
    mutationRate      = simParams.mutationRate;   % U (genome-wide)
    deltaTrait        = simParams.deltaTrait;
    landscapeStdDev   = simParams.landscapeStdDev;
    ellipseParams     = simParams.ellipseParams;
    recombinationRate = simParams.recombinationRate;
    finalFitnessThreshold = 0.99;

    if ~usePreRun && recombinationRate > 0
        warning('Pre-run population data not found for sexual CM. Using monomorphic initialization.');
    end

    % Constant per-individual rates:
    %   beneficial  U_i     = mu_i * f_i * L_i       = U * f_i / 2
    %   deleterious U_i^del = mu_i * (1 - f_i) * L_i = U * (1 - f_i) / 2
    benRate1 = mutationRate * f(1) / 2;
    benRate2 = mutationRate * f(2) / 2;
    delRate1 = mutationRate * (1 - f(1)) / 2;
    delRate2 = mutationRate * (1 - f(2)) / 2;
    resultModularCMConst.beneficialSupply = [benRate1, benRate2];

    for i_pos = 1:length(simParams.initialAngles)
        if usePreRun
            basePopulationMatrix = simParams.populationMatrices{i_pos};
        else
            initialPhenotype = allInitialPhenotypes(i_pos, 1:end);
            % Discretize initial phenotype to delta-lattice and clamp to <= 0
            initialPhenotype(1) = -deltaTrait * round(-initialPhenotype(1) / deltaTrait);
            initialPhenotype(2) = -deltaTrait * round(-initialPhenotype(2) / deltaTrait);
            initialPhenotype = min(initialPhenotype, 0);

            Fitness    = -((initialPhenotype(1)/ellipseParams(1))^2 + (initialPhenotype(2)/ellipseParams(2))^2);
            expFitness = exp(Fitness / (2*landscapeStdDev^2));

            % Columns: [x1, x2, fitness]
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

            meanFitness   = mean(populationMatrix(:, 3));
            meanPhenotype = [mean(populationMatrix(:, 1)), mean(populationMatrix(:, 2))];
            EvolutionLog  = [meanFitness, 1, meanPhenotype];

            simulationEnd  = false;
            reachedFitness = 0;
            t = 1;

            while ~simulationEnd
                %% --- Module 1 mutations ---
                % +delta: constant rate for every individual, with no test on x1.
                ben1_rates = benRate1 * ones(popSize, 1);   % [N x 1]
                totalBen1  = sum(ben1_rates);

                % -delta: constant rate for every individual.
                del1_rates = delRate1 * ones(popSize, 1);
                totalDel1  = sum(del1_rates);

                numBen1 = poissrnd(totalBen1);
                if numBen1 > 0 && totalBen1 > 0
                    mutIDs = datasample(1:popSize, min(numBen1, nnz(ben1_rates>0)), ...
                        'Weights', ben1_rates, 'Replace', false);
                    populationMatrix(mutIDs, 1) = populationMatrix(mutIDs, 1) + deltaTrait;
                end

                numDel1 = poissrnd(totalDel1);
                if numDel1 > 0 && totalDel1 > 0
                    mutIDs = datasample(1:popSize, min(numDel1, nnz(del1_rates>0)), ...
                        'Weights', del1_rates, 'Replace', false);
                    populationMatrix(mutIDs, 1) = populationMatrix(mutIDs, 1) - deltaTrait;
                end

                %% --- Module 2 mutations ---
                ben2_rates = benRate2 * ones(popSize, 1);
                totalBen2  = sum(ben2_rates);

                del2_rates = delRate2 * ones(popSize, 1);
                totalDel2  = sum(del2_rates);

                numBen2 = poissrnd(totalBen2);
                if numBen2 > 0 && totalBen2 > 0
                    mutIDs = datasample(1:popSize, min(numBen2, nnz(ben2_rates>0)), ...
                        'Weights', ben2_rates, 'Replace', false);
                    populationMatrix(mutIDs, 2) = populationMatrix(mutIDs, 2) + deltaTrait;
                end

                numDel2 = poissrnd(totalDel2);
                if numDel2 > 0 && totalDel2 > 0
                    mutIDs = datasample(1:popSize, min(numDel2, nnz(del2_rates>0)), ...
                        'Weights', del2_rates, 'Replace', false);
                    populationMatrix(mutIDs, 2) = populationMatrix(mutIDs, 2) - deltaTrait;
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
                        recMatrix(1:numPairs, 2)     = recMatrix(numPairs+1:end, 2);
                        recMatrix(numPairs+1:end, 2) = temp;

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

        resultModularCMConst.resultTable(i_pos, :)       = temporaryTable;
        resultModularCMConst.terminationStatus(i_pos, :) = terminationStatus;
    end
end
