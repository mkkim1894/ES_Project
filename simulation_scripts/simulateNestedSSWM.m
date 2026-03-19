function [resultNestedSSWM] = simulateNestedSSWM(simParams)
% simulateNestedSSWM - Simulate nested FGM under Strong Selection Weak Mutation.
%
% Description:
%   Runs stochastic simulations of the nested FGM where each module is itself
%   a multi-dimensional FGM. Mutations within a module have effects drawn from
%   the approximate induced distribution on x_i:
%
%       Delta x_i ~ N( -m^2/2 , (m * sqrt(2*|x_i|/n_i))^2 )
%
%   where m = simParams.deltaTrait is the magnitude of the mutation vector
%   in the underlying y-space, and n_i is the dimensionality of module i.
%
%   Beneficial mutation supply is state-dependent. Waiting times are drawn
%   using the total beneficial mutation rate, and mutation effects are drawn
%   from the beneficial tail (Delta x_i > 0).
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
%       .moduleDimension   - Module dimensionalities [n1, n2]
%       .numIteration      - Number of replicate simulations
%
% Outputs:
%   resultNestedSSWM - Structure containing:
%       .resultTable       - Cell array {nAngles x nIterations}, each cell contains
%                            trajectory matrix [fitness, time, trait1, trait2]
%       .terminationStatus - 1: fitness threshold reached
%                            0: no beneficial mutations available
%
% Reference:
%   Kim, M., Ardell, S. M., & Kryazhimskiy, S. (2025).
%   "Module-Selection Balance in the Evolution of Modular Organisms."
%
% See also: simulateNestedCM, simulateModularSSWM
%
% Copyright (c) 2025 Minkyu Kim, Cornell University
% Licensed under MIT License

    resultNestedSSWM.resultTable = cell(length(simParams.initialAngles), simParams.numIteration);
    resultNestedSSWM.terminationStatus = zeros(length(simParams.initialAngles), simParams.numIteration);

    allInitialPhenotypes = simParams.initialPhenotypes;
    popSize = simParams.popSize;
    mutationRate = simParams.mutationRate;     % Per-module mutation opportunity rate
    m = simParams.deltaTrait;                  % Nested FGM mutation-vector magnitude
    landscapeStdDev = simParams.landscapeStdDev;
    ellipseParams = simParams.ellipseParams;
    moduleDimension = simParams.moduleDimension;

    finalFitnessThreshold = 0.99;

    for i_pos = 1:length(simParams.initialAngles)
        initialPhenotypes = allInitialPhenotypes(i_pos, 1:end);

        temporaryTable = cell(1, simParams.numIteration);
        terminationStatus = zeros(1, simParams.numIteration);

        parfor i_repeat = 1:simParams.numIteration
            initialPhenotypeSlice = initialPhenotypes;
            ellipseParamsSlice = ellipseParams;
            moduleDimensionSlice = moduleDimension;

            Fitness = -((initialPhenotypeSlice(1) / ellipseParamsSlice(1))^2 + ...
                        (initialPhenotypeSlice(2) / ellipseParamsSlice(2))^2);
            expFitness = exp(Fitness / (2 * landscapeStdDev^2));

            EvolutionLog = struct('MetaInfo', zeros(1, 4), 'FixationRecords', []);
            EvolutionLog.FixationRecords = [expFitness, 1, initialPhenotypeSlice(1), initialPhenotypeSlice(2)];
            EvolutionLog.MetaInfo(1, 1:4) = [initialPhenotypeSlice(1), initialPhenotypeSlice(2), expFitness, 1];

            simulationEnd = false;
            reachedFitness = 0;
            t = 1;
            changeTime = 1;

            while ~simulationEnd
                x1 = EvolutionLog.MetaInfo(changeTime, 1);
                x2 = EvolutionLog.MetaInfo(changeTime, 2);

                % Approximate mutational effect distribution for each module:
                % Delta x_i ~ N(mu_i, sigma_i^2)
                mu1 = -(m^2) / 2;
                mu2 = -(m^2) / 2;

                sigma1 = m * sqrt(max(0, 2 * abs(x1) / moduleDimensionSlice(1)));
                sigma2 = m * sqrt(max(0, 2 * abs(x2) / moduleDimensionSlice(2)));

                % Probability that a mutation in module i is beneficial: P(Delta x_i > 0)
                pben1 = normalTailAboveZero(mu1, sigma1);
                pben2 = normalTailAboveZero(mu2, sigma2);

                % Total beneficial mutation arrival rate in the population
                % U/2 per module so genome-wide total remains NU
                Ub1 = (mutationRate / 2) * pben1;
                Ub2 = (mutationRate / 2) * pben2;
                totalBeneficialRate = popSize * (Ub1 + Ub2);

                % If no beneficial mutations are effectively available, terminate
                if totalBeneficialRate <= 0 || ~isfinite(totalBeneficialRate)
                    simulationEnd = true;
                    reachedFitness = 0;
                    break;
                end

                %% Draw waiting time until next beneficial mutation appears
                waitingTime = floor(exprnd(1 / totalBeneficialRate));
                t = t + waitingTime;

                %% Choose which module gets the beneficial mutation
                if mybinornd(1, Ub1 / (Ub1 + Ub2)) == 1
                    mutDim = 1;
                    deltaX = sampleBeneficialNestedEffect(mu1, sigma1);
                else
                    mutDim = 2;
                    deltaX = sampleBeneficialNestedEffect(mu2, sigma2);
                end

                traitDisplacement = zeros(1, 2);
                traitDisplacement(mutDim) = deltaX;

                newPhenotype = EvolutionLog.MetaInfo(changeTime, 1:2) + traitDisplacement;
                newFitness = -((newPhenotype(1) / ellipseParamsSlice(1))^2 + ...
                               (newPhenotype(2) / ellipseParamsSlice(2))^2);
                expFitness = exp(newFitness / (2 * landscapeStdDev^2));
                s = log(expFitness) - log(EvolutionLog.MetaInfo(changeTime, 3));

                % Numerical guard in case sampling noise yields s <= 0.
                if s <= 0 || ~isfinite(s)
                    continue;
                end

                Pr_fix = (1 - exp(-2 * s)) / (1 - exp(-2 * popSize * s));
                Fixation = mybinornd(1, Pr_fix);

                if Fixation == 1
                    changeTime = changeTime + 1;
                    EvolutionLog.MetaInfo(changeTime, 1:2) = newPhenotype;
                    EvolutionLog.MetaInfo(changeTime, 3) = expFitness;
                    EvolutionLog.MetaInfo(changeTime, 4) = t;
                    EvolutionLog.FixationRecords = [EvolutionLog.FixationRecords; expFitness, t, newPhenotype];

                    if expFitness >= finalFitnessThreshold
                        simulationEnd = true;
                        reachedFitness = 1;
                    end
                end
            end

            temporaryTable{i_repeat} = EvolutionLog.FixationRecords;
            terminationStatus(i_repeat) = reachedFitness;
        end

        resultNestedSSWM.resultTable(i_pos, :) = temporaryTable;
        resultNestedSSWM.terminationStatus(i_pos, :) = terminationStatus;
    end
end

% -------------------------------------------------------------------------
% P(X > 0) for X ~ N(mu, sigma^2)
% -------------------------------------------------------------------------
function p = normalTailAboveZero(mu, sigma)
    if sigma <= 0
        p = double(mu > 0);
        return;
    end

    z = (0 - mu) / sigma;
    p = 0.5 * erfc(z / sqrt(2));   % 1 - Phi(z)
end

% -------------------------------------------------------------------------
% Sample X ~ N(mu, sigma^2) conditional on X > 0
% Uses inverse-CDF sampling from the truncated normal.
% -------------------------------------------------------------------------
function x = sampleBeneficialNestedEffect(mu, sigma)
    if sigma <= 0
        error('sampleBeneficialNestedEffect:InvalidSigma', ...
              'sigma must be positive to sample a beneficial effect.');
    end

    cdf0 = 0.5 * erfc(-(0 - mu) / (sigma * sqrt(2)));   % Phi((0-mu)/sigma)

    if cdf0 >= 1
        error('sampleBeneficialNestedEffect:NoBeneficialTail', ...
              'Beneficial tail has zero probability.');
    end

    u = cdf0 + rand() * (1 - cdf0);
    x = mu + sigma * sqrt(2) * erfinv(2 * u - 1);

    % Numerical guard
    if x <= 0
        x = eps;
    end
end