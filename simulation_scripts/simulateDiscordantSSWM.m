function [resultDiscordantSSWM] = simulateDiscordantSSWM(simParams)
% simulateDiscordantSSWM - Simulate discordant-module GPM under SSWM regime.
%
% Description:
%   Runs stochastic simulations of adaptive evolution in a two-chromosome
%   discordant-module genotype-phenotype map under the Strong Selection Weak
%   Mutation (SSWM) regime. Structurally identical to simulateModularSSWM
%   with the addition of the M-transformation from latent chromosome states
%   y to functional traits x.
%
%   The clock runs at the beneficial mutation supply only. Both +delta and
%   -delta displacements are proposed with equal probability; the Kimura
%   fixation formula handles deleterious outcomes. After each fixation, y
%   values are rounded to the nearest multiple of delta and clamped to <= 0,
%   enforcing the discrete locus structure of the underlying binary model.
%
%   Trait mapping:
%       x1 = y1*cos(theta1) + y2*cos(theta2)
%       x2 = y1*sin(theta1) + y2*sin(theta2)
%
%   Termination conditions:
%     1. Fitness reaches 0.99 (terminationStatus = 1)
%     2. No beneficial mutations remain (terminationStatus = 0)
%
% Inputs:
%   simParams - Structure containing simulation parameters:
%       .initialPhenotypes - Matrix of initial phenotype coordinates [nPos x 2]
%       .discordantAngles  - [theta1, theta2] chromosome-specific effect angles
%       .popSize           - Population size (N)
%       .mutationRate      - Per-genome mutation rate (U)
%       .deltaTrait        - Mutational step size (delta)
%       .landscapeStdDev   - Fitness landscape width (sigma)
%       .ellipseParams     - Ellipse axes [a1, a2]
%       .geneticTargetSize - Loci per chromosome [K1, K2]
%       .numIteration      - Number of replicate simulations
%
% Outputs:
%   resultDiscordantSSWM - Structure with fields:
%       .resultTable       - Cell array {nPos x nIter}, trajectory [fit,t,x1,x2]
%       .terminationStatus - 1: fitness threshold; 0: no beneficial mutations
%
% Reference:
%   Kim, M., Ardell, S. M., & Kryazhimskiy, S. (2025).
%   "Module-Selection Balance in the Evolution of Modular Organisms."
%
% Copyright (c) 2025 Minkyu Kim, Cornell University
% Licensed under MIT License

    nPos = size(simParams.initialPhenotypes, 1);
    resultDiscordantSSWM.resultTable       = cell(nPos, simParams.numIteration);
    resultDiscordantSSWM.terminationStatus = zeros(nPos, simParams.numIteration);

    if ~isfield(simParams, 'discordantAngles') || numel(simParams.discordantAngles) ~= 2
        error('simulateDiscordantSSWM:MissingAngles', ...
            'simParams.discordantAngles = [theta1, theta2] must be provided.');
    end

    theta1 = simParams.discordantAngles(1);
    theta2 = simParams.discordantAngles(2);

    M = [cos(theta1), cos(theta2);
         sin(theta1), sin(theta2)];

    if abs(det(M)) < 1e-12
        error('simulateDiscordantSSWM:SingularMapping', ...
            'discordantAngles produce a singular mapping matrix.');
    end

    allInitialPhenotypes = simParams.initialPhenotypes;
    popSize           = simParams.popSize;
    mutationRate      = simParams.mutationRate;   % U (genome-wide)
    deltaTrait        = simParams.deltaTrait;
    landscapeStdDev   = simParams.landscapeStdDev;
    ellipseParams     = simParams.ellipseParams;
    geneticTargetSize = simParams.geneticTargetSize;  % [K1, K2]

    % Per-locus rate mu = U/(2*K_k); beneficial supply = mu * b_k = mu * max(-y_k,0)/delta
    invTwoK1 = 1 / (2 * deltaTrait * geneticTargetSize(1));
    invTwoK2 = 1 / (2 * deltaTrait * geneticTargetSize(2));

    for i_pos = 1:nPos
        initialPhenotype = allInitialPhenotypes(i_pos, :)';
        initialY         = M \ initialPhenotype;
        % initialY may not lie on the delta-lattice; it will snap on first fixation.

        temporaryTable    = cell(1, simParams.numIteration);
        terminationStatus = zeros(1, simParams.numIteration);

        parfor i_repeat = 1:simParams.numIteration
            yCurrent = initialY;
            xCurrent = M * yCurrent;

            Fitness    = -((xCurrent(1)/ellipseParams(1))^2 + (xCurrent(2)/ellipseParams(2))^2);
            expFitness = exp(Fitness / (2*landscapeStdDev^2));

            EvolutionLog = struct('MetaInfo', zeros(1,4), 'FixationRecords', []);
            EvolutionLog.FixationRecords  = [expFitness, 1, xCurrent(1), xCurrent(2)];
            EvolutionLog.MetaInfo(1, 1:4) = [yCurrent(1), yCurrent(2), expFitness, 1];

            simulationEnd         = false;
            reachedFitness        = 0;
            t                     = 1;
            changeTime            = 1;
            finalFitnessThreshold = 0.99;

            while ~simulationEnd
                y1cur = EvolutionLog.MetaInfo(changeTime, 1);
                y2cur = EvolutionLog.MetaInfo(changeTime, 2);

                %% Beneficial supply: mu * b_k = U/(2*K_k) * max(-y_k,0)/delta
                beneficialRate1     = mutationRate * max(-y1cur, 0) * invTwoK1;
                beneficialRate2     = mutationRate * max(-y2cur, 0) * invTwoK2;
                totalBeneficialRate = beneficialRate1 + beneficialRate2;

                if totalBeneficialRate <= 0
                    reachedFitness = 0;
                    break;
                end

                waitingTime = floor(exprnd(1 / (popSize * totalBeneficialRate)));
                t = t + waitingTime;

                %% Choose chromosome proportional to beneficial supply
                p_chr1 = beneficialRate1 / totalBeneficialRate;
                if mybinornd(1, p_chr1) == 1
                    mutChr = 1;
                else
                    mutChr = 2;
                end

                %% Propose ±delta with equal probability
                direction = 2 * (rand < 0.5) - 1;
                yNew        = [y1cur; y2cur];
                yNew(mutChr) = yNew(mutChr) + direction * deltaTrait;

                %% Round to nearest delta-multiple and clamp to <= 0
                yNew(mutChr) = -deltaTrait * round(-yNew(mutChr) / deltaTrait);
                yNew(mutChr) = min(yNew(mutChr), 0);

                xNew       = M * yNew;
                newFitness = -((xNew(1)/ellipseParams(1))^2 + (xNew(2)/ellipseParams(2))^2);
                expFitnessNew = exp(newFitness / (2*landscapeStdDev^2));

                s       = log(expFitnessNew) - log(EvolutionLog.MetaInfo(changeTime, 3));
                Pr_fix  = (1 - exp(-2*s)) / (1 - exp(-2*popSize*s));
                Fixation = mybinornd(1, Pr_fix);

                if Fixation == 1
                    changeTime = changeTime + 1;
                    EvolutionLog.MetaInfo(changeTime, 1:2) = yNew';
                    EvolutionLog.MetaInfo(changeTime, 3)   = expFitnessNew;
                    EvolutionLog.MetaInfo(changeTime, 4)   = t;
                    EvolutionLog.FixationRecords = [EvolutionLog.FixationRecords; ...
                        expFitnessNew, t, xNew(1), xNew(2)];

                    if expFitnessNew >= finalFitnessThreshold
                        simulationEnd  = true;
                        reachedFitness = 1;
                    end
                end
            end

            temporaryTable{i_repeat}    = EvolutionLog.FixationRecords;
            terminationStatus(i_repeat) = reachedFitness;
        end

        resultDiscordantSSWM.resultTable(i_pos, :)       = temporaryTable;
        resultDiscordantSSWM.terminationStatus(i_pos, :) = terminationStatus;
    end
end
