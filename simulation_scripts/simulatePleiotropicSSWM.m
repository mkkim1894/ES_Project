function [resultPleiotropicSSWM] = simulatePleiotropicSSWM(simParams, genomeParams)
% simulatePleiotropicSSWM - Simulate pleiotropic FGM under SSWM regime.
%
% Description:
%   Runs stochastic simulations where each mutation affects both traits with
%   a pleiotropic angle. Mutations fix sequentially under SSWM dynamics.
%   Only loci with selection coefficient s > 1/N are considered as beneficial
%   mutation candidates.
%
%   Termination conditions:
%     1. Mean fitness reaches 0.99 (terminationStatus = 1)
%     2. No beneficial mutations available (terminationStatus = 0)
%
% Inputs:
%   simParams    - Structure containing simulation parameters (see initializeSimParams)
%   genomeParams - Structure containing:
%       .genomeTheta    - Vector of pleiotropic angles for each locus
%       .initialGenomes - Matrix of initial allele states [nAngles x nLoci]
%
% Outputs:
%   resultPleiotropicSSWM - Structure containing:
%       .resultTable         - Cell array of trajectory matrices [w, t, x1, x2]
%       .fixedMutationCounts - Count of fixed mutations per replicate
%       .terminationStatus   - 0: no beneficial mutations, 1: fitness threshold
%
% Reference:
%   Kim, M., Ardell, S. M., & Kryazhimskiy, S. (2025).
%   "Module-Selection Balance in the Evolution of Modular Organisms."
%
% See also: simulatePleiotropicCM, initializeGenomeTheta
%
% Copyright (c) 2025 Minkyu Kim, Cornell University
% Licensed under MIT License

    resultPleiotropicSSWM.resultTable = cell(length(simParams.initialAngles), simParams.numIteration);
    resultPleiotropicSSWM.fixedMutationCounts = zeros(length(simParams.initialAngles), simParams.numIteration);
    resultPleiotropicSSWM.terminationStatus = zeros(length(simParams.initialAngles), simParams.numIteration);

    allInitialPhenotypes = simParams.initialPhenotypes;
    popSize = simParams.popSize;
    mutationRate = simParams.mutationRate;
    deltaTrait = simParams.deltaTrait;
    landscapeStdDev = simParams.landscapeStdDev;
    ellipseParams = simParams.ellipseParams;

    numLoci = length(genomeParams.genomeTheta);
    genomeTheta = genomeParams.genomeTheta;
    initialGenomes = genomeParams.initialGenomes;

    % Precompute cos and sin of all locus angles
    cosTheta = cos(genomeTheta);
    sinTheta = sin(genomeTheta);

    % Neutral threshold: mutations with s <= 1/N are effectively neutral
    neutralThreshold = 1 / popSize;

    for i_pos = 1:length(simParams.initialAngles)
        initialPhenotypes = allInitialPhenotypes(i_pos, 1:end);

        temporaryTable = cell(1, simParams.numIteration);
        fixedMutationCounts = zeros(1, simParams.numIteration);
        terminationStatus = zeros(1, simParams.numIteration);

        parfor i_repeat = 1:simParams.numIteration
            % Initialize phenotype and genome for this iteration
            initialPhenotypeSlice = initialPhenotypes;
            initialGenomeSlice = initialGenomes;
            currentGenome = initialGenomeSlice(i_pos, :);

            ellipseParamsSlice = ellipseParams;
            cosThetaSlice = cosTheta;
            sinThetaSlice = sinTheta;

            % Calculate initial fitness
            Fitness = -( (initialPhenotypeSlice(1)/ellipseParamsSlice(1))^2 + ...
                         (initialPhenotypeSlice(2)/ellipseParamsSlice(2))^2 );
            Initial_Fitness = exp(Fitness / (2 * landscapeStdDev^2));

            EvolutionLog = struct('MetaInfo', zeros(1, 4), 'FixationRecords', [], 'Genome', currentGenome);
            EvolutionLog.FixationRecords = [Initial_Fitness, 1, initialPhenotypeSlice(1), initialPhenotypeSlice(2)];
            EvolutionLog.MetaInfo(1, 1:4) = [initialPhenotypeSlice(1), initialPhenotypeSlice(2), Initial_Fitness, 1];

            End = false;
            t = 1;
            ts = 1;
            Final_W = 0.99;
            fixedMutations = 0;
            reachedFitness = 0;
            needRescan = true;

            beneficialIdx = [];
            beneficialS = [];
            numBeneficial = 0;

            while ~End
                currentPhenotype = EvolutionLog.MetaInfo(ts, 1:2);
                currentLogFitness = log(EvolutionLog.MetaInfo(ts, 3));

                %% Scan for beneficial loci (after fixation or at start)
                if needRescan
                    % For allele 0 -> 1: effect is [-delta*cos, -delta*sin], sign = -1
                    % For allele 1 -> 0: effect is [+delta*cos, +delta*sin], sign = +1
                    signs = 2 * currentGenome - 1;
                    newPhenotype1 = currentPhenotype(1) + signs .* deltaTrait .* cosThetaSlice;
                    newPhenotype2 = currentPhenotype(2) + signs .* deltaTrait .* sinThetaSlice;

                    % Compute fitness for all possible single-locus mutations
                    newLogFitness = -( (newPhenotype1 ./ ellipseParamsSlice(1)).^2 + ...
                                       (newPhenotype2 ./ ellipseParamsSlice(2)).^2 ) / (2 * landscapeStdDev^2);
                    allS = newLogFitness - currentLogFitness;

                    % Identify beneficial loci (s > 1/N)
                    beneficialMask = allS > neutralThreshold;
                    beneficialIdx = find(beneficialMask);
                    beneficialS = allS(beneficialMask);
                    numBeneficial = length(beneficialIdx);

                    needRescan = false;
                end

                %% Check termination: no beneficial mutations (drift regime)
                if numBeneficial == 0
                    End = true;
                    reachedFitness = 0;
                    continue;
                end

                %% Draw waiting time based on beneficial mutation supply
                beneficialRate = popSize * mutationRate * (numBeneficial / numLoci);
                r = floor( exprnd(1 / beneficialRate) );
                t = t + r;

                %% Select a random beneficial locus
                pick = randi(numBeneficial);
                mutationIndex = beneficialIdx(pick);
                s = beneficialS(pick);

                %% Apply Kimura fixation probability
                Pr_fix = (1 - exp(-2 * s)) / (1 - exp(-2 * popSize * s));
                Fixation = mybinornd(1, Pr_fix);

                if Fixation == 1
                    % Compute new phenotype BEFORE flipping
                    oldAllele = currentGenome(mutationIndex);
                    if oldAllele == 0
                        newPhenotype = currentPhenotype + [-deltaTrait * cosThetaSlice(mutationIndex), ...
                                                           -deltaTrait * sinThetaSlice(mutationIndex)];
                    else
                        newPhenotype = currentPhenotype + [deltaTrait * cosThetaSlice(mutationIndex), ...
                                                           deltaTrait * sinThetaSlice(mutationIndex)];
                    end

                    % Flip the allele
                    currentGenome(mutationIndex) = 1 - currentGenome(mutationIndex);
                    fixedMutations = fixedMutations + 1;

                    newFitness = -( (newPhenotype(1)/ellipseParamsSlice(1))^2 + ...
                                    (newPhenotype(2)/ellipseParamsSlice(2))^2 );
                    w = exp(newFitness / (2 * landscapeStdDev^2));

                    ts = ts + 1;
                    EvolutionLog.MetaInfo(ts, 1:2) = newPhenotype;
                    EvolutionLog.MetaInfo(ts, 3) = w;
                    EvolutionLog.MetaInfo(ts, 4) = t;
                    EvolutionLog.FixationRecords = [EvolutionLog.FixationRecords; w, t, newPhenotype];
                    EvolutionLog.Genome = currentGenome;

                    % Check termination: fitness threshold reached
                    if w >= Final_W
                        End = true;
                        reachedFitness = 1;
                    else
                        needRescan = true;
                    end
                end
            end

            %% Store results for this replicate
            temporaryTable{i_repeat} = EvolutionLog.FixationRecords;
            fixedMutationCounts(i_repeat) = fixedMutations;
            terminationStatus(i_repeat) = reachedFitness;
        end

        resultPleiotropicSSWM.resultTable(i_pos, :) = temporaryTable;
        resultPleiotropicSSWM.fixedMutationCounts(i_pos, :) = fixedMutationCounts;
        resultPleiotropicSSWM.terminationStatus(i_pos, :) = terminationStatus;
    end
end