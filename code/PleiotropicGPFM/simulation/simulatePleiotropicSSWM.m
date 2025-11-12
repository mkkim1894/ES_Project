function [resultPleiotropicSSWM] = simulatePleiotropicSSWM(simParams, genomeParams)
    % simulatePleiotropicSSWM - Simulates the pleiotropic SSWM model.
    % Inputs:
    %   simParams - Structure containing the simulation parameters.
    %   genomeParams - Structure containing genome parameters, including genomeTheta and initialGenomes.
    %
    % Outputs:
    %   resultPleiotropicSSWM - Structure containing the results of the simulation, including genome states and fixed mutation counts.

    resultPleiotropicSSWM.resultTable = cell(length(simParams.initialAngles), simParams.numIteration);
    resultPleiotropicSSWM.fixedMutationCounts = zeros(length(simParams.initialAngles), simParams.numIteration); % Track fixed mutations

    allInitialPhenotypes = simParams.initialPhenotypes;
    popSize = simParams.popSize;
    mutationRate = simParams.mutationRate;
    deltaTrait = simParams.deltaTrait;
    landscapeStdDev = simParams.landscapeStdDev;
    ellipseParams = simParams.ellipseParams;
    
    numLoci = length(genomeParams.genomeTheta);
    genomeTheta = genomeParams.genomeTheta;
    initialGenomes = genomeParams.initialGenomes;

    for i_pos = 1:length(simParams.initialAngles)
        initialPhenotypes = allInitialPhenotypes(i_pos, 1:end);

        temporaryTable = cell(1, simParams.numIteration);
        fixedMutationCounts = zeros(1, simParams.numIteration); % Initialize mutation counts

        parfor i_repeat = 1:simParams.numIteration
            % Initialize phenotype and genome for this iteration
            initialPhenotypeSlice = initialPhenotypes;
            initialGenomeSlice = initialGenomes;
            currentGenome = initialGenomeSlice(i_pos, :); % Start with the provided initial genome

            ellipseParamsSlice = ellipseParams;
            genomeThetaSlice = genomeTheta;

            % Calculate initial fitness
            Fitness = -( (initialPhenotypeSlice(1)/ellipseParamsSlice(1))^2 + (initialPhenotypeSlice(2)/ellipseParamsSlice(2))^2 );
            Initial_Fitness = exp(Fitness / (2 * landscapeStdDev^2));

            EvolutionLog = struct('MetaInfo', zeros(1, 4), 'FixationRecords', [], 'Genome', currentGenome);
            EvolutionLog.FixationRecords = [Initial_Fitness, 1, initialPhenotypeSlice(1), initialPhenotypeSlice(2)];
            EvolutionLog.MetaInfo(1, 1:4) = [initialPhenotypeSlice(1), initialPhenotypeSlice(2), Initial_Fitness, 1];

            End = false;
            t = 1;
            ts = 1; % timestamps
            Final_W = 0.99; % Ends the simulation at this fitness level.
            fixedMutations = 0; % Initialize fixed mutation counter

            while ~End
                %% Draw a waiting time
                currentPhenotype = EvolutionLog.MetaInfo(ts, 1:2);
                r = floor( exprnd( 1 / (popSize * mutationRate) ) );
                t = t + r;

                %% Mutation selection: Randomly choose a locus to mutate
                mutationIndex = randi(numLoci);
                currentAllele = currentGenome(mutationIndex);

                % Calculate mutational effect
                if currentAllele == 0
                    mutationalEffect = [-deltaTrait * cos(genomeThetaSlice(mutationIndex)), -deltaTrait * sin(genomeThetaSlice(mutationIndex))];
                else
                    mutationalEffect = [deltaTrait * cos(genomeThetaSlice(mutationIndex)), deltaTrait * sin(genomeThetaSlice(mutationIndex))];
                end
                
                newPhenotype = currentPhenotype + mutationalEffect;

                %% Calculate new fitness and probability of fixation
                newFitness = -( (newPhenotype(1)/ellipseParamsSlice(1))^2 + (newPhenotype(2)/ellipseParamsSlice(2))^2 );
                w = exp(newFitness / (2 * landscapeStdDev^2));
                s = log(w) - log(EvolutionLog.MetaInfo(ts, 3));

                Pr_fix = (1 - exp(-2 * s)) / (1 - exp(-2 * popSize * s));
                Fixation = mybinornd(1, Pr_fix);

                if Fixation == 1
                    % Update the genome by flipping the allele at the mutation index
                    currentGenome(mutationIndex) = 1 - currentGenome(mutationIndex);
                    fixedMutations = fixedMutations + 1; % Increment mutation count
                    ts = ts + 1;
                    EvolutionLog.MetaInfo(ts, 1:2) = newPhenotype;
                    EvolutionLog.MetaInfo(ts, 3) = w;
                    EvolutionLog.MetaInfo(ts, 4) = t;
                    EvolutionLog.FixationRecords = [EvolutionLog.FixationRecords; w, t, newPhenotype];
                    EvolutionLog.Genome = currentGenome; % Update genome in log

                    if w >= Final_W
                        End = true;
                    end
                end
            end

            %% Data recording using change points and genome
            temporaryTable{i_repeat} = EvolutionLog.FixationRecords;
            fixedMutationCounts(i_repeat) = fixedMutations; % Store fixed mutations
        end

        resultPleiotropicSSWM.resultTable(i_pos, :) = temporaryTable;
        resultPleiotropicSSWM.fixedMutationCounts(i_pos, :) = fixedMutationCounts; % Store mutation counts
    end
end
