%% SSWM

clearvars
clc

tic

% Declare input parameters
simParams = initializeSimParams('omitParams', {'geneticTargetSize'});
[genomeParams] = initializeGenomeTheta(4000, simParams, 1); % L = 4000

% Run simulation
[resultPleiotropicSSWM] = simulatePleiotropicSSWM(simParams, genomeParams);

% Compute average trajectory
numTimeStamp = 20;
[averageTrajectory] = computeAverageTrajectory(numTimeStamp, simParams, resultPleiotropicSSWM.resultTable);

% Define directory and filename
outputDir = 'results/pleiotropicFGM/SSWM';
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
fname = sprintf('%s/PleiotropicFGM_SSWM_L%d_M%.1e.mat', outputDir, length(genomeParams.genomeTheta), simParams.mutationRate);

% Save results
save(fname, 'simParams', 'genomeParams', 'resultPleiotropicSSWM', 'averageTrajectory');

% Predict analytical trajectories and append to the file
[analyticalTrajectories] = predictPleiotropicSSWM(simParams, averageTrajectory);
save(fname, 'analyticalTrajectories', '-append');

toc

%% CM, asexual

clearvars
clc

tic

% Declare input parameters
simParams = initializeSimParams('numIteration', 200, ...
    'mutationRate', 2*10^-4, 'omitParams', {'geneticTargetSize'});
[genomeParams] = initializeGenomeTheta(400, simParams, 1); % L = 400

% Run simulation
[resultPleiotropicCM] = simulatePleiotropicCM(simParams, genomeParams);

% Compute average trajectory
numTimeStamp = 20;
[averageTrajectory] = computeAverageTrajectory(numTimeStamp, simParams, resultPleiotropicCM.resultTable);

% Define directory and filename
outputDir = 'results/pleiotropicFGM/CM_Asexual';
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
fname = sprintf('%s/PleiotropicFGM_CM_Asexual_L%d_M%.1e.mat', outputDir, length(genomeParams.genomeTheta), simParams.mutationRate);

% Save results
save(fname, 'simParams', 'genomeParams', 'resultPleiotropicCM', 'averageTrajectory');

% Predict analytical trajectories and append to the file
[analyticalTrajectories] = predictPleiotropicSSWM(simParams, averageTrajectory);
save(fname, 'analyticalTrajectories', '-append');

toc

%% CM, sexual

clearvars
clc

tic

% Declare input parameters
simParams = initializeSimParams('numIteration', 200, ...
    'mutationRate', 2*10^-4, 'recombinationRate', 1, 'omitParams', {'geneticTargetSize'});
[genomeParams] = initializeGenomeTheta(400, simParams, 1); % L = 400

% Run simulation
[resultPleiotropicCM] = simulatePleiotropicCM(simParams, genomeParams);

% Compute average trajectory
numTimeStamp = 20;
[averageTrajectory] = computeAverageTrajectory(numTimeStamp, simParams, resultPleiotropicCM.resultTable);

% Define directory and filename
outputDir = 'results/pleiotropicFGM/CM_Sexual';
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
fname = sprintf('%s/PleiotropicFGM_CM_Sexual_L%d_M%.1e_R%.4f.mat', outputDir, length(genomeParams.genomeTheta), simParams.mutationRate, simParams.recombinationRate);

% Save results
save(fname, 'simParams', 'genomeParams', 'resultPleiotropicCM', 'averageTrajectory');

% Predict analytical trajectories and append to the file
[analyticalTrajectories] = predictPleiotropicSSWM(simParams, averageTrajectory);
save(fname, 'analyticalTrajectories', '-append');

toc


%% Create figures
clc
makeFigure2('PleiotropicFGM', 1e-7, 2e-4, 4000, 400)