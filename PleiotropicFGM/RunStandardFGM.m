%% Standard FGM - SSWM
clearvars
clc

tic

% Declare input parameters
simParams = initializeSimParams();

% Run simulation
[resultStandardSSWM] = simulateStandardSSWM(simParams);

% Compute average trajectory
numTimeStamp = 20;
[averageTrajectory] = computeAverageTrajectory(numTimeStamp, simParams, resultStandardSSWM.resultTable);

% Define directory and filename
outputDir = 'results/standardFGM/SSWM';
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
fname = sprintf('%s/StandardFGM_SSWM_M%.1e.mat', outputDir, simParams.mutationRate);

% Save results
save(fname, 'simParams', 'resultStandardSSWM', 'averageTrajectory');

% Predict analytical trajectories and append to the file
[analyticalTrajectories] = predictPleiotropicSSWM(simParams, averageTrajectory);
save(fname, 'analyticalTrajectories', '-append');

toc

%% Standard FGM - CM, asexual
clearvars
tic

% Declare input parameters
simParams = initializeSimParams('numIteration', 200, ...
    'mutationRate', 10^-3);

% Run simulation
[resultStandardCM] = simulateStandardCM(simParams);

% Compute average trajectory
numTimeStamp = 20;
[averageTrajectory] = computeAverageTrajectory(numTimeStamp, simParams, resultStandardCM.resultTable);

% Define directory and filename
outputDir = 'results/standardFGM/CM_Asexual';
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
fname = sprintf('%s/StandardFGM_CM_Asexual_M%.1e.mat', outputDir, simParams.mutationRate);

% Save results
save(fname, 'simParams', 'resultStandardCM', 'averageTrajectory');

% Predict analytical trajectories and append to the file
[analyticalTrajectories] = predictPleiotropicSSWM(simParams, averageTrajectory);
save(fname, 'analyticalTrajectories', '-append');

toc

%% Create Figures
createSimulationFigures('PleiotropicFGM')
