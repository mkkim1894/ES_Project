%% Additional Modular GPFMs Simulations
%%% SSWM
clearvars
clc

tic

% Declare input parameters
simParams = initializeSimParams();

% Run simulation
[resultModularSSWM] = simulateModularSSWM(simParams);

% Compute average trajectory
numTimeStamp = 20;
[averageTrajectory] = computeAverageTrajectory(numTimeStamp, simParams, resultModularSSWM.resultTable);

% Create descriptive file name
outputDir = 'results/SSWM';
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
fname = sprintf('%s/ModularFGM_SSWM_M%.1e.mat', outputDir, simParams.mutationRate);

% Save results
save(fname, 'simParams', 'resultModularSSWM', 'averageTrajectory');

% Predict analytical trajectories and append to the file
[analyticalTrajectories] = predictModularSSWM(simParams, averageTrajectory);
save(fname, 'analyticalTrajectories', '-append');

toc

%%% CM, asexual
clearvars
tic

% Declare input parameters
simParams = initializeSimParams('numIteration', 200, ...
    'mutationRate', 2*10^-4);

[resultModularCM] = simulateModularCM(simParams);

numTimeStamp = 20;
[averageTrajectory] = computeAverageTrajectory(numTimeStamp, simParams, resultModularCM.resultTable);

outputDir = 'results/CM_Asexual';
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
fname = sprintf('%s/ModularFGM_CM_Asexual_M%.1e.mat', outputDir, simParams.mutationRate);

save(fname, 'simParams', 'resultModularCM', 'averageTrajectory');

[analyticalTrajectories] = predictModularCM(simParams);
save(fname, 'analyticalTrajectories', '-append');

toc

%%% CM, sexual
clearvars
tic

% Declare input parameters
simParams = initializeSimParams('numIteration', 200, ...
    'mutationRate', 2*10^-4, 'recombinationRate', 1);
initialParameters = freezeParam(simParams);

numSteps = 200;
simParams.populationMatrices = preRunSimulation(initialParameters, simParams, numSteps, 'rngSeed', 1);
[resultModularCM] = simulateModularCM(simParams);

numTimeStamp = 20;
[averageTrajectory] = computeAverageTrajectory(numTimeStamp, simParams, resultModularCM.resultTable);

outputDir = 'results/CM_Sexual';
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
fname = sprintf('%s/ModularFGM_CM_Sexual_M%.1e_R%.4f.mat', ...
    outputDir, simParams.mutationRate, simParams.recombinationRate);

save(fname, 'simParams', 'resultModularCM', 'averageTrajectory');

[analyticalTrajectories] = predictFullRecomb(simParams, averageTrajectory);
save(fname, 'analyticalTrajectories', '-append');

toc

%% Additional Modular GPFMs Simulations