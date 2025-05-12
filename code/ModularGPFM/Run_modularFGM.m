%% Base
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

%% CM, Intermediate recombinationRates
clearvars
tic

recombinationRateList = [0, 0.0002, 0.001, 0.01, 0.1, 1];
outputDir = 'results/CM_IntermediateRecomb';
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

for r_i = 1:length(recombinationRateList)  
    % Declare input parameters
    simParams = initializeSimParams('numIteration', 200, ...
        'mutationRate', 2*10^-4, 'recombinationRate', recombinationRateList(r_i),...
        'initialAngles', [atan2(1, 6.4), atan2(1, 0.2)]);
    initialParameters = freezeParam(simParams);

    numSteps = 200;
    simParams.populationMatrices = preRunSimulation(initialParameters, simParams, numSteps, 'rngSeed', 1);
    [resultModularCM] = simulateModularCM(simParams);

    numTimeStamp = 20;
    [averageTrajectory] = computeAverageTrajectory(numTimeStamp, simParams, resultModularCM.resultTable);

    fname = sprintf('%s/ModularFGM_CM_Recomb_R%.4f_M%.1e.mat', ...
        outputDir, recombinationRateList(r_i), simParams.mutationRate);

    save(fname, 'simParams', 'resultModularCM', 'averageTrajectory');

    [analyticalTrajectories] = predictFullRecomb(simParams, averageTrajectory);
    save(fname, 'analyticalTrajectories', '-append');
end

toc

%% Generalization
%%% SSWM
clearvars
clc

tic

% Declare input parameters
simParams = initializeSimParams('geneticTargetSize', [10, 40]);

[resultModularSSWM] = simulateModularSSWM(simParams);

numTimeStamp = 20;
[averageTrajectory] = computeAverageTrajectory(numTimeStamp, simParams, resultModularSSWM.resultTable);

outputDir = 'results/Generalization/SSWM';
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
fname = sprintf('%s/ModularFGM_SSWM_M%.1e_Target%d.mat', outputDir, simParams.mutationRate, prod(simParams.geneticTargetSize));

save(fname, 'simParams', 'resultModularSSWM', 'averageTrajectory');

[analyticalTrajectories] = predictModularSSWM(simParams, averageTrajectory);
save(fname, 'analyticalTrajectories', '-append');

toc

%%% CM, asexual
clearvars
tic

simParams = initializeSimParams('numIteration', 200, ...
    'mutationRate', 2*10^-4, 'geneticTargetSize', [10, 40]); % Try this mutation rate

[resultModularCM] = simulateModularCM(simParams);

numTimeStamp = 20;
[averageTrajectory] = computeAverageTrajectory(numTimeStamp, simParams, resultModularCM.resultTable);

outputDir = 'results/Generalization/CM_Asexual';
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
fname = sprintf('%s/ModularFGM_CM_Asexual_M%.1e_Target%d.mat', outputDir, simParams.mutationRate, prod(simParams.geneticTargetSize));

save(fname, 'simParams', 'resultModularCM', 'averageTrajectory');

[analyticalTrajectories] = predictModularCM(simParams);
save(fname, 'analyticalTrajectories', '-append');

toc

%%% CM, sexual
clearvars
tic

simParams = initializeSimParams('numIteration', 200, ...
    'mutationRate', 2*10^-4, 'recombinationRate', 1, ...
    'geneticTargetSize', [10, 40]);
initialParameters = freezeParam(simParams);

numSteps = 200;
simParams.populationMatrices = preRunSimulation(initialParameters, simParams, numSteps, 'rngSeed', 1);
[resultModularCM] = simulateModularCM(simParams);

numTimeStamp = 20;
[averageTrajectory] = computeAverageTrajectory(numTimeStamp, simParams, resultModularCM.resultTable);

outputDir = 'results/Generalization/CM_Sexual';
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
fname = sprintf('%s/ModularFGM_CM_Sexual_M%.1e_R%.4f_Target%d.mat', ...
    outputDir, simParams.mutationRate, simParams.recombinationRate, prod(simParams.geneticTargetSize));

save(fname, 'simParams', 'resultModularCM', 'averageTrajectory');

[analyticalTrajectories] = predictFullRecomb(simParams, averageTrajectory);
save(fname, 'analyticalTrajectories', '-append');

toc

%% NestedFGM
%%% SSWM
clearvars
clc

tic

% Declare input parameters
simParams = initializeSimParams();
simParams.moduleDimension = [10,10];

% Run simulation
[resultNestedSSWM] = simulateNestedSSWM(simParams);

% Compute average trajectory
numTimeStamp = 20;
[averageTrajectory] = computeAverageTrajectory(numTimeStamp, simParams, resultNestedSSWM.resultTable);

% Define directory and filename
outputDir = 'results/Generalization/nestedFGM/SSWM';
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
fname = sprintf('%s/NestedFGM_SSWM_M%.1e.mat', outputDir, simParams.mutationRate);

% Save results
save(fname, 'simParams', 'resultNestedSSWM', 'averageTrajectory');

toc

%%% CM, asexual
clearvars
tic

% Declare input parameters
simParams = initializeSimParams('numIteration', 200, ...
    'mutationRate', 2*10^-4);
simParams.moduleDimension = [10,10];

% Run simulation
[resultNestedCM] = simulateNestedCM(simParams);

% Compute average trajectory
numTimeStamp = 20;
[averageTrajectory] = computeAverageTrajectory(numTimeStamp, simParams, resultNestedCM.resultTable);

% Define directory and filename
outputDir = 'results/Generalization/nestedFGM/CM_Asexual';
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
fname = sprintf('%s/NestedFGM_CM_Asexual_M%.1e.mat', outputDir, simParams.mutationRate);

% Save results
save(fname, 'simParams', 'resultNestedCM', 'averageTrajectory');

toc

%%% CM, sexual
clearvars
tic

% Declare input parameters
simParams = initializeSimParams('numIteration', 200, ...
    'mutationRate', 2*10^-4, ...
    'recombinationRate', 1);
simParams.moduleDimension = [10,10];
initialParameters = freezeParam(simParams);

% Pre-run simulation
numSteps = 200;
simParams.populationMatrices = preRunSimulation(initialParameters, simParams, numSteps, 'rngSeed', 1);

% Run simulation
[resultNestedCM] = simulateNestedCM(simParams);

% Compute average trajectory
numTimeStamp = 20;
[averageTrajectory] = computeAverageTrajectory(numTimeStamp, simParams, resultNestedCM.resultTable);

% Define directory and filename
outputDir = 'results/Generalization/nestedFGM/CM_Sexual';
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
fname = sprintf('%s/NestedFGM_CM_Sexual_M%.1e_R%.4f.mat', outputDir, simParams.mutationRate, simParams.recombinationRate);

% Save results
save(fname, 'simParams', 'resultNestedCM', 'averageTrajectory');

toc

%% Plot Figures
makeFigure3('ModularFGM', 1e-7, 2e-4)
makeFigure4(2e-4)
makeFigure5('ModularFGM', 1e-7, 2e-4)
makeFigure6('NestedFGM', 1e-7, 2e-4);

%% Additional Figures (SI Figs)
