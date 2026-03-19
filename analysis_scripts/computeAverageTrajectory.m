function [averageTrajectory] = computeAverageTrajectory(numTimeStamp, simParams, resultTable)
% computeAverageTrajectory - Calculate mean trajectory from multiple simulation runs.
%
% Description:
%   Computes the average evolutionary trajectory across replicate simulations
%   by sampling at uniformly spaced row indices. For SSWM, rows correspond
%   to fixation events; for CM, rows correspond to generations.
%
% Inputs:
%   numTimeStamp - Number of evenly spaced timestamps to record
%   simParams    - Structure containing simulation parameters (e.g., initial angles)
%   resultTable  - Cell array containing trajectory data from simulations
%
% Outputs:
%   averageTrajectory - Structure containing:
%       .trait1           - Individual trait 1 trajectories
%       .trait2           - Individual trait 2 trajectories
%       .averageTimeStamp - Cell array {nAngles}, each [2 x numTimeStamp] mean trajectory
%       .timeRecord       - Time points for each replicate
%
% Example:
%   ave = computeAverageTrajectory(20, simParams, result.resultTable);
%
% Reference:
%   Kim, M., Ardell, S. M., & Kryazhimskiy, S. (2025).
%   "Module-Selection Balance in the Evolution of Modular Organisms."
%
% Copyright (c) 2025 Minkyu Kim, Cornell University
% Licensed under MIT License

    averageTrajectory = struct();
    
    for i = 1:length(simParams.initialAngles)
        numSimulations = length(resultTable(i,:));
        trait1 = zeros(numSimulations, numTimeStamp);
        trait2 = zeros(numSimulations, numTimeStamp);
        timeRecord = zeros(numSimulations, numTimeStamp);
        
        for j = 1:numSimulations
            phenotypicTrajectory = resultTable{i, j};
            uniformTimeStamp = floor(linspace(1, length(phenotypicTrajectory(:,1)), numTimeStamp));
            timeRecord(j, :) = uniformTimeStamp;
            trait1(j, :) = phenotypicTrajectory(uniformTimeStamp, 3);
            trait2(j, :) = phenotypicTrajectory(uniformTimeStamp, 4);
        end
        
        meanTrait1 = mean(trait1, 1);
        meanTrait2 = mean(trait2, 1);
        averageTimeStamp = [meanTrait1; meanTrait2];
        
        averageTrajectory.trait1{i} = trait1;
        averageTrajectory.trait2{i} = trait2;
        averageTrajectory.averageTimeStamp{i} = averageTimeStamp;
        averageTrajectory.timeRecord{i} = timeRecord;
    end
end