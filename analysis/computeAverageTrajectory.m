function [averageTrajectory] = computeAverageTrajectory(numTimeStamp, simParams, resultTable)
% computeAverageTrajectory - Calculate mean trajectory from multiple simulation runs.
%
% Description:
%   Computes the average evolutionary trajectory across replicate simulations,
%   interpolating to uniform time stamps and calculating log(x2/x1) ratios.
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
%       .averageTimeStamp - Mean trajectory [2 x numTimeStamp]
%       .timeRecord       - Time points for each replicate
%       .logRatio         - log(trait2/trait1) over time
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

    % Initialize the structure to store mean trajectories
    averageTrajectory = struct();
    
    % Loop through each set of initial angles
    for i = 1:length(simParams.initialAngles)
        % Prepare matrices to hold the timestamped data for each simulation run
        numSimulations = length(resultTable(i,:));
        trait1 = zeros(numSimulations, numTimeStamp);
        trait2 = zeros(numSimulations, numTimeStamp);
        timeRecord = zeros(numSimulations, numTimeStamp);
        
        % Loop through each simulation result for the current set of initial angles
        for j = 1:numSimulations
            % Extract the trajectory data from the current simulation
            phenotypicTrajectory = resultTable{i, j};
            
            % Create evenly spaced timestamps for the trajectory
            uniformTimeStamp = floor(linspace(1, length(phenotypicTrajectory(:,1)), numTimeStamp));
            
            % Record trait and time data at each timestamp
            timeRecord(j, :) = uniformTimeStamp;
            trait1(j, :) = phenotypicTrajectory(uniformTimeStamp, 3); % Trait 1 values
            trait2(j, :) = phenotypicTrajectory(uniformTimeStamp, 4); % Trait 2 values
        end
        
        % Compute the mean trajectory at each timestamp
        meanTrait1 = mean(trait1, 1);
        meanTrait2 = mean(trait2, 1);
        averageTimeStamp = [meanTrait1; meanTrait2];
        
        % Compute log ratio log(trait2/trait1), setting invalid values to NaN and issuing a warning
        logRatio = zeros(1, numTimeStamp);
        for k = 1:numTimeStamp
            ratio = meanTrait2(k) / meanTrait1(k);
            if ratio <= 0
                warning('Non-positive ratio encountered at timestamp %d. Setting log(trait2/trait1) to NaN.', k);
                logRatio(k) = NaN; % Set to NaN if ratio is non-positive
            else
                logRatio(k) = log(ratio);
            end
        end
        
        % Store the mean trajectories, log ratios, and individual records in the output structure
        averageTrajectory.trait1{i} = trait1;
        averageTrajectory.trait2{i} = trait2;
        averageTrajectory.averageTimeStamp{i} = averageTimeStamp;
        averageTrajectory.timeRecord{i} = timeRecord;
        averageTrajectory.logRatio{i} = logRatio; % Save log ratio in the structure
    end
end
