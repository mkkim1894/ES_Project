function [averageTrajectory] = computeAverageTrajectory(numTimeStamp, simParams, resultTable)
% computeMeanTrajectory - Calculates the mean trajectory of traits over time
% from multiple simulation runs in a 2D Fisher's Geometric Model (FGM).
%
% Inputs:
%   numTimeStamp - The number of evenly spaced timestamps to record for each trajectory.
%   simParams - A structure containing the simulation parameters (e.g., initial angles).
%   resultTable - A structure containing the result table from the simulations.
%
% Outputs:
%   averageTrajectory - A structure containing the recorded and averaged trajectories for each trait.
%
% Example Usage:
%   [averageTrajectory] = computeMeanTrajectory(20, simParams, resultModularSSWM.resultTable);

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
