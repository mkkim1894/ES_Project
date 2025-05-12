function [analyticalTrajectories] = predictPleiotropicSSWM(simParams, averageTrajectory)
    analyticalTrajectories = cell(1);
    initialAngles = simParams.initialAngles;
    initialPhenotypes = simParams.initialPhenotypes;

    for i_pos = 1:length(initialAngles)
        % Extract initial phenotypic coordinates
        x10 = initialPhenotypes(i_pos, 1);
        x20 = initialPhenotypes(i_pos, 2);

        % Extract final coordinates from averageTrajectory
        x1_final = averageTrajectory.averageTimeStamp{1, i_pos}(1, end);
        x2_final = averageTrajectory.averageTimeStamp{1, i_pos}(2, end);

        % Calculate the distance of the final phenotype from the origin
        finalDistance = sqrt(x1_final^2 + x2_final^2);

        % Generate the trajectory in the trait space
        x1_values = linspace(x10, x1_final, 1000); % Ensure the range goes up to x1_final
        x2_values = x20 * (x1_values / x10).^(simParams.ellipseParams(1)^2 / simParams.ellipseParams(2)^2);

        % Calculate the Euclidean distance for the predicted trajectory
        distances = sqrt(x1_values.^2 + x2_values.^2);

        % Stop the trajectory where the distance matches the finalDistance
        valid_indices = distances >= finalDistance;
        x1_values = x1_values(valid_indices);
        x2_values = x2_values(valid_indices);

        % Store the analytical trajectory
        analyticalTrajectories{i_pos, 1} = [x1_values; x2_values]';
    end
end
