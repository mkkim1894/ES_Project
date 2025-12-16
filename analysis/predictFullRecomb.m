function [analyticalTrajectories] = predictFullRecomb(simParams, averageTrajectory)
% predictFullRecomb - Compute analytical predictions for full recombination regime.
%
% Description:
%   Calculates the expected evolutionary trajectory when recombination is
%   complete (rho = 1), allowing modules to evolve independently.
%
% Inputs:
%   simParams         - Structure containing simulation parameters
%   averageTrajectory - Structure containing averaged simulation trajectories
%
% Outputs:
%   analyticalTrajectories - Cell array of predicted [x1, x2] trajectories
%
% Reference:
%   Kim, M., Ardell, S. M., & Kryazhimskiy, S. (2025).
%   "Module-Selection Balance in the Evolution of Modular Organisms."
%
% See also: predictModularCM, simulateModularCM
%
% Copyright (c) 2025 Minkyu Kim, Cornell University
% Licensed under MIT License

    analyticalTrajectories = cell(1);
    initialAngles = simParams.initialAngles;
    initialPhenotypes = simParams.initialPhenotypes;

    for i_pos = 1:length(initialAngles)
        WT = initialPhenotypes(i_pos,1:end);
        x1_final = averageTrajectory.averageTimeStamp{1,i_pos}(1,end);
        x2_final = averageTrajectory.averageTimeStamp{1,i_pos}(2,end);
           
        A1 = log( (2 * simParams.popSize^2 * simParams.mutationRate)/...
            (simParams.geneticTargetSize(1)*simParams.ellipseParams(1)^2) );  
        A2 = log( (2 * simParams.popSize^2 * simParams.mutationRate)/...
            (simParams.geneticTargetSize(2)*simParams.ellipseParams(2)^2) ); 
 
        alpha1 = (2*simParams.deltaTrait^2)/...
            ( simParams.ellipseParams(1)^2 * (log((2*simParams.geneticTargetSize(1)*simParams.deltaTrait^2)/...
            (simParams.mutationRate*simParams.ellipseParams(1)^2)))^2 );
        
        alpha2 = (2*simParams.deltaTrait^2)/...
            ( simParams.ellipseParams(2)^2 * (log((2*simParams.geneticTargetSize(2)*simParams.deltaTrait^2)/...
            (simParams.mutationRate*simParams.ellipseParams(2)^2)))^2 );
        
        x10 = WT(1);
        x20 = WT(2);
        
        end_norm = norm([x1_final, x2_final]);

        % Initialize the first value of x1
        current_x1 = x10;
        trajectoryRecord = [];
    
        % Loop until the distance to the origin is smaller than the norm
        while true
            % Compute x2 based on the analytical solution
            x2 = -exp(0.5 * ( ((2*log(abs(current_x1)) + A1) / (2*log(abs(x10)) + A1)).^(alpha2 / alpha1) * (2*log(abs(x20)) + A2) - A2 ));
    
            % Calculate the current distance to the origin
            dist_to_origin = norm([current_x1, x2]);

            % Stop if the distance is smaller than the end norm
            if dist_to_origin < end_norm
                break;
            end
    
           % Append current point to the trajectory record
           trajectoryRecord = [trajectoryRecord; current_x1, x2];

           % Update current_x1 by a small step 
           current_x1 = current_x1 + simParams.deltaTrait/10;
        end

        analyticalTrajectories{i_pos, 1} = trajectoryRecord;
    end
end
