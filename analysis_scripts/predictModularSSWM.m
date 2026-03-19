function [analyticalTrajectories] = predictModularSSWM(simParams, averageTrajectory)
% predictModularSSWM - Compute analytical trajectory predictions for modular SSWM regime.
%
% Description:
%   Calculates the expected evolutionary trajectory under SSWM dynamics
%   using the analytical solution for module-specific adaptation rates.
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
% See also: predictModularCM, simulateModularSSWM
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

        alpha1 = 4*simParams.popSize*simParams.mutationRate*simParams.deltaTrait/...
            (simParams.geneticTargetSize(1)*simParams.ellipseParams(1)^2);
        alpha2 = 4*simParams.popSize*simParams.mutationRate*simParams.deltaTrait/...
            (simParams.geneticTargetSize(2)*simParams.ellipseParams(2)^2);

        x10 = WT(1);
        x20 = WT(2);

        % Solve for t when x1 reaches x1_final
        t1 = (1 / (alpha1 * x10)) * (1 - x10 / x1_final);
        
        % Solve for t when x2 reaches x2_final
        t2 = (1 / (alpha2 * x20)) * (1 - x20 / x2_final);

        % Choose the smallest positive t
        t_final = min([t1, t2]);
        if t_final <= 0
            warning('Calculated time is not positive; check the input parameters.');
            continue;
        end

        % Create the analytical trajectory up to t_final
        t = linspace(0, t_final, 1000);
        x1 = x10 ./ (1 - (alpha1 * x10 * t));
        x2 = x20 ./ (1 - (alpha2 * x20 * t));

        analyticalTrajectories{i_pos, 1} = [x1; x2]';
    end
end
