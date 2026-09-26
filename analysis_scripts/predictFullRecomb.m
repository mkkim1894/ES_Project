function [analyticalTrajectories] = predictFullRecomb(simParams, ~)
% predictFullRecomb - Compute analytical predictions for full recombination regime.
%
% Description:
%   Calculates the expected evolutionary trajectory when recombination is
%   complete (rho = 1), allowing modules to evolve independently.
%   The per-locus mutation rate is mu = U/(2*K_i), consistent with the
%   binary locus model used throughout. The second input is accepted for
%   interface compatibility with other predictor functions but is not used.
%
% Inputs:
%   simParams - Structure containing simulation parameters
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

    initialPhenotypes = simParams.initialPhenotypes;
    nPos = size(initialPhenotypes, 1);
    analyticalTrajectories = cell(nPos, 1);

    N     = simParams.popSize;
    U     = simParams.mutationRate;
    L     = simParams.geneticTargetSize;
    a     = simParams.ellipseParams;
    delta = simParams.deltaTrait;
    sigW  = simParams.landscapeStdDev;

    finalFitnessThreshold = 0.99;
    stepSize = delta / 10;
    maxIter = 1e6;

    % Full recombination analytical parameters.
    %
    % UNITS. The manuscript writes F = -sum_i (|x_i|/a_i)^2 while the code writes
    % the equivalent Gaussian form with sigma and the rescaled axes
    % ellipseParams, the two being related by a_i = sqrt(2)*sigma*ellipseParams(i).
    % Equations (gamma, A) in the manuscript are in the MANUSCRIPT convention, so
    % they must be evaluated with aMS, not with ellipseParams. Using
    % ellipseParams here understated every a_i^2 by a factor 2*sigma^2 = 8, which
    % shifted A_i from (5.52, 6.22) to (7.60, 8.29) and gamma2/gamma1 from 1.579
    % to 1.680, moving the drawn path by up to 29% in x_2 and running it roughly
    % 3x closer to the optimum before terminating. The log-ratio panel
    % (predictLogRatioTrajectory) already did the conversion, so panels B and D of
    % Figure 4 disagreed with each other. predictModularSSWM carries the same note.
    %
    % Note that `a` (= ellipseParams) is still the right thing to use in the
    % fitness expression below, because that one is written in the CODE
    % convention, with its own explicit 1/(2*sigW^2).
    aMS = sqrt(2) * sigW * a;

    A1 = log((N^2 * U) / (L(1) * aMS(1)^2));
    A2 = log((N^2 * U) / (L(2) * aMS(2)^2));

    gamma1 = (2 * delta^2) / ...
        (aMS(1)^2 * (log((4 * L(1) * delta^2) / (U * aMS(1)^2)))^2);
    gamma2 = (2 * delta^2) / ...
        (aMS(2)^2 * (log((4 * L(2) * delta^2) / (U * aMS(2)^2)))^2);

    for i_pos = 1:nPos
        WT = initialPhenotypes(i_pos, :);
        x10 = WT(1);
        x20 = WT(2);

        % Sanity checks
        if x10 >= 0 || x20 >= 0
            error('predictFullRecomb:InvalidInitialState', ...
                'Initial phenotypes must satisfy x10 < 0 and x20 < 0.');
        end

        denom1 = 2 * log(-x10) + A1;
        denom2 = 2 * log(-x20) + A2;

        if denom1 <= 0 || denom2 <= 0
            error('predictFullRecomb:InvalidLogDomain', ...
                'Initial condition falls outside the valid log-domain of the analytical solution.');
        end

        current_x1 = x10;
        trajectoryRecord = zeros(0, 2);

        iter = 0;
        while true
            iter = iter + 1;
            if iter > maxIter
                warning('predictFullRecomb:MaxIterReached', ...
                    'Maximum iteration count reached at initial condition %d.', i_pos);
                break;
            end

            if current_x1 >= 0
                break;
            end

            numer1 = 2 * log(-current_x1) + A1;
            if numer1 <= 0
                break;
            end

            % Compute x2 from the analytical full-recombination trajectory
            x2 = -exp(0.5 * ( ...
                (numer1 / denom1)^(gamma2 / gamma1) * denom2 - A2));

            if ~isfinite(x2) || x2 >= 0
                break;
            end

            % Same termination logic as predictModularCM
            logW = -((current_x1 / a(1))^2 + (x2 / a(2))^2) / (2 * sigW^2);
            fitness = exp(logW);

            if fitness >= finalFitnessThreshold
                break;
            end

            trajectoryRecord(end+1, :) = [current_x1, x2]; %#ok<AGROW>
            current_x1 = current_x1 + stepSize;
        end

        analyticalTrajectories{i_pos, 1} = trajectoryRecord;
    end
end