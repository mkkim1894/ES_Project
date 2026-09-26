function [tOut, logROut] = predictLogRatioTrajectory_ConstantSupply(simParams, regime, tMax)
% predictLogRatioTrajectory_ConstantSupply - Theoretical trajectory of
%   log R = log(x_2/x_1) against generations for the constant-supply model.
%
% Description:
%   Constant-supply counterpart of predictLogRatioTrajectory. The trait-space
%   panels of the constant-supply figure already carry analytical curves, but
%   those are time-free paths in the (x_1, x_2) plane. This function supplies
%   the matching time-resolved prediction so that the log-R panels can be
%   overlaid with theory too.
%
%   This overlay is what makes the regime contrast legible. The three regimes
%   do not merely differ in how fast R changes, they differ in whether R
%   approaches anything at all:
%
%     Successive mutations. Each trait decays exponentially at its own rate,
%       x_i(t) = x_i(0) exp(-beta_i t),   beta_i = 2 N U f_i delta^2 / a_i^2,
%     so log R is LINEAR in time,
%       log R(t) = log R_0 - (beta_2 - beta_1) t,
%     and runs off to -Inf or +Inf according to the sign of beta_2 - beta_1.
%     There is no fixed point: the six predicted lines stay parallel and never
%     converge.
%
%     Concurrent mutations, complete linkage. The two modules are coupled
%     through clonal interference, not through their mutational supplies. The
%     coupling alone is enough to produce a fixed point, which sits where the
%     two modules would adapt at the same rate IN ISOLATION,
%       f_DF(s_1, U_1, N) = f_DF(s_2, U_2, N),
%     rather than where their selection coefficients are equal. With U_1 =/= U_2
%     that is a different line from s_1 = s_2, so the balance exists but is not
%     the equal fitness benefits line of the declining-supply model. The
%     predicted curves approach it from both sides.
%
%     Concurrent mutations, free reassortment. Recombination removes the
%     between-module coupling, the two modules integrate independently, and log
%     R again drifts without approaching a limit.
%
%   Curves are truncated where the underlying approximation stops holding
%   (predictModularCM_ConstantSupply and predictFullRecomb_ConstantSupply now
%   drop their trailing non-finite rows), so the spurious late upturn that
%   appears when the Desai-Fisher form is pushed past its domain is not drawn.
%
% Units:
%   The manuscript writes F = -sum_i (|x_i|/a_i)^2, the code writes the
%   equivalent Gaussian form with sigma and the rescaled axes ellipseParams,
%   and the two are related by a_i = sqrt(2)*sigma*ellipseParams(i). Only the
%   successive-mutations branch uses closed-form rate constants, and it is
%   evaluated in manuscript-convention a_i so that the slope against generations
%   is right in absolute terms, not just in ratio. The two concurrent branches
%   delegate to the prediction functions, which work in code units throughout.
%
% Inputs:
%   simParams - Simulation parameter structure as saved with the results. Must
%               carry .beneficialFraction = [f1, f2]. Initial phenotypes are
%               lattice-discretized here, matching the simulations.
%   regime    - 'SSWM', 'CM_Asexual' or 'CM_Sexual'
%   tMax      - Largest generation to draw out to, typically the x-limit of the
%               panel being plotted.
%
% Outputs:
%   tOut, logROut - Cell arrays {nAngles x 1} of time vectors and log R values.
%                   Both are returned empty when no prediction is available.
%
% Reference:
%   Kim, M., Ardell, S. M., & Kryazhimskiy, S. (2025).
%   "Module-Selection Balance in the Evolution of Modular Organisms."
%
% See also: predictLogRatioTrajectory, predictModularSSWM_ConstantSupply,
%           predictModularCM_ConstantSupply, predictFullRecomb_ConstantSupply
%
% Copyright (c) 2025 Minkyu Kim, Cornell University
% Licensed under MIT License

    tOut    = {};
    logROut = {};

    if nargin < 3 || isempty(tMax) || ~isfinite(tMax) || tMax <= 0
        return;
    end
    if ~isfield(simParams, 'beneficialFraction') || isempty(simParams.beneficialFraction)
        return;
    end

    delta = simParams.deltaTrait;

    % Lattice-discretize the initial phenotypes, as the simulations do, so that
    % the theory starts from the same point the simulated populations did.
    X0 = min(-delta * round(-simParams.initialPhenotypes / delta), 0);
    nAngles = size(X0, 1);

    tOut    = cell(nAngles, 1);
    logROut = cell(nAngles, 1);

    % Parameters passed on to the prediction functions start from the same
    % discretized initial condition.
    spd = simParams;
    spd.initialPhenotypes = X0;

    nPts = 400;

    switch upper(regime)

    % ------------------------------------------------------------------
    case 'SSWM'
        N     = simParams.popSize;
        U     = simParams.mutationRate;
        sigma = simParams.landscapeStdDev;
        f     = simParams.beneficialFraction(:)';

        % Manuscript-convention selection parameters
        a = sqrt(2) * sigma * simParams.ellipseParams;

        % beta_i = 2*N*U*f_i*delta^2 / a_i^2, independent of x_i(0) and of L_i
        beta = 2 * N * U .* f .* delta^2 ./ a.^2;

        for j = 1:nAngles
            if X0(j,1) == 0 || X0(j,2) == 0, continue; end
            t = linspace(0, tMax, nPts);
            tOut{j}    = t;
            logROut{j} = log(abs(X0(j,2)/X0(j,1))) - (beta(2) - beta(1)) * t;
        end

    % ------------------------------------------------------------------
    case 'CM_ASEXUAL'
        % Default 'tol' of 10 stops the integration while both modules still
        % satisfy s_i >= 10*U_i, i.e. inside the Desai-Fisher domain.
        [traj, tvec] = predictModularCM_ConstantSupply(spd);
        [tOut, logROut] = odeToLogRatio(traj, tvec, nAngles, tMax);

    % ------------------------------------------------------------------
    case 'CM_SEXUAL'
        % Same validity margin as the linked case.
        [traj, tvec] = predictFullRecomb_ConstantSupply(spd);
        [tOut, logROut] = odeToLogRatio(traj, tvec, nAngles, tMax);

    % ------------------------------------------------------------------
    otherwise
        tOut = {}; logROut = {};
        return;
    end
end

% ======================================================================
function [tOut, logROut] = odeToLogRatio(traj, tvec, nAngles, tMax)
% Convert an ode45 solution into log R against generations, keeping only the
% stretch where both modules are strictly below their optima and the solution
% is finite.

    tOut    = cell(nAngles, 1);
    logROut = cell(nAngles, 1);

    for j = 1:nAngles
        if j > numel(traj) || isempty(traj{j,1}) || isempty(tvec{j,1})
            continue;
        end
        Xj = traj{j,1};
        tj = tvec{j,1};

        ok = all(isfinite(Xj), 2) & Xj(:,1) < 0 & Xj(:,2) < 0 & tj(:) <= tMax;
        if nnz(ok) < 2
            continue;
        end

        tOut{j}    = tj(ok)';
        logROut{j} = log(abs(Xj(ok,2) ./ Xj(ok,1)))';
    end
end
