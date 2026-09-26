function [tOut, logROut] = predictLogRatioTrajectory(simParams, regime, model, tMax)
% predictLogRatioTrajectory - Theoretical trajectory of log R = log(x_2/x_1)
%                             against generations.
%
% Description:
%   The trait-space panels of Figures 2-4 already show analytical trajectories,
%   but those are non-parametric curves in the (x_1, x_2) plane and carry no
%   time axis. This function supplies the missing time-resolved prediction for
%   the module performance ratio R(t) = x_2(t)/x_1(t), so that the log-R panels
%   can be overlaid with theory as well.
%
%   Two things this makes visible:
%     (i)  the late-time rise in the simulated mean log R is a small-|x_i|
%          estimation artifact, not a feature of the model; and
%     (ii) in the modular SSWM regime R(t) approaches the balance value
%          Rbar = a_2^2/a_1^2 only as t -> Inf, so at the end of a finite run the
%          THEORY itself has not reached the orange line either. Drawing it
%          turns an apparent disagreement into visible agreement.
%
% Units. The manuscript writes F = -sum_i (|x_i|/a_i)^2 while the code writes the
%   equivalent Gaussian form with sigma and the rescaled axes ellipseParams. The
%   two are related by a_i = sqrt(2)*sigma*ellipseParams(i), and every formula
%   below is evaluated with the manuscript-convention a_i obtained that way.
%   This is what makes the absolute rate constants (not just their ratios) come
%   out right when plotted against generations.
%
% Inputs:
%   simParams - Simulation parameter structure, as saved with the results.
%               Initial phenotypes are lattice-discretized here, matching the
%               simulations.
%   regime    - 'SSWM', 'CM_Asexual' or 'CM_Sexual'
%   model     - 'Pleiotropic' or 'Modular'
%   tMax      - Largest generation to draw out to (typically the x-limit of the
%               panel being plotted).
%
% Outputs:
%   tOut, logROut - Cell arrays {nAngles x 1} of time vectors and log R values.
%                   Both are returned empty when no time-resolved prediction
%                   exists for the requested combination (see below).
%
% Availability:
%   Pleiotropic, SSWM        - R(t) = R_0 exp(-(beta_2-beta_1) t),
%                              beta_i = N*U*delta^2/a_i^2
%   Pleiotropic, CM (either) - NOT AVAILABLE. The only pleiotropic theory is the
%                              successive-mutations result, used elsewhere as a
%                              time-free SHAPE reference in trait space. Plotting
%                              it against concurrent-mutations generations would
%                              put it on the wrong clock, so nothing is returned.
%   Modular, SSWM            - R(t) = R_0 (1-alpha_1 x_10 t)/(1-alpha_2 x_20 t),
%                              alpha_i = 4*N*mu*delta/a_i^2
%   Modular, CM_Asexual      - from the ode45 solution of predictModularCM,
%                              which now also returns its time vector
%   Modular, CM_Sexual       - closed form
%                              x_i(t) = -exp(-A_i/2 + (log(-x_i0)+A_i/2) e^{-2 gamma_i t}),
%                              valid while 2*log(-x_i) + A_i > 0
%
% Reference:
%   Kim, M., Ardell, S. M., & Kryazhimskiy, S. (2025).
%   "Module-Selection Balance in the Evolution of Modular Organisms."
%
% See also: predictModularSSWM, predictModularCM, predictFullRecomb
%
% Copyright (c) 2025 Minkyu Kim, Cornell University
% Licensed under MIT License

    tOut    = {};
    logROut = {};

    if nargin < 4 || isempty(tMax) || ~isfinite(tMax) || tMax <= 0
        return;
    end

    delta = simParams.deltaTrait;
    N     = simParams.popSize;
    U     = simParams.mutationRate;
    sigma = simParams.landscapeStdDev;

    % Manuscript-convention selection parameters
    a = sqrt(2) * sigma * simParams.ellipseParams;

    % Lattice-discretize the initial phenotypes, as the simulations do
    X0 = min(-delta * round(-simParams.initialPhenotypes / delta), 0);
    nAngles = size(X0, 1);

    tOut    = cell(nAngles, 1);
    logROut = cell(nAngles, 1);

    nPts = 400;

    switch lower(model)

    % ------------------------------------------------------------------
    case 'pleiotropic'
        if ~strcmpi(regime, 'SSWM')
            tOut = {}; logROut = {};   % see Availability note above
            return;
        end
        beta = N * U * delta^2 ./ a.^2;
        for j = 1:nAngles
            if X0(j,1) == 0 || X0(j,2) == 0, continue; end
            t = linspace(0, tMax, nPts);
            tOut{j}    = t;
            logROut{j} = log(abs(X0(j,2)/X0(j,1))) - (beta(2) - beta(1)) * t;
        end

    % ------------------------------------------------------------------
    case 'modular'
        L  = simParams.geneticTargetSize;
        mu = U ./ (2 * L);              % per-site rate, per chromosome

        switch upper(regime)

        case 'SSWM'
            alpha = 4 * N .* mu .* delta ./ a.^2;
            for j = 1:nAngles
                x10 = X0(j,1); x20 = X0(j,2);
                if x10 == 0 || x20 == 0, continue; end
                t   = linspace(0, tMax, nPts);
                num = 1 - alpha(1) * x10 * t;
                den = 1 - alpha(2) * x20 * t;
                ok  = (num > 0) & (den > 0);
                tOut{j}    = t(ok);
                logROut{j} = log(abs((x20/x10) * num(ok) ./ den(ok)));
            end

        case 'CM_ASEXUAL'
            % ode45 solution; predictModularCM returns its time vector as a
            % second output so the trajectory can be plotted against generations.
            spd = simParams;
            spd.initialPhenotypes = X0;
            [traj, tvec] = predictModularCM(spd);
            for j = 1:nAngles
                if isempty(traj{j,1}) || isempty(tvec{j,1}), continue; end
                Xj = traj{j,1};
                tj = tvec{j,1};
                ok = (Xj(:,1) < 0) & (Xj(:,2) < 0) & (tj <= tMax);
                if nnz(ok) < 2, continue; end
                tOut{j}    = tj(ok)';
                logROut{j} = log(abs(Xj(ok,2) ./ Xj(ok,1)))';
            end

        case 'CM_SEXUAL'
            gamma = 2*delta^2 ./ (a.^2 .* log(2*delta^2 ./ (mu .* a.^2)).^2);
            A     = log(2 * N^2 .* mu ./ a.^2);
            for j = 1:nAngles
                x10 = X0(j,1); x20 = X0(j,2);
                if x10 >= 0 || x20 >= 0, continue; end
                t  = linspace(0, tMax, nPts);
                x1 = -exp(-A(1)/2 + (log(-x10) + A(1)/2) * exp(-2*gamma(1)*t));
                x2 = -exp(-A(2)/2 + (log(-x20) + A(2)/2) * exp(-2*gamma(2)*t));
                % The concurrent-mutations approximation holds only while both
                % modules stay far enough from their optima.
                ok = (2*log(-x1) + A(1) > 0) & (2*log(-x2) + A(2) > 0);
                if nnz(ok) < 2, continue; end
                tOut{j}    = t(ok);
                logROut{j} = log(abs(x2(ok) ./ x1(ok)));
            end

        otherwise
            tOut = {}; logROut = {};
            return;
        end

    otherwise
        tOut = {}; logROut = {};
        return;
    end
end
