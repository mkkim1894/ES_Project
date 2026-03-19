function [analyticalTrajectories] = predictModularCM(simParams, varargin)
% predictModularCM - Compute analytical trajectory predictions for modular CM regime.
%
% Description:
%   Solves ODEs for the expected evolutionary trajectory under concurrent
%   mutations, accounting for clonal interference between modules.
%   Uses the closed-form quadratic solution for effective mutation rates
%
% Inputs:
%   simParams - Structure containing simulation parameters
%   varargin  - Optional name-value pairs:
%       'D'   - Rate ratio threshold (default: 100)
%       'tol' - Validity tolerance for s_tilde/U_tilde (default: 10)
%
% Outputs:
%   analyticalTrajectories - Cell array of predicted [x1, x2] trajectories
%
% Reference:
%   Kim, M., Ardell, S. M., & Kryazhimskiy, S. (2025).
%   "Module-Selection Balance in the Evolution of Modular Organisms."
%
% See also: predictModularSSWM, simulateModularCM
%
% Copyright (c) 2025 Minkyu Kim, Cornell University
% Licensed under MIT License

    %% Parse optional arguments
    p = inputParser;
    addParameter(p, 'D', 100, @isnumeric);
    addParameter(p, 'tol', 10, @isnumeric);
    parse(p, varargin{:});

    D   = p.Results.D;
    tol = p.Results.tol;

    initialPhenotypes = simParams.initialPhenotypes;
    nPos = size(initialPhenotypes, 1);

    analyticalTrajectories = cell(nPos, 1);
    finalFitnessThreshold = 0.99;

    for i_pos = 1:nPos
        WT = initialPhenotypes(i_pos, :);

        tspan = [0 10000];
        options = odeset('Events', @(t,x) eventFunction(t, x, simParams, finalFitnessThreshold, D, tol));

        [~, X] = ode45(@(t, x) computeAdaptationRate(t, x, simParams, D, tol), ...
                       tspan, WT, options);

        analyticalTrajectories{i_pos, 1} = X;
    end
end

%--------------------------------------------------------------------------
% Desai-Fisher rate of adaptation
%--------------------------------------------------------------------------
function v = desaiFisher(s, Ub, N)
    if ~isfinite(s) || ~isfinite(Ub) || s <= 0 || Ub <= 0 || N <= 0 || s <= Ub
        v = NaN;
        return;
    end

    logRatio = log(s / Ub);
    if ~isfinite(logRatio) || logRatio <= 0
        v = NaN;
        return;
    end

    v = s^2 * (2 * log(N * s) - logRatio) / logRatio^2;
end

%--------------------------------------------------------------------------
% Compute effective mutation rate U_tilde via the quadratic formula
%--------------------------------------------------------------------------
function Ut = computeEffectiveMutationRate(s_tilde, v_prime, N, tol)
    if ~isfinite(s_tilde) || ~isfinite(v_prime) || s_tilde <= 0 || v_prime <= 0 || N <= 0
        Ut = NaN;
        return;
    end

    A = v_prime;
    B = s_tilde^2;
    C = -2 * s_tilde^2 * log(N * s_tilde);

    discriminant = B^2 - 4 * A * C;
    if ~isfinite(discriminant) || discriminant < 0
        Ut = NaN;
        return;
    end

    x_pos = (-B + sqrt(discriminant)) / (2 * A);
    if ~isfinite(x_pos) || x_pos <= 0
        Ut = NaN;
        return;
    end

    Ut = s_tilde * exp(-x_pos);

    if ~isfinite(Ut) || Ut <= 0 || s_tilde / Ut < tol
        Ut = NaN;
    end
end

%--------------------------------------------------------------------------
% ODE right-hand side
%--------------------------------------------------------------------------
function dxdt = computeAdaptationRate(~, x, simParams, D, tol)
    WT = x;
    delta = simParams.deltaTrait;
    N     = simParams.popSize;
    a1    = simParams.ellipseParams(1);
    a2    = simParams.ellipseParams(2);
    sigW  = simParams.landscapeStdDev;
    K1    = simParams.geneticTargetSize(1);
    K2    = simParams.geneticTargetSize(2);
    U     = simParams.mutationRate;

    % Beneficial mutation rates for each module
    U1 = U * abs(WT(1)) / (2 * delta * K1);
    U2 = U * abs(WT(2)) / (2 * delta * K2);

    % Selection coefficients
    logW0 = -((WT(1) / a1)^2 + (WT(2) / a2)^2) / (2 * sigW^2);
    logW1 = -(((WT(1) + delta) / a1)^2 + (WT(2) / a2)^2) / (2 * sigW^2);
    logW2 = -((WT(1) / a1)^2 + ((WT(2) + delta) / a2)^2) / (2 * sigW^2);

    s1 = logW1 - logW0;
    s2 = logW2 - logW0;

    % Rates of adaptation in isolation
    v1 = desaiFisher(s1, U1, N);
    v2 = desaiFisher(s2, U2, N);

    if any(~isfinite([s1, s2, v1, v2]))
        dxdt = [NaN; NaN];
        return;
    end

    % --- Piecewise approximation ---
    if v1 <= 0 || v2 <= 0
        dxdt = [NaN; NaN];
        return;
    end

    if v2 / v1 > D
        dxdt = [0; (v2 / s2) * delta];

    elseif v1 / v2 > D
        dxdt = [(v1 / s1) * delta; 0];

    else
        s_tilde = (s1^2 + s2^2) / (s1 + s2);

        U1_tilde = computeEffectiveMutationRate(s_tilde, v1, N, tol);
        U2_tilde = computeEffectiveMutationRate(s_tilde, v2, N, tol);

        if isnan(U1_tilde) || isnan(U2_tilde)
            dxdt = [NaN; NaN];
            return;
        end

        U_tilde = U1_tilde + U2_tilde;
        v_total = desaiFisher(s_tilde, U_tilde, N);

        if ~isfinite(v_total) || v_total <= 0
            dxdt = [NaN; NaN];
            return;
        end

        v_12 = (U1_tilde / U_tilde) * v_total;
        v_21 = (U2_tilde / U_tilde) * v_total;

        dxdt = [(v_12 / s1) * delta; (v_21 / s2) * delta];
    end
end

%--------------------------------------------------------------------------
% Event function — stop when fitness threshold reached or NaN encountered
%--------------------------------------------------------------------------
function [value, isterminal, direction] = eventFunction(~, x, simParams, finalFitnessThreshold, D, tol)
    a1   = simParams.ellipseParams(1);
    a2   = simParams.ellipseParams(2);
    sigW = simParams.landscapeStdDev;

    logW0 = -((x(1) / a1)^2 + (x(2) / a2)^2) / (2 * sigW^2);
    fitness = exp(logW0);
    endCondition = fitness - finalFitnessThreshold;

    rij = computeAdaptationRate([], x, simParams, D, tol);
    isNaNCondition = any(isnan(rij));

    value = [endCondition; double(isNaNCondition)];
    isterminal = [1; 1];
    direction  = [0; 0];
end