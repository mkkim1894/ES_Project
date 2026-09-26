function [analyticalTrajectories, timeVectors] = predictModularCM_ConstantSupply(simParams, varargin)
% predictModularCM_ConstantSupply - Analytical trajectory predictions for the
%   modular CM regime (linked modules) with a constant mutational supply.
%
% Description:
%   Identical to predictModularCM except for a single line: the beneficial
%   mutation supply of each module is
%
%       U_i = mu_i * f_i * L_i = U * f_i / 2                      (constant)
%
%   with per-locus rate mu_i = U/(2*L_i) and constant beneficial fraction f_i,
%   rather than
%
%       U_i = U * |x_i| / (2 * delta * L_i)                       (declining)
%
%   as in predictModularCM. Everything downstream - the Desai-Fisher rate of
%   adaptation, the effective mutation rate obtained from the quadratic, the
%   piecewise treatment of the two modules when one dominates by a factor D, and
%   the event-terminated ode45 integration - is unchanged, so any difference in
%   the predicted trajectory is attributable to the supply alone.
%
%   With the declining supply, U_i and s_i both scale with |x_i|, which is what
%   allows the two modules to reach a fixed performance ratio (module-selection
%   balance). With a constant supply only s_i scales with |x_i|, so the module
%   under stronger selection (relative to its own beneficial fraction) keeps
%   adapting faster, and log(x_2/x_1) drifts without bound.
%
% Inputs:
%   simParams - Structure containing simulation parameters, including
%               .beneficialFraction = [f1, f2]
%   varargin  - Optional name-value pairs:
%       'D'   - Rate ratio threshold (default: 100)
%       'tol' - Validity tolerance, applied both to s_tilde/U_tilde inside
%               computeEffectiveMutationRate and, via the event function, to
%               each module's own s_i/U_i (default: 10)
%
%   DOMAIN OF VALIDITY. With a declining supply, s_i and U_i shrink together and
%   the Desai-Fisher requirement s_i/U_i >> 1 holds along the whole trajectory.
%   With a constant supply it does not: s_i/U_i falls steadily and the
%   approximation expires part-way in. Left unguarded, the integration runs to
%   the hard boundary where desaiFisher returns NaN, the affected module's rate
%   snaps to zero, and the frozen module makes log(x_2/x_1) reverse. The event
%   function below stops each trajectory while both modules still satisfy
%   s_i >= tol*U_i, so nothing is drawn outside the regime the formula describes.
%   The convergence this model shows is established well inside that boundary:
%   at s_2/U_2 = 10 the spread of log R across the six initial conditions has
%   already fallen from 3.41 to 0.50.
%
% Outputs:
%   analyticalTrajectories - Cell array of predicted [x1, x2] trajectories
%   timeVectors            - Cell array of the matching ode45 time vectors, in
%                            generations. Needed to plot the predicted R(t)
%                            against generations rather than only as a
%                            time-free curve in the (x1, x2) plane.
%
% Reference:
%   Kim, M., Ardell, S. M., & Kryazhimskiy, S. (2025).
%   "Module-Selection Balance in the Evolution of Modular Organisms."
%
% See also: predictModularCM, simulateModularCM_ConstantSupply
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

    if ~isfield(simParams, 'beneficialFraction') || isempty(simParams.beneficialFraction)
        error('predictModularCM_ConstantSupply:missingFraction', ...
              'simParams.beneficialFraction = [f1, f2] must be supplied.');
    end

    initialPhenotypes = simParams.initialPhenotypes;
    nPos = size(initialPhenotypes, 1);

    analyticalTrajectories = cell(nPos, 1);
    timeVectors            = cell(nPos, 1);
    finalFitnessThreshold = 0.99;

    for i_pos = 1:nPos
        WT = initialPhenotypes(i_pos, :);

        tspan = [0 10000];
        options = odeset('Events', @(t,x) eventFunction(t, x, simParams, ...
                                          finalFitnessThreshold, D, tol));

        [T, X] = ode45(@(t, x) computeAdaptationRate(t, x, simParams, D, tol), ...
                       tspan, WT, options);

        % The integration terminates on the NaN event, which leaves a trailing
        % non-finite row. Drop it here so that every consumer of this output
        % sees only points where the approximation was still valid.
        keep = all(isfinite(X), 2);
        analyticalTrajectories{i_pos, 1} = X(keep, :);
        timeVectors{i_pos, 1}            = T(keep);
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
    WT    = x;
    delta = simParams.deltaTrait;
    N     = simParams.popSize;
    a1    = simParams.ellipseParams(1);
    a2    = simParams.ellipseParams(2);
    sigW  = simParams.landscapeStdDev;
    U     = simParams.mutationRate;
    f     = simParams.beneficialFraction;

    % ---- The one substantive difference from predictModularCM ----
    % Constant supply per module: U_i = mu_i*f_i*L_i = U*f_i/2, independent of
    % x_i everywhere in trait space. It is never switched off; a module stops at
    % its optimum because s_i <= 0 there, which desaiFisher already returns NaN
    % for, not because the supply is altered.
    U1 = U * f(1) / 2;
    U2 = U * f(2) / 2;
    % --------------------------------------------------------------

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
% Event function - stop at the fitness threshold, at the edge of the
% Desai-Fisher domain, or where the rate law breaks down
%--------------------------------------------------------------------------
function [value, isterminal, direction] = eventFunction(~, x, simParams, ...
                                                        finalFitnessThreshold, D, tol)
    a1   = simParams.ellipseParams(1);
    a2   = simParams.ellipseParams(2);
    sigW = simParams.landscapeStdDev;
    d    = simParams.deltaTrait;
    U    = simParams.mutationRate;
    f    = simParams.beneficialFraction;

    logW0 = -((x(1) / a1)^2 + (x(2) / a2)^2) / (2 * sigW^2);
    fitness = exp(logW0);
    endCondition = fitness - finalFitnessThreshold;

    % Continuous validity margin: each module's selection coefficient against
    % its own constant supply. Stopping on this rather than on the NaN indicator
    % keeps the solver away from the discontinuity instead of making it
    % integrate across one, and it is a proper bracketable crossing rather than
    % a 0/1 flag that sits at zero for the whole run.
    s1 = -(2*x(1)*d + d^2) / (2 * sigW^2 * a1^2);
    s2 = -(2*x(2)*d + d^2) / (2 * sigW^2 * a2^2);
    validity = min(s1 - tol*U*f(1)/2, s2 - tol*U*f(2)/2);

    % Retained as a backstop for parameter combinations where the piecewise rate
    % law fails before the validity margin is reached.
    rij = computeAdaptationRate([], x, simParams, D, tol);
    isNaNCondition = any(isnan(rij));

    value      = [endCondition; validity; double(isNaNCondition) - 0.5];
    isterminal = [1; 1; 1];
    direction  = [0; -1; 1];
end
