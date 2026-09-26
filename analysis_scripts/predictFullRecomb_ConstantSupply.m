function [analyticalTrajectories, timeVectors] = predictFullRecomb_ConstantSupply(simParams, ~, varargin)
% predictFullRecomb_ConstantSupply - Analytical predictions for the full
%   recombination regime with a constant mutational supply.
%
% Description:
%   Under complete reassortment of the two module chromosomes (rho = 1), the two
%   modules do not interfere with each other, so each adapts at its own
%   Desai-Fisher rate given its own beneficial supply and selection coefficient:
%
%       dx_i/dt = ( v_i / s_i ) * delta,      v_i = v_DF( s_i, U_i, N ),
%
%   with the constant beneficial supply
%
%       U_i = mu_i * f_i * L_i = U * f_i / 2,     mu_i = U / (2*L_i),
%
%   and s_i = -(2*x_i*delta + delta^2) / (2*sigma^2*a_i^2). Within-module clonal
%   interference is retained through v_DF; only the between-module interference is
%   removed by recombination.
%
%   DOMAIN OF VALIDITY - why this integration has to be stopped early.
%   The Desai-Fisher form requires s_i / U_i >> 1. With a DECLINING supply both
%   s_i and U_i shrink together as the module improves, so the ratio stays
%   roughly constant and the approximation holds all the way in; equation
%   "recomb R" carries its own explicit cut, 2*log(-x_i) + A_i > 0. With a
%   CONSTANT supply U_i does not shrink, so s_i / U_i falls steadily and the
%   approximation expires part-way along the trajectory. At the default
%   parameters module 2 reaches s_2 = U_2 at x_2 = -0.054.
%
%   Integrating past that point is not merely inaccurate, it is actively
%   misleading. desaiFisher returns NaN there and the rate for that module snaps
%   to zero, which (i) puts a hard discontinuity in the right-hand side that the
%   adaptive stepper has to hammer its way across, leaving a visible kink and a
%   dense cluster of output points, and (ii) freezes x_2 while x_1 keeps
%   improving, so log(x_2/x_1) REVERSES and climbs. Measured at the default
%   parameters, 49-74%% of the curve for the four most x_1-lagging initial
%   conditions lay past that boundary, log R for the first condition swung back
%   from -3.31 to -1.62, and all four converged onto the same final value purely
%   because they had all frozen at x_2 = -0.054 and were stopped by the W >= 0.99
%   rule at the same x_1. Plotted, that reads as a module-selection balance under
%   recombination. There is none: inside the domain of validity the spread of
%   log R across initial conditions is 3.15, against 3.41 at the start.
%
%   The integration is therefore terminated by a continuous event on
%   min_i (s_i - margin*U_i), with margin defaulting to 10 - the same tolerance
%   predictModularCM_ConstantSupply already applies to s_tilde/U_tilde.
%
%   predictFullRecomb solves the corresponding declining-supply problem in closed
%   form, using the fact that U_i and s_i both scale with |x_i| so that the
%   logarithm log(s_i/U_i) is nearly constant along the trajectory. That
%   simplification is not available here, because with a constant U_i the ratio
%   s_i/U_i varies by orders of magnitude along the trajectory. The two ODEs are
%   therefore integrated numerically instead; they are uncoupled, so this is
%   inexpensive and involves no additional approximation.
%
%   Because dx_i/dt depends only on x_i, the two modules approach their optima
%   independently and at different rates, and log(x_2/x_1) drifts monotonically
%   without approaching any attractor.
%
% Inputs:
%   simParams - Structure containing simulation parameters, including
%               .beneficialFraction = [f1, f2]. The second input is accepted for
%               interface compatibility with predictFullRecomb and is not used.
%   varargin  - Optional name-value pair:
%       'validityMargin' - Stop when min_i s_i/U_i falls to this value.
%                          Default 10. Set to 1 to integrate to the hard
%                          Desai-Fisher boundary, which is not recommended.
%
% Outputs:
%   analyticalTrajectories - Cell array of predicted [x1, x2] trajectories
%   timeVectors            - Cell array of the matching ode45 time vectors, in
%                            generations, so that the predicted R(t) can be
%                            plotted against generations.
%
% Reference:
%   Kim, M., Ardell, S. M., & Kryazhimskiy, S. (2025).
%   "Module-Selection Balance in the Evolution of Modular Organisms."
%
% See also: predictFullRecomb, predictModularCM_ConstantSupply,
%           simulateModularCM_ConstantSupply
%
% Copyright (c) 2025 Minkyu Kim, Cornell University
% Licensed under MIT License

    q = inputParser;
    addParameter(q, 'validityMargin', 10, @(v) isnumeric(v) && isscalar(v) && v >= 1);
    parse(q, varargin{:});
    validityMargin = q.Results.validityMargin;

    if ~isfield(simParams, 'beneficialFraction') || isempty(simParams.beneficialFraction)
        error('predictFullRecomb_ConstantSupply:missingFraction', ...
              'simParams.beneficialFraction = [f1, f2] must be supplied.');
    end

    initialPhenotypes = simParams.initialPhenotypes;
    nPos = size(initialPhenotypes, 1);
    analyticalTrajectories = cell(nPos, 1);
    timeVectors            = cell(nPos, 1);

    finalFitnessThreshold = 0.99;

    for i_pos = 1:nPos
        WT = initialPhenotypes(i_pos, :);

        if WT(1) >= 0 || WT(2) >= 0
            error('predictFullRecomb_ConstantSupply:InvalidInitialState', ...
                'Initial phenotypes must satisfy x10 < 0 and x20 < 0.');
        end

        tspan   = [0 1e7];
        options = odeset('Events', @(t,x) eventFunction(t, x, simParams, ...
                                          finalFitnessThreshold, validityMargin), ...
                         'RelTol', 1e-6, 'AbsTol', 1e-8);

        [T, X] = ode45(@(t, x) independentModuleRates(t, x, simParams), tspan, WT(:), options);

        % The event above already stops the integration inside the domain of
        % validity; this only removes any trailing row the solver may leave.
        keep = all(isfinite(X), 2) & X(:,1) < 0 & X(:,2) < 0;
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
% Uncoupled module dynamics under full reassortment
%--------------------------------------------------------------------------
function dxdt = independentModuleRates(~, x, simParams)
    delta = simParams.deltaTrait;
    N     = simParams.popSize;
    a     = simParams.ellipseParams;
    sigW  = simParams.landscapeStdDev;
    U     = simParams.mutationRate;
    f     = simParams.beneficialFraction;

    dxdt = [0; 0];

    for i = 1:2
        xi = x(i);

        % Constant supply for this module, independent of x_i everywhere in trait
        % space. There is deliberately no test on x_i here.
        Ui = U * f(i) / 2;

        % Selection coefficient of a +delta step in module i alone. This is what
        % stops the module at its optimum: s_i <= 0 once x_i >= -delta/2, and the
        % module then ceases to advance regardless of how much supply it has.
        si = -(2*xi*delta + delta^2) / (2 * sigW^2 * a(i)^2);

        vi = desaiFisher(si, Ui, N);
        if ~isfinite(vi) || vi <= 0 || si <= 0
            dxdt(i) = 0;
        else
            dxdt(i) = (vi / si) * delta;
        end
    end
end

%--------------------------------------------------------------------------
% Event function - stop at the fitness threshold or at the edge of the
% Desai-Fisher domain, whichever comes first
%--------------------------------------------------------------------------
function [value, isterminal, direction] = eventFunction(~, x, simParams, ...
                                                        finalFitnessThreshold, validityMargin)
    a    = simParams.ellipseParams;
    sigW = simParams.landscapeStdDev;
    d    = simParams.deltaTrait;
    U    = simParams.mutationRate;
    f    = simParams.beneficialFraction;

    logW0   = -((x(1) / a(1))^2 + (x(2) / a(2))^2) / (2 * sigW^2);
    fitness = exp(logW0);

    % Selection coefficient of a +delta step in each module, against that
    % module's constant beneficial supply. This is a CONTINUOUS function of x,
    % which matters: the previous version signalled stalling with
    % double(all(dxdt == 0)), a value that is identically 0 while the population
    % is moving. An event function that sits exactly on zero is degenerate - the
    % solver cannot bracket a crossing - and it was a large part of why the
    % integration crawled along at steps of ~0.1 generations and emitted a
    % visibly piecewise-linear curve.
    s1 = -(2*x(1)*d + d^2) / (2 * sigW^2 * a(1)^2);
    s2 = -(2*x(2)*d + d^2) / (2 * sigW^2 * a(2)^2);
    validity = min(s1 - validityMargin*U*f(1)/2, ...
                   s2 - validityMargin*U*f(2)/2);

    value      = [fitness - finalFitnessThreshold; validity];
    isterminal = [1; 1];
    direction  = [0; -1];
end
