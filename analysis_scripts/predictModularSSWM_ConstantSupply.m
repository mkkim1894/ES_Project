function [analyticalTrajectories, decayRates] = predictModularSSWM_ConstantSupply(simParams, averageTrajectory)
% predictModularSSWM_ConstantSupply - Analytical trajectories for the modular SSWM
%   model with a constant supply of module-improving mutations.
%
% Description:
%   Derivation. Module i is below its optimum at x_i < 0 and
%       log W = -[ (x_1/a_1)^2 + (x_2/a_2)^2 ] / (2*sigma^2),
%   so a +delta step in module i has selection coefficient
%       s_i = -(2*x_i*delta + delta^2) / (2*sigma^2*a_i^2) ~ |x_i|*delta / (sigma^2*a_i^2)
%   and fixes with probability P_fix ~ 2*s_i when N*s_i >> 1.
%
%   With a CONSTANT beneficial fraction f_i, the beneficial supply per genome per
%   generation is U_i = mu_i * f_i * L_i = U*f_i/2 with mu_i = U/(2*L_i), so
%   beneficial mutations fix at rate N*U_i*2*s_i and each moves the trait by delta:
%
%       dx_i/dt = N * U_i * 2*s_i * delta = -beta_i * x_i,
%       beta_i  = 2*N*U_i*delta^2 / (sigma^2*a_i^2)
%               = N*U*f_i*delta^2 / (sigma^2*a_i^2).
%
%   Hence each trait decays EXPONENTIALLY,
%       x_i(t) = x_i(0) * exp(-beta_i * t),
%   the phenotypic trajectory is the power-law curve
%       x_2 = x_2(0) * ( x_1 / x_1(0) )^(beta_2/beta_1),
%   and the log ratio of module performances is linear in time,
%       log( x_2(t)/x_1(t) ) = log( x_2(0)/x_1(0) ) - (beta_2 - beta_1) * t,
%   diverging to -Inf or +Inf according to the sign of beta_2 - beta_1, with
%       beta_2 / beta_1 = (f_2/f_1) * (a_1^2/a_2^2).
%   There is no module-selection balance: no attractor, and permanent memory of
%   the initial condition.
%
%   The supply is constant everywhere in trait space and is never switched off, so
%   the exponential above holds for as long as s_i > 0, i.e. while |x_i| >> delta.
%   A module then stops advancing not because it runs out of mutations but because
%   s_i -> 0, which is the alternative mechanism: with constant
%   evolvability, the RATE of evolution decays as the optimum is approached because
%   the selection gradient decays. This closed form is not expected to hold within
%   about one delta of the optimum, where the discreteness of the lattice and the
%   vanishing of s_i take over.
%
%   Contrast with predictModularSSWM, where the supply declines as
%   b_i = |x_i|/delta, the flux is quadratic in |x_i|, both traits decay as a power
%   law x_i ~ -1/(alpha_i*t), and the ratio converges to the attractor
%   x_2/x_1 -> alpha_1/alpha_2 = (L_2*a_2^2)/(L_1*a_1^2), at which s_1 = s_2.
%
%   Knife-edge case. If f_2/f_1 = a_2^2/a_1^2 exactly, then beta_1 = beta_2 and the
%   ratio is constant in time. This is NOT module-selection balance: each initial
%   condition simply retains its own ratio forever, so the trajectories remain a
%   family of distinct rays instead of collapsing onto a single one.
%
% Inputs:
%   simParams         - Structure containing simulation parameters, including
%                       .beneficialFraction = [f1, f2]. The initial phenotypes are
%                       assumed already lattice-discretized (see the Run_ script).
%   averageTrajectory - Structure containing averaged simulation trajectories,
%                       used only to decide how far along the curve to draw.
%
% Outputs:
%   analyticalTrajectories - Cell array {nAngles x 1} of predicted [x1, x2] curves
%   decayRates             - Structure containing:
%       .beta          - [nAngles x 2] decay rates beta_i (per generation)
%       .logRatioSlope - [nAngles x 1] slope of log(x2/x1) vs t, = -(beta_2-beta_1)
%       .logRatio0     - [nAngles x 1] intercept log(x_2(0)/x_1(0))
%       .tFinal        - [nAngles x 1] time horizon used for each curve
%
% Reference:
%   Kim, M., Ardell, S. M., & Kryazhimskiy, S. (2025).
%   "Module-Selection Balance in the Evolution of Modular Organisms."
%
% See also: predictModularSSWM, simulateModularSSWM_ConstantSupply
%
% Copyright (c) 2025 Minkyu Kim, Cornell University
% Licensed under MIT License

    initialAngles     = simParams.initialAngles;
    initialPhenotypes = simParams.initialPhenotypes;
    nAngles           = length(initialAngles);

    f = simParams.beneficialFraction(:)';

    analyticalTrajectories = cell(nAngles, 1);

    decayRates.beta          = zeros(nAngles, 2);
    decayRates.logRatioSlope = zeros(nAngles, 1);
    decayRates.logRatio0     = zeros(nAngles, 1);
    decayRates.tFinal        = zeros(nAngles, 1);

    N     = simParams.popSize;
    U     = simParams.mutationRate;
    delta = simParams.deltaTrait;
    sigma = simParams.landscapeStdDev;
    a     = simParams.ellipseParams;

    % beta_i = N*U*f_i*delta^2 / (sigma^2*a_i^2); independent of L_i and of x_i(0)
    beta1 = N * U * f(1) * delta^2 / (sigma^2 * a(1)^2);
    beta2 = N * U * f(2) * delta^2 / (sigma^2 * a(2)^2);

    for i_pos = 1:nAngles
        WT  = initialPhenotypes(i_pos, 1:end);
        x10 = WT(1);
        x20 = WT(2);

        decayRates.beta(i_pos, 1:2)     = [beta1, beta2];
        decayRates.logRatioSlope(i_pos) = -(beta2 - beta1);
        if x10 ~= 0 && x20 ~= 0
            decayRates.logRatio0(i_pos) = log(abs(x20) / abs(x10));
        else
            decayRates.logRatio0(i_pos) = NaN;
        end

        % Draw the curve out to wherever the simulated average trajectory ended.
        x1_final = averageTrajectory.averageTimeStamp{1,i_pos}(1,end);
        x2_final = averageTrajectory.averageTimeStamp{1,i_pos}(2,end);

        t1 = Inf;
        t2 = Inf;
        if beta1 > 0 && x1_final ~= 0
            t1 = log(abs(x10) / abs(x1_final)) / beta1;
        end
        if beta2 > 0 && x2_final ~= 0
            t2 = log(abs(x20) / abs(x2_final)) / beta2;
        end

        t_final = min([t1, t2]);
        if ~isfinite(t_final) || t_final <= 0
            warning('predictModularSSWM_ConstantSupply:badHorizon', ...
                    'Calculated time horizon is not positive for angle %d; skipping.', i_pos);
            decayRates.tFinal(i_pos) = NaN;
            continue;
        end
        decayRates.tFinal(i_pos) = t_final;

        t  = linspace(0, t_final, 1000);
        x1 = x10 * exp(-beta1 * t);
        x2 = x20 * exp(-beta2 * t);

        analyticalTrajectories{i_pos, 1} = [x1; x2]';
    end
end
