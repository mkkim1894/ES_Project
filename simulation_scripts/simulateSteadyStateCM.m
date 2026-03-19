function [WFsim, WFsim_vi, Parameter, Prediction, Prediction_U_weights, Prediction_s_weights] = simulateSteadyStateCM(N, varargin)
% simulateSteadyStateCM - Validate adaptation rate predictions under concurrent mutations.
%
% Description:
%   Runs Wright-Fisher simulations to measure adaptation rates (v1, v2) for two
%   modules under concurrent mutations, and compares against analytical predictions
%   with different weighting schemes. Tests the theory across a parameter space
%   of selection coefficients and mutation rates.
%
% Inputs:
%   N        - Population size (default: 10^4)
%   varargin - Optional name-value pairs:
%       'generations'  - Total generations to simulate (default: 15000)
%       't_burnin'     - Burn-in generations before measuring (default: 10000)
%       'n_v'          - Number of adaptation rate values to test (default: 6)
%       'n_a'          - Number of s/U ratio values to test (default: 4)
%       'savePath'     - Path to save results (default: current directory)
%
% Outputs:
%   WFsim                 - Structure with simulation results (v1_2, v2_1, F1_2, F2_1)
%   WFsim_vi              - Structure with single-trait validation results
%   Parameter             - Structure with parameter space (s1, s2, U1, U2)
%   Prediction            - Analytical predictions with equal weights
%   Prediction_U_weights  - Analytical predictions with U-based weights
%   Prediction_s_weights  - Analytical predictions with s-based weights
%
% Example:
%   [WFsim, ~, Parameter, Pred] = simulateSteadyStateCM(10^4);
%
% Reference:
%   Kim, M., Ardell, S. M., & Kryazhimskiy, S. (2025).
%   "Module-Selection Balance in the Evolution of Modular Organisms."
%
% See also: makeFigureS_SteadyStateCM, predictModularCM
%
% Copyright (c) 2025 Minkyu Kim, Cornell University
% Licensed under MIT License

%% Parse inputs
p = inputParser;
addRequired(p, 'N', @isnumeric);
addParameter(p, 'generations', 15000, @isnumeric);
addParameter(p, 't_burnin', 10000, @isnumeric);
addParameter(p, 'n_v', 6, @isnumeric);
addParameter(p, 'n_a', 4, @isnumeric);
addParameter(p, 'savePath', '', @ischar);
parse(p, N, varargin{:});

generations = p.Results.generations;
t_burnin = p.Results.t_burnin;
n_v = p.Results.n_v;
n_a = p.Results.n_a;
savePath = p.Results.savePath;

%% Parameter Setup
Parameter.N = N;

Parameter.v1 = exp(linspace(log(5*10^-6), log(1*10^-2), n_v));
Parameter.v2 = flip(exp(linspace(log(5*10^-6), log(1*10^-2), n_v)))';

% Create the full meshgrid
[X, Y] = meshgrid(Parameter.v1, Parameter.v2);

% Focus on the regime where v1 < v2
M = Y >= X;

X = M.*X;
v1 = X(X > 0);

Y = M.*Y;
v2 = Y(Y > 0);

% Set a, where a is s = a*U. a must be >> 1
a = round(exp(linspace(log(1*10^1), log(1*10^2), n_a)));
Parameter.a = a;

%% Create Parameter Space
fprintf('Creating parameter space...\n');
tic
for i = 1:numel(v1)
    target_v1 = v1(i);
    target_v2 = v2(i);
    for j = 1:numel(a)
        a0 = a(j);

        % Solve f_DF(a0*U, U, N) = target_v numerically for U
        U1 = solveForU(target_v1, a0, N);
        s1 = a0 * U1;

        U2 = solveForU(target_v2, a0, N);
        s2 = a0 * U2;

        Parameter.s1U1{i,j} = [s1; U1];
        Parameter.s2U2{i,j} = [s2; U2];
    end
end

% Check accuracy
for i = 1:numel(v1)
    for j = 1:numel(Parameter.a)
        s1U1 = Parameter.s1U1{i,j};
        s2U2 = Parameter.s2U2{i,j};

        Parameter.v1_check(i,j) = desaiFisher(s1U1(1), s1U1(2), N);
        Parameter.v2_check(i,j) = desaiFisher(s2U2(1), s2U2(2), N);
    end
end
fprintf('Parameter space created in %.2f seconds\n', toc);

%% WF Simulation
fprintf('Running Wright-Fisher simulations...\n');
tic

[ID_1, ID_2] = meshgrid(1:numel(a), 1:numel(a));

WFsim.generations = generations;
WFsim.t_burnin = t_burnin;

WFsim.F1_2 = cell(numel(v1), 1);
WFsim.F2_1 = cell(numel(v1), 1);
WFsim.v1_2 = cell(numel(v1), 1);
WFsim.v2_1 = cell(numel(v1), 1);

% Set broadcast variables
WFsim_F1_2_all = WFsim.F1_2;
WFsim_F2_1_all = WFsim.F2_1;
WFsim_v1_all = WFsim.v1_2;
WFsim_v2_all = WFsim.v2_1;

% Transform 2D cell arrays into 1D cell arrays
v1_param_all = Parameter.s1U1(:,1:end);
v2_param_all = Parameter.s2U2(:,1:end);
v1_param_all_1D = arrayfun(@(i) v1_param_all(i, :), 1:size(v1_param_all, 1), 'UniformOutput', false);
v2_param_all_1D = arrayfun(@(i) v2_param_all(i, :), 1:size(v2_param_all, 1), 'UniformOutput', false);

parfor i = 1:numel(v1)
    v1_param = v1_param_all_1D{i};
    v2_param = v2_param_all_1D{i};

    Rec_F1_2 = zeros(numel(ID_1), generations);
    Rec_F2_1 = zeros(numel(ID_2), generations);
    Rec_v1  = zeros(1,numel(ID_1));
    Rec_v2  = zeros(1,numel(ID_2));

    for j = 1:numel(ID_1)
        s1 = v1_param{ID_1(j)}(1);
        U1 = v1_param{ID_1(j)}(2);
        s2 = v2_param{ID_2(j)}(1);
        U2 = v2_param{ID_2(j)}(2);

        pop = [N 0 0];
        for t = 1:generations
            min_fitness = min( (pop(:,2)*s1) + (pop(:,3)*s2) );
            w = exp( ((pop(:,2)*s1) + (pop(:,3)*s2)) - min_fitness );
            fitness = pop(:,1).* w;
            PROB = fitness/sum(fitness);
            parents_counts = poissrnd(N * PROB);

            U = U1+U2;
            tot_mutations = binornd(parents_counts, U);
            num_mutations1 = binornd(tot_mutations, U1/U);
            parents_counts = parents_counts - num_mutations1;
            num_mutations2 = tot_mutations - num_mutations1;
            parents_counts = parents_counts - num_mutations2;

            mutation1_indices = (num_mutations1 > 0);
            mutation2_indices = (num_mutations2 > 0);

            pop_new = [parents_counts, pop(:, 2:3)];
            pop_new1 = zeros(0, 3);
            pop_new2 = zeros(0, 3);

            if any(mutation1_indices)
               pop_new1 = [num_mutations1(mutation1_indices), pop(mutation1_indices, 2) + 1, pop(mutation1_indices, 3)];
            end
            if any(mutation2_indices)
               pop_new2 = [num_mutations2(mutation2_indices), pop(mutation2_indices, 2), pop(mutation2_indices, 3) + 1];
            end
            pop = [pop_new; pop_new1; pop_new2];

            pop(pop(:, 1) == 0, :) = [];

            Rec_F1_2(j, t) = sum(pop(:, 1) .* pop(:, 2) * s1 ) / sum(pop(:,1));
            Rec_F2_1(j, t) = sum(pop(:, 1) .* pop(:, 3) * s2 ) / sum(pop(:,1));

            if t == generations
                Rec_v1(j) = (Rec_F1_2(j, generations) - Rec_F1_2(j, t_burnin)) / (generations-t_burnin+1);
                Rec_v2(j) = (Rec_F2_1(j, generations) - Rec_F2_1(j, t_burnin)) / (generations-t_burnin+1);
            end
        end
    end

    WFsim_F1_2_all{i} = Rec_F1_2;
    WFsim_F2_1_all{i} = Rec_F2_1;
    WFsim_v1_all{i} = Rec_v1';
    WFsim_v2_all{i} = Rec_v2';
end

WFsim.F1_2 = WFsim_F1_2_all;
WFsim.F2_1 = WFsim_F2_1_all;
WFsim.v1_2 = WFsim_v1_all;
WFsim.v2_1 = WFsim_v2_all;

fprintf('WF simulations completed in %.2f seconds\n', toc);

%% WF Simulation for v' estimate (single trait)
fprintf('Running single-trait validation simulations...\n');
tic

WFsim_vi.generations = generations;
WFsim_vi.t_burnin = t_burnin;

vectorStrings = Parameter.s1U1;
v_prime_param_all(1,1:n_a) = vectorStrings(1,1:n_a);
v_prime_param_all(2,1:n_a) = vectorStrings(7,1:n_a);
v_prime_param_all(3,1:n_a) = vectorStrings(12,1:n_a);
v_prime_param_all(4,1:n_a) = vectorStrings(16,1:n_a);
v_prime_param_all(5,1:n_a) = vectorStrings(19,1:n_a);
v_prime_param_all(6,1:n_a) = vectorStrings(21,1:n_a);
v_prime_param_all = v_prime_param_all';

WFsim_vi.F1 = cell(numel(v_prime_param_all), 1);
WFsim_vi.v1 = cell(numel(v_prime_param_all), 1);

WFsim_F1_all = WFsim_vi.F1;
WFsim_v1_all = WFsim_vi.v1;

parfor i = 1:numel(v_prime_param_all)
    Rec_F1 = zeros(1, generations);
    Rec_v1  = zeros(1);

    s1 = v_prime_param_all{i}(1);
    U1 = v_prime_param_all{i}(2);

    pop = [N 0 0];
    for t = 1:generations
        min_fitness = min( (pop(:,2)*s1) );
        w = exp( (pop(:,2)*s1) - min_fitness );
        fitness = pop(:,1).* w;
        PROB = fitness/sum(fitness);
        parents_counts = poissrnd(N * PROB);

        num_mutations1 = binornd(parents_counts, U1);
        parents_counts = parents_counts - num_mutations1;
        mutation1_indices = (num_mutations1 > 0);
        pop_new = [parents_counts, pop(:, 2:3)];
        pop_new1 = zeros(0, 3);

        if any(mutation1_indices)
            pop_new1 = [num_mutations1(mutation1_indices), pop(mutation1_indices, 2) + 1, pop(mutation1_indices, 3)];
        end
        pop = [pop_new; pop_new1];
        pop(pop(:, 1) == 0, :) = [];

        Rec_F1(1, t) = sum(pop(:, 1) .* pop(:, 2) * s1 ) /sum(pop(:,1));

        if t == generations
            Rec_v1 = (Rec_F1(1, generations) - Rec_F1(1, t_burnin)) / (generations-t_burnin+1);
        end
    end

    WFsim_F1_all{i} = Rec_F1;
    WFsim_v1_all{i} = Rec_v1';
end

WFsim_vi.F1 = WFsim_F1_all;
WFsim_vi.v1 = WFsim_v1_all;

fprintf('Single-trait validation completed in %.2f seconds\n', toc);

%% Analytical Predictions (closed-form quadratic — no symbolic toolbox)
fprintf('Computing analytical predictions...\n');
tic

tol = 10;

% Equal weights
v1_2 = zeros(numel(WFsim.F1_2), numel(ID_1), 1);
v2_1 = zeros(numel(WFsim.F2_1), numel(ID_2), 1);

for i = 1:numel(WFsim.F1_2)
    for j = 1:numel(ID_2)
        v2_val = Parameter.v2_check(i, ID_2(j));
        s1 = Parameter.s1U1{i, ID_1(j)}(1);

        v1_val = Parameter.v1_check(i, ID_1(j));
        s2 = Parameter.s2U2{i, ID_2(j)}(1);

        s_tilde = (s1 + s2) / 2;     % equal weights: w_i = 1/2

        U1_tilde = computeEffectiveMutationRate(s_tilde, v1_val, N, tol);
        U2_tilde = computeEffectiveMutationRate(s_tilde, v2_val, N, tol);

        if ~isnan(U1_tilde) && ~isnan(U2_tilde)
            U_tilde = U1_tilde + U2_tilde;
            v_total = desaiFisher(s_tilde, U_tilde, N);
            v1_2(i,j) = (U1_tilde / U_tilde) * v_total;
            v2_1(i,j) = (U2_tilde / U_tilde) * v_total;
        else
            v1_2(i,j) = NaN;
            v2_1(i,j) = NaN;
        end
    end
end
Prediction.v1_2 = v1_2;
Prediction.v2_1 = v2_1;

% U weights
v1_2 = zeros(numel(WFsim.F1_2), numel(ID_1), 1);
v2_1 = zeros(numel(WFsim.F2_1), numel(ID_2), 1);

for i = 1:numel(WFsim.F1_2)
    for j = 1:numel(ID_2)
        v1_val = Parameter.v1_check(i, ID_1(j));
        s1 = Parameter.s1U1{i, ID_1(j)}(1);
        U1 = Parameter.s1U1{i, ID_1(j)}(2);

        v2_val = Parameter.v2_check(i, ID_2(j));
        s2 = Parameter.s2U2{i, ID_2(j)}(1);
        U2 = Parameter.s2U2{i, ID_2(j)}(2);

        U_tot = U1 + U2;
        s_tilde = (U1/U_tot)*s1 + (U2/U_tot)*s2;    % U weights: w_i = U_i/U

        U1_tilde = computeEffectiveMutationRate(s_tilde, v1_val, N, tol);
        U2_tilde = computeEffectiveMutationRate(s_tilde, v2_val, N, tol);

        if ~isnan(U1_tilde) && ~isnan(U2_tilde)
            U_tilde = U1_tilde + U2_tilde;
            v_total = desaiFisher(s_tilde, U_tilde, N);
            v1_2(i,j) = (U1_tilde / U_tilde) * v_total;
            v2_1(i,j) = (U2_tilde / U_tilde) * v_total;
        else
            v1_2(i,j) = NaN;
            v2_1(i,j) = NaN;
        end
    end
end
Prediction_U_weights.v1_2 = v1_2;
Prediction_U_weights.v2_1 = v2_1;

% s weights
v1_2 = zeros(numel(WFsim.F1_2), numel(ID_1), 1);
v2_1 = zeros(numel(WFsim.F2_1), numel(ID_2), 1);

for i = 1:numel(WFsim.F1_2)
    for j = 1:numel(ID_2)
        v1_val = Parameter.v1_check(i, ID_1(j));
        s1 = Parameter.s1U1{i, ID_1(j)}(1);

        v2_val = Parameter.v2_check(i, ID_2(j));
        s2 = Parameter.s2U2{i, ID_2(j)}(1);

        s_tot = s1 + s2;
        s_tilde = (s1^2 + s2^2) / s_tot;             % s weights: w_i = s_i/(s1+s2)

        U1_tilde = computeEffectiveMutationRate(s_tilde, v1_val, N, tol);
        U2_tilde = computeEffectiveMutationRate(s_tilde, v2_val, N, tol);

        if ~isnan(U1_tilde) && ~isnan(U2_tilde)
            U_tilde = U1_tilde + U2_tilde;
            v_total = desaiFisher(s_tilde, U_tilde, N);
            v1_2(i,j) = (U1_tilde / U_tilde) * v_total;
            v2_1(i,j) = (U2_tilde / U_tilde) * v_total;
        else
            v1_2(i,j) = NaN;
            v2_1(i,j) = NaN;
        end
    end
end
Prediction_s_weights.v1_2 = v1_2;
Prediction_s_weights.v2_1 = v2_1;

fprintf('Analytical predictions completed in %.2f seconds\n', toc);

%% Save results
if ~isempty(savePath)
    fname = fullfile(savePath, sprintf('SteadyStateCM_N_%d.mat', N));
    save(fname, 'WFsim', 'WFsim_vi', 'Parameter', 'Prediction', 'Prediction_U_weights', 'Prediction_s_weights');
    fprintf('Results saved to %s\n', fname);
end

end

%% ====================== Helper functions ================================

function v = desaiFisher(s, Ub, N)
% Desai-Fisher rate of adaptation
    logRatio = log(s / Ub);
    v = s^2 * (2*log(N*s) - logRatio) / logRatio^2;
end

function Ut = computeEffectiveMutationRate(s_tilde, v_prime, N, tol)
% Solve f_DF(s_tilde, U_tilde, N) = v_prime for U_tilde via quadratic formula.
%
%   Quadratic in x = log(s_tilde / U_tilde):
%       v_prime * x^2 + s_tilde^2 * x - 2 * s_tilde^2 * log(N * s_tilde) = 0
%
%   Take the positive root so that s_tilde / U_tilde >> 1.
    A = v_prime;
    B = s_tilde^2;
    C = -2 * s_tilde^2 * log(N * s_tilde);

    discriminant = B^2 - 4*A*C;
    x_pos = (-B + sqrt(discriminant)) / (2*A);

    Ut = s_tilde * exp(-x_pos);

    if s_tilde / Ut < tol
        Ut = NaN;
    end
end

function U_sol = solveForU(target_v, a0, N)
% Numerically solve f_DF(a0*U, U, N) = target_v for U.
%
%   f_DF(a0*U, U, N) = (a0*U)^2 * (2*log(N*a0*U) - log(a0)) / log(a0)^2
%
%   Uses fzero on log(U) for numerical stability.
    log_a0 = log(a0);

    objFun = @(logU) a0^2 * exp(2*logU) * (2*log(N * a0 * exp(logU)) - log_a0) / log_a0^2 - target_v;

    % Bracket: U in [1e-15, 1]
    logU_sol = fzero(objFun, [-35, 0]);
    U_sol = exp(logU_sol);
end