function [analyticalTrajectories] = predictCM_ThresholdD(summary, varargin)
% predictCM_ThresholdD - Compute analytical trajectories with rate-ratio threshold D.
%
% Description:
%   Solves ODEs for expected evolutionary trajectory under concurrent mutations,
%   using a threshold parameter D to determine when one module dominates adaptation.
%   When v2/v1 > D, only module 2 adapts; when v1/v2 > D, only module 1 adapts.
%
% Inputs:
%   summary - Structure containing:
%       .InitialAngle    - Vector of initial angles
%       .SelectionBias   - Selection pressure [a1, a2]
%       .MutationBias    - Mutation bias [rho1, rho2]
%       .step            - Mutational step sizes [k1, k2]
%       .sigW            - Fitness landscape width
%       .d_ref           - Reference distance for fitness
%       .W_ref           - Reference fitness
%       .N               - Population size
%       .U               - Mutation rate
%       .D               - Rate ratio threshold (e.g., 100)
%   varargin - Optional:
%       'timestamp'      - Start from nth timestamp of existing trajectory
%       'Mean_TimeStamp' - Existing trajectory to start from
%
% Outputs:
%   analyticalTrajectories - Cell array {nAngles x 2}:
%       Column 1: [x1, x2] trajectory matrix
%       Column 2: Time points
%
% Example:
%   summary.D = 100;
%   trajectories = predictCM_ThresholdD(summary);
%
% Reference:
%   Kim, M., Ardell, S. M., & Kryazhimskiy, S. (2025).
%   "Module-Selection Balance in the Evolution of Modular Organisms."
%
% See also: simulateModularCM, predictModularCM
%
% Copyright (c) 2025 Minkyu Kim, Cornell University
% Licensed under MIT License

%% Parse optional arguments
timestamp = [];
Mean_TimeStamp = [];
if nargin > 1
    for i = 1:2:length(varargin)
        switch lower(varargin{i})
            case 'timestamp'
                timestamp = varargin{i+1};
            case 'mean_timestamp'
                Mean_TimeStamp = varargin{i+1};
        end
    end
end

%% Initialize
analyticalTrajectories = cell(length(summary.InitialAngle), 2);
InitialAngle = summary.InitialAngle;
sigW = summary.sigW;
tol = 10;
endW = 0.99;

%% Solve for each initial condition
for i_pos = 1:length(InitialAngle)
    
    % Determine initial point
    if ~isempty(timestamp) && ~isempty(Mean_TimeStamp)
        WT0 = Mean_TimeStamp(:, timestamp);
    else
        theta = InitialAngle(i_pos); 
        d = findInitialPoint(theta, summary.SelectionBias, summary.W_ref, sigW);
        d = d(d > 0);
        WT0 = d .* [-cos(theta); -sin(theta)];
    end

    % Solve ODEs
    tspan = [0 10000];
    options = odeset('Events', @(t,x) eventFunction(t, x, summary, sigW, endW, tol));
    [t, X] = ode45(@(t, x) computeAdaptationRate(t, x, summary, tol), tspan, WT0, options);

    analyticalTrajectories{i_pos, 1} = X;
    analyticalTrajectories{i_pos, 2} = t;
end

end

%% ODE function
function dxdt = computeAdaptationRate(~, x, summary, tol)

N = summary.N;
U = summary.U;
D = summary.D;

k1 = summary.step(1);
k2 = summary.step(2);
sigW = summary.sigW;
rho = abs(summary.d_ref(1)) + abs(summary.d_ref(2));
rho1 = summary.MutationBias(1) / rho;
rho2 = summary.MutationBias(2) / rho;

WT = x;

U1 = (abs(WT(1)) * rho1) * U;
U2 = (abs(WT(2)) * rho2) * U;

X = WT(1) + k1;
Y = WT(2) + k2;

d0 = sqrt(summary.SelectionBias(1) * WT(1).^2 + summary.SelectionBias(2) * WT(2).^2);
d1 = sqrt(summary.SelectionBias(1) * X.^2 + summary.SelectionBias(2) * WT(2).^2);
d2 = sqrt(summary.SelectionBias(1) * WT(1).^2 + summary.SelectionBias(2) * Y.^2);

W0 = exp(-d0.^2 / 2 / sigW^2);
W1 = exp(-d1.^2 / 2 / sigW^2);
W2 = exp(-d2.^2 / 2 / sigW^2);

s1 = log(W1) - log(W0);
s2 = log(W2) - log(W0);

syms s0 U0 N0
AdaptRate = ((s0^2) * ((2 * log(N0 * s0)) - log(s0 / U0))) / (log(s0 / U0))^2;

v1 = double(subs(AdaptRate, {s0, U0, N0}, {s1, U1, N}));
v2 = double(subs(AdaptRate, {s0, U0, N0}, {s2, U2, N}));

%% Apply threshold D logic
if v2 / v1 > D
    % Module 2 dominates
    rij = [0; (v2 / s2) * k2];
elseif v1 / v2 > D
    % Module 1 dominates
    rij = [(v1 / s1) * k1; 0];
else
    % Both modules contribute
    s_tilde = (s1 / (s1 + s2)) * s1 + (s2 / (s1 + s2)) * s2;
    
    U2_tilde = double(solve(subs(AdaptRate, {s0, N0}, {s_tilde, N}) == v2, U0));
    U2_tilde = U2_tilde(s_tilde ./ U2_tilde > tol);
    
    U1_tilde = double(solve(subs(AdaptRate, {s0, N0}, {s_tilde, N}) == v1, U0));       
    U1_tilde = U1_tilde(s_tilde ./ U1_tilde > tol);   

    if numel(U1_tilde) == 1 && numel(U2_tilde) == 1         
        temp_v1_2 = double((U1_tilde / (U1_tilde + U2_tilde)) * subs(AdaptRate, {s0, U0, N0}, {s_tilde, U1_tilde + U2_tilde, N}));
        temp_v2_1 = double((U2_tilde / (U1_tilde + U2_tilde)) * subs(AdaptRate, {s0, U0, N0}, {s_tilde, U1_tilde + U2_tilde, N}));
    else    
        temp_v1_2 = NaN;
        temp_v2_1 = NaN;  
    end     
    v_12 = temp_v1_2;    
    v_21 = temp_v2_1;

    rij = [(v_12 / s1) * k1; (v_21 / s2) * k2];
end

dxdt = rij;
end

%% Event function to stop integration
function [value, isterminal, direction] = eventFunction(t, x, summary, sigW, endW, tol)
persistent last_t last_x last_rij
if isempty(last_t)
    last_t = NaN;
    last_x = NaN(size(x));
    last_rij = NaN(size(x));
end

% Recompute if needed
if t ~= last_t || any(x ~= last_x)
    last_rij = computeAdaptationRate(t, x, summary, tol);
    last_t = t;
    last_x = x;
end

d0 = sqrt(summary.SelectionBias(1) * x(1).^2 + summary.SelectionBias(2) * x(2).^2);
W0 = exp(-d0.^2 / 2 / sigW^2);
endWCondition = W0 - endW;

isNaNCondition = any(isnan(last_rij));

value = [endWCondition; isNaNCondition];
isterminal = [1; 1];
direction = [0; 0];
end

%% Helper function
function d = findInitialPoint(theta, ShapeParameter, Default_fitness, sigW)
if nargin < 3
    Default_fitness = 0.000123409804086679;
    sigW = 1;
end

a1 = ShapeParameter(1);
a2 = ShapeParameter(2); 

syms r
eqn = Default_fitness == exp(-sqrt(a1*(r*cos(theta))^2 + a2*(r*sin(theta))^2)^2/2/sigW^2); 
sol = solve(eqn, r);
d = double(sol);
end
