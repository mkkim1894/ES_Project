function [Theory] = Prediction_CM_D_check( summary, Mean_TimeStamp, timestamp )
%% Set structure "Final" that stores all final meta-information

Theory.Record = cell(1);

%% << Setup Parameters and variables >>

InitialAngle = summary.InitialAngle;

for i_pos = 1:length(InitialAngle)
sigW = summary.sigW;

if nargin > 2
    % If the optional argument is provided, calculate analytical
    % trajectories from the nth timestamp location of simulations
    WT0 = Mean_TimeStamp(:,timestamp); 
else
    theta = InitialAngle(i_pos); 
    d = Find_InitialPoint( theta, summary.SelectionBias, summary.W_ref, sigW );
    d = d(d>0);

    WT0 = d.*[-cos(theta); -sin(theta)];
end

% Solve the ODEs using ode45
tspan = [0 10000];
tol = 10;
endW = 0.99;
options = odeset('Events', @(t,x) event_function_CM(t, x, summary, sigW, endW, tol));
[t, X] = ode45(@(t, x) rij_function_CM(t, x, summary, tol), tspan, WT0, options);

Theory.Record{i_pos, 1} = X;
Theory.Record{i_pos, 2} = t;
end
        
end

%--------------------------------------------------------------------------

function dxdt = rij_function_CM(~, x, summary, tol)
    N = summary.N;
    U = summary.U;
    
    k1 = summary.step(1);
    k2 = summary.step(2);
    sigW = summary.sigW;
    rho = abs(summary.d_ref(1)) + abs(summary.d_ref(2));
    rho1 = summary.MutationBias(1)/rho;
    rho2 = summary.MutationBias(2)/rho;
    
    WT = x;

    U1 = (abs(WT(1))*rho1) * U;
    U2 = (abs(WT(2))*rho2) * U;

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
    AdaptRate = ( (s0^2) * ( (2 * log(N0 * s0)) - log(s0 / U0) ) ) / ( log(s0 / U0) )^2;

    v1 = double(subs(AdaptRate, {s0, U0, N0}, {s1, U1, N}));
    v2 = double(subs(AdaptRate, {s0, U0, N0}, {s2, U2, N}));

    if v2 / v1 > summary.D
        rij = [0; (v2 / s2) * k2];
    elseif v1 / v2 > summary.D
        rij = [(v1 / s1) * k1; 0];
    else
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
    
    % The derivatives are represented by the variable rij here
    dxdt = rij;
    
end

%--------------------------------------------------------------------------

function [value, isterminal, direction] = event_function_CM(t, x, summary, sigW, endW, tol)
    persistent last_t last_x last_rij
    if isempty(last_t)
        last_t = NaN;
        last_x = NaN(size(x));
        last_rij = NaN(size(x));
    end
    
    % Check if we need to recompute rij
    if t ~= last_t || any(x ~= last_x)
        last_rij = rij_function_CM(t, x, summary, tol);
        last_t = t;
        last_x = x;
    end
    
    d0 = sqrt(summary.SelectionBias(1) * x(1).^2 + summary.SelectionBias(2) * x(2).^2);
    W0 = exp(-d0.^2 / 2 / sigW^2);
    endWCondition = W0 - endW;
    
    % Additional condition to check if any value of rij is NaN
    isNaNCondition = any(isnan(last_rij));
    
    % Combine the conditions
    value = [endWCondition; isNaNCondition]; % The conditions to be checked
    isterminal = [1; 1];  % Stop the integration if either condition is met
    direction = [0; 0];   % Zero-crossing in any direction will trigger the event
end

%--------------------------------------------------------------------------

function d = Find_InitialPoint( theta, ShapeParameter, Default_fitness, sigW )

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