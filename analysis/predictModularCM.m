function [analyticalTrajectories] = predictModularCM( simParams )
% predictModularCM - Compute analytical trajectory predictions for modular CM regime.
%
% Description:
%   Solves ODEs for the expected evolutionary trajectory under concurrent
%   mutations, accounting for clonal interference between modules.
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
% See also: predictModularSSWM, simulateModularCM
%
% Copyright (c) 2025 Minkyu Kim, Cornell University
% Licensed under MIT License

analyticalTrajectories = cell(1);
initialAngles = simParams.initialAngles;    
initialPhenotypes = simParams.initialPhenotypes;

for i_pos = 1:length(initialAngles)
    WT = initialPhenotypes(i_pos,1:end);
    
    % Solve the ODEs using ode45
    tspan = [0 10000];
    tol = 10;
    finalFitnessThreshold = 0.99;
    options = odeset('Events', @(t,x) eventFunction(t, x, simParams, finalFitnessThreshold, tol));
    [~, X] = ode45(@(t, x) computeAdaptationRate(t, x, simParams, tol), tspan, WT, options);

    analyticalTrajectories{i_pos, 1} = X;
end
end

%--------------------------------------------------------------------------

function dxdt = computeAdaptationRate(~, x, simParams, tol)    
    WT = x;
 
    U1 = simParams.mutationRate * abs(WT(1)) * ( 1/(simParams.deltaTrait * simParams.geneticTargetSize(1)) );
    U2 = simParams.mutationRate * abs(WT(2)) * ( 1/(simParams.deltaTrait * simParams.geneticTargetSize(2)) );
    
    traitDelta1 = WT(1) + simParams.deltaTrait;
    traitDelta2 = WT(2) + simParams.deltaTrait;
            
    fitnessDelta0 = -( ( WT(1)/simParams.ellipseParams(1) )^2 + ( WT(2)/simParams.ellipseParams(2) )^2 ); 
    fitnessDelta1 = -( ( traitDelta1/simParams.ellipseParams(1) )^2 + ( WT(2)/simParams.ellipseParams(2) )^2 );
    fitnessDelta2 = -( ( WT(1)/simParams.ellipseParams(1) )^2 + ( traitDelta2/simParams.ellipseParams(2) )^2 );
    
    expFitnessDelta0 = exp( fitnessDelta0./(2*simParams.landscapeStdDev^2) ); 
    expFitnessDelta1 = exp( fitnessDelta1./(2*simParams.landscapeStdDev^2) ); 
    expFitnessDelta2 = exp( fitnessDelta2./(2*simParams.landscapeStdDev^2) ); 
    
    s1 = log(expFitnessDelta1) - log(expFitnessDelta0);
    s2 = log(expFitnessDelta2) - log(expFitnessDelta0);

    syms s0 U0 N0
    AdaptRate = ( (s0^2) * ( (2 * log(N0 * s0)) - log(s0 / U0) ) ) / ( log(s0 / U0) )^2;

    v1 = double(subs(AdaptRate, {s0, U0, N0}, {s1, U1, simParams.popSize}));
    v2 = double(subs(AdaptRate, {s0, U0, N0}, {s2, U2, simParams.popSize}));

    if v2 / v1 > 100
        rij = [0; (v2 / s2) * simParams.deltaTrait];
    elseif v1 / v2 > 100
        rij = [(v1 / s1) * simParams.deltaTrait; 0];
    else
        s_new = (s1 / (s1 + s2)) * s1 + (s2 / (s1 + s2)) * s2;
        
        U2_new = double(solve(subs(AdaptRate, {s0, N0}, {s_new, simParams.popSize}) == v2, U0));
        U2_new = U2_new( (s_new ./ U2_new) > tol);
        
        U1_new = double(solve(subs(AdaptRate, {s0, N0}, {s_new, simParams.popSize}) == v1, U0));       
        U1_new = U1_new( (s_new ./ U1_new) > tol);   

        if numel(U1_new) == 1 && numel(U2_new) == 1         
            temp_v12 = double((U1_new / (U1_new + U2_new)) * subs(AdaptRate, {s0, U0, N0}, {s_new, U1_new + U2_new, simParams.popSize}));
            temp_v21 = double((U2_new / (U1_new + U2_new)) * subs(AdaptRate, {s0, U0, N0}, {s_new, U1_new + U2_new, simParams.popSize}));
        else    
            temp_v12 = NaN;
            temp_v21 = NaN;  
        end     
        v_12 = temp_v12;    
        v_21 = temp_v21;

        rij = [(v_12 / s1) * simParams.deltaTrait; (v_21 / s2) * simParams.deltaTrait];
    end
    
    % The derivatives are represented by the variable rij here
    dxdt = rij;
end

%--------------------------------------------------------------------------

function [value, isterminal, direction] = eventFunction(t, x, simParams, finalFitnessThreshold, tol)
    persistent last_t last_x last_rij
    if isempty(last_t)
        last_t = NaN;
        last_x = NaN(size(x));
        last_rij = NaN(size(x));
    end
    
    % Check if we need to recompute rij
    if t ~= last_t || any(x ~= last_x)
        last_rij = computeAdaptationRate(t, x, simParams, tol);
        last_t = t;
        last_x = x;
    end
    
    Fitness = -( (x(1)/simParams.ellipseParams(1))^2 + (x(2)/simParams.ellipseParams(2))^2 );
    expFitness = exp( Fitness./(2*simParams.landscapeStdDev^2) ); 
    endCondition = expFitness - finalFitnessThreshold;
    
    % Additional condition to check if any value of rij is NaN
    isNaNCondition = any(isnan(last_rij));
    
    % Combine the conditions
    value = [endCondition; isNaNCondition]; % The conditions to be checked
    isterminal = [1; 1];  % Stop the integration if either condition is met
    direction = [0; 0];   % Zero-crossing in any direction will trigger the event
end
