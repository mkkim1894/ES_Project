function initialParameters = freezeParam(simParams)
% freezeParam - Precompute mutation rates and selection coefficients.
%
% Description:
%   Calculates the mutation supply rates (U1, U2) and selection coefficients
%   (s1, s2) for each initial phenotype, used for pre-run equilibration.
%
% Inputs:
%   simParams - Structure containing simulation parameters
%
% Outputs:
%   initialParameters - Matrix [nAngles x 4] where each row contains [U1, U2, s1, s2]
%
% Reference:
%   Kim, M., Ardell, S. M., & Kryazhimskiy, S. (2025).
%   "Module-Selection Balance in the Evolution of Modular Organisms."
%
% See also: preRunSimulation, initializeSimParams
%
% Copyright (c) 2025 Minkyu Kim, Cornell University
% Licensed under MIT License
    
    % Initialize the output matrix
    initialParameters = zeros(length(simParams.initialAngles), 4);
    
    % Loop through each initial phenotype
    for i_pos = 1:length(simParams.initialAngles)
        WT = simParams.initialPhenotypes(i_pos, :);
        
        % Compute mutation rates
        U1 = simParams.mutationRate * ( abs(WT(1)) * 1 / (simParams.deltaTrait * simParams.geneticTargetSize(1)) );
        U2 = simParams.mutationRate * ( abs(WT(2)) * 1 / (simParams.deltaTrait * simParams.geneticTargetSize(2)) );
        
        % Compute fitness differences for mutational effects
        traitDelta1 = WT(1) + simParams.deltaTrait;
        traitDelta2 = WT(2) + simParams.deltaTrait;
        
        fitnessDelta0 = -( ( WT(1) / simParams.ellipseParams(1) )^2 + ( WT(2) / simParams.ellipseParams(2) )^2 ); 
        fitnessDelta1 = -( ( traitDelta1 / simParams.ellipseParams(1) )^2 + ( WT(2) / simParams.ellipseParams(2) )^2 );
        fitnessDelta2 = -( ( WT(1) / simParams.ellipseParams(1) )^2 + ( traitDelta2 / simParams.ellipseParams(2) )^2 );
        
        % Compute selection coefficients
        expFitnessDelta0 = exp( fitnessDelta0 / (2 * simParams.landscapeStdDev^2) ); 
        expFitnessDelta1 = exp( fitnessDelta1 / (2 * simParams.landscapeStdDev^2) ); 
        expFitnessDelta2 = exp( fitnessDelta2 / (2 * simParams.landscapeStdDev^2) ); 

        s1 = log(expFitnessDelta1) - log(expFitnessDelta0);
        s2 = log(expFitnessDelta2) - log(expFitnessDelta0);
        
        % Store parameters in the output matrix
        initialParameters(i_pos, :) = [U1, U2, s1, s2];
    end
end
