function initialParameters = freezeParam(simParams)
% freezeParam - Precomputes mutation rates (U1, U2) and selection coefficients (s1, s2)
% for each initial phenotype based on simulation parameters.
%
% Inputs:
%   simParams - Structure containing simulation parameters.
%
% Outputs:
%   initialParameters - A matrix where each row corresponds to an initial phenotype
%                       and contains [U1, U2, s1, s2].
    
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
