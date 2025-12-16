% Function: findInitialPhenotype
% Description:
%   This function calculates the initial phenotypic coordinates in a 2D 
%   Fisher's Geometric Model (FGM) given the initial angles, ellipse 
%   parameters representing selection pressures, initial fitness, and the 
%   standard deviation of the fitness landscape.
%
% Inputs:
%   initialAngles - A vector containing initial angles (in radians) for the 
%                   phenotypic space, defining the direction of each starting point.
%   ellipseParams - A 1x2 vector [a1, a2] representing the parameters of 
%                   the elliptical fitness landscape (axes lengths).
%   initialFitness - A scalar representing the initial fitness value that the 
%                    phenotypic coordinates should correspond to.
%   landscapeStdDev - A scalar representing the standard deviation of the 
%                     fitness landscape.
%
% Outputs:
%   initialPhenotypes - An Nx2 matrix where N is the number of initial angles. 
%                      Each row contains the (x, y) coordinates of the initial 
%                      phenotype for the corresponding angle.
%
% Example Usage:
%   initialPhenotypes = findInitialPhenotypes([atan2(1, 6.4), atan2(1, 3.2)], [sqrt(2), 1], 0.25, 2);
%
% Reference:
%   Kim, M., Ardell, S. M., & Kryazhimskiy, S. (2025).
%   "Module-Selection Balance in the Evolution of Modular Organisms."
%
% Copyright (c) 2025 Minkyu Kim, Cornell University
% Licensed under MIT License

function initialPhenotypes = findInitialPhenotypes(initialAngles, ellipseParams, initialFitness, landscapeStdDev)    
    % Preallocate matrix for initial phenotype coordinates
    initialPhenotypes = zeros(length(initialAngles), 2);
    
    % Loop through each angle to calculate the corresponding initial phenotype
    for i = 1:length(initialAngles)
        theta = initialAngles(i);  % Current angle in radians
        
        % Symbolic variable for solving radius
        syms r
        % Define the fitness equation to solve for r
        eqn = initialFitness == exp(-sqrt((r * cos(theta) / ellipseParams(1))^2 + (r * sin(theta) / ellipseParams(2))^2)^2 / (2 * landscapeStdDev^2));
        
        % Solve for the positive radius value
        sol = solve(eqn, r);
        Radius = double(sol);
        Radius = Radius(Radius > 0);  % Only keep positive solutions
        
        % Store the calculated coordinates in the output matrix
        initialPhenotypes(i, 1:end) = Radius * [-cos(theta); -sin(theta)];
    end
end
