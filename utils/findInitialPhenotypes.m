function initialPhenotypes = findInitialPhenotypes(initialAngles, ellipseParams, initialFitness, landscapeStdDev) 
% findInitialPhenotypes - Find initial phenotype coordinates on a fitness contour.
%
% Description:
%   For each initial angle, computes the phenotype coordinates (x1, x2) that
%   lie on the fitness contour W = initialFitness in the direction of that angle.
%   Uses the symbolic solver to find the radius along each direction.
%
% Inputs:
%   initialAngles   - Vector of angles (radians) defining phenotypic directions
%   ellipseParams   - [a1, a2] ellipse axes of the fitness landscape
%   initialFitness  - Target fitness value W0 for the initial contour
%   landscapeStdDev - Fitness landscape width parameter (sigma)
%
% Outputs:
%   initialPhenotypes - Matrix [nAngles x 2], each row is (x1, x2) on the W0 contour
%
% Example:
%   initialPhenotypes = findInitialPhenotypes([atan2(1, 6.4), atan2(1, 3.2)], [1, 1/sqrt(2)], 0.25, 2);
%
% Reference:
%   Kim, M., Ardell, S. M., & Kryazhimskiy, S. (2025).
%   "Module-Selection Balance in the Evolution of Modular Organisms."
%
% See also: initializeSimParams
%
% Copyright (c) 2025 Minkyu Kim, Cornell University
% Licensed under MIT License

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
