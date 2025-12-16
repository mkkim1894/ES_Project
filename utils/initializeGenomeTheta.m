function [genomeParams] = initializeGenomeTheta(L, simParams, seed)
% initializeGenomeTheta - Initialize genome parameters for pleiotropic FGM.
%
% Description:
%   Generates pleiotropic angle assignments and initial genome states for
%   the pleiotropic Fisher's Geometric Model. Each locus is assigned a random
%   angle that determines how mutations at that locus affect phenotype.
%
% Inputs:
%   L         - Number of loci in the genome
%   simParams - Structure containing simulation parameters including:
%       .initialPhenotypes - Matrix of initial phenotype coordinates
%       .deltaTrait        - Mutational step size
%   seed      - Random seed for reproducibility
%
% Outputs:
%   genomeParams - Structure containing:
%       .genomeTheta   - Vector [1 x L] of pleiotropic angles
%       .initialGenomes - Matrix [nAngles x L] of initial allele states
%
% Example:
%   genomeParams = initializeGenomeTheta(1000, simParams, 42);
%
% Reference:
%   Kim, M., Ardell, S. M., & Kryazhimskiy, S. (2025).
%   "Module-Selection Balance in the Evolution of Modular Organisms."
%
% See also: simulatePleiotropicSSWM, simulatePleiotropicCM
%
% Copyright (c) 2025 Minkyu Kim, Cornell University
% Licensed under MIT License

    % Set the seed for reproducibility
    if nargin > 2
        rng(seed);
    end

    % Extract parameters from simParams
    targetPhenotypes = simParams.initialPhenotypes;
    delta = simParams.deltaTrait;

    % Number of target phenotypes (assumed to be 6)
    numPhenotypes = size(targetPhenotypes, 1);

    % Initialize outputs
    genomeTheta = 2 * pi * rand(1, L); % Sample theta values uniformly from [0, 2*pi]
    initialGenomes = zeros(numPhenotypes, L); % Start with all loci set to 0
    currentPhenotypes = zeros(numPhenotypes, 2); % Start at [0, 0]

    % Iterate through each target phenotype
    for i = 1:numPhenotypes
        for ell = 1:L
            % Calculate the effect of flipping locus ell
            mutationalEffect = [-delta * cos(genomeTheta(ell)), -delta * sin(genomeTheta(ell))];
            newPhenotype = currentPhenotypes(i, :) + mutationalEffect;

            % Calculate distances to the target phenotype
            distanceBefore = norm(currentPhenotypes(i, :) - targetPhenotypes(i, :));
            distanceAfter = norm(newPhenotype - targetPhenotypes(i, :));

            % Accept the flip if it reduces the distance to the target
            if distanceAfter < distanceBefore
                initialGenomes(i, ell) = 1; % Flip the bit
                currentPhenotypes(i, :) = newPhenotype; % Update the current phenotype
            end
        end
    end

    % Store outputs in genomeParams structure
    genomeParams.genomeTheta = genomeTheta;
    genomeParams.initialGenomes = initialGenomes;
    genomeParams.currentPhenotypes = currentPhenotypes;
end
