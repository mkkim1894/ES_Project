function [genomeParams] = initializeGenomeTheta(L, simParams, seed)
% initializeGenomeTheta - Initialize genome for pleiotropic FGM simulations.
%
% Description:
%   Assigns random pleiotropic angles to L loci, then constructs an initial
%   genome for each target phenotype via a greedy pass: locus ell is set to
%   allele 1 if doing so reduces the Euclidean distance to the target phenotype.
%   Starting from the all-zeros genome (phenotype [0,0]), this produces a genome
%   whose encoded phenotype approximates the target.
%
% Inputs:
%   L         - Number of loci in the genome
%   simParams - Structure containing simulation parameters:
%       .initialPhenotypes - Matrix of target phenotype coordinates [nAngles x 2]
%       .deltaTrait        - Mutational step size (delta)
%   seed      - Random seed for reproducibility
%
% Outputs:
%   genomeParams - Structure containing:
%       .genomeTheta       - Vector [1 x L] of pleiotropic angles, uniform on [0, 2*pi)
%       .initialGenomes    - Matrix [nAngles x L] of initial allele states
%       .currentPhenotypes - Matrix [nAngles x 2] of phenotypes encoded by initialGenomes
%
% Example:
%   genomeParams = initializeGenomeTheta(400, simParams, 1);
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

    % Number of target phenotypes
    numPhenotypes = size(targetPhenotypes, 1);

    % Initialize outputs
    genomeTheta = 2 * pi * rand(1, L); % Sample theta values uniformly from [0, 2*pi]
    initialGenomes = zeros(numPhenotypes, L); % Start with all loci set to 0
    currentPhenotypes = zeros(numPhenotypes, 2); % Allele-0 genome produces zero phenotypic displacement

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
