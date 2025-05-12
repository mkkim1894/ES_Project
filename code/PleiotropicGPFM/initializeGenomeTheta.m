function [genomeParams] = initializeGenomeTheta(L, simParams, seed)
% initializeGenomeTheta - Finds a genotype by iterating through loci and
% flipping bits if they move the current phenotype closer to the target phenotype.
% Generates a theta vector for each site in the genome.
%
% Inputs:
%   L - Number of loci (e.g., 1000 or 10000).
%   simParams - Structure containing the simulation parameters, including:
%       initialPhenotypes - 6x2 matrix of initial phenotype coordinates.
%       deltaTrait - Scalar effect size for each locus.
%   seed - Integer value to set the seed for the random number generator for reproducibility.
%
% Outputs:
%   genomeParams - Structure containing the generated genomeTheta, initialGenomes, and currentPhenotypes.
%
% Example Usage:
%   genomeParams = initializeGenomeTheta(1000, simParams, 42);

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
