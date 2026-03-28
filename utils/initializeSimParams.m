function simParams = initializeSimParams(varargin)
% initializeSimParams - Initializes simulation parameters for a 2D Fisher's Geometric Model (FGM).
% 
% Description:
%   Sets default simulation parameters and allows customization through
%   optional name-value pairs. Returns a structure containing all parameters
%   and generates initial phenotypes from the provided angles.
% 
% Inputs:
%   varargin - Optional name-value pairs to override default parameter values or opt-out parameters.
%              Example: 'numIteration', 500, 'popSize', 5e4, 'ellipseRatio', sqrt(2)
%              To opt-out, use: 'omitParams', {'ellipseRatio', 'mutationRate'}
% 
% Outputs:
%   simParams - A structure containing all simulation parameters, including 
%               initial phenotypes.
% 
% Example Usage:
%   simParams = initializeSimParams('numIteration', 500, 'popSize', 5e4, ...
%                                   'discordantAngles', [pi/16, 7*pi/16]);
%
% Reference:
%   Kim, M., Ardell, S. M., & Kryazhimskiy, S. (2025).
%   "Module-Selection Balance in the Evolution of Modular Organisms."
%
% See also: findInitialPhenotypes, simulateModularSSWM
%
% Copyright (c) 2025 Minkyu Kim, Cornell University
% Licensed under MIT License

    % Initialize input parser
    p = inputParser;

    % Define default parameter values
    addParameter(p, 'numIteration', 1000); % Number of iterations
    addParameter(p, 'landscapeStdDev', 2); % Standard deviation of the fitness landscape
    addParameter(p, 'deltaTrait', 0.1); % Step size for mutational changes in traits
    addParameter(p, 'ellipseRatio', sqrt(2)); % Desired aspect ratio for ellipse parameters
    addParameter(p, 'geneticTargetSize', [10, 10]); % Number of loci per module [L1, L2]
    addParameter(p, 'initialFitness', 0.25); % Initial fitness value
    addParameter(p, 'initialAngles', [atan2(1, 6.4), atan2(1, 3.2), ...
        atan2(1, 1.6), atan2(1, 0.8), atan2(1, 0.4), atan2(1, 0.2)]); % Initial angles in phenotype space
    addParameter(p, 'popSize', 10^4); % Population size
    addParameter(p, 'mutationRate', 10^-7); % Mutation rate
    addParameter(p, 'recombinationRate', 0); % Recombination rate
    addParameter(p, 'discordantAngles', []); % [theta1, theta2] for discordant-module GPM
    addParameter(p, 'omitParams', {}); % List of parameters to exclude

    % Parse input arguments
    parse(p, varargin{:});
    simParams = p.Results;

    % Remove the omitted parameters
    omitParams = simParams.omitParams;
    simParams = rmfield(simParams, intersect(fieldnames(simParams), omitParams));
    
    % Generate normalized ellipse parameters based on the provided ellipseRatio if it's not omitted
    if ~isfield(simParams, 'ellipseRatio')
        warning('Ellipse ratio is not provided or omitted. Default [1, 1] ellipse parameters are used.');
        simParams.ellipseParams = [1, 1];
    else
        if simParams.ellipseRatio >= 1
            simParams.ellipseParams = [1, 1/simParams.ellipseRatio];
        else
            simParams.ellipseParams = [simParams.ellipseRatio, 1];
        end
    end

    % Generate initial phenotypes based on input or default parameters
    if isfield(simParams, 'initialAngles')
        simParams.initialPhenotypes = findInitialPhenotypes(simParams.initialAngles, ...
            simParams.ellipseParams, simParams.initialFitness, simParams.landscapeStdDev);
    else
        warning('Initial angles omitted. Initial phenotypes are not generated.');
    end

    % Remove 'omitParams' itself — it is an internal bookkeeping field, not a simulation parameter
    simParams = rmfield(simParams, 'omitParams');
end