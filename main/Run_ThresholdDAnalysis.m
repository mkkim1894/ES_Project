function Run_ThresholdDAnalysis(varargin)
% Run_ThresholdDAnalysis - Driver script for threshold D sensitivity analysis.
%
% Description:
%   Computes analytical evolutionary trajectories for different rate-ratio
%   threshold D values. The threshold D determines when adaptation is dominated
%   by a single module (when v2/v1 > D or v1/v2 > D).
%
% Inputs:
%   varargin - Optional name-value pairs:
%       'D_values'  - Vector of threshold values (default: [10000, 1000, 100, 10])
%       'outputDir' - Directory for output (default: current)
%       'outputFile'- Output filename (default: 'ThresholdD_Analysis.mat')
%
% Outputs:
%   Saves .mat file containing:
%       analyticalTrajectories - Cell array of trajectories for each D value
%       D_values              - Threshold values used
%       summary               - Parameter structure
%
% Example:
%   Run_ThresholdDAnalysis();
%   Run_ThresholdDAnalysis('D_values', [1000, 100, 10], 'outputDir', 'results/');
%
% Reference:
%   Kim, M., Ardell, S. M., & Kryazhimskiy, S. (2025).
%   "Module-Selection Balance in the Evolution of Modular Organisms."
%
% See also: predictCM_ThresholdD, makeFigureS_ThresholdDTrajectories
%
% Copyright (c) 2025 Minkyu Kim, Cornell University
% Licensed under MIT License

%% Parse inputs
p = inputParser;
addParameter(p, 'D_values', [10000, 1000, 100, 10], @isnumeric);
addParameter(p, 'outputDir', '.', @ischar);
addParameter(p, 'outputFile', 'ThresholdD_Analysis.mat', @ischar);
parse(p, varargin{:});

D_values = p.Results.D_values;
outputDir = p.Results.outputDir;
outputFile = p.Results.outputFile;

if ~isfolder(outputDir)
    mkdir(outputDir);
end

%% Parameter Setup
fprintf('Setting up parameters...\n');

summary.step = [0.1, 0.1];                  % Step size for evolutionary process
summary.sigW = 2;                           % Standard deviation of fitness landscape
summary.SelectionBias = [1, 2];             % Selection bias in each direction
summary.MutationBias = [1, 1];              % Mutation bias in each direction
summary.InitialAngle = [atan2(1, 6.4), ...  % Initial angles
                        atan2(1, 3.2), ...
                        atan2(1, 1.6), ...
                        atan2(1, 0.8), ...
                        atan2(1, 0.4), ...
                        atan2(1, 0.2)];
summary.d_ref = -2 * [1; 1];                % Reference distance
summary.W_ref = exp(-vecnorm(summary.d_ref)^2 / (2 * summary.sigW^2));

% Evolutionary regime parameters
summary.N = 10^4;                           % Population size
summary.U = 10^-3;                          % Mutation rate

%% Main loop over D values
fprintf('Computing analytical trajectories for %d threshold values...\n', length(D_values));
tic;

analyticalTrajectories = cell(length(summary.InitialAngle), length(D_values));

for i_d = 1:length(D_values)
    summary.D = D_values(i_d);
    fprintf('  Processing D = %d...\n', summary.D);
    
    Theory = predictCM_ThresholdD(summary);
    
    % Store trajectories for this D value
    for i_pos = 1:length(summary.InitialAngle)
        analyticalTrajectories{i_pos, i_d} = Theory{i_pos, 1};
    end
end

elapsed = toc;
fprintf('Completed in %.2f seconds\n', elapsed);

%% Save results
outputPath = fullfile(outputDir, outputFile);
save(outputPath, 'analyticalTrajectories', 'D_values', 'summary');
fprintf('Results saved to %s\n', outputPath);

%% Display summary
fprintf('\n=== Analysis Summary ===\n');
fprintf('Population size N: %d\n', summary.N);
fprintf('Mutation rate U: %.1e\n', summary.U);
fprintf('D values tested: ');
fprintf('%d ', D_values);
fprintf('\n');
fprintf('Initial angles: %d conditions\n', length(summary.InitialAngle));
fprintf('Output file: %s\n', outputPath);

end
