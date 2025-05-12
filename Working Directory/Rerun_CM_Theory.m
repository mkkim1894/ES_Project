clc;
clearvars;
tic;

%% Parameter Setup
summary.step = [0.1, 0.1];                  % Step size for evolutionary process
summary.sigW = 2;                           % Standard deviation of bivariate Gaussian fitness landscape
summary.SelectionBias = [1, 2];             % Selection bias in each direction
summary.MutationBias = [1, 1];              % Mutation bias in each direction
summary.InitialAngle = [atan2(1, 6.4), ...  % Initial angles for various mutation directions
                        atan2(1, 3.2), ...
                        atan2(1, 1.6), ...
                        atan2(1, 0.8), ...
                        atan2(1, 0.4), ...
                        atan2(1, 0.2)];
summary.d_ref = -2 * [1; 1];                % Reference distance for fitness calculation
summary.W_ref = exp(-vecnorm(summary.d_ref)^2 / (2 * summary.sigW^2));  % Reference fitness calculation

% Parameters defining evolutionary regimes
summary.N = 10^4;                           % Population size
summary.U = 10^-3;                          % Mutation rate

% List of D values for different evolutionary regimes
Dlist = [10000, 1000, 100, 10];             % D values to iterate over
fname = 'Theory_D_check.mat';               % Output file name

%% Main Loop for D values
for i_n = 1:length(Dlist)
    summary.D = Dlist(i_n);                 % Set the current D value
    Theory = Prediction_CM_D_check(summary); % Call the prediction function
    
    % Store results from Theory.Record for current D
    Theory_D_check(1:length(Theory.Record(:, 1)), i_n) = Theory.Record(:, 1);
end

% Save results to file
save(fname, 'Theory_D_check');

toc; % Display elapsed time
