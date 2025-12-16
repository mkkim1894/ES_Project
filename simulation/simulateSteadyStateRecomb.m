function [results, summary] = simulateSteadyStateRecomb(traitMode, varargin)
% simulateSteadyStateRecomb - Analyze phenotypic variance at mutation-selection-recombination balance.
%
% Description:
%   Runs Wright-Fisher simulations to measure trait variance and covariance
%   at steady state under different recombination regimes. Compares single-trait
%   models against the full two-trait model.
%
% Inputs:
%   traitMode - Which traits to simulate:
%       'both'    - Two-trait model (default)
%       'trait1'  - Only trait 1 active
%       'trait2'  - Only trait 2 active
%   varargin  - Optional name-value pairs:
%       'numRepeat'   - Number of replicates (default: 100)
%       'generations' - Total generations (default: 15000)
%       't_burnin'    - Burn-in generations (default: 10000)
%       'savePath'    - Path to save results (default: '')
%
% Outputs:
%   results - Cell array of simulation results for each initial condition
%   summary - Structure containing:
%       .var1       - Variance of trait 1 across replicates
%       .var2       - Variance of trait 2 across replicates
%       .cov        - Covariance between traits
%       .parameter  - Parameter values [s1, s2, U1, U2, N, R]
%
% Example:
%   [results, summary] = simulateSteadyStateRecomb('both');
%   [results, summary] = simulateSteadyStateRecomb('trait1', 'numRepeat', 50);
%
% Reference:
%   Kim, M., Ardell, S. M., & Kryazhimskiy, S. (2025).
%   "Module-Selection Balance in the Evolution of Modular Organisms."
%
% See also: makeFigureS_SteadyStateRecomb, simulateSteadyStateCM
%
% Copyright (c) 2025 Minkyu Kim, Cornell University
% Licensed under MIT License

%% Parse inputs
p = inputParser;
addRequired(p, 'traitMode', @(x) ismember(x, {'both', 'trait1', 'trait2'}));
addParameter(p, 'numRepeat', 100, @isnumeric);
addParameter(p, 'generations', 15000, @isnumeric);
addParameter(p, 't_burnin', 10000, @isnumeric);
addParameter(p, 'savePath', '', @ischar);
parse(p, traitMode, varargin{:});

numRepeat = p.Results.numRepeat;
generations = p.Results.generations;
t_burnin = p.Results.t_burnin;
savePath = p.Results.savePath;

% Convert traitMode to Trait vector
switch traitMode
    case 'both'
        Trait = [1, 1];
    case 'trait1'
        Trait = [1, 0];
    case 'trait2'
        Trait = [0, 1];
end

%% Parameter Setup
summary.step = 0.1;
summary.sigW = 2;
summary.SelectionBias = [1, 2]; 
summary.InitialAngle = [atan2(1, 6.4), atan2(1, 3.2), atan2(1, 1.6), ...
                        atan2(1, 0.8), atan2(1, 0.4), atan2(1, 0.2)];
summary.d_ref_List = {-1*[1; 1], -2*[1; 1]};

summary.N = 10^4; 
summary.U = 10^-3;
summary.R = tern(strcmp(traitMode, 'both'), 1, 0);  % Full recomb for two-trait, none for single

summary.generations = generations;
summary.t_burnin = t_burnin;

%% Run simulations
fprintf('Running steady-state recombination simulations (mode: %s)...\n', traitMode);
tic

results = cell(1, length(summary.InitialAngle));

for idx = 1:2
    summary.d_ref = summary.d_ref_List{idx};
    summary.W_ref = exp(-(vecnorm(summary.d_ref))^2/2/summary.sigW^2);
    
    for i_pos = 1:length(summary.InitialAngle)
        res = runSteadyStateSimulation(numRepeat, i_pos, summary, Trait);
        j = (idx-1)*length(summary.InitialAngle) + i_pos;
        summary.parameter(j,:) = res.parameter;
        summary.var1(j, 1:numRepeat) = mean(res.var1, 2);
        summary.var2(j, 1:numRepeat) = mean(res.var2, 2);
        summary.cov(j, 1:numRepeat) = mean(res.Covariance, 2);
        
        results{i_pos} = res;
    end
end

fprintf('Simulations completed in %.2f seconds\n', toc);

%% Save results
if ~isempty(savePath)
    switch traitMode
        case 'both'
            fname = fullfile(savePath, 'SteadyStateRecomb_TwoTraits.mat');
        case 'trait1'
            fname = fullfile(savePath, 'SteadyStateRecomb_Trait1.mat');
        case 'trait2'
            fname = fullfile(savePath, 'SteadyStateRecomb_Trait2.mat');
    end
    save(fname, 'results', 'summary');
    fprintf('Results saved to %s\n', fname);
end

end

%% Helper Functions

function res = runSteadyStateSimulation(numRepeat, i_pos, summary, Trait)
% runSteadyStateSimulation - Core simulation loop for steady-state analysis

k = summary.step;
sigW = summary.sigW;

theta = summary.InitialAngle(i_pos); 
d = findInitialPoint(theta, summary.SelectionBias, summary.W_ref, sigW);
d = d(d > 0);
WT = d .* [-cos(theta); -sin(theta)];

rho = abs(summary.d_ref(1)) + abs(summary.d_ref(2));

a = summary.SelectionBias(1);
b = summary.SelectionBias(2);

N = summary.N;
U = summary.U;
R = summary.R;

%% Compute selection coefficients
F0 = -sqrt(a*WT(1)^2 + b*WT(2)^2)^2 / 2 / sigW^2;

if Trait(1) == 1 && Trait(2) == 1
    s1 = (-sqrt(a*(WT(1) + k)^2 + b*WT(2)^2).^2/2/sigW^2) - F0;
    s2 = (-sqrt(a*WT(1)^2 + b*(WT(2) + k)^2).^2/2/sigW^2) - F0;
    U1 = (abs(WT(1))/rho)*U;
    U2 = (abs(WT(2))/rho)*U;
elseif Trait(1) == 1 && Trait(2) == 0
    s1 = (-sqrt(a*(WT(1) + k)^2 + b*WT(2)^2).^2/2/sigW^2) - F0;
    s2 = 0;
    U1 = (abs(WT(1))/rho)*U;
    U2 = 0;
elseif Trait(1) == 0 && Trait(2) == 1
    s1 = 0;
    s2 = (-sqrt(a*WT(1)^2 + b*(WT(2) + k)^2).^2/2/sigW^2) - F0;
    U1 = 0;
    U2 = (abs(WT(2))/rho)*U;
else
    s1 = 0; s2 = 0; U1 = 0; U2 = 0;
end

generations = summary.generations;
t_burnin = summary.t_burnin;

var1 = zeros(numRepeat, generations);
var2 = zeros(numRepeat, generations);
Covariance = zeros(numRepeat, generations);

parfor i_rep = 1:numRepeat
    % Initialize population
    Popmat = zeros(6, N);
    Popmat(1,:) = WT(1);
    Popmat(2,:) = WT(2);
    Popmat(3,:) = 0;
    Popmat(4,:) = 0;
    Popmat(5,:) = U1;
    Popmat(6,:) = U2;

    for t = 1:generations
        %% Mutation
        if U1 ~= 0
            nMut_X = poissrnd(sum(Popmat(5,:)));
            MutID_X = datasample(1:N, nMut_X, 'Weights', Popmat(5,:), 'Replace', false);
            Popmat(1, MutID_X) = Popmat(1, MutID_X) + k;
            Popmat(3, MutID_X) = Popmat(3, MutID_X) + 1;
        end
        
        if U2 ~= 0
            nMut_Y = poissrnd(sum(Popmat(6,:)));
            MutID_Y = datasample(1:N, nMut_Y, 'Weights', Popmat(6,:), 'Replace', false);
            Popmat(2, MutID_Y) = Popmat(2, MutID_Y) + k;
            Popmat(4, MutID_Y) = Popmat(4, MutID_Y) + 1;
        end

        %% Recombination
        if R ~= 0
            if R == 1
                nPair = N/2;
            else
                nPair = mybinornd(N, R/2);
                if 2*nPair > N 
                    nPair = N/2;
                end
            end
            
            RecID = randperm(N, 2*nPair);
            temp = Popmat(2, RecID(1:nPair));
            Popmat(2, RecID(1:nPair)) = Popmat(2, RecID(nPair+1:end));
            Popmat(2, RecID(nPair+1:end)) = temp;
            
            temp = Popmat(4, RecID(1:nPair));
            Popmat(4, RecID(1:nPair)) = Popmat(4, RecID(nPair+1:end));
            Popmat(4, RecID(nPair+1:end)) = temp;
        end
        
        %% Selection 
        F = Popmat(3,:)*s1 + Popmat(4,:)*s2;
        min_fitness = min(F);
        w = exp(F - min_fitness);
        PROB = w / sum(w);
        
        offspring = mnrnd(N, PROB);
        
        newPopmat = zeros(6, N);
        idx = 1;
        for i = 1:N
            if offspring(i) > 0
                newPopmat(:, idx:idx+offspring(i)-1) = repmat(Popmat(:,i), 1, offspring(i));
                idx = idx + offspring(i);
            end
        end
        Popmat = newPopmat;
        
        %% Record variance/covariance
        CovMat = cov(Popmat(1,:), Popmat(2,:));
        var1(i_rep, t) = CovMat(1,1);
        var2(i_rep, t) = CovMat(2,2);
        Covariance(i_rep, t) = CovMat(2,1);
    end
end

res.parameter = [s1; s2; U1; U2; N; R];
res.var1 = var1(:, t_burnin+1:generations);
res.var2 = var2(:, t_burnin+1:generations);
res.Covariance = Covariance(:, t_burnin+1:generations);

end

function d = findInitialPoint(theta, ShapeParameter, Default_fitness, sigW)
% findInitialPoint - Find initial distance for given fitness level

if nargin < 3
    Default_fitness = 0.000123409804086679;
    sigW = 1;
end

a1 = ShapeParameter(1);
a2 = ShapeParameter(2); 

syms r
eqn = Default_fitness == exp(-sqrt(a1*(r*cos(theta))^2 + a2*(r*sin(theta))^2)^2/2/sigW^2); 
sol = solve(eqn, r);
d = double(sol);
end

function res = mybinornd(N, p)
% mybinornd - Fast binomial random draw
res = sum(rand(1, N) < p);
end

function y = tern(cond, a, b)
% tern - Ternary operator
if cond, y = a; else, y = b; end
end
