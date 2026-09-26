function runModularND(varargin)
% runModularND  Modular GPFM with n modules in the concurrent mutations regime.
%
% Simulates replicate populations from each of several initial conditions
% lying on a common fitness contour and differing in how the initial deficit
% is distributed across modules. Output is written to
% results/nModule/Modular_nModule_*.mat.
%
% Name-value pairs
%   'n'                  Number of modules. Default 10.
%   'popSize'            Population size N. Default 1e4.
%   'numGenerations'     Generations per run. Default 1e4.
%   'numReplicates'      Replicates per initial condition. Default 30.
%   'conditions'         nCond-by-2 matrix of initial conditions, each row
%                        [k, c]. Default as below.
%   'geneticTargetSize'  Loci per module L_i. Default [], in which case L_i
%                        accommodates the largest deficit required by any
%                        initial condition.
%
% Initial conditions
%   Each condition distributes the deficit over a two-level profile in which
%   k modules carry share 1 and the remaining n - k carry share c < 1. This
%   represents an organism in which k functions are substantially compromised
%   and the remainder are mildly suboptimal. Small k concentrates adaptation
%   in few modules, so that the effective number of module targets n_g is
%   initially low; c > 0 leaves the remaining modules with beneficial
%   mutations available, so that they can become targets later. Defaults for
%   n = 10:
%
%       [ 1, 0.15]   concentrated
%       [ 3, 0.25]   intermediate
%       [10, 1.00]   uniform (control)
%
% Requires the Statistics and Machine Learning Toolbox.
%
% Reference
%   Kim, M., Ardell, S. M., & Kryazhimskiy, S. (2025). Module-Selection
%   Balance in the Evolution of Modular Organisms.

% ----------------------------- Parse inputs -----------------------------
p = inputParser;
addParameter(p, 'n',                 10,  @isscalar);
addParameter(p, 'popSize',           1e4, @isscalar);
addParameter(p, 'numGenerations',    1e4, @isscalar);
addParameter(p, 'numReplicates',     30,  @isscalar);
addParameter(p, 'conditions',        [],  @(v) isempty(v) || size(v,2) == 2);
addParameter(p, 'geneticTargetSize', [],  @(v) isempty(v) || isscalar(v));
addParameter(p, 'outRoot',           '',  @ischar);   % '' = <project>/results
parse(p, varargin{:});

n              = p.Results.n;
popSize        = p.Results.popSize;
numGenerations = p.Results.numGenerations;
numReplicates  = p.Results.numReplicates;
conditions     = p.Results.conditions;
userTargetSize = p.Results.geneticTargetSize;
outRoot        = p.Results.outRoot;
if isempty(outRoot)
    outRoot = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'results');
end

if isempty(conditions), conditions = defaultConditions(n); end

% ------------------------------ Parameters ------------------------------
deltaTrait      = 0.1;
landscapeStdDev = 2;
mu              = 1e-5;
initialFitness  = 0.25;

ellipseParams = ones(1, n);     % a_i, equal across modules

twoSigSq = 2 * landscapeStdDev^2;
invD     = 1 / deltaTrait;
fitCol   = n + 1;
nCond    = size(conditions, 1);
D0       = -twoSigSq * log(initialFitness);

% -------------------- Initial conditions and genome size ----------------
Xc = zeros(nCond, n);
for i = 1:nCond
    u = twoLevelDirection(n, conditions(i,1), conditions(i,2));
    r = sqrt( D0 / sum((u ./ ellipseParams).^2) );
    Xc(i, :) = -r * u;
end

maxDeficit = max(abs(Xc(:)));
if isempty(userTargetSize)
    targetSizeScalar = max(40, ceil(1.15 * maxDeficit / deltaTrait));
else
    targetSizeScalar = userTargetSize;
end

targetSize = targetSizeScalar * ones(1, n);
totalLoci  = sum(targetSize);
capacity   = targetSize * deltaTrait;

simParams = struct( ...
    'model',             'Modular', ...
    'n',                 n, ...
    'popSize',           popSize, ...
    'deltaTrait',        deltaTrait, ...
    'landscapeStdDev',   landscapeStdDev, ...
    'ellipseParams',     ellipseParams, ...
    'geneticTargetSize', targetSize, ...
    'mutationRate',      mu, ...
    'initialFitness',    initialFitness, ...
    'conditions',        conditions, ...
    'numReplicates',     numReplicates, ...
    'numGenerations',    numGenerations);

fprintf('\n=============== Modular GPFM, n = %d ===============\n', n);
fprintf('N = %g, delta = %.2f, sigma = %g, W0 = %.2f\n', ...
        popSize, deltaTrait, landscapeStdDev, initialFitness);
fprintf('L_i = %d (auto-sized to max deficit %.2f), total loci = %d\n', ...
        targetSizeScalar, maxDeficit, totalLoci);
fprintf('mu = %.3g, MSB ceiling = %.5f\n', mu, exp(-mu * totalLoci));
fprintf('%d conditions x %d replicates, %g generations per run\n\n', ...
        nCond, numReplicates, numGenerations);

% Round onto the delta-lattice and check capacity.
initialPhenotypes = zeros(nCond, n);
nx0 = zeros(nCond, 1);  ng0 = zeros(nCond, 1);  nLive = zeros(nCond, 1);

for i = 1:nCond
    x = -deltaTrait * round(-Xc(i, :) / deltaTrait);
    x = min(x, 0);

    if any(x < -capacity - eps(1))
        bad = find(x < -capacity - eps(1), 1);
        error('runModularND:CapacityExceeded', ...
            'Condition %d puts module %d beyond its capacity L*delta = %.2f.', ...
            i, bad, capacity(bad));
    end
    if all(abs(x) < deltaTrait)
        error('runModularND:EmptyCondition', 'Condition %d has no deficit.', i);
    end

    initialPhenotypes(i, :) = x;
    nx0(i)   = 1 / sum( (abs(x) / sum(abs(x))).^2 );
    ng0(i)   = 1 / sum( (x.^2 / sum(x.^2)).^2 );   % from the initial rates r_i ~ x_i^2
    nLive(i) = nnz(abs(x) >= deltaTrait);
end

condLabels = cell(1, nCond);
for i = 1:nCond
    condLabels{i} = sprintf('k=%d,c=%.2f', conditions(i,1), conditions(i,2));
end
fprintf('condition  : %s\n', sprintf('%14s', condLabels{:}));
fprintf('n_x(0)     : %s\n', sprintf('%9.2f', nx0));
fprintf('n_g(0)     : %s\n', sprintf('%9.2f', ng0));
fprintf('live       : %s\n', sprintf('%9d',   nLive));
fprintf('tail loci  : %s\n', sprintf('%9d',   round(min(abs(initialPhenotypes), [], 2)/deltaTrait)));
fprintf('fix / run  : %s\n', sprintf('%9d',   round(sum(abs(initialPhenotypes), 2)/deltaTrait)));
fprintf('N*U_b(0)   : %s\n\n', ...
        sprintf('%9.1f', popSize * mu * sum(abs(initialPhenotypes), 2) / deltaTrait));

% ------------------------------ Simulation ------------------------------
resultTable = cell(nCond, numReplicates);

for i = 1:nCond
    tCond = tic;
    x0 = initialPhenotypes(i, :);
    W0 = exp( -sum((x0 ./ ellipseParams).^2) / twoSigSq );

    repLogs = cell(1, numReplicates);

    % Replicates are independent. Seeding within the loop body makes the
    % result independent of execution order.
    parfor rep = 1:numReplicates
        rng(1e4*i + rep);

        populationMatrix = zeros(popSize, fitCol);
        populationMatrix(:, 1:n)    = repmat(x0, popSize, 1);
        populationMatrix(:, fitCol) = W0;

        EvolutionLog = zeros(numGenerations, n + 2);   % [W, t, x_1 ... x_n]
        EvolutionLog(1, :) = [W0, 1, x0];

        t = 1;
        while t < numGenerations
            % Mutation. In module m an individual carries b_m = |x_m|/delta
            % derived alleles, each reverting at rate mu, and L_m - b_m
            % intact loci, each breaking at rate mu.
            for m = 1:n
                xcol = populationMatrix(:, m);

                benRates = mu * max(-xcol, 0) * invD;
                delRates = mu * max(capacity(m) + xcol, 0) * invD;

                ids = weightedSample(benRates, poissrnd(sum(benRates)));
                if ~isempty(ids)
                    populationMatrix(ids, m) = populationMatrix(ids, m) + deltaTrait;
                end

                ids = weightedSample(delRates, poissrnd(sum(delRates)));
                if ~isempty(ids)
                    populationMatrix(ids, m) = populationMatrix(ids, m) - deltaTrait;
                end
            end

            populationMatrix(:, fitCol) = ...
                exp( -sum((populationMatrix(:, 1:n) ./ ellipseParams).^2, 2) / twoSigSq );

            % Selection: multinomial offspring numbers, probabilities
            % proportional to Wrightian fitness.
            probs         = populationMatrix(:, fitCol) / sum(populationMatrix(:, fitCol));
            numOffspring  = mnrnd(popSize, probs(:)');
            parentIndices = repelem(1:popSize, numOffspring);
            populationMatrix = populationMatrix(parentIndices(1:popSize), :);

            t = t + 1;
            EvolutionLog(t, :) = [mean(populationMatrix(:, fitCol)), t, ...
                                  mean(populationMatrix(:, 1:n), 1)];
        end

        repLogs{rep} = EvolutionLog(1:t, :);
    end

    resultTable(i, :) = repLogs;

    Wend  = cellfun(@(L) L(end,1), repLogs);
    nxEnd = cellfun(@(L) 1/sum((abs(L(end,3:end))/sum(abs(L(end,3:end)))).^2), repLogs);
    fprintf(['condition %d (k=%d, c=%.2f): W %.4f -> %.4f +/- %.4f, ' ...
             'n_x %.2f -> %.2f +/- %.2f, %.1f s\n'], ...
            i, conditions(i,1), conditions(i,2), W0, mean(Wend), std(Wend), ...
            nx0(i), mean(nxEnd), std(nxEnd), toc(tCond));
end

% -------------------------------- Save ----------------------------------
outDir = fullfile(outRoot, 'nModule');
if ~isfolder(outDir), mkdir(outDir); end

outFile = fullfile(outDir, sprintf('Modular_nModule_n%d_N%.0e.mat', n, popSize));
save(outFile, 'resultTable', 'simParams', ...
     'initialPhenotypes', 'nx0', 'ng0', 'nLive', '-v7.3');

fprintf('\nSaved %s\n', outFile);
end

% ======================================================================
% Helper subfunctions
% ======================================================================

function C = defaultConditions(n)
% Two-level initial conditions [k, c]; see the function header.
    C = [ 1, 0.15;      % concentrated
          3, 0.25;      % intermediate
          n, 1.00 ];    % uniform (control)
end

function u = twoLevelDirection(n, k, c)
% Unit vector in trait space for a two-level deficit profile. Duplicated in
% runPleiotropicND so that both models use identical trait-space targets;
% makeFigure_nModule verifies that the saved conditions agree.
    k = min(max(round(k), 1), n);
    q = [ones(1, k), c * ones(1, n - k)];
    u = q / norm(q);
end

function ids = weightedSample(w, howMany)
% Draw howMany distinct indices with probability proportional to w.
    ids = zeros(1, 0);
    if howMany <= 0, return; end

    % Only indices with positive weight are eligible. Restricting to these
    % is required as well as efficient: a module at its optimum gives its
    % carriers zero beneficial rate, and a cumulative sum over a vector
    % containing zeros has repeated entries, which is not a valid set of
    % bin edges.
    nz = find(w > 0);
    if isempty(nz), return; end

    howMany = min(howMany, numel(nz));

    wn         = w(nz);
    edges      = [0; cumsum(wn(:))];
    edges      = edges / edges(end);
    edges(end) = 1;

    if any(diff(edges) <= 0)
        % Rounding has left the edges non-monotonic; use datasample rather
        % than sampling from an incorrect distribution.
        ids = datasample(nz(:)', howMany, 'Weights', wn, 'Replace', false);
        return
    end

    % Inverse-CDF draw with de-duplication. Mutation counts are small
    % relative to popSize, so repeated passes are rarely needed.
    picked = zeros(1, 0);
    for attempt = 1:10
        need = howMany - numel(picked);
        if need <= 0, break; end
        draw   = discretize(rand(1, ceil(need * 1.5) + 1), edges);
        picked = unique([picked, draw(~isnan(draw))], 'stable');
    end

    ids = nz(picked(1:min(howMany, numel(picked))));
    ids = ids(:)';
end
