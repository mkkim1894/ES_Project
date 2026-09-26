function runPleiotropicND(varargin)
% runPleiotropicND  Pleiotropic GPFM with n traits in the concurrent
% mutations regime.
%
% Simulates replicate populations from each of several initial conditions
% lying on a common fitness contour and differing in how the initial deficit
% is distributed across traits. Initial conditions are constructed as in
% runModularND, so that the two genotype-phenotype maps begin from
% approximately the same points in trait space. Output is written to
% results/nModule/Pleiotropic_nModule_*.mat.
%
% Name-value pairs
%   'n'                  Number of traits. Default 10.
%   'popSize'            Population size N. Default 1e4.
%   'numGenerations'     Generations per run. Default 1e4.
%   'numReplicates'      Replicates per initial condition. Default 30.
%   'conditions'         nCond-by-2 matrix of initial conditions, each row
%                        [k, c]. Must match the value used for runModularND;
%                        see that function for the two-level profile.
%   'geneticTargetSize'  Loci per module L_i, giving a genome size of n*L_i.
%                        Default [], sized as in runModularND so that both
%                        models use the same genome size.
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
resyncInterval  = 1000;
directionSeed   = 2;

ellipseParams = ones(1, n);     % a_i, equal across traits

twoSigSq = 2 * landscapeStdDev^2;
nCond    = size(conditions, 1);
D0       = -twoSigSq * log(initialFitness);

% -------------------- Initial conditions and genome size ----------------
targetPhenotypes = zeros(nCond, n);
for i = 1:nCond
    u = twoLevelDirection(n, conditions(i,1), conditions(i,2));
    r = sqrt( D0 / sum((u ./ ellipseParams).^2) );
    targetPhenotypes(i, :) = -r * u;
end

maxDeficit = max(abs(targetPhenotypes(:)));
if isempty(userTargetSize)
    targetSizeScalar = max(40, ceil(1.15 * maxDeficit / deltaTrait));
else
    targetSizeScalar = userTargetSize;
end

totalLoci      = n * targetSizeScalar;
genomeWideRate = mu * totalLoci;    % U, per individual per generation

simParams = struct( ...
    'model',            'Pleiotropic', ...
    'n',                n, ...
    'popSize',          popSize, ...
    'deltaTrait',       deltaTrait, ...
    'landscapeStdDev',  landscapeStdDev, ...
    'ellipseParams',    ellipseParams, ...
    'totalLoci',        totalLoci, ...
    'mutationRate',     mu, ...
    'initialFitness',   initialFitness, ...
    'conditions',       conditions, ...
    'numReplicates',    numReplicates, ...
    'numGenerations',   numGenerations, ...
    'resyncInterval',   resyncInterval, ...
    'directionSeed',    directionSeed);

fprintf('\n============= Pleiotropic GPFM, n = %d =============\n', n);
fprintf('N = %g, delta = %.2f, sigma = %g, W0 = %.2f\n', ...
        popSize, deltaTrait, landscapeStdDev, initialFitness);
fprintf('total loci = %d (auto-sized to max deficit %.2f), mu = %.3g\n', ...
        totalLoci, maxDeficit, mu);
fprintf('U = %.4f, MSB ceiling = %.5f\n', genomeWideRate, exp(-mu * totalLoci));
fprintf('%d conditions x %d replicates, %g generations per run\n\n', ...
        nCond, numReplicates, numGenerations);

% ---------------------------- Locus directions --------------------------
% Each locus carries a fixed direction drawn from the uniform distribution
% on the unit sphere in n dimensions, generalizing the uniformly
% distributed angle of the two-trait model. Directions are drawn once from
% a private stream and held fixed across all conditions and replicates.
s = RandStream('threefry', 'Seed', directionSeed);
v = randn(s, totalLoci, n);
alleleDirections = v ./ vecnorm(v, 2, 2);
alleleEffects    = -deltaTrait * alleleDirections;

% -------------------------- Initial genomes -----------------------------
% Initial genotypes are SAMPLED at the target phenotype rather than built
% greedily, so that this model is initialized the same way as the two-module
% simulations (utils/initializeGenomeSampled).
%
% WHY. The greedy construction switches on whichever locus is best aligned
% with the remaining displacement, so the genotype it returns is not typical
% of the genotypes sitting at that phenotype: it uses few, well-aligned loci.
% In two traits that measurably distorted the supply of beneficial mutations.
% In ten traits it does not - measured at these targets, the number of
% beneficial mutations and their mean selection coefficient are the same to
% within a few per cent for genotypes carrying anywhere from 100 to 280 of the
% 400 loci - but the two model families should be initialized the same way
% regardless, and the Methods says they are.
%
% HOW. Greedy still runs, but only to choose a feasible number of 1s. The
% genotype is then redrawn at random with that many 1s and moved onto the
% target by swaps: exchange a 1 with a 0, keep the exchange when it reduces
% the distance. Swaps preserve the number of 1s, so what comes out is a
% genotype chosen without regard to alignment.
%
% WHY NOT THE WALK-AND-ANNEAL SAMPLER. Its first stage walks outward until it
% crosses the target fitness contour. In ten dimensions 400 isotropic steps of
% delta cover an RMS distance of only delta*sqrt(400) = 2.0, against a contour
% at 3.33, so the walk never arrives and every draw would have to be restarted.
MAX_SWAPS = 2.5e5;
SWAP_TOL  = 0.05;
initSeed  = 7;
rInit     = RandStream('threefry', 'Seed', initSeed);

initialPhenotypes = zeros(nCond, n);
initialGenomes    = false(nCond, totalLoci);
matchError        = zeros(nCond, 1);
nx0 = zeros(nCond, 1);  ng0 = zeros(nCond, 1);

for i = 1:nCond
    target = targetPhenotypes(i, :);

    % pass 1: greedy, used only to fix the number of 1s
    gGreedy = false(1, totalLoci);
    xGreedy = zeros(1, n);
    for ell = 1:totalLoci
        xNew = xGreedy + alleleEffects(ell, :);
        if norm(xNew - target) < norm(xGreedy - target)
            gGreedy(ell) = true;
            xGreedy      = xNew;
        end
    end
    nOnes   = sum(gGreedy);
    dGreedy = norm(xGreedy - target);

    % pass 2: random genotype with the same number of 1s, landed by swaps
    genome            = false(1, totalLoci);
    genome(randperm(rInit, totalLoci, nOnes)) = true;
    x        = sum(alleleEffects(genome, :), 1);
    d        = norm(x - target);
    onesIdx  = find(genome);
    zerosIdx = find(~genome);

    for it = 1:MAX_SWAPS
        if d <= SWAP_TOL, break; end
        p1 = randi(rInit, numel(onesIdx));
        p0 = randi(rInit, numel(zerosIdx));
        i1 = onesIdx(p1);  i0 = zerosIdx(p0);
        xNew = x - alleleEffects(i1, :) + alleleEffects(i0, :);
        dNew = norm(xNew - target);
        if dNew < d
            genome(i1) = false;  genome(i0) = true;
            onesIdx(p1) = i0;    zerosIdx(p0) = i1;
            x = xNew;  d = dNew;
        end
    end

    % Never land worse than greedy did.
    if d > dGreedy
        fprintf(['  cond %d: sampled genotype landed at %.3f vs greedy %.3f - ' ...
                 'keeping greedy\n'], i, d, dGreedy);
        genome = gGreedy;  x = xGreedy;  d = dGreedy;
    end

    initialPhenotypes(i, :) = x;
    initialGenomes(i, :)    = genome;
    matchError(i)           = d;
    nx0(i) = 1 / sum( (abs(x) / sum(abs(x))).^2 );
    ng0(i) = 1 / sum( (x.^2 / sum(x.^2)).^2 );
end

condLabels = cell(1, nCond);
for i = 1:nCond
    condLabels{i} = sprintf('k=%d,c=%.2f', conditions(i,1), conditions(i,2));
end
fprintf('condition  : %s\n', sprintf('%14s', condLabels{:}));
fprintf('n_x(0)     : %s\n', sprintf('%9.2f', nx0));
fprintf('n_g(0)     : %s\n', sprintf('%9.2f', ng0));
fprintf('match err  : %s   (delta = %.2f)\n', sprintf('%9.2f', matchError), deltaTrait);
fprintf('derived    : %s   of %d loci\n\n', ...
        sprintf('%9d', sum(initialGenomes, 2)), totalLoci);

if any(matchError > 5 * deltaTrait)
    warning('runPleiotropicND:PoorGenomeMatch', ...
        ['Match error exceeds 5*delta for at least one condition; the two ' ...
         'GPFMs may not begin from comparable states.']);
end

% ------------------------------ Simulation ------------------------------
resultTable = cell(nCond, numReplicates);

for i = 1:nCond
    tCond = tic;
    g0 = initialGenomes(i, :);

    repLogs = cell(1, numReplicates);

    % Replicates are independent. Seeding within the loop body makes the
    % result independent of execution order.
    parfor rep = 1:numReplicates
        rng(2e4*i + rep);

        genomeMatrix      = repmat(g0, popSize, 1);
        currentPhenotypes = double(genomeMatrix) * alleleEffects;
        currentFitness    = exp( -sum((currentPhenotypes ./ ellipseParams).^2, 2) / twoSigSq );

        EvolutionLog = zeros(numGenerations, n + 2);   % [W, t, x_1 ... x_n]
        EvolutionLog(1, :) = [mean(currentFitness), 1, mean(currentPhenotypes, 1)];

        t = 1;
        while t < numGenerations
            t = t + 1;

            % Mutation. Events per generation are Poisson with mean N*U;
            % each flips one allele and displaces that individual's
            % phenotype by the corresponding effect vector.
            totalMutations = poissrnd(popSize * genomeWideRate);
            if totalMutations > 0
                mutInd = randi(popSize,   totalMutations, 1);
                mutLoc = randi(totalLoci, totalMutations, 1);

                for idx = 1:totalMutations
                    ind = mutInd(idx);
                    loc = mutLoc(idx);
                    if genomeMatrix(ind, loc)
                        currentPhenotypes(ind, :) = currentPhenotypes(ind, :) - alleleEffects(loc, :);
                        genomeMatrix(ind, loc) = false;
                    else
                        currentPhenotypes(ind, :) = currentPhenotypes(ind, :) + alleleEffects(loc, :);
                        genomeMatrix(ind, loc) = true;
                    end
                end
            end

            currentFitness = exp( -sum((currentPhenotypes ./ ellipseParams).^2, 2) / twoSigSq );

            % Selection: multinomial offspring numbers, probabilities
            % proportional to Wrightian fitness.
            probs         = currentFitness / sum(currentFitness);
            numOffspring  = mnrnd(popSize, probs(:)');
            parentIndices = repelem(1:popSize, numOffspring);
            parentIndices = parentIndices(1:popSize);

            genomeMatrix      = genomeMatrix(parentIndices, :);
            currentPhenotypes = currentPhenotypes(parentIndices, :);
            currentFitness    = currentFitness(parentIndices);

            % Phenotypes are updated incrementally at mutation; the exact
            % genotype-phenotype product is re-evaluated periodically to
            % bound accumulated rounding.
            if mod(t, resyncInterval) == 0
                currentPhenotypes = double(genomeMatrix) * alleleEffects;
                currentFitness    = exp( -sum((currentPhenotypes ./ ellipseParams).^2, 2) / twoSigSq );
            end

            EvolutionLog(t, :) = [mean(currentFitness), t, mean(currentPhenotypes, 1)];
        end

        repLogs{rep} = EvolutionLog(1:t, :);
    end

    resultTable(i, :) = repLogs;

    Wend  = cellfun(@(L) L(end,1), repLogs);
    nxEnd = cellfun(@(L) 1/sum((abs(L(end,3:end))/sum(abs(L(end,3:end)))).^2), repLogs);
    fprintf(['condition %d (k=%d, c=%.2f): W %.4f -> %.4f +/- %.4f, ' ...
             'n_x %.2f -> %.2f +/- %.2f, %.1f s\n'], ...
            i, conditions(i,1), conditions(i,2), ...
            mean(cellfun(@(L) L(1,1), repLogs)), mean(Wend), std(Wend), ...
            nx0(i), mean(nxEnd), std(nxEnd), toc(tCond));
end

% -------------------------------- Save ----------------------------------
outDir = fullfile(outRoot, 'nModule');
if ~isfolder(outDir), mkdir(outDir); end

outFile = fullfile(outDir, sprintf('Pleiotropic_nModule_n%d_N%.0e.mat', n, popSize));
save(outFile, 'resultTable', 'simParams', ...
     'initialPhenotypes', 'targetPhenotypes', 'initialGenomes', ...
     'alleleDirections', 'matchError', 'nx0', 'ng0', '-v7.3');

fprintf('\nSaved %s\n', outFile);
end

% ======================================================================
% Helper subfunctions
% ======================================================================

function C = defaultConditions(n)
% Two-level initial conditions [k, c]. Must match runModularND.
    C = [ 1, 0.15;      % concentrated
          3, 0.25;      % intermediate
          n, 1.00 ];    % uniform (control)
end

function u = twoLevelDirection(n, k, c)
% Unit vector in trait space for a two-level deficit profile. Duplicated in
% runModularND so that both models use identical trait-space targets;
% makeFigure_nModule verifies that the saved conditions agree.
    k = min(max(round(k), 1), n);
    q = [ones(1, k), c * ones(1, n - k)];
    u = q / norm(q);
end
