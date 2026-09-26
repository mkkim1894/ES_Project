function D = computeDPEAtPhenotype(varargin)
% computeDPEAtPhenotype - the distribution of phenotypic effects of the
%   mutations available at a phenotype, for each of several initial conditions.
%
%   No simulation is involved. At each initial condition this draws genotypes
%   from the maximum-entropy ensemble at that phenotype, enumerates every
%   single-locus mutation of every drawn genotype, and returns the directions
%   of the beneficial ones together with their fixation probabilities. It is
%   therefore the mutational supply the population sees before any evolution
%   happens.
%
%   The numerics live here rather than in the figure function so that the
%   quantity and the panel that shows it stay separable: the top row of
%   makeFigure_RestrictedTheta draws what this returns.
%
% Name-value pairs
%   'theta'      the locus angles, 1 x L. Given, these are used as they stand
%                and 'L' and 'seed' do not apply to them; pass
%                genomeParams.genomeTheta to describe the genome a simulation
%                actually used rather than another draw from the same cone.
%   'L'          number of loci, when 'theta' is not given. Default 400.
%   'delta'      phenotypic effect size of a mutation. Default 0.1.
%   'N'          population size, for the Kimura fixation probability. 1e4.
%   'a'          ellipseParams. Default [1, 1/sqrt(2)].
%   'sigma'      landscapeStdDev. Default 2.
%   'W0'         fitness of the contour the initial conditions sit on. 0.25.
%   'thetaRange' the mutational cone. Default [0, pi/2].
%   'X0'         nCond x 2 phenotypes to evaluate at, given directly. Takes
%                precedence over 'R0list'; pass simParams.initialPhenotypes to
%                evaluate at the initial conditions a simulation actually used.
%   'R0list'     initial x_2/x_1 ratios, one per condition, placed on the W0
%                contour. Used only when 'X0' is empty. [0.16, 0.625, 5].
%   'nGeno'      genotypes pooled per condition. Default 2000.
%   'tol'        tolerance window on the target phenotype. Default 0.02.
%   'batch'      draws per rejection batch. Default 20000.
%   'maxDraws'   give up after this many draws. Default 4e6.
%   'seed'       seeds both the locus angles and the draws. Default 1.
%
% Output, a struct with one entry per field and one cell or row per condition
%   ang{i}     directions of the available beneficial mutations, degrees
%   pfx{i}     their fixation probabilities, same length as ang{i}
%   optDeg(i)  direction from the phenotype to the optimum, degrees
%   gradDeg(i) direction of the fitness gradient there, degrees
%   bMean(i)   mean number of beneficial mutations available per genotype
%   X0         nCond x 2 initial phenotypes on the W0 contour
%   R0list     the ratios, echoed
%   params     the parsed parameters, so a cached result can be checked
%
% See also: solveMaxEntTilt, drawDPEPanel, makeFigure_RestrictedTheta
%
% Copyright (c) 2025 Minkyu Kim, Cornell University
% Licensed under MIT License

p = inputParser;
addParameter(p, 'theta',      []);
addParameter(p, 'L',          400);
addParameter(p, 'delta',      0.1);
addParameter(p, 'N',          1e4);
addParameter(p, 'a',          [1, 1/sqrt(2)]);
addParameter(p, 'sigma',      2);
addParameter(p, 'W0',         0.25);
addParameter(p, 'thetaRange', [0, pi/2]);
addParameter(p, 'X0',         []);
addParameter(p, 'R0list',     [0.16, 0.625, 5]);
addParameter(p, 'nGeno',      2000);
addParameter(p, 'tol',        0.02);
addParameter(p, 'batch',      20000);
addParameter(p, 'maxDraws',   4e6);
addParameter(p, 'seed',       1);
addParameter(p, 'verbose',    true);
parse(p, varargin{:});
P = p.Results;

if isempty(P.X0)
    % initial phenotypes on the W0 contour, the same construction as
    % findInitialPhenotypes but solved in closed form rather than symbolically
    F0 = log(P.W0) * 2 * P.sigma^2;
    x1 = -sqrt(-F0 ./ (1./P.a(1)^2 + P.R0list.^2./P.a(2)^2));
    X0 = [x1(:), (P.R0list(:) .* x1(:))];
else
    X0 = P.X0;
    if size(X0, 2) ~= 2
        error('computeDPEAtPhenotype:BadX0', '''X0'' must be nCond x 2.');
    end
end
R0list = (X0(:,2) ./ X0(:,1))';
nC     = size(X0, 1);

rng(P.seed, 'twister');
if isempty(P.theta)
    theta = P.thetaRange(1) + diff(P.thetaRange) * rand(1, P.L);
else
    theta = reshape(P.theta, 1, []);
    P.L   = numel(theta);
end
cosT  = cos(theta);  sinT = sin(theta);
e     = [cosT(:), sinT(:)];
logW  = @(u, v) -((u./P.a(1)).^2 + (v./P.a(2)).^2) ./ (2*P.sigma^2);

D.ang = cell(nC,1);  D.pfx = cell(nC,1);
D.optDeg = zeros(nC,1);  D.gradDeg = zeros(nC,1);  D.bMean = zeros(nC,1);

if P.verbose
    fprintf('computeDPEAtPhenotype: pooling %d genotypes at each of %d conditions\n', ...
            P.nGeno, nC);
end

for i = 1:nC
    xT = X0(i, :);
    y0 = -xT(:) / P.delta;
    [lam, pl] = solveMaxEntTilt(e, y0, sprintf('R_0 = %g', R0list(i)));

    % draw from the tilted ensemble, keep those inside the tolerance window,
    % importance weight back to uniform over the window, then resample nGeno
    tolY = P.tol / P.delta;
    Gk = false(0, P.L);  Wk = zeros(0,1);  nDraws = 0;
    while size(Gk,1) < P.nGeno && nDraws < P.maxDraws
        B  = min(P.batch, P.maxDraws - nDraws);
        G  = rand(B, P.L) < pl';
        dy = double(G) * e - y0';
        hit = hypot(dy(:,1), dy(:,2)) <= tolY;
        nDraws = nDraws + B;
        if any(hit)
            Gk = [Gk; G(hit,:)];                       %#ok<AGROW>
            Wk = [Wk; -(dy(hit,:) * lam)];             %#ok<AGROW>
        end
    end
    if isempty(Gk)
        error('computeDPEAtPhenotype:NoAccept', ...
              'No draw landed within %g of condition %d (R_0 = %g).', ...
              P.tol, i, R0list(i));
    end
    % resample with probability proportional to the importance weight, by
    % inverse CDF. randsample would need the Statistics Toolbox, which the rest
    % of this repository avoids.
    w  = exp(Wk - max(Wk));  w = w / sum(w);
    cw = [0; cumsum(w(:))];  cw(end) = 1;
    u  = min(rand(P.nGeno, 1), 1 - eps);
    pick = discretize(u, cw);
    Gk = double(Gk(pick, :));

    % phenotypic effects of every single-locus mutation, for every genotype
    signs = 2*Gk - 1;                                  % nG x L
    d1 = signs .* (P.delta * cosT);
    d2 = signs .* (P.delta * sinT);
    xq = [-P.delta * (Gk * cosT'), -P.delta * (Gk * sinT')];   % nG x 2
    s  = logW(xq(:,1) + d1, xq(:,2) + d2) - logW(xq(:,1), xq(:,2));
    ben = s > 0;

    D.ang{i} = rad2deg(atan2(d2(ben), d1(ben)));
    sb       = s(ben);
    D.pfx{i} = (1 - exp(-2*sb)) ./ (1 - exp(-2*P.N*sb));
    D.bMean(i) = mean(sum(ben, 2));

    D.optDeg(i) = rad2deg(atan2(-xT(2), -xT(1)));
    gr = [-xT(1)/P.a(1)^2, -xT(2)/P.a(2)^2];
    D.gradDeg(i) = rad2deg(atan2(gr(2), gr(1)));
end

D.X0     = X0;
D.R0list = R0list;
D.params = P;
end
