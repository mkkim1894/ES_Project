function [lam, p] = solveMaxEntTilt(e, y0, condLabel)
% solveMaxEntTilt - find the tilt that makes a target phenotype typical.
%
%   Solves the 2-by-2 (or n-by-n) system
%
%       sum_l sigmoid(<lambda, e_l>) e_l = y_0
%
%   which is equation (3) of the maximum-entropy construction . Under the
%   tilted ensemble P_lambda(g) ~ exp(<lambda, y(g)>) the loci are independent
%   Bernoulli with p_l = sigmoid(<lambda, e_l>), and this choice of lambda puts
%   the mean of the ensemble on y_0, so the target stops being a rare event.
%
%   The Jacobian sum_l p_l(1-p_l) e_l e_l' is positive definite, so the solution
%   is unique wherever it exists and Newton converges; a backtracking line
%   search keeps it from overshooting.
%
% Inputs
%   e          L x d matrix of unit site vectors, one row per locus
%   y0         d x 1 target, in units of delta (y = -x/delta)
%   condLabel  optional string used in error messages
%
% Outputs
%   lam        d x 1 tilt vector
%   p          L x 1 Bernoulli probabilities sigmoid(<lambda, e_l>)
%
% Errors if y0 lies outside the reachable set {sum_l c_l e_l : c_l in [0,1]},
% in which case no genotype has that phenotype and no tilt exists.
%
% See also: initializeGenomeMaxEnt, computeDPEAtPhenotype
%
% Copyright (c) 2025 Minkyu Kim, Cornell University
% Licensed under MIT License

    if nargin < 3, condLabel = ''; end
    y0 = y0(:);
    assertReachable(e, y0, condLabel);

    lam = zeros(size(e, 2), 1);
    for it = 1:200
        p = 1 ./ (1 + exp(-(e * lam)));
        F = (p' * e)' - y0;
        if max(abs(F)) < 1e-12, return; end
        w = p .* (1 - p);
        J = e' * (w .* e);
        if rcond(J) < 1e-14
            error('solveMaxEntTilt:Singular', ...
                  ['%sTilt Jacobian is singular; the target is probably on the ' ...
                   'boundary of the reachable set.'], tag(condLabel));
        end
        step = -(J \ F);
        t = 1; f0 = norm(F);
        for bt = 1:40
            pn = 1 ./ (1 + exp(-(e * (lam + t*step))));
            if norm((pn' * e)' - y0) < f0, break; end
            t = t / 2;
        end
        lam = lam + t*step;
    end
    error('solveMaxEntTilt:NoConvergence', ...
          '%sNewton did not converge on lambda in 200 iterations.', tag(condLabel));
end

% ---------------------------------------------------------------------------
function assertReachable(e, y0, condLabel)
% y0 must lie inside the zonotope {sum_l c_l e_l : c_l in [0,1]}. Checked by the
% support function along a sweep of directions (exact in 2D up to the sweep
% resolution; a necessary condition in higher dimensions).
    d = size(e, 2);
    if d == 2
        a = linspace(0, 2*pi, 721)';
        U = [cos(a), sin(a)];
    else
        U = randn(4000, d);
        U = U ./ vecnorm(U, 2, 2);
    end
    support = sum(max(0, e * U'), 1)';
    need    = U * y0;
    slack   = min(support - need);
    if slack <= 0
        error('solveMaxEntTilt:Unreachable', ...
            ['%sNo genotype has this phenotype. The target needs %.2f units of ' ...
             'displacement in a direction where the loci can supply at most ' ...
             '%.2f. Move the target inside the cone, lower W0, or raise L.'], ...
            tag(condLabel), max(need), max(support));
    end
    if slack < 0.5
        warning('solveMaxEntTilt:NearBoundary', ...
            ['%sTarget sits within %.2f of the edge of the reachable set. Almost ' ...
             'every locus that can point that way must carry allele 1, so the ' ...
             'ensemble is nearly a single genotype.'], tag(condLabel), slack);
    end
end

function s = tag(condLabel)
    if isempty(condLabel), s = ''; else, s = [condLabel ': ']; end
end
