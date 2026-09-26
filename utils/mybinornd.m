function r = mybinornd(n, p)
% mybinornd - Binomial random draw, without the Statistics Toolbox.
%
%   REIMPLEMENTATION. The original mybinornd used by the published runs is not
%   in the repository. This is a drop-in with the same contract: r ~ Binomial(n, p),
%   elementwise if n or p is an array. It is statistically equivalent, but it
%   consumes the random stream differently, so a run using this file will not
%   reproduce a published run seed-for-seed. Results are unchanged in
%   distribution. If the original file turns up, prefer it.
%
%   The overwhelmingly common call in this codebase is mybinornd(1, p) - a single
%   Bernoulli draw, used for whether a beneficial mutation escapes drift - which
%   is special-cased.
%
% Inputs:
%   n - number of trials (scalar or array of non-negative integers)
%   p - success probability (scalar or array in [0, 1])
%
% Outputs:
%   r - binomial draws, size max(size(n), size(p))
%
% Copyright (c) 2025 Minkyu Kim, Cornell University
% Licensed under MIT License

    if isscalar(n) && isscalar(p) && n == 1
        r = double(rand() < p);                 % the hot path
        return;
    end

    if isscalar(n) && ~isscalar(p), n = repmat(n, size(p)); end
    if isscalar(p) && ~isscalar(n), p = repmat(p, size(n)); end
    if ~isequal(size(n), size(p))
        error('mybinornd:SizeMismatch', 'n and p must be the same size, or scalar.');
    end

    p = min(max(p, 0), 1);
    r = zeros(size(n));
    for k = 1:numel(n)
        nk = n(k);
        if nk <= 0 || p(k) <= 0
            r(k) = 0;
        elseif p(k) >= 1
            r(k) = nk;
        else
            r(k) = sum(rand(1, nk) < p(k));
        end
    end
end
