function bandData = computeTrajectoryPerpendicularBand(resultTableRow, numTimeStamp)
% computeTrajectoryPerpendicularBand - Error band for a 2D phase-plane
% trajectory, built by projecting per-replicate deviation from the mean
% path onto the axis perpendicular to the path's local direction of
% travel, rather than showing raw 2D scatter or separate x1/y1 error bars.
%
% Rationale: at each point along the mean trajectory, a replicate's
% deviation from that point is a 2D vector. Decomposing it into a
% component ALONG the local tangent direction (mostly reflects timing/
% speed differences between replicates -- how far ahead or behind
% schedule that replicate happens to be) and a component PERPENDICULAR to
% it (reflects the replicate's path actually deviating in shape from the
% mean path) isolates the part of the variance that's meaningful to draw
% as a band around the curve. Along-path variance is not drawn: showing
% it would just smear the band along the direction of travel rather than
% thickening it, and is usually not what "uncertainty in the trajectory
% shape" is meant to capture.
%
% Sampling convention matches computeAverageTrajectory.m exactly (each
% replicate sampled at floor(linspace(1,length(traj),numTimeStamp)) row
% indices -- i.e. uniformly spaced BY ROW FRACTION, not by generation),
% so the resulting mean curve here should match the existing
% averageTrajectory/ave data already used for the phase-plane plot.
% This function is intentionally independent of computeAverageTrajectory.m
% (rather than modifying its return signature) so existing callers of
% that function are untouched; it recomputes the same per-replicate
% samples from the raw resultTable in order to also retain what
% computeAverageTrajectory.m discards -- the per-replicate values needed
% to compute any spread at all.
%
% Local tangent direction is estimated from the MEAN path (central
% difference at interior points, forward/backward difference at the two
% endpoints), not from each replicate's own path -- deviations are
% measured relative to a single consistent reference direction at each
% point, not a different axis per replicate.
%
% Inputs:
%   resultTableRow - Cell array of trajectories for ONE starting
%                    configuration (one angle/start), across replicates,
%                    each cell a matrix with columns [fitness, t, x1, x2].
%                    Only columns 3:4 (x1,x2) are used -- this is
%                    specific to 2D phase-plane panels.
%   numTimeStamp   - Number of row-fraction sample points (should match
%                    whatever value produced the existing averageTrajectory)
%
% Outputs (struct):
%   .meanX1, .meanX2       - [1 x numTimeStamp] mean trajectory (should
%                             match the existing ave/averageTrajectory data)
%   .perpStd               - [1 x numTimeStamp] std of perpendicular
%                             deviation at each point
%   .upperX1, .upperX2     - [1 x numTimeStamp] band edge, +1 perpStd
%   .lowerX1, .lowerX2     - [1 x numTimeStamp] band edge, -1 perpStd
%                             (upper/lower here mean "one side/other side
%                             of the path," not literally up/down)
%
% Usage (drawing the ribbon with fill(), same pattern as the existing
% log-ratio shaded bands):
%   fill([bandData.upperX1, fliplr(bandData.lowerX1)], ...
%        [bandData.upperX2, fliplr(bandData.lowerX2)], ...
%        color, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
%
% See also: computeAverageTrajectory

    numSims = numel(resultTableRow);
    allX1 = zeros(numSims, numTimeStamp);
    allX2 = zeros(numSims, numTimeStamp);

    for k = 1:numSims
        traj = resultTableRow{k};
        idx = floor(linspace(1, size(traj,1), numTimeStamp));
        allX1(k, :) = traj(idx, 3)';
        allX2(k, :) = traj(idx, 4)';
    end

    meanX1 = mean(allX1, 1);
    meanX2 = mean(allX2, 1);

    % Local tangent direction of the MEAN path.
    dirX = zeros(1, numTimeStamp);
    dirY = zeros(1, numTimeStamp);
    dirX(1) = meanX1(2) - meanX1(1);
    dirY(1) = meanX2(2) - meanX2(1);
    dirX(end) = meanX1(end) - meanX1(end-1);
    dirY(end) = meanX2(end) - meanX2(end-1);
    for t = 2:numTimeStamp-1
        dirX(t) = meanX1(t+1) - meanX1(t-1);
        dirY(t) = meanX2(t+1) - meanX2(t-1);
    end

    dirNorm = sqrt(dirX.^2 + dirY.^2);
    dirNorm(dirNorm == 0) = 1;   % guard against a locally stalled mean path
    dirX = dirX ./ dirNorm;
    dirY = dirY ./ dirNorm;

    % Perpendicular direction: +90 degree rotation, consistent orientation
    % along the whole path (not re-chosen per point), so the resulting
    % band edges trace a smooth ribbon rather than crossing themselves.
    perpX = -dirY;
    perpY = dirX;

    perpDevAll = zeros(numSims, numTimeStamp);
    for k = 1:numSims
        dx = allX1(k, :) - meanX1;
        dy = allX2(k, :) - meanX2;
        perpDevAll(k, :) = dx .* perpX + dy .* perpY;
    end

    perpStd = std(perpDevAll, 0, 1);

    bandData.meanX1 = meanX1;
    bandData.meanX2 = meanX2;
    bandData.perpStd = perpStd;
    bandData.upperX1 = meanX1 + perpStd .* perpX;
    bandData.upperX2 = meanX2 + perpStd .* perpY;
    bandData.lowerX1 = meanX1 - perpStd .* perpX;
    bandData.lowerX2 = meanX2 - perpStd .* perpY;
end