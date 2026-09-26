function st = drawDPEPanel(D, i, varargin)
% drawDPEPanel - one panel of the distribution of phenotypic effects.
%
%   Draws into the current axes the distribution over theta of the beneficial
%   mutations available at condition i of a struct from computeDPEAtPhenotype:
%
%     grey filled   the distribution, pooled over the max-entropy ensemble at
%                   that phenotype
%     blue stairs   the same distribution weighted by fixation probability
%     orange line   the canonical FGM, which over this window is flat
%     dotted line   the direction to the optimum
%     solid line    the fitness gradient
%
%   The two black lines are directions in the trait space at this phenotype;
%   the orange line is a different model's prediction. They are drawn in
%   different colours because they are different kinds of object.
%
%   The orange line is a SHAPE reference, not a density comparison: the
%   canonical FGM spreads its beneficial mutations evenly over the half-plane
%   facing the gradient, of which this 0-90 degree window shows half, so the
%   line is that flat distribution renormalised to the window.
%
%   The blue curve lies close to the grey one, which is itself the point: the
%   shape is set by the mutational supply, not by selection.
%
%   Lives apart from the figure that draws it so that the panel and the
%   numerics behind it stay in one place.
%
% Name-value pairs
%   'nBins'         histogram bins across the window. Default 30.
%   'showFixed'     overlay the fixation-weighted curve. Default true.
%   'showIsotropic' the flat canonical-FGM reference. Default true.
%   'showKey'       draw the line-sample key in this panel. Default false.
%   'showSkew'      print the skew of the distribution. Default true.
%   'tickSize'      axis font size. Default 8.
%   'labelSize'     axis label font size. Default 9.
%   'showYLabel'    label the y axis. Default true.
%
% Output
%   st  struct of the numbers the panel is built on: mean, meanFixed, sd, skew
%
% See also: computeDPEAtPhenotype, makeFigure_RestrictedTheta
%
% Copyright (c) 2025 Minkyu Kim, Cornell University
% Licensed under MIT License

p = inputParser;
addParameter(p, 'nBins',         30);
addParameter(p, 'showFixed',     true);
addParameter(p, 'showIsotropic', true);
addParameter(p, 'showKey',       false);
addParameter(p, 'showSkew',      true);
addParameter(p, 'tickSize',      8);
addParameter(p, 'labelSize',     9);
addParameter(p, 'showYLabel',    true);
parse(p, varargin{:});
Q = p.Results;

GREY   = [0.68 0.68 0.68];  BLUE  = [42 120 214]/255;  BLACK = [0 0 0];
ORANGE = [235 104 52]/255;

lo   = rad2deg(D.params.thetaRange(1));
hi   = rad2deg(D.params.thetaRange(2));
span = hi - lo;
edges = linspace(lo, hi, Q.nBins + 1);

ang = D.ang{i};  pfx = D.pfx{i};

hold on
histogram(ang, edges, 'Normalization', 'pdf', 'FaceColor', GREY, ...
          'EdgeColor', 'none');

if Q.showFixed
    % The fixation-weighted curve has to be binned by hand: histogram() takes
    % no Weights option, and neither does histcounts.
    bin  = discretize(ang, edges);
    keep = ~isnan(bin);
    hw   = accumarray(bin(keep), pfx(keep), [numel(edges)-1, 1]);
    hw   = hw / (sum(hw) * (edges(2) - edges(1)));    % normalise to a density
    histogram('BinEdges', edges, 'BinCounts', hw, 'DisplayStyle', 'stairs', ...
              'EdgeColor', BLUE, 'LineWidth', 1.5);
end
if Q.showIsotropic
    yline(1/span, '-', 'Color', ORANGE, 'LineWidth', 1.4);
end
xline(D.optDeg(i),  ':', 'Color', BLACK, 'LineWidth', 1.3);
xline(D.gradDeg(i), '-', 'Color', BLACK, 'LineWidth', 1.3);

xlim([lo hi]); box off
set(gca, 'FontSize', Q.tickSize, 'TickDir', 'out', ...
         'XTick', round(linspace(lo, hi, 4)));
xlabel('\theta (deg)', 'FontSize', Q.labelSize);
if Q.showYLabel, ylabel('density', 'FontSize', Q.labelSize); end

st.mean      = mean(ang);
st.sd        = std(ang);
st.skew      = mean((ang - st.mean).^3) / st.sd^3;
st.meanFixed = sum(pfx .* ang) / sum(pfx);

% Skew is the quantity this panel exists to show, so it goes on the panel.
% Both it and the key sit on whichever side the distribution is not.
if st.mean > (lo + hi)/2, xt = 0.04; ha = 'left'; else, xt = 0.96; ha = 'right'; end
if Q.showSkew
    text(xt, 0.95, sprintf('skew %+.1f', st.skew), 'Units', 'normalized', ...
         'HorizontalAlignment', ha, 'FontSize', Q.tickSize + 0.5);
end

% Short line samples rather than coloured words, so the three reference lines
% cannot be confused with one another.
if Q.showKey
    yl = ylim;
    if st.mean > (lo + hi)/2, x0k = lo + 0.06*span; else, x0k = lo + 0.50*span; end
    keyLS  = {'-', ':', '-'};
    keyCol = {BLACK, BLACK, ORANGE};
    keyTxt = {'fitness gradient', 'to optimum', 'isotropic'};
    if ~Q.showIsotropic
        keyLS(3) = [];  keyCol(3) = [];  keyTxt(3) = [];
    end
    if Q.showFixed
        keyLS  = [{'-'}, keyLS];
        keyCol = [{BLUE}, keyCol];
        keyTxt = [{'fixed'}, keyTxt];
    end
    for q = 1:numel(keyTxt)
        yk = yl(1) + (0.80 - 0.105*(q-1)) * diff(yl);
        plot([x0k, x0k + 0.10*span], [yk yk], keyLS{q}, ...
             'Color', keyCol{q}, 'LineWidth', 1.3);
        text(x0k + 0.13*span, yk, keyTxt{q}, 'FontSize', Q.tickSize - 0.5, ...
             'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
    end
end
end
