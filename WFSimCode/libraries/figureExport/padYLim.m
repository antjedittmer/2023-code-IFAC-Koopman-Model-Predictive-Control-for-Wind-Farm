function padYLim(h, frac)
%PADYLIM Y axis limits not completely tight: lines at the data extremes
%stay visible next to the axis
%   padYLim(h) widens the y limits of all axes of figure h so that the
%   lines (line, stair) within the x limits keep a margin of 5 % of the
%   data range to the axis; limits set with a larger margin are kept (only
%   widened), and data outside the current limits (clipped on purpose,
%   e.g. the start-up) are not taken into account. Log axes: margin in
%   decades. Call it after all limits are set.
%   padYLim(h, frac) uses the margin frac (fraction of the data range).

if nargin < 2 || isempty(frac)
    frac = 0.05;
end
for ax = findall(h, 'Type', 'axes')'
    isLog = strcmp(ax.YScale, 'log');
    xl = ax.XLim;
    y = [];
    for hl = [findall(ax, 'Type', 'line'); findall(ax, 'Type', 'stair')]'
        x = hl.XData; yd = hl.YData;
        if numel(x) == numel(yd)
            yd = yd(x >= xl(1) & x <= xl(2));
        end
        y = [y, yd(:)']; %#ok<AGROW>
    end
    y = y(isfinite(y));
    yl = ax.YLim;
    if isLog
        y = log10(y(y > 0));
        yl = log10(yl);
    end
    y = y(y >= yl(1) & y <= yl(2)); % data clipped on purpose (e.g. start-up) stay outside
    if isempty(y)
        continue
    end
    yMin = min(y); yMax = max(y);
    d = yMax - yMin;
    if d == 0 % constant: margin relative to the value
        d = max(abs(yMax), 1);
    end
    yl = [min(yl(1), yMin - frac*d), max(yl(2), yMax + frac*d)];
    if isLog
        yl = 10.^yl;
    end
    ax.YLim = yl;
end
end
