function setFigureSize(h, factor, posScreen)
%SETFIGURESIZE Fixed figure size for display and export, independent of the screen
%   setFigureSize(h, factor) sets the export (paper) size of figure h to
%   factor = [fw fh] times the standard figure of 560 x 420 px, i.e.
%   14.82 x 11.11 cm at 96 dpi, e.g. [1 1.6] for the time series with 5
%   subplots, [1 1.3] for 4 subplots. printFigure exports with this size,
%   so that the output does not depend on the window size, which Windows
%   shrinks to the screen (as in the qLMPC code).
%   setFigureSize(h, factor, posScreen) also sets the window position
%   posScreen ([left bottom width height] in px).

sizeCm = [14.82 11.11].*factor;
if nargin > 2 && ~isempty(posScreen)
    set(h, 'Position', posScreen);
end
set(h, 'PaperUnits', 'centimeters', 'PaperPositionMode', 'manual', ...
    'PaperSize', sizeCm, 'PaperPosition', [0 0 sizeCm]);
end
