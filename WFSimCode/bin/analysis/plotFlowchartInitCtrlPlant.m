function notCalled = plotFlowchartInitCtrlPlant(variant, profInfo, dirFig)
%PLOTFLOWCHARTINITCTRLPLANT WFSim program structure diagram (thesis figure)
%   plotFlowchartInitCtrlPlant(variant) draws the WFSim function call tree
%   of a closed-loop simulation in three boxes: initialization (black),
%   wind farm model (blue) and controller (red), and saves it as
%   FlowchartInitCtrlPlant_<variant> (png, eps) in dirFig (default
%   figuresBlockDiagram_IFAC23 next to this file):
%   - variant 1: original WFSim, nine turbines, NMPC with YALMIP
%   - variant 2: modified WFSim, two turbines, Koopman MPC (QP). Modified
%     functions in the colour of their box, functions added for the
%     Koopman MPC in bold, unchanged functions grey.
%   The functions and their nesting are those of a MATLAB profiler run of
%   WFSim_demo, reordered by their dependencies (order of the calls).
%
%   notCalled = plotFlowchartInitCtrlPlant(variant, profInfo) also checks the
%   diagram against the profiler data profInfo = profile('info') of a
%   WFSim_demo run (profile on; WFSim_demo; profInfo = profile('info')):
%   notCalled lists the functions of the diagram that the run did not call.
%   For variant 1 only the initialization and the farm model can be checked
%   with the current WFSim_demo (two turbines, no YALMIP).

if nargin < 1 || isempty(variant)
    variant = 2;
end
if nargin < 2
    profInfo = [];
end
if nargin < 3 || isempty(dirFig)
    dirFig = fullfile(fileparts(mfilename('fullpath')), 'figuresBlockDiagram_IFAC23');
end

% {text, indentation level, box (1 init, 2 farm model, 3 controller),
%  style ('' unchanged: grey, 'm' modified: box colour, 'n' new: box colour, bold)}
farmModel = {
    'WFSim_timestepping',   1, 2, ''
    'Make_Ax_b',            2, 2, ''
    'SpatialDiscr_Hybrid',  3, 2, ''
    'Turbulence',           4, 2, ''
    'ConstructLmu',         5, 2, ''
    'diskfilter',           6, 2, ''
    'Dynamical',            3, 2, ''
    'Actuator',             3, 2, 'm'
    'BoundaryConditions',   3, 2, 'm'
    'MakingSparseMatrix (for A_x)', 3, 2, ''
    'MakingSparseMatrix (for A_y)', 3, 2, ''
    'Computesol',           2, 2, ''
    'MapSolution',          2, 2, ''
    'Updateboundaries',     3, 2, ''};
if variant == 1
    farmModel(:,4) = {''};                                                 % original: all grey
    entries = [{
        'WFSim_demo',                           0, 1, ''
        'layoutSet_sowfa_9turb_apc_alm_turbl',  1, 1, ''
        'controlSet_sowfa_9turb_apc_alm_turbl', 1, 1, ''
        'solverSet_default',                    1, 1, ''
        'WFSim_addpaths',                       1, 1, ''
        'InitWFSim',                            1, 1, ''
        'meshing',                              2, 1, ''
        'Compute_B1_B2_bc',                     2, 1, ''}
        farmModel
        {'NMPCcontroller',                      1, 3, ''
        'MPCinit',                              2, 3, ''
        'yalmip',                               2, 3, ''
        'wfmodel',                              3, 3, ''
        'matsys',                               3, 3, ''
        'optimize',                             3, 3, ''
        'wfmodel',                              3, 3, ''
        'matsys',                               3, 3, ''
        'optimize',                             3, 3, ''}];
    yTop = 44;
else
    entries = [{
        'WFSim_demo',                           0, 1, 'm'
        'WFSim_addpaths',                       1, 1, 'm'
        'layoutSet_sowfa_2turb_alm_turbl',      1, 1, 'm'
        'controlSet_sowfa_2turb_alm_turbl',     1, 1, 'm'
        'solverSet_default',                    1, 1, ''
        'InitWFSim',                            1, 1, ''
        'meshing',                              2, 1, ''
        'Compute_B1_B2_bc',                     2, 1, ''
        'MPCinitK',                             1, 1, 'n'
        'koopmanstateextensionWFSim',           2, 1, 'n'}
        farmModel
        {'qLPVMPCcontroller_dU_quadprog',       1, 3, 'n'
        'MPCinit_ctrl',                         2, 3, 'n'
        'wfmodel',                              2, 3, 'm'
        'Matsys',                               2, 3, 'm'
        'quadprog_sol',                         2, 3, 'n'
        'H_G_matrix_dU',                        3, 3, 'n'}];
    yTop = 51;
end

% Layout in pixels of the png (150 dpi, y downwards), as the thesis figure
wFig = 875; hFig = 919;
x0 = 113; dx = 34.3;                                                       % first indentation, per level
y0 = 70; dy = 24.3;                                                        % first line, line spacing
xBox = [95, 630];
xLabel = 404;
boxCol = {[0 0 0], [0 0 1], [1 0 0]};
boxName = {'Initialization', 'Farm Model', 'Controller'};
grey = 0.5*[1 1 1];
fs = 10; fsLabel = 16;

fig = figure('Color', 'w', 'Units', 'pixels');
fig.Position(3:4) = [wFig, hFig]*96/150;
ax = axes(fig, 'Position', [0 0 1 1], 'XLim', [0 wFig], 'YLim', [0 hFig], 'YDir', 'reverse');
axis(ax, 'off'); hold(ax, 'on');
nE = size(entries, 1);
y = y0 + (0:nE-1)*dy;
box = [entries{:,3}];
xRight = zeros(nE, 1);                                                     % right end of the texts
for k = 1:nE
    [str, lvl, b, st] = entries{k,:};
    cl = grey; fw = 'normal';
    if ~isempty(st), cl = boxCol{b}; end
    if strcmp(st, 'n'), fw = 'bold'; end
    str = regexprep(strrep(str, '_', '\_'), 'A\\_(x|y)', 'A_$1');          % names with underscores, A_x/A_y subscript
    h = text(ax, x0 + lvl*dx, y(k), str, 'FontSize', fs, 'Color', cl, 'FontWeight', fw, ...
        'VerticalAlignment', 'middle', 'Interpreter', 'tex');
    xRight(k) = sum(h.Extent([1 3]));
end
% Boxes: top line of a box is the bottom line of the previous one; box
% name next to the last line of the box that leaves room for it
yEdge = [yTop, arrayfun(@(b) y(find(box == b, 1, 'last')) + 0.45*dy, 1:3)];
for b = 1:3
    rectangle(ax, 'Position', [xBox(1), yEdge(b), diff(xBox), yEdge(b+1) - yEdge(b)], ...
        'EdgeColor', boxCol{b}, 'LineWidth', 0.5);
    kL = find(box(:) == b & xRight < xLabel - 5, 1, 'last');
    text(ax, xLabel, y(kL) + 0.45*dy - 2, boxName{b}, 'FontSize', fsLabel, 'Color', boxCol{b}, ...
        'VerticalAlignment', 'bottom');
end

if ~isfolder(dirFig)
    mkdir(dirFig);
end
fileName = fullfile(dirFig, sprintf('FlowchartInitCtrlPlant_%d', variant));
print(fig, fileName, '-dpng', '-r150');
print(fig, fileName, '-depsc');

% Check against the profiler run
notCalled = {};
if ~isempty(profInfo)
    called = regexprep({profInfo.FunctionTable.FunctionName}, '>.*$', ''); % file names (no subfunctions)
    called = regexprep(called, '^.*[\\/]', '');
    names = unique(regexprep(entries(:,1), ' .*$', ''), 'stable');
    if variant == 1
        names = names(ismember(names, entries([entries{:,3}] < 3, 1)) & ~contains(names, '9turb'));
    end
    notCalled = names(~ismember(lower(names), lower(called)));
    if isempty(notCalled)
        fprintf('Flowchart %d: all functions called in the profiler run\n', variant);
    else
        fprintf('Flowchart %d: not called in the profiler run: %s\n', variant, strjoin(notCalled, ', '));
    end
end
end
