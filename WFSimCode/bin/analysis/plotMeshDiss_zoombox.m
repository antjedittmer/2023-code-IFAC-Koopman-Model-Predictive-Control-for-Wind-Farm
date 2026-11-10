function gridmesh_zoombox(T2,useLatex)

if ~nargin
    T2 = 9;
end

% --- LATEX SWITCH ---
if nargin <2
useLatex = false; % Set to true for LaTeX font, false for standard MATLAB font
% Note to self: D:\ditt_aj\Projekte\260205_WFadaptiveKMPC\adaptivekoopmanmpcwindfarm\WFSimCode
end

if useLatex
    interpStr = 'latex';
    Vstr = '$V_{\infty}$';
    fs = 12;
else
    interpStr = 'tex';
    Vstr = 'V_{\infty}';
end


dirFigureOutput = 'dirFigWFSimLayout';

if exist(dirFigureOutput,'dir')~=7
    mkdir(dirFigureOutput);
end


if T2 == 2
    filename =  'WorkspaceMesh.mat';
    fsT = 12;
    zoomT1 = 0.75;
    j0 = 1;
    my = 0.25; % Vertical offset for main labels
    myZ = 0.4;
else
    filename = 'WorkspaceMesh9.mat';
    fsT = 12;
    zoomT1 = 0.8;
    j0 = 2;
    my = 0.45; % Fixed the syntax error here
    myZ = 0.3;
end


load(filename); %#ok<LOAD>

figure;

% --- Colors and Mesh Data ---
cl = [unique(lines,'rows','stable'); [0,0,0]; 0.5 *[1,1,1]];
bl = [0,0,0.5];
gr = 0.7 *[1,1,1];
Z1 = -2*ones(size(ldxx));
Z2 = -1*ones(size(ldxx2));

% --- Main Plot ---
mesh(ldyy,ldxx,Z1,'FaceAlpha',0,'EdgeColor',bl,'LineWidth',0.2,'DisplayName','Primary mesh (wind speeds)')
hold on;
mesh(ldyy2,ldxx2,Z2,'FaceAlpha',0,'EdgeColor',gr,'LineWidth',0.5,'DisplayName','Secondary mesh (pressure)')

set(findall(gcf,'-property','FontSize'),'FontSize',12)

% --- Loop: Plot Main Turbines & Labels ---
for j = 1:length(Wp.turbine.Cry)
    % Used 'my' variable here so your T2 settings apply

    text(Wp.turbine.Cry(j) + Wp.turbine.Cry(1)*.02, ...
        Wp.turbine.Crx(j) + Wp.turbine.Crx(1) * my, ...
        ['WT',num2str(j)], ...
        'Fontsize', fsT, ...
        'BackgroundColor', 0.9*ones(3,1), ...
        'Margin', 0.01);
    hold on;
    plot([Wp.turbine.Cry(j)-Wp.turbine.Drotor/2; Wp.turbine.Cry(j)+Wp.turbine.Drotor/2],...
        [Wp.turbine.Crx(j);Wp.turbine.Crx(j)],...
        'LineWidth',4.0,'DisplayName',['Turbine WT' num2str(j)],'color', cl(j,:))
end

axis equal; xlim([-0.1*Wp.mesh.Ly 1.2*Wp.mesh.Ly]); ylim([-0.1*Wp.mesh.Lx 1.1*Wp.mesh.Lx]);

% Applied interpStr to labels
xlabel('Lateral direction y (m)','FontSize',13, 'Interpreter', interpStr);
ylabel('Longitudinal direction x (m)','FontSize',13, 'Interpreter', interpStr);

view(0,90);
% Applied interpStr to legend
legend('-dynamicLegend', 'Interpreter', interpStr);

view(0,90); legend('-dynamicLegend');
legend('boxoff');

ax = gca; ax.GridColor = [1 1 1]; ax.MinorGridColor = [1 1 1];

annotation(gcf,'arrow',[0.184 0.184], [0.125 0.179]);


if T2 == 2, txV = 0.19; else txV = 0.23; end
text(txV, 0.045, 'V_{\infty}', 'Units','normalized','FontSize',13,'Color','k');

mainAx = gca; % keep a handle to the main axes for the zoom-box callouts below

if T2 ~= 2
pos0 = mainAx.Position;
mainAx.Position = [0.12 pos0(2:4)];
end

% =========================================================
% LOOP FOR ZOOM INSETS
% =========================================================

offset_x = -0.01;
offset_y = 0.0;
insetPositions = [0.55+offset_x 0.08+offset_y  0.425 0.25;
    0.55+offset_x 0.425+offset_y 0.425 0.25];

numZooms = 1;
if T2 == 2
    numZooms = 2;
end

for idx = 1:numZooms

    zoom_radius = zoomT1 * Wp.turbine.Drotor;
    x_center = Wp.turbine.Crx(idx);
    y_center = Wp.turbine.Cry(idx);

    x_zoom = [x_center - zoom_radius*0.4, x_center + zoom_radius*myZ];
    y_zoom = [y_center - zoom_radius, y_center + zoom_radius];

    % --- NEW LOGIC: TEXT LIMITS ---
    % If T2=2, text limits are the same as the view (shows everything).
    % If T2=9, text limits are restricted to just the Rotor Diameter.
    if T2 == 2
        y_text_min = y_zoom(1);
        y_text_max = y_zoom(2);
    else
        % Strictly under the turbine blades (plus small buffer)
        y_text_min = y_center - Wp.turbine.Drotor/2 - 5;
        y_text_max = y_center + Wp.turbine.Drotor/2 + 5;
    end

    insetAx = axes('Position', insetPositions(idx, :));

    mesh(ldyy, ldxx, Z1, 'FaceAlpha',0, 'EdgeColor',bl,'LineWidth',0.75); hold on;
    mesh(ldyy2, ldxx2, Z2, 'FaceAlpha',0, 'EdgeColor',gr,'LineWidth',1.2);

    set(findall(insetAx,'-property','FontSize'), 'FontSize', 12);

    [ny, nx] = size(ldxx);


    for i = 1:ny
        for j = j0:2:nx
            % 1. Check if point is inside the VIEW (Zoom window)
            inView = (ldxx(i,j) >= x_zoom(1) && ldxx(i,j) <= x_zoom(2)) && ...
                (ldyy(i,j) >= y_zoom(1) && ldyy(i,j) <= y_zoom(2));

            % 2. Check if point is inside the TEXT LIMIT (Tighter for T2=9)
            inTextRegion = (ldyy(i,j) >= y_text_min && ldyy(i,j) <= y_text_max);

            if inView && inTextRegion && i~= 18 && i ~= 17

                text(ldyy(i,j), ldxx(i,j)+ 6, 0, sprintf('(%d,%d)', i, j), ...
                    'FontSize', 11, 'Color', bl, 'HorizontalAlignment', 'center', ...
                    'BackgroundColor', 'w', 'Margin', 0.01);
            elseif inView && inTextRegion && i~= 18 && i == 17
                text(ldyy(i,j), ldxx(i,j)+ 10, 0, sprintf('(%d,%d)', i, j), ...
                    'FontSize', 11, 'Color', bl, 'HorizontalAlignment', 'center', ...
                    'BackgroundColor', 'w', 'Margin', 0.01);

            end
        end
    end


    plot([Wp.turbine.Cry(idx)-Wp.turbine.Drotor/2; Wp.turbine.Cry(idx)+Wp.turbine.Drotor/2], ...
        [Wp.turbine.Crx(idx); Wp.turbine.Crx(idx)], ...
        'color',cl(idx,:), 'LineWidth',4.0);

    axis equal; xlim(y_zoom); ylim(x_zoom); view(0,90);
    box on
    ax = gca; ax.GridColor = [1 1 1]; ax.MinorGridColor = [1 1 1];

    % =====================================================
    % ZOOM-BOX CALLOUT: rectangle in main plot + connector lines
    % (replaces the removed inset title)
    % =====================================================
    axes(mainAx); % switch context back to the main axes to draw the box

    rectangle('Position', [y_zoom(1), x_zoom(1), diff(y_zoom), diff(x_zoom)], ...
        'EdgeColor','k', 'LineWidth', 0.75, 'LineStyle', '-'); % cl(idx,:)

    % Box corners on the side facing the insets (right side, since
    % insets sit at x >= 0.55 in the figure, right of the main plot)
    boxTR = [y_zoom(2), x_zoom(2)]; % top-right corner, main-axes data coords
    boxBR = [y_zoom(2), x_zoom(1)]; % bottom-right corner

    % Matching corners on the inset (its left side faces the main plot)
    insetXlim = xlim(insetAx);
    insetYlim = ylim(insetAx);

    % Pull the target points slightly inward from the exact top/bottom
    % limits so the lines don't sweep straight through the y-tick labels
    % that sit right at those limits (e.g. '380' at the bottom edge).
    cornerInsetFrac = 0.06; % fraction of the inset's y-range -- TUNE THIS
    yInsetRange = diff(insetYlim);
    insetTL = [insetXlim(1), insetYlim(2) - cornerInsetFrac*yInsetRange]; % near top-left
    insetBL = [insetXlim(1), insetYlim(1) + cornerInsetFrac*yInsetRange]; % near bottom-left

    figBoxTR   = dsxy2figxy(mainAx,  boxTR(1),   boxTR(2));
    figBoxBR   = dsxy2figxy(mainAx,  boxBR(1),   boxBR(2));
    figInsetTL = dsxy2figxy(insetAx, insetTL(1), insetTL(2));
    figInsetBL = dsxy2figxy(insetAx, insetBL(1), insetBL(2));

    stopFrac = 0.7; % fraction of the way from start to original end point
    startFrac = 0;
    startFrac1 = 0;
    if idx == 2
        startFrac1 = 0.2;  % skip the first 10% of the line, leaving a gap at the box
    end


    lineStartTR = figBoxTR + startFrac1 * (figInsetTL - figBoxTR);
    lineEndTR   = figBoxTR + (figInsetTL - figBoxTR);

    lineStartBR = figBoxBR + startFrac * (figInsetBL - figBoxBR);
    lineEndBR   = figBoxBR + stopFrac  * (figInsetBL - figBoxBR);


    % lineEndTR = figBoxTR + (figInsetTL - figBoxTR);
    % lineEndBR = figBoxBR + stopFrac * (figInsetBL - figBoxBR);

    annotation('line', [lineStartTR(1) lineEndTR(1)], [lineStartTR(2) lineEndTR(2)], ...
        'Color', 'k', 'LineStyle', '--', 'LineWidth', 0.75);
    annotation('line', [lineStartBR(1) lineEndBR(1)], [lineStartBR(2) lineEndBR(2)], ...
        'Color', 'k', 'LineStyle', '--', 'LineWidth', 0.75);
    axes(insetAx); % restore current axes so any following code targets the inset
end

if T2 ~= 2, name = 'gridmesh9'; else, name = 'gridmesh2'; end



filenameFigure = fullfile(dirFigureOutput,[name,'_',interpStr]);
print(gcf,filenameFigure,'-dpng');
print(gcf,filenameFigure,'-depsc');
end


% =========================================================
% Local helper function
% =========================================================
function figPos = dsxy2figxy(ax, x, y)
% Convert a data-space coordinate (x,y) in axes ax into figure-normalized
% coordinates, so it can be used with annotation(), which only accepts
% figure-normalized positions.
%
% This version explicitly accounts for the "letterboxing" that can occur
% with 'axis equal': the axes' Position rectangle may be larger than the
% actual rendered plot box, with blank margin centered on one axis. We
% compute the true plotted pixel box (assuming equal data-to-pixel scale
% in x and y, matching 'axis equal') rather than assuming Position maps
% 1:1 onto xlim/ylim. If there is no letterboxing, this reduces to the
% same result as a direct linear mapping, so it's safe to use generally.

origAxUnits = ax.Units;
ax.Units = 'pixels';
axPixPos = ax.Position; % [left bottom width height] in pixels, relative to figure
ax.Units = origAxUnits;

xl = xlim(ax);
yl = ylim(ax);
xRange = diff(xl);
yRange = diff(yl);

axPixWidth  = axPixPos(3);
axPixHeight = axPixPos(4);

% Limiting scale (pixels per data unit) -- the smaller of the two
% determines the actual plotted box size under equal data aspect.
scaleX = axPixWidth  / xRange;
scaleY = axPixHeight / yRange;
scale  = min(scaleX, scaleY);

plotWidthPix  = xRange * scale;
plotHeightPix = yRange * scale;

marginXPix = (axPixWidth  - plotWidthPix)  / 2; % centered letterbox margin
marginYPix = (axPixHeight - plotHeightPix) / 2;

% Position of (x,y) within the actual plotted box, in pixels
xPixInPlot = (x - xl(1)) * scale;
yPixInPlot = (y - yl(1)) * scale;

% Absolute pixel position within the figure
xPixAbs = axPixPos(1) + marginXPix + xPixInPlot;
yPixAbs = axPixPos(2) + marginYPix + yPixInPlot;

% Convert to figure-normalized coordinates
fig = ancestor(ax, 'figure');
origFigUnits = fig.Units;
fig.Units = 'pixels';
figPixPos = fig.Position; % [left bottom width height]
fig.Units = origFigUnits;

figPos(1) = xPixAbs / figPixPos(3);
figPos(2) = yPixAbs / figPixPos(4);
end
