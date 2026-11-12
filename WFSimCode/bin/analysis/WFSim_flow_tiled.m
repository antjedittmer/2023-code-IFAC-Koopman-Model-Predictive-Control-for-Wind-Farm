function hfig = WFSim_flow_tiled(Wp, solCell, deltaCell, sigmaYCell, xDistCell, hfig, opts)
%WFSIM_FLOW_TILED Longitudinal wind fields of several WFSim solutions, flow left to right
%   WFSim_flow_tiled(Wp, solCell, deltaCell, sigmaYCell, xDistCell, hfig)
%   plots one tile per WFSim solution in solCell (tiled layout, one below
%   the other) with x horizontal (flow from left to right), one colorbar
%   for all tiles and the Bastankhah Gaussian wake (centerline delta
%   dashed, expansion sigma_y dash-dotted) in bold blue lines with white
%   outline. Default colormap oslo (black - dark blue - light blue -
%   white, brightness increasing evenly, also readable in grey-scale
%   print; upper limit 1.05*uRef as in the original WFSim figures), rotors
%   as black lines, labels WT1 blue and WT2 red. The free
%   stream V_inf is shown in the lower left corner, the yaw angle of WT1
%   in the lower right corner.
%   deltaCell, sigmaYCell, xDistCell: centerline, wake width and
%   downstream distance from the rotor of WT1 (runWFSimAnimationPlots);
%   empty deltaCell: no wake lines.
%   opts (optional struct):
%   - titles: cell with one title per tile instead of the yaw angle
%   - uRef: free stream speed (default Wp.site.u_Inf)
%   - relative: true plots u - V_inf with the inflow speed V_inf of each
%     tile, e.g. for snapshots with different inflow speeds (default false)
%   - cLim: colour limits; the same cLim in several figures gives the same
%     colour for the same speed
%   - cMapStyle: 'oslo' (default; black - dark blue - light blue - white,
%     approximation of Crameri's oslo), 'hot' (black - red - yellow -
%     white), 'cold' (hot with swapped colour channels: black - blue -
%     cyan - white), 'sky' (black - blue-grey - sky - white) or
%     'diverging' (as the MAECOsim figures: blue - white - red with white
%     at the free stream, blue deficit, red speed-up, each over its own
%     range); all but diverging with upper limit 1.05*uRef. Wake lines
%     blue with light grey outline, labels in tight light grey boxes
%   - rotorColor: colour of the rotor lines (default black)

if nargin < 6 || isempty(hfig)
    hfig = figure;
end
if nargin < 7
    opts = struct();
end
fs = 12;
lwWake = 2;                                                                % bold Bastankhah lines
clWake = [0 0 1];
clT = {[0 0 1], [1 0 0], [0 0 0]};                                         % WT1 blue, WT2 red (as in the power plots)
cMapStyle = 'oslo';
if isfield(opts, 'cMapStyle'), cMapStyle = opts.cMapStyle; end
relative = isfield(opts, 'relative') && opts.relative;

Dr  = Wp.turbine.Drotor;
Crx = Wp.turbine.Crx; Cry = Wp.turbine.Cry;
N   = Wp.turbine.N;
uRef = Wp.site.u_Inf;
if isfield(opts, 'uRef'), uRef = opts.uRef; end
x = Wp.mesh.ldxx2(:,1)';
y = Wp.mesh.ldyy(1,:);
idxX = 5:numel(x)-4;                                                       % without the boundary cells
idxY = 5:numel(y)-4;

% Fields: u, or u - V_inf of each tile (inflow boundary)
nT = numel(solCell);
F = cell(nT, 1);
uIn = zeros(nT, 1);
for iT = 1:nT
    uIn(iT) = solCell{iT}.u(1, ceil(end/2));
    F{iT} = solCell{iT}.u(idxX,idxY);
    if relative
        F{iT} = F{iT} - uIn(iT);
    end
end
uWhite = uRef;                                                             % free stream: white
if relative, uWhite = 0; end
fMin = min(cellfun(@(f) min(f, [], 'all'), F));
fMax = max(cellfun(@(f) max(f, [], 'all'), F));

clBox = [0.9 0.9 0.9];                                                     % light grey: line outline and label boxes
if any(strcmp(cMapStyle, {'hot', 'cold', 'oslo'}))
    % hot: MATLAB hot (black - red - yellow - white); cold: hot with the
    % colour channels swapped (black - blue - cyan - white); oslo: black -
    % dark blue - light blue - white (approximation of Crameri's oslo,
    % brightness increasing evenly). Upper limit 1.05*uRef as in the
    % original WFSim animation figures
    cLim = [fMin, uWhite + 0.05*uRef];
    if isfield(opts, 'cLim'), cLim = opts.cLim; end
    if strcmp(cMapStyle, 'oslo')
        cNodes = [0.004 0.004 0.005; 0.067 0.165 0.282; 0.133 0.329 0.553; ...
            0.333 0.525 0.796; 0.643 0.698 0.788; 1 1 1];
    else
        cNodes = hot(256);
        if strcmp(cMapStyle, 'cold')
            cNodes = cNodes(:, [3 2 1]);
        end
    end
    uNodes = linspace(cLim(1), cLim(2), size(cNodes, 1));
elseif strcmp(cMapStyle, 'sky')
    % Black (wake core) - blue-grey - sky (dark to light) - white (free
    % stream): brightness increases with the wind speed; nodes in m/s
    cLim = [fMin, uWhite + 0.05*uRef];
    if isfield(opts, 'cLim'), cLim = opts.cLim; end
    cSky = flipud(sky(256));
    uNodes = uWhite + [cLim(1) - uWhite, -0.5*uRef, -0.38*uRef, -0.05*uRef, cLim(2) - uWhite];
    uNodes(2:end-1) = max(uNodes(2:end-1), cLim(1) + (1:3)*1e-3);       % increasing also for low minima
    cNodes = [0 0 0; 0.33 0.39 0.48; cSky(1,:); cSky(end,:); 1 1 1];
else
    % As MAECOsim: dark blue - light blue - white (free stream) - light
    % red - red (ColorBrewer RdBu); blue over the deficit, red over the
    % speed-up, each half over its own range
    cLim = [min(fMin, uWhite - 0.01), max(fMax, uWhite + 0.01)];
    if isfield(opts, 'cLim'), cLim = opts.cLim; end
    % red range at least a quarter of the blue range: small speed-ups are
    % not shown as strong as the deficit in the wake
    cLim(2) = max(cLim(2), uWhite + 0.25*(uWhite - cLim(1)));
    cBlue = [0.020 0.188 0.380; 0.129 0.400 0.674; 0.263 0.576 0.765; 0.573 0.773 0.871; 0.820 0.898 0.941];
    cRed  = [0.992 0.859 0.780; 0.957 0.647 0.510; 0.839 0.376 0.302; 0.698 0.094 0.169];
    uNodes = [linspace(cLim(1), uWhite, 6), linspace(uWhite, cLim(2), 5)];
    uNodes(7) = [];                                                        % white once at the free stream
    cNodes = [cBlue; 1 1 1; cRed];
end
t = linspace(cLim(1), cLim(2), 256)';
cMap = [interp1(uNodes, cNodes(:,1), t), interp1(uNodes, cNodes(:,2), t), interp1(uNodes, cNodes(:,3), t)];

tl = tiledlayout(hfig, nT, 1, 'TileSpacing', 'tight', 'Padding', 'tight');
ax = gobjects(nT, 1);
hLab = gobjects(0);                                                        % wake labels (boxes at the end)
for iT = 1:nT
    sol = solCell{iT};
    ax(iT) = nexttile(tl);
    % contour levels: 50 below and 30 above the free stream (smooth also
    % for the small range of the speed-up)
    if cLim(1) < uWhite && uWhite < cLim(2)
        levels = unique([linspace(cLim(1), uWhite, 50), linspace(uWhite, cLim(2), 30)]);
    else
        levels = linspace(cLim(1), cLim(2), 80);
    end
    contourf(ax(iT), x(idxX), y(idxY), F{iT}', levels, 'LineColor', 'none');
    hold(ax(iT), 'on');
    colormap(ax(iT), cMap); clim(ax(iT), cLim);

    % Turbines: rotor in the turbine colour (WT1 blue, WT2 red) with white
    % outline (visible in the wake)
    phi = sol.turbInput.phi;
    for kk = 1:N
        dx = 0.5*Dr*sind(phi(kk)); dy = 0.5*Dr*cosd(phi(kk));
        xr = Crx(kk) + [-dx, dx]; yr = Cry(kk) + [-dy, dy];
        plot(ax(iT), xr, yr, 'w', 'LineWidth', 4.5);
        clRotor = [0 0 0];                                                 % black rotors, labels in the turbine colours
        if isfield(opts, 'rotorColor'), clRotor = opts.rotorColor; end
        plot(ax(iT), xr, yr, 'Color', clRotor, 'LineWidth', 2.5);
        text(ax(iT), Crx(kk) - 40, Cry(kk) + 0.5*Dr + 50, sprintf('WT%d', kk), 'FontSize', fs, ...
            'Color', clT{kk}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    end

    % Bastankhah wake of WT1: centerline delta and expansion +/- sigma_y,
    % bold blue on a solid, slightly transparent white line (readable on
    % the blue wake; a dashed white line does not match the blue dashes,
    % as the dash length scales with the line width)
    if ~isempty(deltaCell)
        delta = deltaCell{iT}(:)'; sigmaY = sigmaYCell{iT}(:)';
        xW = Crx(1) + xDistCell{iT}(:)'; yC = Cry(1) - delta;
        yLines = {yC, yC + sigmaY, yC - sigmaY};
        styles = {'--', '-.', '-.'};
        for iL = 1:3
            plot(ax(iT), xW, yLines{iL}, '-', 'Color', [clBox 0.8], 'LineWidth', lwWake + 1.5);
            plot(ax(iT), xW, yLines{iL}, styles{iL}, 'Color', clWake, 'LineWidth', lwWake);
        end
        iS = find(~isnan(sigmaY), 1);
        iM = round(0.12*numel(xW));                                        % delta below the centerline, close to WT1
        % labels in light grey boxes, far enough below the lines that the
        % boxes do not touch them
        hLab(end+1) = text(ax(iT), xW(iM), yC(iM) - 28, '\delta', 'FontSize', 16, 'Color', clWake, ...
            'FontWeight', 'bold', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center'); %#ok<AGROW>
        iSig = iS + round(0.3*(numel(xW) - iS));                           % sigma_y below the lower expansion line
        hLab(end+1) = text(ax(iT), xW(iSig), yC(iSig) - sigmaY(iSig) - 26, '\sigma_y', 'FontSize', fs+2, ...
            'Color', clWake, 'FontWeight', 'bold', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center'); %#ok<AGROW>
    end

    % Free stream at the inflow boundary: arrow bottom left, away from the turbine labels
    text(ax(iT), x(idxX(1)) + 30, y(idxY(1)) + 15, sprintf('V_\\infty = %.0f m/s \\rightarrow', uIn(iT)), ...
        'FontSize', 14, 'VerticalAlignment', 'bottom');                    % subscript inside the axes
    if ~isfield(opts, 'titles')                                            % without titles: yaw angle of WT1, lower right corner
        text(ax(iT), x(idxX(end)) - 30, y(idxY(1)) + 15, sprintf('\\gamma_1 = %d%s', round(phi(1)), char(176)), ...
            'FontSize', 14, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
    end

    axis(ax(iT), 'equal', 'tight');
    if isfield(opts, 'titles')
        title(ax(iT), opts.titles{iT}, 'FontSize', fs, 'FontWeight', 'normal');
    end
    ylabel(ax(iT), 'y (m)');
    set(ax(iT), 'FontSize', fs, 'TickLabelInterpreter', 'tex', 'Layer', 'top', 'Box', 'on');
    if iT < nT
        ax(iT).XTickLabel = {};
    end
end
xlabel(ax(end), 'x (m)');
cb = colorbar(ax(end));
cb.Layout.Tile = 'east';
if relative
    cb.Label.String = 'Long. wind u - V_\infty (m/s)';
else
    cb.Label.String = 'Long. wind u (m/s)';
end
cb.FontSize = fs;

% Tight light grey boxes behind the wake labels: the text background box
% includes the full line height of the font, so the boxes are drawn from
% the text extent (after the final layout) without the empty space above
% and below the glyphs
drawnow;
for iL = 1:numel(hLab)
    e = hLab(iL).Extent;                                                   % [x y w h] in data units
    padX = 0.08*e(4);
    xB = [e(1) - padX, e(1) + e(3) + padX];
    yB = e(2) + e(4)*[0.12, 0.80];
    hP = patch(hLab(iL).Parent, xB([1 2 2 1]), yB([1 1 2 2]), clBox, 'EdgeColor', 'none');
    uistack(hLab(iL), 'top');
    hP.PickableParts = 'none';
end
end
