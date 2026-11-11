function hfig = WFSim_flow_tiled(Wp, solCell, deltaCell, sigmaYCell, xDistCell, hfig)
%WFSIM_FLOW_TILED Longitudinal wind fields of several yaw angles, flow left to right
%   WFSim_flow_tiled(Wp, solCell, deltaCell, sigmaYCell, xDistCell, hfig)
%   plots one tile per WFSim solution in solCell (tiled layout, one below
%   the other) with x horizontal (flow from left to right), the colormap
%   sky extended by black (low wind speed) and white (free stream), one
%   colorbar for all tiles and the Bastankhah Gaussian wake (centerline
%   delta dashed, expansion sigma_y dash-dotted) in bold blue lines.
%   deltaCell, sigmaYCell, xDistCell: centerline, wake width and
%   downstream distance from the rotor of WT1 (runWFSimAnimationPlots).

if nargin < 6 || isempty(hfig)
    hfig = figure;
end
fs = 12;
lwWake = 2;                                                                % bold Bastankhah lines
clWake = [0 0 1];
clT = {[0 0 1], [1 0 0], [0 0 0]};                                         % WT1 blue, WT2 red (as in the power plots)

Dr  = Wp.turbine.Drotor;
Crx = Wp.turbine.Crx; Cry = Wp.turbine.Cry;
N   = Wp.turbine.N;
uInf = Wp.site.u_Inf;
x = Wp.mesh.ldxx2(:,1)';
y = Wp.mesh.ldyy(1,:);
idxX = 5:numel(x)-4;                                                       % without the boundary cells
idxY = 5:numel(y)-4;

% Colormap: black (wake core) - dark red - sky (dark to light) - white
% (free stream): brightness increases with the wind speed, the hue changes
% inside the wake (more detail than a single blue gradient), and the blue
% wake lines stay readable on dark red. Nodes in m/s, so that the same
% speed has the same colour in all figures.
cLim = [min(cellfun(@(s) min(s.u(idxX,idxY), [], 'all'), solCell)), uInf*1.05];
cSky = flipud(sky(256));
uNodes = [cLim(1), 0.5*uInf, 0.62*uInf, 0.95*uInf, cLim(2)];
uNodes(2:end-1) = max(uNodes(2:end-1), cLim(1) + (1:3)*1e-3);           % increasing also for low minima
cNodes = [0 0 0; 0.55 0 0; cSky(1,:); cSky(end,:); 1 1 1];
t = linspace(cLim(1), cLim(2), 256)';
cMap = [interp1(uNodes, cNodes(:,1), t), interp1(uNodes, cNodes(:,2), t), interp1(uNodes, cNodes(:,3), t)];

nT = numel(solCell);
tl = tiledlayout(hfig, nT, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = gobjects(nT, 1);
for iT = 1:nT
    sol = solCell{iT};
    ax(iT) = nexttile(tl);
    contourf(ax(iT), x(idxX), y(idxY), sol.u(idxX,idxY)', linspace(cLim(1), cLim(2), 60), ...
        'LineColor', 'none');
    hold(ax(iT), 'on');
    colormap(ax(iT), cMap); clim(ax(iT), cLim);

    % Turbines: rotor in the turbine colour (WT1 blue, WT2 red) with white
    % outline (visible in the wake)
    phi = sol.turbInput.phi;
    for kk = 1:N
        dx = 0.5*Dr*sind(phi(kk)); dy = 0.5*Dr*cosd(phi(kk));
        xr = Crx(kk) + [-dx, dx]; yr = Cry(kk) + [-dy, dy];
        plot(ax(iT), xr, yr, 'w', 'LineWidth', 4.5);
        plot(ax(iT), xr, yr, 'Color', clT{kk}, 'LineWidth', 2.5);
        text(ax(iT), Crx(kk) - 40, Cry(kk) + 0.5*Dr + 50, sprintf('WT%d', kk), 'FontSize', fs, ...
            'Color', clT{kk}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    end

    % Bastankhah wake of WT1: centerline delta and expansion +/- sigma_y
    delta = deltaCell{iT}(:)'; sigmaY = sigmaYCell{iT}(:)';
    xW = Crx(1) + xDistCell{iT}(:)'; yC = Cry(1) - delta;
    plot(ax(iT), xW, yC, '--', 'Color', clWake, 'LineWidth', lwWake);
    plot(ax(iT), xW, yC + sigmaY, '-.', 'Color', clWake, 'LineWidth', lwWake);
    plot(ax(iT), xW, yC - sigmaY, '-.', 'Color', clWake, 'LineWidth', lwWake);
    iS = find(~isnan(sigmaY), 1);
    iM = round(0.3*numel(xW));                                             % delta below the centerline, before the sigma_y lines
    text(ax(iT), xW(iM), yC(iM) - 10, '\delta', 'FontSize', fs+2, 'Color', clWake, ...
        'FontWeight', 'bold', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center');
    text(ax(iT), xW(iS) - 15, yC(iS) + sigmaY(iS), '\sigma_y', 'FontSize', fs+2, ...
        'Color', clWake, 'FontWeight', 'bold', 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right');

    % Free stream: arrow at the inflow (bottom left, away from the turbine labels)
    text(ax(iT), x(idxX(1)) + 30, y(idxY(1)) + 50, sprintf('V_\\infty = %.0f m/s \\rightarrow', uInf), ...
        'FontSize', 14, 'VerticalAlignment', 'middle');

    axis(ax(iT), 'equal', 'tight');
    title(ax(iT), sprintf('\\gamma_1 = %d%s', round(phi(1)), char(176)), 'FontSize', fs, ...
        'FontWeight', 'normal');
    ylabel(ax(iT), 'y (m)');
    set(ax(iT), 'FontSize', fs, 'TickLabelInterpreter', 'tex', 'Layer', 'top', 'Box', 'on');
    if iT < nT
        ax(iT).XTickLabel = {};
    end
end
xlabel(ax(end), 'x (m)');
cb = colorbar(ax(end));
cb.Layout.Tile = 'east';
cb.Label.String = 'Long. wind u (m/s)';
cb.FontSize = fs;
end
