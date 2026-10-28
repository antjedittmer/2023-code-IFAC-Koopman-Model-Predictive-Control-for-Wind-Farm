function fig = plotWFSimTimeseriesTex(pd, dirFig, fileBase)
% plotWFSimTimeseriesTex closed-loop time series as plotWFSimTimeseriesLatex2
% with the TeX interpreter, the operating mode in the title and the time in
% s on the x axis (thesis figures). Instead of legends, the title of the
% first subplot names the signals in their line colors (Farm signals: P_ref, P_WF,
% P_WF,est (WRC-AIC with the input-only model, pd.P_est), P_greedy, V_inf;
% Turbine signals: V_r1, V_r2 and their estimates (AIC); the turbine
% colors also apply to C_T' and yaw), so no data are covered; line styles
% go into the caption. The data pd are saved per run by
% WFSim_demo (<figure name>_plotData.mat), so the figure can be recreated
% without simulation: tmp = load(file); plotWFSimTimeseriesTex(tmp.pd, dirFig, fileBase)
%
% Inputs
% - pd: struct with time, Power, CT_prime, Phi, V, V_est, u1, uInf, Pref,
%   Pgreedy, JQ (MW), JR and the settings controller, OptYaw, KoopAIC, KoopFeat;
%   optional: modeStr (title), update (update flags of the adaptive AIC,
%   shown in an additional subplot), showGreedy, showYaw (default true)
% - dirFig, fileBase (optional): save as <fileBase>_Tex (.fig, .png, .eps)

time = pd.time;
N = size(pd.CT_prime,1);
fs = 12;
lw = 1;
cl = lines(7);

if pd.OptYaw == 0
    modeStr = 'Baseline AIC';
elseif pd.OptYaw == 2 || pd.KoopAIC
    modeStr = 'AIC with static yaw';
else
    modeStr = 'WRC-AIC';
end
if isfield(pd,'modeStr'), modeStr = pd.modeStr; end
titleStr = sprintf('%s: TE = %.1f kW, AA = %.1e', modeStr, pd.JQ*1e3, pd.JR);
showUpdate = isfield(pd,'update') && any(pd.update);
showGreedy = ~isfield(pd,'showGreedy') || pd.showGreedy;
showYaw = ~isfield(pd,'showYaw') || pd.showYaw;
nSub = 4 + showYaw + showUpdate;
iSub = cumsum([1, 1, 1, showYaw, 1, showUpdate]); % subplot index of P, V, CT, yaw, Ur, update
% Estimates: Vr1, Vr2 for the AIC (Ur model, or input-only power model
% without wind), farm power for the WRC-AIC with the input-only power model
isFeat = isfield(pd,'KoopFeat') && pd.KoopFeat > 0 && isfield(pd,'P_est');
showEstP = isFeat && pd.OptYaw == 1;
showEst = pd.KoopAIC == 1 || (isFeat && pd.OptYaw ~= 1);
% VAF of the shown estimates over the plotted time in the first title line
vaf = @(y,yh) 100*(1 - var(y - yh)/var(y));
if showEstP
    titleStr = sprintf('%s, VAF_{P} = %.1f %%', titleStr, vaf(sum(pd.Power(:,time)), pd.P_est(time)));
elseif showEst
    titleStr = sprintf('%s, VAF_{V} = %.1f/%.1f %%', titleStr, ...
        vaf(pd.V(1,time), pd.V_est(1,time)), vaf(pd.V(2,time), pd.V_est(2,time)));
end

fig = figure;
posDefault = get(fig, 'position');
set(fig, 'position', [posDefault(1),posDefault(2)-150,posDefault(3),posDefault(4)*(1.6 + 0.3*(nSub-5))]);

ax = gobjects(nSub,1);
ax(1) = subplot(nSub,1,1);
stairs(time,pd.Pref(time)/1e6,'k','linewidth',lw); hold on; grid on;
stairs(time,sum(pd.Power(:,time))/1e6,'-.','linewidth',lw+0.2,'color',cl(4,:));
farmNames = {'P_{ref}','P_{WF}'};
farmColors = [0 0 0; cl(4,:)];
if showEstP % estimated farm power (dark purple)
    clPest = [0.30 0.05 0.35];
    stairs(time,pd.P_est(time)/1e6,':','linewidth',lw+0.4,'color',clPest);
    farmNames{end+1} = 'P_{WF,est}';
    farmColors(end+1,:) = clPest;
end
if showGreedy
    stairs(time,pd.Pgreedy*ones(size(time))/1e6,':','linewidth',lw+0.2,'color',cl(5,:));
    farmNames{end+1} = 'P_{greedy}';
    farmColors(end+1,:) = cl(5,:);
end
farmNames{end+1} = 'V_\infty';
farmColors(end+1,:) = [0 0 0];
turbNames = {'V_{r1}','V_{r2}'};
turbColors = cl(1:2,:);
if showEst % estimated effective wind speeds (dark blue, dark red)
    turbNames = [turbNames, {'V_{r1,est}','V_{r2,est}'}];
    turbColors = [turbColors; 0 0 0.7; 0.7 0 0];
end
strLeg = ['Farm sig.: ', colorTitle(farmNames, farmColors), ...
    '; Turbine sig.: ', colorTitle(turbNames, turbColors)];
ylabel('P_{WF} (MW)');
title({titleStr, strLeg},'FontWeight','normal');

ax(2) = subplot(nSub,1,iSub(2));
stairs(time,pd.uInf(time),'k','linewidth',lw); hold on; grid on;
stairs(time,pd.u1(time),'color',cl(1,:),'linewidth',lw);
ylabel('V (m/s)');

ax(3) = subplot(nSub,1,iSub(3));
for idx = 1:N
    stairs(time,pd.CT_prime(idx,time),'color',cl(idx,:),'linewidth',lw); hold on;
end
grid on; ylim([-0.1 2.1]);
ylabel('C_T'' (-)');

if showYaw
    ax(4) = subplot(nSub,1,iSub(4));
    for idx = 1:N
        stairs(time,pd.Phi(idx,time),'color',cl(idx,:),'linewidth',lw); hold on;
    end
    grid on;
    ylabel('\gamma (deg)');
end

ax(5) = subplot(nSub,1,iSub(5));
stairs(time,pd.V(1,time),'color',cl(1,:),'linewidth',lw); hold on; grid on;
if showEst
    stairs(time,pd.V_est(1,time),'--','linewidth',lw,'color',[0,0,0.7]);
end
stairs(time,pd.V(2,time),'color',cl(2,:),'linewidth',lw);
if showEst
    stairs(time,pd.V_est(2,time),'--','linewidth',lw,'color',[0.7,0,0]);
end
ylabel('V_r (m/s)');
if showUpdate
    ax(6) = subplot(nSub,1,iSub(6));
    stairs(time,pd.update(time),'k','linewidth',lw); grid on;
    ylabel('Update (-)');
end
xlabel('Time (s)');

ax = ax(isgraphics(ax));
for idx = 1:numel(ax)
    axis(ax(idx),'tight');
    set(ax(idx),'FontSize',fs,'TickLabelInterpreter','tex');
end
ax(1).TitleFontSizeMultiplier = 11/fs; % title 11 pt (default 1.1*fs), fits the figure width
axP = ax(1); axCT = ax(3); axUr = ax(end - showUpdate);
ylim(axCT,[-0.1 2.1]);
tP = time(time >= 60); % without the start-up
pAll = [pd.Pref(tP)'; sum(pd.Power(:,tP))]/1e6;
if showGreedy, pAll(end+1,:) = pd.Pgreedy/1e6; end
if showEstP, pAll(end+1,:) = pd.P_est(tP)/1e6; end
pAll = pAll(:,~any(isnan(pAll),1));
ylim(axP,[min(pAll(:))-0.05, max(pAll(:))+0.05]);
set(axP,'YTick',0:0.4:10);
if showYaw && max(abs(pd.Phi(:,time)),[],'all') < 1 % yaw (numerically) zero
    ylim(ax(4),[-1 1]);
end
vAll = [pd.uInf(time), pd.u1(time)];
ylim(ax(2),[min(vAll)-0.3, max(vAll)+0.3]);
yl = ylim(axUr); ylim(axUr,yl + [-0.3 0.3]);
if showUpdate, ylim(ax(end),[-0.1 1.1]); end

if nargin > 2
    saveas(fig,fullfile(dirFig,[fileBase,'_Tex']),'fig');
    printFixedWidth(fig, fullfile(dirFig,[fileBase,'_Tex']));
end
end

function str = colorTitle(names, colors)
% colorTitle signal names in their line colors (TeX) as title instead of a legend
parts = cell(size(names));
for idx = 1:numel(names)
    parts{idx} = sprintf('{\\color[rgb]{%.3f,%.3f,%.3f}%s}', colors(idx,:), names{idx});
end
str = strjoin(parts, ', ');
end

function printFixedWidth(h, fName)
% printFixedWidth saves figure h as vector .eps and 300 dpi .png with the
% width of the standard MATLAB figure (560 px = 14.82 cm), the height from
% the aspect ratio of the window: the output is laid out for this paper
% size, not for the window, and the eps bounding box is the full paper
% size (-loose), so that all thesis figures have the same width and font
% size in LaTeX (e.g. two side by side with width=0.495	extwidth)
pos = get(h, 'Position');
sizeCm = 14.82*[1, pos(4)/pos(3)];
set(h, 'PaperUnits', 'centimeters', 'PaperPositionMode', 'manual', ...
    'PaperSize', sizeCm, 'PaperPosition', [0 0 sizeCm]);
print(h, fName, '-depsc', '-painters', '-loose');
print(h, fName, '-dpng', '-r300');
end
