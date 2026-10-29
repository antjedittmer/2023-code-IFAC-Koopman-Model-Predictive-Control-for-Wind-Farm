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
%   shown in an additional subplot, also if all zero: non-adaptive run),
%   showGreedy, showYaw (default true)
% - dirFig, fileBase (optional): save as <fileBase>_Tex (.fig, .png, .eps)

time = pd.time;
N = size(pd.CT_prime,1);
fs = 12;
lw = 1;
cl = lines(7);

if pd.OptYaw == 0
    modeStr = 'Baseline AIC';
elseif pd.OptYaw == 2 || pd.KoopAIC
    modeStr = 'AIC, opt. yaw';
else
    modeStr = 'WRC-AIC';
end
if isfield(pd,'modeStr'), modeStr = pd.modeStr; end
titleStr = sprintf('%s: TE = %.1f kW, 1000\\cdotAA = %.1f', modeStr, pd.JQ*1e3, 1e3*pd.JR);
showUpdate = isfield(pd,'update'); % adaptive AIC runs, also without update (all zero)
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
% height: 5 subplots 1.6, 4 subplots 1.3 times the standard figure
hFac = 1.6 + 0.3*(nSub-5);
if nSub == 4, hFac = 1.3; end
set(fig, 'position', [posDefault(1),posDefault(2)-150,posDefault(3),posDefault(4)*hFac]);
setFigureSize(fig, [1 hFac]); % export size independent of the screen

ax = gobjects(nSub,1);
% tiled layout: less white space between and around the axes than subplot
tl = tiledlayout(fig, nSub, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
ax(1) = nexttile(tl,1);
stairs(time,pd.Pref(time)/1e6,'k','linewidth',lw); hold on; grid on;
stairs(time,sum(pd.Power(:,time))/1e6,'-.','linewidth',lw,'color',cl(4,:)); % P_WF dash-dotted, 1
farmNames = {'P_{ref}','P_{WF}'};
farmColors = [0 0 0; cl(4,:)];
if showEstP % estimated farm power (pink, dotted, 1.5 as the estimated wind speeds)
    clPest = [0.93 0.35 0.65];
    stairs(time,pd.P_est(time)/1e6,':','linewidth',lw+0.5,'color',clPest);
    farmNames{end+1} = 'P_{WF,est}';
    farmColors(end+1,:) = clPest;
end
if showGreedy
    stairs(time,pd.Pgreedy*ones(size(time))/1e6,'--','linewidth',lw,'color',cl(5,:)); % P_greedy dashed green
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
strLeg = ['Farm: ', colorTitle(farmNames, farmColors), ...
    '; Turbine: ', colorTitle(turbNames, turbColors)];
ylabel('Power P_{WF} (MW)'); % name and symbol as in chapter 4
title({titleStr, strLeg});

ax(2) = nexttile(tl,iSub(2));
stairs(time,pd.uInf(time),'k','linewidth',lw); hold on; grid on;
stairs(time,pd.u1(time),'color',cl(1,:),'linewidth',lw);
ylabel('Wind V (m/s)');

ax(3) = nexttile(tl,iSub(3));
for idx = 1:N
    stairs(time,pd.CT_prime(idx,time),'color',cl(idx,:),'linewidth',lw); hold on;
end
grid on; ylim([-0.1 2.1]);
ylabel('Thrust C_T'' (-)');

if showYaw
    ax(4) = nexttile(tl,iSub(4));
    for idx = 1:N
        stairs(time,pd.Phi(idx,time),'color',cl(idx,:),'linewidth',lw); hold on;
    end
    grid on;
    ylabel(['Yaw \gamma (',char(176),')']); % degree sign as in the qLMPC plots
end

ax(5) = nexttile(tl,iSub(5));
stairs(time,pd.V(1,time),'color',cl(1,:),'linewidth',lw); hold on; grid on;
if showEst
    stairs(time,pd.V_est(1,time),'--','linewidth',lw+0.5,'color',[0,0,0.7]);
end
stairs(time,pd.V(2,time),'color',cl(2,:),'linewidth',lw);
if showEst
    stairs(time,pd.V_est(2,time),'--','linewidth',lw+0.5,'color',[0.7,0,0]);
end
ylabel('Wind V_r (m/s)');
if showUpdate
    ax(6) = nexttile(tl,iSub(6));
    stairs(time,pd.update(time),'k','linewidth',lw); grid on;
    ylabel('Update (-)');
end
xlabel('Time (s)');

ax = ax(isgraphics(ax));
for idx = 1:numel(ax)
    axis(ax(idx),'tight');
    set(ax(idx),'FontSize',fs,'TickLabelInterpreter','tex');
end
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
setFigureFont(fig, fs); % all text 12 pt, bold titles

if nargin > 2
    saveas(fig,fullfile(dirFig,[fileBase,'_Tex']),'fig');
    printFigure(fig, fullfile(dirFig,[fileBase,'_Tex'])); % fixed size (setFigureSize)
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
