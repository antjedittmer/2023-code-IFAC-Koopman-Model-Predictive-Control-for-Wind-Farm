function fig = plotWFSimTimeseriesTex(pd, dirFig, fileBase)
% plotWFSimTimeseriesTex closed-loop time series as plotWFSimTimeseriesLatex2
% with the TeX interpreter, the operating mode in the title and the time in
% s on the x axis (thesis figures). The data pd are saved per run by
% WFSim_demo (<figure name>_plotData.mat), so the figure can be recreated
% without simulation: tmp = load(file); plotWFSimTimeseriesTex(tmp.pd, dirFig, fileBase)
%
% Inputs
% - pd: struct with time, Power, CT_prime, Phi, V, V_est, u1, uInf, Pref,
%   Pgreedy, JQ (MW), JR and the settings controller, OptYaw, KoopAIC, KoopFeat
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
titleStr = sprintf('%s: TE = %.1f kW, AA = %.2e', modeStr, pd.JQ*1e3, pd.JR);
showEst = pd.KoopAIC == 1; % estimated Ur used by the controller

fig = figure;
posDefault = get(fig, 'position');
set(fig, 'position', [posDefault(1),posDefault(2)-75,posDefault(3),posDefault(4)*1.35]);

ax = gobjects(5,1);
ax(1) = subplot(5,1,1);
stairs(time,pd.Pref(time)/1e6,'k','linewidth',lw); hold on; grid on;
stairs(time,sum(pd.Power(:,time))/1e6,'--','linewidth',lw+0.2,'color',cl(4,:));
stairs(time,pd.Pgreedy*ones(size(time))/1e6,':','linewidth',lw+0.2,'color',cl(5,:));
legend('P_{ref}','P_T','P_{greedy}','Orientation','horizontal','Location','South');
ylabel('\Sigma P_i (MW)');
title(titleStr,'FontWeight','normal');

ax(2) = subplot(5,1,2);
stairs(time,pd.uInf(time),'k','linewidth',lw); hold on; grid on;
stairs(time,pd.u1(time),'color',cl(1,:),'linewidth',lw);
legend('V_\infty','V_1','Location','NorthWest','Orientation','horizontal');
ylabel('V (m/s)');

ax(3) = subplot(5,1,3);
for idx = 1:N
    stairs(time,pd.CT_prime(idx,time),'color',cl(idx,:),'linewidth',lw); hold on;
end
grid on; ylim([-0.1 2.1]);
legend(compose('WT%d',1:N),'Orientation','horizontal','Location','South');
ylabel('C_T'' (-)');

ax(4) = subplot(5,1,4);
for idx = 1:N
    stairs(time,pd.Phi(idx,time),'color',cl(idx,:),'linewidth',lw); hold on;
end
grid on;
legend(compose('WT%d',1:N),'Orientation','horizontal','Location','SouthWest');
ylabel('\gamma (deg)');

ax(5) = subplot(5,1,5);
stairs(time,pd.V(1,time),'color',cl(1,:),'linewidth',lw); hold on; grid on;
if showEst
    stairs(time,pd.V_est(1,time),'--','linewidth',lw,'color',[0,0,0.7]);
end
stairs(time,pd.V(2,time),'color',cl(2,:),'linewidth',lw);
if showEst
    stairs(time,pd.V_est(2,time),'--','linewidth',lw,'color',[0.7,0,0]);
    legend({'WT1','WT1 est.','WT2','WT2 est.'},'Orientation','horizontal','Location','NorthEast');
else
    legend({'WT1','WT2'},'Orientation','horizontal','Location','NorthEast');
end
ylabel('U_r (m/s)');
xlabel('Time (s)');

for idx = 1:5
    axis(ax(idx),'tight');
    set(ax(idx),'FontSize',fs,'TickLabelInterpreter','tex');
    legend(ax(idx),'boxoff');
end
ylim(ax(3),[-0.1 2.1]);
pAll = [pd.Pref(time)'; sum(pd.Power(:,time))]/1e6;
ylim(ax(1),[min(pAll(:))-0.05, max(pAll(:))+0.05]);
set(ax(1),'YTick',0:0.4:10);
if max(abs(pd.Phi(:,time)),[],'all') < 1 % yaw (numerically) zero
    ylim(ax(4),[-1 1]);
end
vAll = [pd.uInf(time), pd.u1(time)];
ylim(ax(2),[min(vAll)-0.3, max(vAll)+0.3]);
yl = ylim(ax(5)); ylim(ax(5),[yl(1), 1.3*yl(2)]);

if nargin > 2
    saveas(fig,fullfile(dirFig,[fileBase,'_Tex']),'fig');
    print(fig,fullfile(dirFig,[fileBase,'_Tex']),'-dpng');
    print(fig,fullfile(dirFig,[fileBase,'_Tex']),'-depsc');
end
end
