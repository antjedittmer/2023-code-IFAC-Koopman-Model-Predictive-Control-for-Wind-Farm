function fig1 = plotEDMDinputsEffWindUpdate(ysim_val,Inputs_val,Deterministic_val,FITje_val,dirFig,Vin,nx,strVal,Pvec,displayOff,Vinf,tVec,interp,marks)
% plotEDMDinputsEffWind plots the simulated outputs of the eDMD estimations
% and prints them into the directory given in dirFig.
% Optional tVec: time (s) of the samples, x axis 'Time (s)' (default: sample index)
% Optional interp: 'latex' (default) or 'tex' (labels converted: no $, units
% in (), no legends: the first title names the measured signals WT1, WT2,
% V_inf and the estimated signals V_r1,est, V_r2,est in their colors)
% Optional marks: table with Time_s and Color (plotAICnoV1Sigma), colored
% lines at these times in the update subplot (tVec required)

if nargin < 6
    Vin = 0;
end
if nargin < 7
    nx = 0;
end

if nargin < 8
    strVal = 'Val.';
end

if nargin < 9
    displayOff = 0;
end

if nargin < 11
    Vinf = [];
end

if nargin < 12
    tVec = [];
end

if nargin < 13
    interp = 'latex';
end

if nargin < 14
    marks = [];
end

if nx > 0
    strK =  sprintf(', $n_K$ = %d', nx);
else
    strK = '';
end

fig1 = figure;
posDefault = get(gcf, 'position');
set(gcf, 'position', [posDefault(1:2),posDefault(3),posDefault(4)*1.2]);
cl = colormap('lines');
lw = 1;
fs = 10;
if isempty(tVec)
    kplot = 60: length(ysim_val);
    xk = 1:length(kplot);
else % time axis: the data start after the start-up (t0), show all
    kplot = 1: length(ysim_val);
    xk = tVec(kplot);
end

noSub = size(Inputs_val,1) + 2;
subplot(noSub,1,1)
plot(xk,Inputs_val(1,kplot),'LineWidth',lw,'Color', cl(1,:)); %1
hold on;
plot(xk,Inputs_val(2,kplot),'LineWidth',lw,'Color', cl(2,:)) %3
grid on; axis tight; ylim([0 2+0.2])
ylabel('$C_T$ [-]','interpreter','latex','Fontsize',fs);

if displayOff
    title(sprintf('%s: VAF Ur1: %d, Ur2: %d',...
    strVal,round(FITje_val)),'interpreter','latex','Fontsize',fs);
else
title(sprintf('%s Set%s:  VAF Ur1: %d, Ur2: %d',...
    strVal,strK,round(FITje_val)),'interpreter','latex','Fontsize',fs);
end


strLeg = sprintf('WT%d,',1:2);
legend(regexp(strLeg(1:end-1), ',', 'split'),'interpreter','latex','Fontsize',fs, 'Orientation','horizontal', 'location','southeast');
legend('boxoff')
set(gca,'fontsize', fs)

if noSub > 3
    subplot(noSub,1,2)
    if isempty(Vinf)
    plot(xk,Inputs_val(3,kplot),'LineWidth',lw,'Color', cl(1,:)); %1
    else
         plot(xk,Vinf(kplot),'k','LineWidth',lw); hold on; plot(xk,Inputs_val(3,kplot),'LineWidth',lw,'Color', cl(1,:));
    end
    grid on; axis tight
    if Vin && ~isempty(Vinf)
        ylabel('$V$ [m/s]','interpreter','latex','Fontsize',fs);
        legend('$V_{\inf}$','$V_1$','interpreter','latex','Fontsize',fs,'orientation', 'horizontal','Location','NorthWest');
        legend('boxoff')
    elseif Vin
        ylabel('Input $V_1$ [m/s]','interpreter','latex','Fontsize',fs);
    else
        ylabel('Input $\Phi$ [deg]','interpreter','latex','Fontsize',fs);
    end
    set(gca,'fontsize', fs)
end

subplot(noSub,1,noSub-2)
plot(xk,Deterministic_val(1,kplot),'LineWidth',lw,'Color', cl(1,:)); %1
hold on;
plot(xk,ysim_val(1,kplot),'Color',[0,0,0.6],'LineWidth',lw,'LineStyle','--') %3
grid on; axis tight;
pos = axis; axis([pos(1:3), pos(4)+0.3]);
%xlabel('k [-]','interpreter','latex')
ylabel('$U_{r1}$ [m/s]','interpreter','latex')
lgh = legend({'Real','Est.'},'Location','southeast','interpreter','latex','fontsize',fs,'orientation','horizontal')
%  0.6155    0.4802    0.2702    0.0397
lgh.Position = [lgh.Position(1) + 0.02,lgh.Position(2:4)];
legend('boxoff')
set(gca,'fontsize', fs)

subplot(noSub,1,noSub-1)
plot(xk,Deterministic_val(2,kplot),'LineWidth',lw,'Color', cl(2,:)); %1
hold on;
plot(xk,ysim_val(2,kplot),'Color',[0.6, 0, 0.0],'LineWidth',lw,'LineStyle','--') %3
grid on; axis tight
pos = axis; axis([pos(1:3), pos(4)+0.3]);
ylabel('$U_{r2}$ [m/s]','interpreter','latex','fontsize', fs)
lgh = legend({'Real','Est.'},'Location','southeast','interpreter','latex','fontsize',fs,'orientation','horizontal')
lgh.Position = [0.2595    0.3585    0.2702    0.0397];
legend('boxoff')
set(gca,'fontsize', fs)

subplot(noSub,1,noSub)
if isempty(tVec)
    plot(min((Pvec-1),1),'k','LineWidth',lw); %1
else
    plot(xk,min((Pvec(kplot)-1),1),'k','LineWidth',lw);
end
if ~isempty(marks) % updates after which the model is shown elsewhere
    hold on;
    for idx = find(marks.Time_s(:)' >= xk(1) & marks.Time_s(:)' <= xk(end))
        % first sample (initial model): shifted to be visible next to the axis
        tMark = max(marks.Time_s(idx), xk(1) + 0.004*(xk(end)-xk(1)));
        plot(tMark*[1 1], [0 1], 'Color', marks.Color(idx,:), 'LineWidth', 2);
    end
end

grid on; axis tight
pos = axis; axis([pos(1:3), pos(4)+0.3]);
if isempty(tVec)
    xlabel('k [-]','interpreter','latex')
else
    xlabel('Time (s)','interpreter','latex')
end
ylabel('$Update$ ','interpreter','latex')
axis tight;

%legend('P','Location','northeast','interpreter','latex','fontsize',fs,'orientation','horizontal')
%legend('boxoff')
set(gca,'fontsize', fs)

if strcmpi(interp,'tex') % convert all LaTeX strings to TeX, color key instead of legends
    for hTxt = findall(fig1,'Type','text')'
        set(hTxt,'String',latex2tex(get(hTxt,'String')),'Interpreter','tex');
    end
    delete(findall(fig1,'Type','legend'));
    axs = findall(fig1,'Type','axes');
    set(axs,'TickLabelInterpreter','tex');
    axTop = axs(end); % first subplot
    axTop.Title.String = {axTop.Title.String, ...
        ['Measured sig.: ', colorKey({'WT1','WT2','V_\infty'}, [cl(1,:); cl(2,:); 0 0 0]), ...
        '; Estimated sig.: ', colorKey({'V_{r1,est}','V_{r2,est}'}, [0 0 0.6; 0.6 0 0])]};
    axTop.Title.FontWeight = 'normal';
    axs(1).YLabel.String = 'Update (-)';
end

%sys_red{1,1}.Xint = xo;
strVal1 = matlab.lang.makeValidName(strVal);
if strcmpi(interp,'tex') % thesis: fixed width
    printFixedWidth(fig1, fullfile(dirFig,['ValP',strVal1]));
else
    print(gcf,fullfile(dirFig,['ValP',strVal1]), '-dpng');
    print(gcf,fullfile(dirFig,['ValP',strVal1]), '-depsc');
end
end

function str = latex2tex(str)
% latex2tex LaTeX label to TeX: no $, \inf -> \infty, units [..] -> (..),
% U_r -> V_r
str = strrep(str,'$','');
str = strrep(str,'\inf}','\infty}');
str = regexprep(str,'\[([^\]]*)\]','($1)');
str = regexprep(str,'\<Ur(\d)','V_{r$1}'); % VAF Ur1 -> V_{r1}
str = regexprep(str,', r(\d):',', V_{r$1}:'); % VAF ..., r2 -> V_{r2}
str = regexprep(str,'\<U(_?\{?r)','V$1'); % effective wind speed V_r as in the thesis
end

function str = colorKey(names, colors)
% colorKey signal names in their line colors (TeX), used instead of a legend
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
