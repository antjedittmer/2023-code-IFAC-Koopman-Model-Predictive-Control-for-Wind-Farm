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

addpath(fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'WFSimCode','libraries','figureExport')); % setFigureSize, printFigure
noSub = size(Inputs_val,1) + 2;
fig1 = figure;
tl = tiledlayout(fig1, noSub, 1, 'TileSpacing', 'compact', 'Padding', 'compact'); % less white space than subplot
posDefault = get(gcf, 'position');
hFac = 1.3; % 4 subplots: 1.3 times the standard figure
if noSub == 5 % as the closed-loop plots (plotWFSimTimeseriesTex)
    hFac = 1.6;
end
set(gcf, 'position', [posDefault(1),posDefault(2)-150,posDefault(3),posDefault(4)*hFac]);
setFigureSize(fig1, [1 hFac]); % export size independent of the screen
cl = colormap('lines');
lw = 1;
fs = 12;
if isempty(tVec)
    kplot = 60: length(ysim_val);
    xk = 1:length(kplot);
else % time axis: the data start after the start-up (t0), show all
    kplot = 1: length(ysim_val);
    xk = tVec(kplot);
end

nexttile(tl,1)
plot(xk,Inputs_val(1,kplot),'LineWidth',lw,'Color', cl(1,:)); %1
hold on;
plot(xk,Inputs_val(2,kplot),'LineWidth',lw,'Color', cl(2,:)) %3
grid on; axis tight; ylim([0 2+0.2])
ylabel('Thrust $C_T^\prime$ [-]','interpreter','latex','Fontsize',fs); % name and symbol as in chapter 4

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
    nexttile(tl,2)
    if isempty(Vinf)
    plot(xk,Inputs_val(3,kplot),'LineWidth',lw,'Color', cl(1,:)); %1
    else
         plot(xk,Vinf(kplot),'k','LineWidth',lw); hold on; plot(xk,Inputs_val(3,kplot),'LineWidth',lw,'Color', cl(1,:));
    end
    grid on; axis tight
    if Vin && ~isempty(Vinf)
        ylabel('Wind $V$ [m/s]','interpreter','latex','Fontsize',fs);
        legend('$V_{\inf}$','$V_1$','interpreter','latex','Fontsize',fs,'orientation', 'horizontal','Location','NorthWest');
        legend('boxoff')
    elseif Vin
        ylabel('Wind $V_1$ [m/s]','interpreter','latex','Fontsize',fs);
    else
        ylabel('Yaw $\gamma$ [deg]','interpreter','latex','Fontsize',fs);
    end
    set(gca,'fontsize', fs)
end

nexttile(tl,noSub-2)
plot(xk,Deterministic_val(1,kplot),'LineWidth',lw,'Color', cl(1,:)); %1
hold on;
plot(xk,ysim_val(1,kplot),'Color',[0,0,0.6],'LineWidth',lw,'LineStyle','--') %3
grid on; axis tight;
pos = axis; axis([pos(1:3), pos(4)+0.3]);
%xlabel('k [-]','interpreter','latex')
ylabel('Wind $U_{r1}$ [m/s]','interpreter','latex')
lgh = legend({'Real','Est.'},'Location','southeast','interpreter','latex','fontsize',fs,'orientation','horizontal');
%  0.6155    0.4802    0.2702    0.0397
lgh.Position = [lgh.Position(1) + 0.02,lgh.Position(2:4)];
legend('boxoff')
set(gca,'fontsize', fs)

nexttile(tl,noSub-1)
plot(xk,Deterministic_val(2,kplot),'LineWidth',lw,'Color', cl(2,:)); %1
hold on;
plot(xk,ysim_val(2,kplot),'Color',[0.6, 0, 0.0],'LineWidth',lw,'LineStyle','--') %3
grid on; axis tight
pos = axis; axis([pos(1:3), pos(4)+0.3]);
ylabel('Wind $U_{r2}$ [m/s]','interpreter','latex','fontsize', fs)
lgh = legend({'Real','Est.'},'Location','southeast','interpreter','latex','fontsize',fs,'orientation','horizontal');
lgh.Position = [0.2595    0.3585    0.2702    0.0397];
legend('boxoff')
set(gca,'fontsize', fs)

nexttile(tl,noSub)
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
        plot(tMark*[1 1], [0 1], 'Color', marks.Color(idx,:), 'LineWidth', 2, 'LineStyle', '-'); % markers, not signals: 2
        % time stamp at the bottom, right of the line (left if TextSide < 0,
        % with more space, the updates before are dense), the magnitude of
        % TextSide scales the distance
        side = 1;
        if ismember('TextSide', marks.Properties.VariableNames), side = marks.TextSide(idx); end
        hAlign = {'right','left'};
        gap = 0.008 + 0.012*(side < 0);
        text(tMark + side*gap*(xk(end)-xk(1)), 0.05, sprintf('%d s', marks.Time_s(idx)), ...
            'Color', marks.Color(idx,:), 'FontSize', 12, 'FontWeight', 'bold', ...
            'VerticalAlignment', 'bottom', 'HorizontalAlignment', hAlign{(side > 0) + 1});
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
if ~isempty(tVec) % zero line visible above the axis, as the closed loop
    ylim([-0.1 1.1]);
end

%legend('P','Location','northeast','interpreter','latex','fontsize',fs,'orientation','horizontal');
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
    axs(1).YLabel.String = 'Update (-)';
end
setFigureFont(fig1, fs); % all text 12 pt, bold titles

%sys_red{1,1}.Xint = xo;
strVal1 = matlab.lang.makeValidName(strVal);
if strcmpi(interp,'tex') % thesis figures: y axes not completely tight
    padYLim(fig1);
end
printFigure(fig1, fullfile(dirFig,['ValP',strVal1])); % fixed size (setFigureSize)
end

function str = latex2tex(str)
% latex2tex LaTeX label to TeX: no $, \inf -> \infty, units [..] -> (..),
% U_r -> V_r, C_T^\prime -> C_T'
if iscell(str), str = cellfun(@latex2tex, str, 'UniformOutput', false); return; end
str = strrep(str,'$','');
str = strrep(str,'^\prime','''');
str = strrep(str,'\inf}','\infty}');
str = regexprep(str,'\[([^\]]*)\]','($1)');
str = strrep(str,'(deg)',['(',char(176),')']); % degree sign as in the qLMPC plots
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
