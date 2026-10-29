function fig1 = plotEDMDinputsEffWind(ysim_val,Inputs_val,Deterministic_val,FITje_val,dirFig,Vin,nx,strVal,displayOff,Vinf,strOut,tVec,interp)
% plotEDMDinputsEffWind plots the simulated outputs of the eDMD estimations
% and prints them into the directory given in dirFig.
% Optional tVec: time (s) of the samples, x axis 'Time (s)' (default: sample index)
% Optional interp: 'latex' (default) or 'tex' (labels converted: no $, units
% in (), no legends: the first title names the measured (WT1, WT2) and the
% estimated signals in their colors)

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

if nargin < 10
    Vinf = [];
end

if nargin < 11
    strOut = 'U';
end

if nargin < 12
    tVec = [];
end

if nargin < 13
    interp = 'latex';
end

if nx > 0
    strK =  sprintf(', $n_K$ = %d', nx);
else
    strK = '';
end


if size(ysim_val,2) < size(ysim_val,1)
    ysim_val = ysim_val';
end


addpath(fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'WFSimCode','libraries','figureExport')); % setFigureSize, printFigure
noSubIn = size(Inputs_val,1);
noSubOut =  min(size(ysim_val));
noSub = noSubOut + noSubIn -1;
fig1 = figure;
tl = tiledlayout(fig1, noSub, 1, 'TileSpacing', 'compact', 'Padding', 'compact'); % less white space than subplot
hFac = 1; % height: 5 subplots 1.6 (as the closed-loop plots), 4 subplots 1.3, 3 subplots 1.15 times the standard figure
if noSub == 5, hFac = 1.6; elseif noSub == 4, hFac = 1.3; elseif noSub == 3, hFac = 1.15; end
if hFac > 1
    posDefault = get(gcf, 'position');
    set(gcf, 'position', [posDefault(1),posDefault(2)-150,posDefault(3),posDefault(4)*hFac]);
end
setFigureSize(fig1, [1 hFac]); % export size independent of the screen
cl = lines;
lw = 1;
fs = 12;
kplot = 60: length(ysim_val);
if isempty(tVec)
    xk = 1:length(kplot);
else
    xk = tVec(kplot);
end

nexttile(tl,1)
plot(xk,Inputs_val(1,kplot),'LineWidth',lw,'Color', cl(1,:)); %1
hold on;
plot(xk,Inputs_val(2,kplot),'LineWidth',lw,'Color', cl(2,:)) %3
grid on; axis tight; ylim([0 2+0.2])
ylabel('Thrust $C_T^\prime$ [-]','interpreter','latex','Fontsize',fs); % name and symbol as in chapter 4

if displayOff && length(FITje_val) >1
    title(sprintf('%s:  VAF %sr1: %d, r2: %d',...
        strVal,strOut,round(FITje_val(1:2))),'interpreter','latex','Fontsize',fs);
elseif length(FITje_val)>1
    title(sprintf('%s Set%s:  VAF %sr1: %d, r2: %d',...
        strVal,strK,strOut,round(FITje_val(1:2))),'interpreter','latex','Fontsize',fs);
elseif displayOff
    title(sprintf('%s:  VAF %sT: %d',...
        strVal,strOut,round(FITje_val)),'interpreter','latex','Fontsize',fs);
else
    title(sprintf('%s Set%s:  VAF %sT: %d',...
        strVal,strK,strOut,round(FITje_val)),'interpreter','latex','Fontsize',fs);
end


strLeg = sprintf('WT%d,',1:2);
legend(regexp(strLeg(1:end-1), ',', 'split'),'interpreter','latex','Fontsize',fs, 'Orientation','horizontal', 'location','southeast');
%legend('boxoff')
set(gca,'fontsize', fs)

if noSubIn > 2
    nexttile(tl,2)
    if isempty(Vinf)
        plot(xk,Inputs_val(3,kplot),'LineWidth',lw,'Color', cl(1,:)); %1
    else
        plot(xk,Vinf(kplot),'k','LineWidth',lw); hold on; plot(xk,Inputs_val(3,kplot),'LineWidth',lw,'Color', cl(1,:));
    end
    grid on; axis tight
    if ~isempty(Vinf) % V_inf line not on the axis limit
        ylim(ylim + [-0.2 0.2]);
    end
    if Vin && ~isempty(Vinf)
        ylabel('Wind $V$ [m/s]','interpreter','latex','Fontsize',fs);
        legend('$V_{\inf}$','$V_1$','interpreter','latex','Fontsize',fs,'orientation', 'horizontal','Location','NorthWest');
        %legend('boxoff')
    elseif (Vin && noSub == 4) || (~Vin && noSub == 3) 
         ylabel('Yaw $\gamma$ [deg]','interpreter','latex','Fontsize',fs);       
    else
        ylabel('Wind $V_1$ [m/s]','interpreter','latex','Fontsize',fs);
    end
    set(gca,'fontsize', fs)
    
    if noSubIn > 3 %&& noSub ~= 4
        nexttile(tl,3)
        plot(xk,Inputs_val(4,kplot),'LineWidth',lw,'Color', cl(1,:)); %1
        ylabel('Wind $V_1$ [m/s]','interpreter','latex','Fontsize',fs);
        grid on; axis tight
    end
    
end

if noSubOut == 1
    strAdd = {'_T$ [MW]'};
else
    strAdd = {'_{r1}$ [m/s]','_{r2}$ [m/s]'};
end

for idx = 1: noSubOut
    idxT = 2-mod(idx,2);
    nexttile(tl,noSub - (noSubOut-idx))
    plot(xk,Deterministic_val(idx,kplot),'LineWidth',lw,'Color', cl(idxT,:)); %1
    hold on;
    plot(xk,ysim_val(idx,kplot),'Color',[0,0,0],'LineWidth',lw,'LineStyle','--') %3
    grid on; axis tight;
    pos = axis; axis([pos(1:3), pos(4)+0.3]);

    %xlabel('k [-]','interpreter','latex')
    sigName = 'Wind';
    if noSubOut == 1, sigName = 'Power'; end
    ylabel([sigName,' $',strOut,strAdd{idxT}],'interpreter','latex')
    legend({'Real','Est.'},'Location','northwest','interpreter','latex','fontsize',fs,'orientation','horizontal')
    %legend('boxoff')
    set(gca,'fontsize', fs)
    
end
if ~isempty(tVec)
    xlabel('Time (s)','interpreter','latex','Fontsize',fs);
end

if strcmpi(interp,'tex') % convert all LaTeX strings to TeX
    for hTxt = findall(fig1,'Type','text')'
        set(hTxt,'String',latex2tex(get(hTxt,'String')),'Interpreter','tex');
    end
    delete(findall(fig1,'Type','legend'));
    axs = findall(fig1,'Type','axes');
    set(axs,'TickLabelInterpreter','tex');
    axTop = axs(end); % first subplot
    if noSubOut == 1 % estimated output(s), black dashed
        estNames = {[strOut,'_{T,est}']};
    else
        estNames = {'V_{r1,est}','V_{r2,est}'};
    end
    measNames = {'WT1','WT2'}; measColors = [cl(1,:); cl(2,:)];
    if Vin && ~isempty(Vinf) % free stream wind speed (black) in the wind subplot
        measNames{end+1} = 'V_\infty'; measColors(end+1,:) = [0 0 0];
    end
    axTop.Title.String = {axTop.Title.String, ...
        ['Measured sig.: ', colorKey(measNames, measColors), ...
        '; Estimated sig.: ', colorKey(estNames, zeros(numel(estNames),3))]};
end
setFigureFont(fig1, fs); % all text 12 pt, bold titles
for ax = findall(fig1,'Type','axes')' % time in s without x10^4 (long identification sets): no extra row
    ax.XAxis.Exponent = 0;
end

%sys_red{1,1}.Xint = xo;
strVal1 = matlab.lang.makeValidName(strVal);
printFigure(fig1, fullfile(dirFig,['Val',strVal1])); % fixed size (setFigureSize)
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
str = regexprep(str,'VAF PT:','VAF P_T:'); % farm power as in the y label
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
