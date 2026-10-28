function [fig, sel] = plotAICnoV1Sigma(dirFig, filenameVal, t0, ad, tSel)
% plotAICnoV1Sigma singular values of the Koopman effective wind speed model
% without V1 (identifyAICnoV1: inputs CT1, CT2, outputs Vr1, Vr2) during the
% online update, as Fig. 3 of the CDC22 paper: the model is updated in the
% open-loop validation with the wind speed ramp from 7 to 9 m/s
% (simAdaptiveKoopman, settings as validateAICnoV1) and compared with the
% models identified offline at constant 7 m/s (initial model) and 9 m/s.
% The singular values are those of the lifted linear model
% phi(k+1) = A phi(k) + B w(k), [Vr1; Vr2] = [I 0] phi(k), sample time 1 s.
%
% Inputs
% - dirFig: folder of AdaptiveAIC_OL_Sigma (.eps, .png), default
%   DataInOutWfSim/AdaptiveAIC
% - filenameVal, t0, ad: validation data and update settings (validateAICnoV1)
% - tSel: times (s) of the shown models, at most 4 (default: initial, after a
%   quarter, half and all of the updates)
%
% Outputs
% - fig: figure handle
% - sel: table of the shown models (time, number of updates so far, line
%   color), input marks of validateAICnoV1

mainDir = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(mainDir,'WFSimCode','bin','Kmpc')); % koopmanAdaptStep
dirOut = fullfile(mainDir,'DataInOutWfSim','AdaptiveAIC');
if nargin < 1 || isempty(dirFig), dirFig = dirOut; end
if nargin < 2 || isempty(filenameVal), filenameVal = 'sowfa_2turb_ct_levels_val'; end
if nargin < 3 || isempty(t0), t0 = 240; end
if nargin < 4 || isempty(ad), ad = struct('h',5,'nu',0.15,'resetState',true,'trigChannels',1); end
if nargin < 5, tSel = []; end

% Offline models at 7 m/s (initial model of the update) and 9 m/s
file7 = fullfile(dirOut,'KoopAICnoV1.mat');
if ~isfile(file7), identifyAICnoV1(); end
file9 = fullfile(dirOut,'KoopAICnoV1_Vinf9.mat');
if ~isfile(file9), identifyAICnoV1('Vinf9dot0_sowfa_2turb_ct_levels', [], 'KoopAICnoV1_Vinf9.mat'); end
tmp = load(file7,'model'); model7 = tmp.model;
tmp = load(file9,'model'); model9 = tmp.model;

% Validation data with the ramp (as validateAICnoV1)
f = dir(fullfile(mainDir,'DataT2OLWFSim','Vinf7dot0To9dot0_OL_Ct',['*',filenameVal,'.mat']));
S = load(fullfile(f(1).folder,f(1).name),'Ur1','Ur2','Ct1','Ct2');
k = t0:numel(S.Ur1);
Ur = [S.Ur1(k); S.Ur2(k)];
W = [S.Ct1(k); S.Ct2(k)];

% Shown models: at the times tSel (s) or initial, after a quarter, half and
% all of the updates (a single update hardly changes the model, the initial
% correlation matrices are sums over the whole identification data)
[~, upd] = simAdaptiveKoopman(model7, Ur, W, true, ad);
if ~isempty(tSel)
    saveAt = tSel - k(1) + 1;
else
    iUpd = find(upd);
    saveAt = unique([1, iUpd(ceil(numel(iUpd)*[0.25 0.5 1]))]);
end
[~, ~, ~, ~, models] = simAdaptiveKoopman(model7, Ur, W, true, ad, saveAt);
% Colors: initial model gray, updated models in colors not used for WT1,
% WT2 (lines 1, 2), the same colors mark the update times in the
% validation plot (validateAICnoV1)
cl = lines(5);
cl = [0.5 0.5 0.5; cl(3:2+numel(models)-1,:)];
sel = table(k([models.k])', arrayfun(@(m) sum(upd(1:m.k)), models)', cl, ...
    'VariableNames', {'Time_s','NoUpdates','Color'});

% Singular values
omega = logspace(-5, log10(pi), 600);
sv = @(A,B) sigma(ss(A, B, eye(2,size(A,1)), zeros(2,size(B,2)), 1), omega);

fig = figure;
loglog(omega, sv(model9.A, model9.B)', 'k--', 'LineWidth', 2); hold on; % below the updated models
ls = {'-','-.','--',':'};
lw = [0.75 0.75 0.75 1];
names = cell(1,numel(models));
for idx = 1:numel(models)
    loglog(omega, sv(models(idx).A, models(idx).B)', 'Color', cl(idx,:), ...
        'LineStyle', ls{idx}, 'LineWidth', lw(idx));
    names{idx} = sprintf('%d s', k(models(idx).k));
end
grid on; axis tight;
xlabel('Frequency (rad/s)');
ylabel('Singular values (m/s)');
set(gca, 'FontSize', 11, 'TickLabelInterpreter', 'tex');
title({'Model C_{T1}, C_{T2} \rightarrow V_{r1}, V_{r2}, initial: offline ID 7 m/s', ...
    ['Model at t = ', colorKey(names, cl), '; offline ID ', colorKey({'9 m/s'}, [0 0 0])]}, ...
    'FontWeight', 'normal', 'Interpreter', 'tex');
printFixedWidth(fig, fullfile(dirFig, 'AdaptiveAIC_OL_Sigma'));
disp(sel)
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
% the aspect ratio of the window (as plotEDMDinputsEffWindUpdate)
pos = get(h, 'Position');
sizeCm = 14.82*[1, pos(4)/pos(3)];
set(h, 'PaperUnits', 'centimeters', 'PaperPositionMode', 'manual', ...
    'PaperSize', sizeCm, 'PaperPosition', [0 0 sizeCm]);
print(h, fName, '-depsc', '-painters', '-loose');
print(h, fName, '-dpng', '-r300');
end
