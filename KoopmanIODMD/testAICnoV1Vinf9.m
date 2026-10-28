function [R, fig] = testAICnoV1Vinf9(dirFig, t0, ad)
% testAICnoV1Vinf9 test of the Koopman effective wind speed models without
% V1 (inputs CT1, CT2, outputs Vr1, Vr2) on unseen data at constant
% V_inf = 9 m/s (validation CT levels, seed 4, not used for any model):
%   - offline ID 7 m/s: initial model of the online update
%   - offline ID 9 m/s: identified on the 9 m/s training data (seed of
%     sowfa_2turb_ct_levels)
%   - learned: 7 m/s model after the online update on the ramp from 7 to
%     9 m/s (validateAICnoV1, plotAICnoV1Sigma), frozen (an update on the
%     test data would not trigger: the Vr1 error stays below nu)
% Estimation as in the MPC (simAdaptiveKoopman, free run from the measured
% initial state).
%
% Inputs
% - dirFig: folder of AdaptiveAIC_Test9 (.eps, .png), default
%   DataInOutWfSim/AdaptiveAIC
% - t0: first sample (s), default 240 (wake established, as validateAICnoV1)
% - ad: settings of the online update (validateAICnoV1)
%
% Outputs
% - R: table of the VAF and RMSE of Vr1, Vr2
% - fig: figure handle

mainDir = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(mainDir,'WFSimCode','bin','Kmpc')); % koopmanAdaptStep
addpath(fullfile(mainDir,'WFSimCode','libraries','figureExport')); % setFigureSize, printFigure
koopmanDir = fileparts(mfilename('fullpath')); % vaf
addpath(genpath(fullfile(koopmanDir,'Functions')));
dirOut = fullfile(mainDir,'DataInOutWfSim','AdaptiveAIC');
if nargin < 1 || isempty(dirFig), dirFig = dirOut; end
if nargin < 2 || isempty(t0), t0 = 240; end
if nargin < 3, ad = struct('h',5,'nu',0.15,'resetState',true,'trigChannels',1); end

% Offline models
file7 = fullfile(dirOut,'KoopAICnoV1.mat');
if ~isfile(file7), identifyAICnoV1(); end
file9 = fullfile(dirOut,'KoopAICnoV1_Vinf9.mat');
if ~isfile(file9), identifyAICnoV1('Vinf9dot0_sowfa_2turb_ct_levels', [], 'KoopAICnoV1_Vinf9.mat'); end
tmp = load(file7,'model'); model7 = tmp.model;
tmp = load(file9,'model'); model9 = tmp.model;

% Learned model: online update on the ramp validation data (as validateAICnoV1)
f = dir(fullfile(mainDir,'DataT2OLWFSim','Vinf7dot0To9dot0_OL_Ct','*sowfa_2turb_ct_levels_val.mat'));
S = load(fullfile(f(1).folder,f(1).name),'Ur1','Ur2','Ct1','Ct2');
k = 240:numel(S.Ur1);
[~, ~, ~, modelL] = simAdaptiveKoopman(model7, [S.Ur1(k); S.Ur2(k)], [S.Ct1(k); S.Ct2(k)], true, ad);

% Test data at 9 m/s
S = load(fullfile(mainDir,'DataT2OLWFSim','Vinf9dot0_OL_Ct', ...
    'Vinf9dot0_sowfa_2turb_ct_levels_val.mat'),'Ur1','Ur2','Ct1','Ct2');
k = t0:numel(S.Ur1);
Ur = [S.Ur1(k); S.Ur2(k)];
W = [S.Ct1(k); S.Ct2(k)];

% The VAF ignores offsets (offline ID 7 m/s: about 1.5 m/s too low), hence
% also the RMSE
cases = {'Offline ID 7 m/s', model7
         'Offline ID 9 m/s', model9
         'Learned', modelL};
nC = size(cases,1);
UrEst = cell(nC,1);
R = table();
for idx = 1:nC
    UrEst{idx} = simAdaptiveKoopman(cases{idx,2}, Ur, W, false);
    rmse = sqrt(mean((Ur - UrEst{idx}).^2, 2));
    R = [R; table(string(cases{idx,1}), vaf(Ur(1,:), UrEst{idx}(1,:)), ...
        vaf(Ur(2,:), UrEst{idx}(2,:)), rmse(1), rmse(2), ...
        'VariableNames', {'Model','VAF_Ur1','VAF_Ur2','RMSE_Ur1','RMSE_Ur2'})]; %#ok<AGROW>
end

% Plot: measured and estimated Vr1, Vr2
fig = figure;
cl = [0 0 0; lines(nC)];
iPlot = 1:nC;
for idxV = 1:2
    subplot(2,1,idxV);
    plot(k, Ur(idxV,:), 'Color', cl(1,:), 'LineWidth', 1.2); hold on;
    for idx = iPlot
        plot(k, UrEst{idx}(idxV,:), 'Color', cl(idx+1,:), 'LineWidth', 1);
    end
    grid on; axis tight;
    ylabel(sprintf('V_{r%d} (m/s)', idxV));
    set(gca, 'FontSize', 11);
    if idxV == 2, xlabel('Time (s)', 'FontSize', 12); end
    rmseStr = sprintf('%.2f/', R{iPlot, sprintf('RMSE_Ur%d', idxV)});
    rmseStr = ['RMSE = ', rmseStr(1:end-1), ' m/s'];
    if idxV == 1
        title({['Test data 9 m/s (unseen), V_{r1}: ', rmseStr], ...
            ['Measured, ', colorKey(cellstr(R.Model(iPlot)), cl(iPlot+1,:))]}, ...
            'FontWeight', 'normal', 'Interpreter', 'tex');
    else
        title(['V_{r2}: ', rmseStr], 'FontWeight', 'normal', 'Interpreter', 'tex');
    end
end
setFigureSize(fig, [1 1]); % standard figure, independent of the screen
printFigure(fig, fullfile(dirFig, 'AdaptiveAIC_Test9'));
disp(R)
end

function str = colorKey(names, colors)
% colorKey signal names in their line colors (TeX), used instead of a legend
parts = cell(size(names));
for idx = 1:numel(names)
    parts{idx} = sprintf('{\\color[rgb]{%.3f,%.3f,%.3f}%s}', colors(idx,:), names{idx});
end
str = strjoin(parts, ', ');
end
