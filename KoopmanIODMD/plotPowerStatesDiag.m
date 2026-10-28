function R = plotPowerStatesDiag(modelFile, filenameVal, t0)
% plotPowerStatesDiag diagnostic plots of the open-loop identification of
% the input-only Koopman power model (KoopFeatModel<k>.mat, runKoopStatesCL):
% per validation set one figure with the rotor effective wind speeds Ur1,
% Ur2 and one with the farm power, WFSim and the linear model simulated
% from the measured initial state (as testPowerStates), VAF in the legend.
% Plots with plotEDMDinputsEffWind (as MainWfSimVinYawPowerOut), saved in
% DataInOutWfSim/PowerInputsOnly/Diag_<model>_Ur and Diag_<model>_P
%
% Inputs
% - modelFile: default KoopFeatModel3.mat ('Ur n=6')
% - filenameVal: validation sets (cell)
% - t0: first sample (s) as in the identification

mainDir = fileparts(fileparts(mfilename('fullpath')));
dirOut = fullfile(mainDir,'DataInOutWfSim','PowerInputsOnly');
if nargin < 1 || isempty(modelFile), modelFile = 'KoopFeatModel3.mat'; end
if nargin < 2 || isempty(filenameVal)
    filenameVal = {'Vinf8dot0_sowfa_2turb_yaw_ct_levels_val','Vinf8dot0_sowfa_2turb_yaw_steps', ...
        'Vinf8dot0_sowfa_2turb_yaw_steps_Ct_comb'};
end
if nargin < 3, t0 = 240; end
addpath(fullfile(mainDir,'WFSimCode','bin','Kmpc')); % powerFeatures, koopFeatState

tmp = load(fullfile(dirOut,modelFile),'model');
model = tmp.model;
sys = ss(model.A, model.B, model.C, model.D, 1);
[~, modelName] = fileparts(modelFile);
dirUr = fullfile(dirOut, sprintf('Diag_%s_Ur', modelName));
dirP = fullfile(dirOut, sprintf('Diag_%s_P', modelName));
if ~isfolder(dirUr), mkdir(dirUr); end
if ~isfolder(dirP), mkdir(dirP); end
R = table();
for idxV = 1:numel(filenameVal)
    evalc('[~, D] = testPowerInputsOnly(model.filenameId, filenameVal{idxV}, 8, t0);');
    close(findobj('Type','figure','Name','PowerInputsOnly'));
    S = D.val;

    % States from the second sample on (initial state: measured Ur, history)
    urS = struct('M',S.lifted(1:2,1),'UrPrev',S.lifted(1:2,1),'k',1);
    X = zeros(model.n, size(S.lifted,2));
    for k = 1:size(S.lifted,2)
        [X(:,k), urS] = koopFeatState(model, S.lifted(1:2,k), urS);
    end
    W = powerFeatures(S.lifted(end-2:end,:), model.featureSet, model.tau);
    W = W(:,2:end); X = X(:,2:end);
    P = S.P(2:end)/1e6;
    Ur = S.lifted(1:2,2:end);
    [y, ~, xSim] = lsim(sys, W', [], X(:,1));
    y = y'; UrSim = xSim(:,1:2)';

    vafUr = [vaf(Ur(1,:), UrSim(1,:)), vaf(Ur(2,:), UrSim(2,:))];
    vafP = vaf(P, y);
    setName = strrep(filenameVal{idxV},'Vinf8dot0_sowfa_2turb_','');
    R = [R; table(string(setName), vafUr(1), vafUr(2), vafP, ...
        'VariableNames', {'ValSet','VAF_Ur1','VAF_Ur2','VAF_P'})]; %#ok<AGROW>

    % Diagnostic plots as MainWfSimVinYawPowerOut (plotEDMDinputsEffWind):
    % inputs ct1, ct2, phi1 and WFSim (Real) vs. model (Est.)
    % Title: set, states and (added by plotEDMDinputsEffWind) the VAF
    strVal = sprintf('Val. %s, $n_K$ = %d (%s)', strrep(setName,'_','-'), model.n, stateStr(model.stateSet));
    inp = S.lifted(end-2:end,2:end);
    t = t0 + (1:numel(P)); % time (s), samples from t0 + 1 on (dt = 1 s)
    fig = plotEDMDinputsEffWind(UrSim, inp, Ur, vafUr, dirUr, 1, model.n, strVal, 1, [], 'U', t, 'tex');
    set(fig,'Name',sprintf('Ur %s',setName));
    fig = plotEDMDinputsEffWind(y, inp, P, vafP, dirP, 0, model.n, strVal, 1, [], 'P', t, 'tex'); % Vin = 0: 3rd input labelled gamma
    set(fig,'Name',sprintf('P %s',setName));
end
disp(R)
end

function str = stateStr(stateSet)
% Koopman states of the state set (testPowerStates) as LaTeX string
nL = sscanf(stateSet,'Ur n=%d');
names = {'$U_r$','$U_r^2$','$U_r^3$','$\Delta U_r$'};
str = strjoin(names(1:nL/2), ', ');
if endsWith(stateSet,'+d')
    str = [str, ', $U_r(k-1)$'];
end
end

