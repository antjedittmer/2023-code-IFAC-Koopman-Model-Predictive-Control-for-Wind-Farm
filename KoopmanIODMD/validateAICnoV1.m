function R = validateAICnoV1(filenameVal, t0, ad)
% validateAICnoV1 open-loop validation of the Koopman effective wind speed
% model without V1 (identifyAICnoV1, inputs CT1, CT2 only, no yaw) on the
% data set with the wind speed ramp from 7 to 9 m/s: non-adaptive and
% online updated model (simAdaptiveKoopman), plotted as in the thesis
% (plotEDMDinputsEffWindUpdate: CT, V_inf and V1, Ur1, Ur2, update signal).
% V1 is only shown, it is not used by the model.
%
% Inputs
% - filenameVal: open-loop data with the ramp (default: CT levels, ramp at 5000 s)
% - t0: first sample (s)
% - ad: settings of the online update (koopmanAdaptStep, simAdaptiveKoopman),
%   default h = 5, nu = 0.15 m/s, resetState = true, trigger on Vr1 only
%   (trigChannels = 1, as in the controller)

mainDir = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(mainDir,'WFSimCode','bin','Kmpc')); % koopmanAdaptStep
koopmanDir = fileparts(mfilename('fullpath')); % vaf, plotEDMDinputsEffWindUpdate
addpath(fullfile(koopmanDir,'2.DYNAMIC_MODE_DECOMPOSITION'), genpath(fullfile(koopmanDir,'Functions')));
if nargin < 1 || isempty(filenameVal)
    filenameVal = 'sowfa_2turb_ct_levels_val';
end
if nargin < 2, t0 = 240; end
if nargin < 3, ad = struct('h',5,'nu',0.15,'resetState',true,'trigChannels',1); end
dirOut = fullfile(mainDir,'DataInOutWfSim','AdaptiveAIC');
modelFile = fullfile(dirOut,'KoopAICnoV1.mat');
if ~isfile(modelFile)
    identifyAICnoV1();
end
tmp = load(modelFile,'model');
model = tmp.model;

% Validation data (folder of WFSim_demo for the wind time series)
f = dir(fullfile(mainDir,'DataT2OLWFSim','Vinf7dot0To9dot0_OL_Ct',['*',filenameVal,'.mat']));
S = load(fullfile(f(1).folder,f(1).name),'Ur1','Ur2','Ct1','Ct2','u');
k = t0:numel(S.Ur1);
Ur = [S.Ur1(k); S.Ur2(k)];
W = [S.Ct1(k); S.Ct2(k)];

% V_inf and V1 (u in the grid cell in front of WT1) for the plot only
Wp.mesh = struct('Lx',1882.1,'Ly',800.0,'Nx',50,'Ny',25);
nTx = floor(0.4e3/Wp.mesh.Lx*Wp.mesh.Nx);
nTy = round(400/Wp.mesh.Ly*Wp.mesh.Ny);
utemp = reshape(S.u, Wp.mesh.Nx, Wp.mesh.Ny, []);
V1 = squeeze(utemp(nTx,nTy,k))';
Vinf = squeeze(utemp(1,nTy,k))';

R = table();
for adaptive = [false true]
    [UrEst, upd] = simAdaptiveKoopman(model, Ur, W, adaptive, ad);
    vafUr = [vaf(Ur(1,:), UrEst(1,:)), vaf(Ur(2,:), UrEst(2,:))];
    R = [R; table(adaptive, vafUr(1), vafUr(2), sum(upd), ...
        'VariableNames', {'Adaptive','VAF_Ur1','VAF_Ur2','NoUpdates'})]; %#ok<AGROW>
    if adaptive
        strVal = 'Val. adaptive';
    else
        strVal = 'Val. non-adaptive';
    end
    fig = plotEDMDinputsEffWindUpdate(UrEst, [W; V1], Ur, max(vafUr,0), dirOut, 1, model.n, ...
        strVal, 1 + upd, 1, Vinf, k, 'tex'); % time axis (s), TeX
    set(fig,'Name',strVal);
end
disp(R)
end
