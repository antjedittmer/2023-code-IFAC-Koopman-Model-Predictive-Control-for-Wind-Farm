function [R, files] = plotAICnoV1IdVal(t0, tEndVal)
% plotAICnoV1IdVal open-loop plots of the Koopman effective wind speed
% model without V1 (identifyAICnoV1, 7 m/s, non-adaptive) on its
% identification set (CT levels at V_inf = 7 m/s) and on the validation
% set up to the wind speed ramp (also 7 m/s, other CT sequence), plotted as
% the diagnostic plots (plotEDMDinputsEffWind: CT, V_inf and V1, Vr1, Vr2).
% The validation after the ramp is shown by validateAICnoV1.
%
% Inputs
% - t0: first sample (s), default 240 (wake established)
% - tEndVal: last sample (s) of the validation set, default 5000 (ramp)
%
% Outputs
% - R: VAF of Vr1, Vr2 on both sets
% - files: png/eps files without extension (DataInOutWfSim/AdaptiveAIC),
%   {identification, validation}

if nargin < 1 || isempty(t0), t0 = 240; end
if nargin < 2 || isempty(tEndVal), tEndVal = 5000; end
mainDir = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(mainDir,'WFSimCode','bin','Kmpc')); % simAdaptiveKoopman
koopmanDir = fileparts(mfilename('fullpath')); % vaf, plotEDMDinputsEffWind
addpath(fullfile(koopmanDir,'2.DYNAMIC_MODE_DECOMPOSITION'), genpath(fullfile(koopmanDir,'Functions')));
dirOut = fullfile(mainDir,'DataInOutWfSim','AdaptiveAIC');
tmp = load(fullfile(dirOut,'KoopAICnoV1.mat'),'model');
model = tmp.model;

% {data file, last sample, title}
sets = {fullfile('Vinf7dot0_OL_Ct',model.filenameId), inf, 'Id. 7 m/s'
    fullfile('Vinf7dot0To9dot0_OL_Ct','Vinf9dot0_sowfa_2turb_ct_levels_val'), tEndVal, 'Val. 7 m/s'};
Wp.mesh = struct('Lx',1882.1,'Ly',800.0,'Nx',50,'Ny',25);
nTx = floor(0.4e3/Wp.mesh.Lx*Wp.mesh.Nx);
nTy = round(400/Wp.mesh.Ly*Wp.mesh.Ny);
R = table();
files = cell(1,size(sets,1));
for idx = 1:size(sets,1)
    S = load(fullfile(mainDir,'DataT2OLWFSim',sets{idx,1}),'Ur1','Ur2','Ct1','Ct2','u');
    k = t0:min(numel(S.Ur1), sets{idx,2});
    Ur = [S.Ur1(k); S.Ur2(k)];
    W = [S.Ct1(k); S.Ct2(k)];
    utemp = reshape(S.u, Wp.mesh.Nx, Wp.mesh.Ny, []); % V_inf and V1 for the plot only
    V1 = squeeze(utemp(nTx,nTy,k))';
    Vinf = squeeze(utemp(1,nTy,k))';

    UrEst = simAdaptiveKoopman(model, Ur, W, false); % non-adaptive
    vafUr = [vaf(Ur(1,:), UrEst(1,:)), vaf(Ur(2,:), UrEst(2,:))];
    R = [R; table(string(sets{idx,3}), vafUr(1), vafUr(2), ...
        'VariableNames', {'Set','VAF_Ur1','VAF_Ur2'})]; %#ok<AGROW>
    fig = plotEDMDinputsEffWind(UrEst, [W; V1], Ur, max(vafUr,0), dirOut, 1, model.n, ...
        sets{idx,3}, 1, Vinf, 'U', k, 'tex');
    set(fig,'Name',sets{idx,3});
    files{idx} = fullfile(dirOut,['Val',matlab.lang.makeValidName(sets{idx,3})]); % as plotEDMDinputsEffWind
end
disp(R)
end
