function [R, D] = testPowerInputsOnly(filenameId,filenameVal,noStates,t0)
% testPowerInputsOnly compares Koopman power models for constant free
% stream wind: the model of MainWfSimVinYawPowerOut with the measured wind
% speed V1 (u in the grid cell of turbine 1) as known input, and models with
% the control inputs ct1, ct2, phi1 only. Data, lifting and identification
% as in MainWfSimVinYawPowerOut (eDMD_RT_Vin_P, yawmode 1).
%
% Models (output: farm power PT1 + PT2)
%   'V1':   inputs ct1, ct2, phi1, V1 (as in MainWfSimVinYawPowerOut)
%   'u':    inputs ct1, ct2, phi1 (no wind speed, no offset)
%   'u+1':  inputs ct1, ct2, phi1 and a constant one (offset); for constant
%           wind the free stream wind speed is part of this offset
%
% Each model is simulated over the validation set from its first sample
%   'nonlin': lifting of the simulated Ur1, Ur2 in each step (eDMD_RT_Vin_P)
%   'lin':    linear state space model sys_red (as used in the MPC)
%   'lin,ss': as 'lin', initial state the steady state for the first inputs
%             (no measurement at all, model 'u+1' only)

if ~nargin
    filenameId = 'Vinf8dot0_sowfa_2turb_yaw_noise_step';
    filenameVal = 'Vinf8dot0_sowfa_2turb_yaw_steps';
    noStates = 6;
    t0 = 240;
end
yawmode = 1;
poly = 0;
n = noStates;

KoopmanDir = fileparts(mfilename('fullpath'));
addpath(fullfile(KoopmanDir,'2.DYNAMIC_MODE_DECOMPOSITION'), genpath(fullfile(KoopmanDir,'Functions')));
DataIn = fullfile(fileparts(KoopmanDir),'DataT2OLWFSim');
dirFig = fullfile(fileparts(KoopmanDir),'DataInOutWfSim','PowerInputsOnly');
if ~isfolder(dirFig)
    mkdir(dirFig);
end

%% Load data (as MainWfSimVinYawPowerOut, useVal = 1)
vars = {'Ct1','Ct2','PT1','PT2','Ur1','Ur2','phi1','u'};
if ischar(filenameId) % several identification sets (cell): prepared from t0 each and appended
    filenameId = {filenameId};
end
Ids = cell(size(filenameId));
for idxId = 1:numel(filenameId)
    sepStr = strfind(filenameId{idxId},'_');
    Ids{idxId} = load(fullfile(DataIn,[filenameId{idxId}(1:sepStr-1),'_OL_Ct'],filenameId{idxId}),vars{:});
end
sepStr = strfind(filenameVal,'_');
Val = load(fullfile(DataIn,[filenameVal(1:sepStr-1),'_OL_Ct'],filenameVal),vars{:});

Wp.mesh = struct('Lx',1882.1,'Ly',800.0,'Nx',50,'Ny',25);
Crx1 = 0.4000e3; Cry1 = 400;
nTx = floor(Crx1/Wp.mesh.Lx*Wp.mesh.Nx);
nTy = round(Cry1/Wp.mesh.Ly*Wp.mesh.Ny);

D = struct();
[D.id, Vinf] = prepareSet(Ids{1}, t0, nTx, nTy, Wp, true(size(Ids{1}.Ct1)), n);
for idxId = 2:numel(Ids)
    Sadd = prepareSet(Ids{idxId}, t0, nTx, nTy, Wp, true(size(Ids{idxId}.Ct1)), n);
    for fn = {'V1','P','P12','lifted'}
        D.id.(fn{1}) = [D.id.(fn{1}), Sadd.(fn{1})];
    end
end
D.val = prepareSet(Val, t0, nTx, nTy, Wp, Val.phi1 <= 30, n);
fprintf('Free stream wind: %.2f m/s (id.), V1 from %.2f to %.2f m/s\n', Vinf, min(D.id.V1), max(D.id.V1));

%% Identify and simulate the models
models = {'V1','u','u+1'};
R = table();
for idxM = 1:numel(models)
    S = {D.id, D.val};
    for k = 1:2
        switch models{idxM}
            case 'V1',  S{k}.extra = S{k}.V1;
            case 'u',   S{k}.extra = zeros(0,numel(S{k}.P));
            case 'u+1', S{k}.extra = ones(1,numel(S{k}.P));
        end
        S{k}.Z = [S{k}.P; S{k}.lifted; S{k}.extra];
    end
    stateName = ['PT;', S{1}.stateName, repmat(';extra', 1, size(S{1}.extra,1))];
    [sys_red,FITid,~,~,~,FITval,~,~,~,~,~,ysim_val] = ...
        eDMD_RT_Vin_P(S{1}.Z,stateName,S{2}.Z,poly,n,yawmode);

    % Linear model (MPC): inputs [ct1;ct2;phi1;extra], x0 from the data
    uVal = S{2}.Z(n+2:end,1:end-1);
    pVal = S{2}.P(1:end-1);
    x0 = S{2}.lifted(1:n,1);
    yLin = lsim(sys_red, uVal', [], x0)';
    row = {models{idxM}, FITid, FITval, vaf(pVal(1:end-1), ysim_val), rmse(pVal(1:end-1), ysim_val), ...
        vaf(pVal, yLin), rmse(pVal, yLin), NaN, NaN, max(abs(eig(sys_red.A)))};
    yLinSS = nan(size(yLin));
    if strcmp(models{idxM},'u+1')
        xss = (eye(n) - sys_red.A) \ (sys_red.B*uVal(:,1));
        yLinSS = lsim(sys_red, uVal', [], xss)';
        row(8:9) = {vaf(pVal, yLinSS), rmse(pVal, yLinSS)};
    end
    R = [R; cell2table(row, 'VariableNames', {'Model','VAFid','VAFval','VAFval_nonlin', ...
        'RMSEval_nonlin_MW','VAFval_lin','RMSEval_lin_MW','VAFval_linSS','RMSEval_linSS_MW','maxAbsEig'})]; %#ok<AGROW>
    Y.(matlab.lang.makeValidName(models{idxM})) = struct('nonlin', ysim_val, 'lin', yLin, 'linSS', yLinSS, 'sys', sys_red);
end
disp(R)

%% Plot validation
t = (0:numel(D.val.P)-2);
cl = lines(4);
fig = figure('Name','PowerInputsOnly');
ax = subplot(3,1,1);
plot(ax, t, D.val.lifted(n+1:n+2,1:end-1)'); hold(ax,'on');
yyaxis(ax,'right'); plot(ax, t, D.val.lifted(n+3,1:end-1),'k'); ylabel(ax,'\phi_1 (deg)');
yyaxis(ax,'left'); ylabel(ax,'C_T'''); legend(ax,'C_{T,1}''','C_{T,2}''','\phi_1'); grid(ax,'on'); axis(ax,'tight');
title(ax, sprintf('Validation %s, V_\\infty = %.1f m/s', strrep(filenameVal,'_','\_'), Vinf));
ax = subplot(3,1,2);
plot(ax, t, D.val.P(1:end-1)/1e6, 'k', 'LineWidth', 1.5); hold(ax,'on');
plot(ax, t, Y.V1.lin/1e6, '--', 'Color', cl(1,:));
plot(ax, t, Y.u.lin/1e6, ':', 'Color', cl(2,:), 'LineWidth', 1.2);
plot(ax, t, Y.u_1.lin/1e6, '-', 'Color', cl(3,:));
plot(ax, t, Y.u_1.linSS/1e6, '-.', 'Color', cl(4,:));
ylabel(ax,'P_T (MW)'); grid(ax,'on'); axis(ax,'tight');
legend(ax,'WFSim', ...
    sprintf('V1 input (%.0f%%)', R.VAFval_lin(1)), sprintf('inputs only (%.0f%%)', R.VAFval_lin(2)), ...
    sprintf('inputs + offset (%.0f%%)', R.VAFval_lin(3)), sprintf('inputs + offset, x_0 = x_{ss} (%.0f%%)', R.VAFval_linSS(3)), ...
    'Location','best');
title(ax,'Linear model (as in the MPC), VAF in brackets');
ax = subplot(3,1,3);
plot(ax, t, (Y.V1.lin - D.val.P(1:end-1))/1e6, '--', 'Color', cl(1,:)); hold(ax,'on');
plot(ax, t, (Y.u_1.lin - D.val.P(1:end-1))/1e6, '-', 'Color', cl(3,:));
plot(ax, t, (Y.u_1.linSS - D.val.P(1:end-1))/1e6, '-.', 'Color', cl(4,:));
ylabel(ax,'error (MW)'); xlabel(ax,'time (s)'); grid(ax,'on'); axis(ax,'tight');
legend(ax,'V1 input','inputs + offset','inputs + offset, x_0 = x_{ss}','Location','best');
set(fig,'Position',[100 100 700 700]);
set(fig,'PaperUnits','centimeters','PaperPositionMode','manual','PaperSize',[17 17],'PaperPosition',[0 0 17 17]);
print(fig, fullfile(dirFig,'PowerInputsOnly_val'), '-dpng', '-r200');
print(fig, fullfile(dirFig,'PowerInputsOnly_val'), '-depsc', '-painters', '-loose');
save(fullfile(dirFig,'PowerInputsOnly.mat'), 'R', 'Y', 'D');
end

function [S, Vinf] = prepareSet(T, t0, nTx, nTy, Wp, idxKeep, n)
% prepareSet lifted states (K_psi) and power of one data set from t0 on
kk = length(T.u)/Wp.mesh.Ny;
utemp = reshape(T.u, Wp.mesh.Nx, Wp.mesh.Ny, kk);
V1 = squeeze(utemp(nTx,nTy,:))';
Vinf = utemp(1,1,t0);
f = @(x) x(idxKeep);
Ur = [f(T.Ur1); f(T.Ur2)];
inp = [f(T.Ct1); f(T.Ct2); f(T.phi1)];
V1 = f(V1);
idx = t0:numel(Ur(1,:));
Ur = Ur(:,idx); inp = inp(:,idx);
S.V1 = V1(idx);
S.P12 = [f(T.PT1); f(T.PT2)]; S.P12 = S.P12(:,idx); % turbine powers (W)
S.P = sum(S.P12, 1);
structPrev = struct('Ur1_prev1',Ur(1,1),'Ur2_prev1',Ur(2,1),'dUr1_prev1',0,'dUr2_prev1',0, ...
    'M1',Ur(1,1),'M2',Ur(2,1),'k',1);
S.lifted = nan(n+3, numel(idx));
for k = 1:numel(idx)
    [S.lifted(:,k),S.stateName,structPrev] = K_psi([Ur(:,k); inp(:,k)], 0, n, structPrev, 1);
end
end

function e = rmse(y, yhat)
e = sqrt(mean((y(:) - yhat(:)).^2))/1e6;
end
