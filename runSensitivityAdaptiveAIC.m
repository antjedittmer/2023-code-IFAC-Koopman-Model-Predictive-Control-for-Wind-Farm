% runSensitivityAdaptiveAIC sensitivity of the adaptive AIC without the
% wind speed V1 (effective wind speed model with the inputs CT1, CT2) to
% the settings of the online update (koopmanAdaptStep): trigger threshold
% nu and window/batch length h. For each combination:
%   - Open loop on the validation data with the wind ramp from 7 to 9 m/s
%     at t = 5000 s (as validateAICnoV1): number of updates, VAF, RMSE
%     after the ramp, and the learned model (frozen after the run) on the
%     unseen test data at constant 9 m/s (as testAICnoV1Vinf9)
%   - Closed loop with the qLMPC, V_inf ramp from 7 to 9 m/s at k = 2000
%     (as mainGeneratePlotThesisAdaptiveAIC): tracking error TE, actuator
%     activity AA, number of updates, VAF of the estimates
% The non-adaptive model is added as reference.
%
% Outputs in figSensitivityAdaptiveAIC:
%   - OL_nu<100*nu>_h<h>.png/.eps, CL_nu<100*nu>_h<h>.png/.eps (e.g.
%     OL_nu015_h05: nu = 0.15 m/s, h = 5) and OL_nonAdaptive,
%     CL_nonAdaptive, so that all figures can be compared together
%   - sensitivityAdaptiveAIC.mat/.csv: all results
%   - tab_sensitivityAdaptiveAIC.tex: LaTeX table (booktabs, adjustbox)
% The open-loop data and the Koopman models are created first with
% mainGeneratePlotThesisAdaptiveAIC if they do not exist. Closed-loop runs
% whose plot data exist are not repeated (WFSim_demo), so a second call
% only recreates the figures and the table.
clc; clear; close all;
mainDir = fileparts(mfilename('fullpath'));
cd(mainDir);
dirOL = fullfile(mainDir,'DataT2OLWFSim');
dirAd = fullfile(mainDir,'DataInOutWfSim','AdaptiveAIC');
needed = {fullfile(dirOL,'Vinf9dot0_OL_Ct','Vinf9dot0_sowfa_2turb_ct_levels_val.mat'), ...
    fullfile(dirAd,'KoopAICnoV1.mat'), fullfile(dirAd,'KoopAICnoV1_Vinf9.mat')};
if ~all(cellfun(@isfile, needed)) || ...
        isempty(dir(fullfile(dirOL,'Vinf7dot0To9dot0_OL_Ct','*sowfa_2turb_ct_levels_val.mat')))
    mainGeneratePlotThesisAdaptiveAIC; % open-loop data, models (clears the workspace)
    close all;
    mainDir = fileparts(mfilename('fullpath'));
    cd(mainDir);
end
addpath('WFSimCode');
addpath('KoopmanIODMD');
addpath(fullfile('WFSimCode','bin','Kmpc'));
koopmanDir = fullfile(mainDir,'KoopmanIODMD');
addpath(fullfile(koopmanDir,'2.DYNAMIC_MODE_DECOMPOSITION'), genpath(fullfile(koopmanDir,'Functions')));
warning('off','MATLAB:nearlySingularMatrix'); % online updates with little excitation

nuVec = [0.3 0.15 0.1 0.05 0.02]; % trigger threshold (m/s)
hVec = [5 10 20]; % window/batch length (samples)
figDir = fullfile(mainDir,'figSensitivityAdaptiveAIC');
if ~isfolder(figDir), mkdir(figDir); end
tmpDir = fullfile(figDir,'tmp');
if ~isfolder(tmpDir), mkdir(tmpDir); end

%% Data and models
dirAd = fullfile(mainDir,'DataInOutWfSim','AdaptiveAIC');
tmp = load(fullfile(dirAd,'KoopAICnoV1.mat'),'model'); model = tmp.model;
t0 = 240; % wake established (validateAICnoV1)
f = dir(fullfile(mainDir,'DataT2OLWFSim','Vinf7dot0To9dot0_OL_Ct','*sowfa_2turb_ct_levels_val.mat'));
S = load(fullfile(f(1).folder,f(1).name),'Ur1','Ur2','Ct1','Ct2','u');
k = t0:numel(S.Ur1);
Ur = [S.Ur1(k); S.Ur2(k)];
W = [S.Ct1(k); S.Ct2(k)];
afterStep = k >= 5000;
% V_inf and V1 (u in the grid cell in front of WT1) for the plot only
Wp.mesh = struct('Lx',1882.1,'Ly',800.0,'Nx',50,'Ny',25);
nTx = floor(0.4e3/Wp.mesh.Lx*Wp.mesh.Nx);
nTy = round(400/Wp.mesh.Ly*Wp.mesh.Ny);
utemp = reshape(S.u, Wp.mesh.Nx, Wp.mesh.Ny, []);
V1 = squeeze(utemp(nTx,nTy,k))';
Vinf = squeeze(utemp(1,nTy,k))';
% Unseen test data at constant 9 m/s
S9 = load(fullfile(mainDir,'DataT2OLWFSim','Vinf9dot0_OL_Ct', ...
    'Vinf9dot0_sowfa_2turb_ct_levels_val.mat'),'Ur1','Ur2','Ct1','Ct2');
k9 = t0:numel(S9.Ur1);
Ur9 = [S9.Ur1(k9); S9.Ur2(k9)];
W9 = [S9.Ct1(k9); S9.Ct2(k9)];
rmse = @(y,yh) sqrt(mean((y - yh).^2, 2))';

%% Closed-loop settings (as mainGeneratePlotThesisAdaptiveAIC)
R = 1e-6; refstairs = 20; measured = 0; KoopmanStates = 6; PolyLiftingFunction = 0;
controller = 2; ControlSetStr = 'sowfa_2turb_aic_5000'; vinfStr = ''; kTimestep = 2;
OptYaw = 0; KoopAIC = 1; Rvec = [1,1,0.001]; RU = 100; KoopFeat = 0; AICnoV1 = 1;
VinfCL = min(7 + 0.1*max((1:5001) - 2000, 0), 9);
dirCL = fullfile(mainDir,'DataT2OLWFSim','Vinf7dot0To9dot0_CL_K06_P0');

%% Runs: non-adaptive (nu = inf) and all combinations of nu and h
[hGrid, nuGrid] = ndgrid(hVec, nuVec); % nu as given, h ascending
cases = [inf, NaN; nuGrid(:), hGrid(:)];
res = table();
for idx = 1:size(cases,1)
    nu = cases(idx,1); h = cases(idx,2);
    adaptive = ~isinf(nu);
    if adaptive
        ad = struct('h',h,'nu',nu,'resetState',true,'trigChannels',1);
        tag = sprintf('nu%03d_h%02d', round(100*nu), h);
        strVal = sprintf('$\\nu$ = %.2f m/s, $h$ = %d', nu, h);
    else
        ad = [];
        tag = 'nonAdaptive';
        strVal = 'Non-adaptive';
    end

    % Open loop: ramp data, then the learned model on the 9 m/s test data
    [UrEst, upd, ~, modelEnd] = simAdaptiveKoopman(model, Ur, W, adaptive, ad);
    vafOL = [vaf(Ur(1,:), UrEst(1,:)), vaf(Ur(2,:), UrEst(2,:))];
    rmseOL = rmse(Ur(:,afterStep), UrEst(:,afterStep));
    UrEst9 = simAdaptiveKoopman(modelEnd, Ur9, W9, false);
    rmse9 = rmse(Ur9, UrEst9);
    fig = plotEDMDinputsEffWindUpdate(UrEst, [W; V1], Ur, max(vafOL,0), tmpDir, 1, model.n, ...
        sprintf('%s (%d upd.)', strVal, sum(upd)), 1 + upd, 1, Vinf, k, 'tex', []);
    moveFigure(tmpDir, figDir, ['OL_',tag]);
    close(fig);

    % Closed loop (plot data are reused if they exist)
    clear qLPVMPCcontroller_dU_quadprog % reset the persistent buffers of the controller
    rng(0); % excitation signal of the adaptive scheme
    WFSim_demo(R,refstairs,measured,KoopmanStates,PolyLiftingFunction,controller, ...
        ControlSetStr,VinfCL,vinfStr,kTimestep,adaptive,OptYaw,KoopAIC,Rvec,RU,KoopFeat,AICnoV1,ad);
    baseCL = 'Stairs20R1dot0e_neg06_AICnoV1';
    if adaptive, baseCL = sprintf('%s_nu%03d_h%d', baseCL, round(100*nu), h); end
    baseCL = sprintf('%s_Update%d', baseCL, adaptive);
    tmp = load(fullfile(dirCL,[baseCL,'_plotData.mat']),'pd'); pd = tmp.pd;
    for ext = {'.png','.eps'}
        copyfile(fullfile(dirCL,[baseCL,'_Tex',ext{1}]), fullfile(figDir,['CL_',tag,ext{1}]));
    end
    close all;
    t = pd.time;
    updCL = pd.update(t) > 0;

    res = [res; table(nu, h, ...
        sum(upd), vafOL(1), vafOL(2), rmseOL(1), rmseOL(2), max(abs(eig(modelEnd.A))), ...
        vaf(Ur9(1,:),UrEst9(1,:)), vaf(Ur9(2,:),UrEst9(2,:)), rmse9(1), rmse9(2), ...
        sum(updCL), max([t(updCL), NaN]), pd.JQ*1e3, pd.JR, ...
        vaf(pd.V(1,t),pd.V_est(1,t)), vaf(pd.V(2,t),pd.V_est(2,t)), ...
        'VariableNames', {'nu','h', ...
        'OL_NoUpd','OL_VAF_Vr1','OL_VAF_Vr2','OL_RMSE_Vr1_afterStep','OL_RMSE_Vr2_afterStep','OL_maxEigA', ...
        'Test9_VAF_Vr1','Test9_VAF_Vr2','Test9_RMSE_Vr1','Test9_RMSE_Vr2', ...
        'CL_NoUpd','CL_tLastUpd','CL_TE_kW','CL_AA','CL_VAF_Vr1','CL_VAF_Vr2'})]; %#ok<AGROW>
    fprintf('%s: OL %d upd., CL TE = %.1f kW, AA = %.4g\n', tag, sum(upd), pd.JQ*1e3, pd.JR);
end
rmdir(tmpDir,'s');

%% Results: mat, csv and LaTeX table
save(fullfile(figDir,'sensitivityAdaptiveAIC.mat'),'res');
writetable(res, fullfile(figDir,'sensitivityAdaptiveAIC.csv'));
writeLatexTable(res, fullfile(figDir,'tab_sensitivityAdaptiveAIC.tex'));
format short g
disp(res)

function moveFigure(srcDir, dstDir, name)
% moveFigure moves the .png and .eps of the only figure in srcDir to
% dstDir/name (plotEDMDinputsEffWindUpdate names the files after the title)
for ext = {'.png','.eps'}
    f = dir(fullfile(srcDir,['*',ext{1}]));
    movefile(fullfile(f(1).folder,f(1).name), fullfile(dstDir,[name,ext{1}]));
end
end

function writeLatexTable(res, fileName)
% writeLatexTable sensitivity table (booktabs, adjustbox): VAF <= 0 as
% $\leq 0$, a diverging estimate of the frozen learned model as div.
fid = fopen(fileName,'w');
c = onCleanup(@() fclose(fid));
L = {['% Created by runSensitivityAdaptiveAIC on ', char(datetime('now','Format','yyyy-MM-dd'))]
    '% Adaptive AIC without V1: trigger threshold nu and window/batch length h of the online update'
    '\begin{table}[htbp]'
    '    \centering'
    ['    \caption[Sensitivity of the adaptive AIC to the update settings]{Sensitivity of the adaptive AIC to the ', ...
    'trigger threshold $\nu$ and the window and batch length $h$ of the online update. Open loop: validation data ', ...
    'with the wind ramp from 7 to 9\,m/s at $t = 5000$\,s (from 240\,s), RMSE after the ramp; the model learned ', ...
    'on these data is frozen and tested on unseen data at constant 9\,m/s. Closed loop: qLMPC with the ramp from 7 to ', ...
    '9\,m/s at $k = 2000$; TE tracking error, AA actuator activity. $N_\text{upd}$: number of model updates.}']
    '    \label{tab:sensitivityAdaptiveAIC}'
    '    \begin{adjustbox}{max width=\textwidth}'
    '    \begin{tabular}{rr rrrrr rr rrrr}'
    '        \toprule'
    ['        & & \multicolumn{5}{c}{Open loop, wind ramp} & \multicolumn{2}{c}{Test 9\,m/s} & ', ...
    '\multicolumn{4}{c}{Closed loop} \\']
    '        \cmidrule(lr){3-7} \cmidrule(lr){8-9} \cmidrule(lr){10-13}'
    ['        $\nu$ (m/s) & $h$ & $N_\text{upd}$ & VAF $V_{r1}$ (\%) & VAF $V_{r2}$ (\%) & ', ...
    'RMSE $V_{r1}$ (m/s) & RMSE $V_{r2}$ (m/s) & RMSE $V_{r1}$ (m/s) & RMSE $V_{r2}$ (m/s) & ', ...
    '$N_\text{upd}$ & TE (kW) & AA ($10^{-3}$) & VAF $V_{r2}$ (\%) \\']
    '        \midrule'};
for idx = 1:height(res)
    r = res(idx,:);
    if isinf(r.nu)
        first = '\multicolumn{2}{l}{non-adaptive}';
    else
        first = [sprintf('%.2f', r.nu), ' & ', sprintf('%d', r.h)];
    end
    cells = {first, sprintf('%d', r.OL_NoUpd), fmtVaf(r.OL_VAF_Vr1), fmtVaf(r.OL_VAF_Vr2), ...
        fmtNum(r.OL_RMSE_Vr1_afterStep,'%.3f'), fmtNum(r.OL_RMSE_Vr2_afterStep,'%.2f'), ...
        fmtNum(r.Test9_RMSE_Vr1,'%.3f'), fmtNum(r.Test9_RMSE_Vr2,'%.2f'), ...
        sprintf('%d', r.CL_NoUpd), sprintf('%.1f', r.CL_TE_kW), sprintf('%.2f', 1e3*r.CL_AA), fmtVaf(r.CL_VAF_Vr2)};
    line = ['        ', strjoin(cells, ' & '), ' \\'];
    if idx == 1 || (idx < height(res) && res.nu(idx+1) ~= r.nu)
        line = [line, ' \midrule']; %#ok<AGROW>
    end
    L{end+1,1} = line; %#ok<AGROW>
end
L = [L; {'        \bottomrule'; '    \end{tabular}'; '    \end{adjustbox}'; '\end{table}'}];
fprintf(fid, '%s\n', L{:});
end

function s = fmtVaf(v)
% fmtVaf VAF in %, <= 0 as $\leq 0$
if v <= 0
    s = '$\leq 0$';
else
    s = sprintf('%.1f', v);
end
end

function s = fmtNum(v, fmt)
% fmtNum number, NaN/Inf (diverging estimate) as div.
if ~isfinite(v)
    s = 'div.';
else
    s = sprintf(fmt, v);
end
end
