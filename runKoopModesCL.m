% Closed loop with the Koopman MPC (KoopAIC = 0) in three operating modes:
% pure AIC (OptYaw = 0), AIC with the static yaw of findMaxYaw (OptYaw = 2)
% and AIC + WRC with the yaw optimized by the MPC (OptYaw = 1), each with
% the Koopman power model with V1 (KoopFeat = 0) and the input-only model
% with input features (KoopFeat = 1). Other settings as step 3 of mainGeneratePlots.
addpath('WFSimCode'); addpath('KoopmanIODMD');
R = 1e-6; refstairs = -1; measured = 0; controller = 2; Vinf = 8;
ControlSetStr = 'sowfa_2turb_alm_turbl'; PolyLiftingFunction = 0; KoopmanStates = 6;
vinfStr = ''; kTimestep = 2; sol_update = 0; KoopAIC = 0; Rvec = [1,1,0.001]; RU = 100;
modes = [0 2 1]; modeNames = {'pure AIC','AIC + static yaw','AIC + WRC'};
res = table();
for KoopFeat = [0 1]
    for idxM = 1:numel(modes)
        [~,JR,JQ] = WFSim_demo(R,refstairs,measured,KoopmanStates,PolyLiftingFunction,...
            controller,ControlSetStr,Vinf,vinfStr,kTimestep,sol_update,modes(idxM),KoopAIC,Rvec,RU,KoopFeat);
        close all
        res = [res; table(KoopFeat, string(modeNames{idxM}), JQ*1e3, JR, ...
            'VariableNames', {'KoopFeat','Mode','TE_kW','AA'})]; %#ok<AGROW>
        disp(res(end,:))
    end
end
disp(res)
save(fullfile('DataInOutWfSim','PowerInputsOnly','KoopModesCL.mat'), 'res');
