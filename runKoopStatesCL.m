% Input-only Koopman power models with different lifted Ur states in the
% Koopman MPC (KoopAIC = 0): identification and validation (VAF) with
% testPowerStates, closed loop in the three operating modes of
% runKoopModesCL (pure AIC, AIC + static yaw, AIC + WRC).
% The models are saved as KoopFeatModel<k>.mat, k = 2..6 (WFSim_demo with KoopFeat = k).
addpath('WFSimCode'); addpath('KoopmanIODMD'); addpath(fullfile('WFSimCode','bin','Kmpc'));
dirOut = fullfile('DataInOutWfSim','PowerInputsOnly');

%% Identification and validation
idN = {'Vinf8dot0_sowfa_2turb_yaw_noise_step','Vinf8dot0_sowfa_2turb_yaw_ct_levels'};
valN = {'Vinf8dot0_sowfa_2turb_yaw_ct_levels_val','Vinf8dot0_sowfa_2turb_yaw_steps','Vinf8dot0_sowfa_2turb_yaw_steps_Ct_comb'};
stateSets = {'Ur n=4','Ur n=6','Ur n=8','Ur n=4+d','Ur n=6+d'};
[Tval, Mdl] = testPowerStates(idN, valN, {'NL+1+M+x'}, 240, 100, stateSets);
close all
for idxS = 1:numel(Mdl)
    model = Mdl(idxS);
    model.t0 = 240;
    model.filenameId = idN;
    save(fullfile(dirOut,sprintf('KoopFeatModel%d.mat',idxS+1)), 'model');
end

%% Closed loop (settings as runKoopModesCL)
R = 1e-6; refstairs = -1; measured = 0; controller = 2; Vinf = 8;
ControlSetStr = 'sowfa_2turb_alm_turbl'; PolyLiftingFunction = 0; KoopmanStates = 6;
vinfStr = ''; kTimestep = 2; sol_update = 0; KoopAIC = 0; Rvec = [1,1,0.001]; RU = 100;
modes = [0 2 1]; modeNames = {'pure AIC','AIC + static yaw','AIC + WRC'};
res = table();
for idxS = 1:numel(stateSets)
    for idxM = 1:numel(modes)
        [~,JR,JQ] = WFSim_demo(R,refstairs,measured,KoopmanStates,PolyLiftingFunction,...
            controller,ControlSetStr,Vinf,vinfStr,kTimestep,sol_update,modes(idxM),KoopAIC,Rvec,RU,idxS+1);
        close all
        res = [res; table(string(stateSets{idxS}), string(modeNames{idxM}), JQ*1e3, JR, ...
            'VariableNames', {'States','Mode','TE_kW','AA'})]; %#ok<AGROW>
        disp(res(end,:))
        save(fullfile(dirOut,'KoopStatesCL.mat'), 'res', 'Tval');
    end
end
disp(Tval)
disp(res)
save(fullfile(dirOut,'KoopStatesCL.mat'), 'res', 'Tval');
