% Closed-loop test case of the adaptive AIC (thesis Section 5.2) with the
% Koopman effective wind speed model without V1 (inputs CT1, CT2 only, no
% yaw): 5000 s, V_inf ramp from 7 to 9 m/s with 0.1 m/s per sample at
% k = 2000, P_ref step from 1.72 MW to 2.58 MW at k = 2500 (refstairs = 20).
% Non-adaptive (sol_update = 0) and adaptive (sol_update = 1) qLMPC with
% the estimated effective wind speeds (KoopAIC = 1, OptYaw = 0).
addpath('WFSimCode'); addpath('KoopmanIODMD');
if ~isfile(fullfile('DataInOutWfSim','AdaptiveAIC','KoopAICnoV1.mat'))
    identifyAICnoV1();
end
R = 1e-6; refstairs = 20; measured = 0; KoopmanStates = 6; PolyLiftingFunction = 0;
controller = 2; ControlSetStr = 'sowfa_2turb_aic_5000'; vinfStr = ''; kTimestep = 2;
OptYaw = 0; KoopAIC = 1; Rvec = [1,1,0.001]; RU = 100; KoopFeat = 0; AICnoV1 = 1;
Vinf = min(7 + 0.1*max((1:5001) - 2000, 0), 9);

res = table();
for sol_update = [0 1]
    clear qLPVMPCcontroller_dU_quadprog % reset the persistent buffers of the controller
    rng(0); % excitation signal of the adaptive scheme
    try
        [~,JR,JQ] = WFSim_demo(R,refstairs,measured,KoopmanStates,PolyLiftingFunction,controller, ...
            ControlSetStr,Vinf,vinfStr,kTimestep,sol_update,OptYaw,KoopAIC,Rvec,RU,KoopFeat,AICnoV1);
    catch ME % the non-adaptive estimate without V1 diverges after the wind ramp
        fprintf('sol_update = %d aborted: %s\n', sol_update, ME.message);
        JR = NaN; JQ = NaN;
    end
    res = [res; table(sol_update, JQ*1e3, JR, 'VariableNames', {'Adaptive','TE_kW','AA'})]; %#ok<AGROW>
    disp(res(end,:))
end
disp(res)
