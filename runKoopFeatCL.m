% Closed loop with the Koopman MPC as the third run of step 3 in
% mainGeneratePlots (KoopAIC = 0, OptYaw = 1), power prediction with the
% input-only Koopman power model with input features (KoopFeat = 1)
addpath('WFSimCode'); addpath('KoopmanIODMD');
R = 1e-6; refstairs = -1; measured = 0; controller = 2; Vinf = 8;
ControlSetStr = 'sowfa_2turb_alm_turbl'; PolyLiftingFunction = 0; KoopmanStates = 6;
vinfStr = ''; kTimestep = 2; sol_update = 0; OptYaw = 1; KoopAIC = 0;
Rvec = [1,1,0.001]; RU = 100; KoopFeat = 1;
[~,JR,JQ] = WFSim_demo(R,refstairs,measured,KoopmanStates,PolyLiftingFunction,...
    controller,ControlSetStr,Vinf,vinfStr,kTimestep,sol_update,OptYaw,KoopAIC,Rvec,RU,KoopFeat);
fprintf('KoopFeat closed loop: JR = %g, JQ = %g\n', JR, JQ);
