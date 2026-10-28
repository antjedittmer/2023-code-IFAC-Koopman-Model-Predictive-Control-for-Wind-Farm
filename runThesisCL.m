% Closed-loop runs for the thesis: the three published cases of
% mainGeneratePlots (step 3) with the measured wind V1, and the same
% operating modes with the input-only Koopman power model 'Ur n=6'
% (KoopFeatModel3.mat, runKoopStatesCL) without wind in the power model.
% Each run saves <figure name>_plotData.mat; repeated calls only replot
% (plotWFSimTimeseriesTex, figures <figure name>_Tex).
addpath('WFSimCode'); addpath('KoopmanIODMD');
R = 1e-6; refstairs = -1; measured = 0; controller = 2; Vinf = 8;
ControlSetStr = 'sowfa_2turb_alm_turbl'; PolyLiftingFunction = 0; KoopmanStates = 6;
vinfStr = ''; kTimestep = 2; sol_update = 0; Rvec = [1,1,0.001]; RU = 100;

% {name, OptYaw, KoopAIC, KoopFeat}
cases = {'Baseline AIC (paper)',          0, 1, 0
         'AIC with static yaw (paper)',   1, 1, 0
         'WRC-AIC (paper)',               1, 0, 0
         'Baseline AIC, no wind',         0, 0, 3
         'AIC with static yaw, no wind',  2, 0, 3
         'WRC-AIC, no wind',              1, 0, 3};
res = table();
for idxC = 1:size(cases,1)
    [~,JR,JQ] = WFSim_demo(R,refstairs,measured,KoopmanStates,PolyLiftingFunction,...
        controller,ControlSetStr,Vinf,vinfStr,kTimestep,sol_update,cases{idxC,2},cases{idxC,3},Rvec,RU,cases{idxC,4});
    close all
    res = [res; table(string(cases{idxC,1}), JQ*1e3, JR, 'VariableNames', {'Case','TE_kW','AA'})]; %#ok<AGROW>
    disp(res(end,:))
end
disp(res)
