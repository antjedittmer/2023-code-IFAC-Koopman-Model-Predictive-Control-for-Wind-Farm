% mainGeneratePlotThesisAdaptiveAIC creates the thesis plots of the
% adaptive AIC without the wind speed V1 (thesis Section 5.2 test case,
% effective wind speed model with the inputs CT1, CT2 only, no yaw):
%   - Open loop: validation of the non-adaptive and the online updated
%     model on the data set with the wind ramp (validateAICnoV1) and the
%     singular values of the updated models compared with the model
%     identified offline at 9 m/s (plotAICnoV1Sigma); test of the offline
%     and the learned models on unseen data at 9 m/s (testAICnoV1Vinf9)
%   - Closed loop: non-adaptive and adaptive qLMPC, 5000 s, V_inf ramp
%     from 7 to 9 m/s with 0.1 m/s per sample at k = 2000, P_ref step from
%     1.72 MW to 2.58 MW at k = 2500 (plotWFSimTimeseriesTex)
% Steps whose results exist are not repeated (open-loop data, Koopman
% model, closed-loop plot data <figure name>_plotData.mat), so a second
% call only replots. The figures are copied into figThesisAdaptiveAIC.
clc; clear; close all;
addpath('WFSimCode');
addpath('KoopmanIODMD');
addpath(fullfile('WFSimCode','bin','Kmpc'));
warning('off','MATLAB:nearlySingularMatrix'); % online updates with little excitation
mainDir = pwd;

%% 1. Open-loop data with WFSim (as runAdaptiveAICOL)
% Identification: random CT_prime levels at constant V_inf = 7 m/s;
% validation: random CT_prime levels, ramp from 7 to 9 m/s at t = 5000 s
R = 1e-6; refstairs = 0; measured = 1; KoopmanStates = 6; PolyLiftingFunction = 0;
controller = 0; vinfstr = '';
if ~isfile(fullfile(mainDir,'DataT2OLWFSim','Vinf7dot0_OL_Ct','Vinf7dot0_sowfa_2turb_ct_levels.mat'))
    WFSim_demo(R,refstairs,measured,KoopmanStates,PolyLiftingFunction,controller, ...
        'sowfa_2turb_ct_levels',7,vinfstr);
    close all
end
if ~isfile(fullfile(mainDir,'DataT2OLWFSim','Vinf9dot0_OL_Ct','Vinf9dot0_sowfa_2turb_ct_levels.mat'))
    % offline model at 9 m/s for the comparison of the singular values
    WFSim_demo(R,refstairs,measured,KoopmanStates,PolyLiftingFunction,controller, ...
        'sowfa_2turb_ct_levels',9,vinfstr);
    close all
end
if ~isfile(fullfile(mainDir,'DataT2OLWFSim','Vinf9dot0_OL_Ct','Vinf9dot0_sowfa_2turb_ct_levels_val.mat'))
    % unseen test data at 9 m/s (validation CT levels)
    WFSim_demo(R,refstairs,measured,KoopmanStates,PolyLiftingFunction,controller, ...
        'sowfa_2turb_ct_levels_val',9,vinfstr);
    close all
end
if isempty(dir(fullfile(mainDir,'DataT2OLWFSim','Vinf7dot0To9dot0_OL_Ct','*sowfa_2turb_ct_levels_val.mat')))
    VinfRamp = min(7 + 0.1*max((1:10001) - 5000, 0), 9);
    WFSim_demo(R,refstairs,measured,KoopmanStates,PolyLiftingFunction,controller, ...
        'sowfa_2turb_ct_levels_val',VinfRamp,vinfstr);
    close all
end

%% 2. Koopman effective wind speed model without V1 (inputs CT1, CT2)
if ~isfile(fullfile(mainDir,'DataInOutWfSim','AdaptiveAIC','KoopAICnoV1.mat'))
    identifyAICnoV1();
end
if ~isfile(fullfile(mainDir,'DataInOutWfSim','AdaptiveAIC','KoopAICnoV1_Vinf9.mat'))
    identifyAICnoV1('Vinf9dot0_sowfa_2turb_ct_levels', [], 'KoopAICnoV1_Vinf9.mat');
end

%% 3. Open-loop validation: non-adaptive and online updated model,
% singular values of the updated models vs. the offline model at 9 m/s
cd KoopmanIODMD
[~, sigmaSel] = plotAICnoV1Sigma([], [], [], [], [240 5100 6000 8000]); % models at these times (s)
VAFol = validateAICnoV1([], [], [], sigmaSel); % updates of the shown models in color
test9 = testAICnoV1Vinf9();
cd(mainDir)

%% 4. Closed loop: non-adaptive and adaptive AIC without V1
R = 1e-6; refstairs = 20; measured = 0; KoopmanStates = 6; PolyLiftingFunction = 0;
controller = 2; ControlSetStr = 'sowfa_2turb_aic_5000'; vinfStr = ''; kTimestep = 2;
OptYaw = 0; KoopAIC = 1; Rvec = [1,1,0.001]; RU = 100; KoopFeat = 0; AICnoV1 = 1;
Vinf = min(7 + 0.1*max((1:5001) - 2000, 0), 9);
resCL = table();
for sol_update = [0 1]
    clear qLPVMPCcontroller_dU_quadprog % reset the persistent buffers of the controller
    rng(0); % excitation signal of the adaptive scheme
    [~,JR,JQ] = WFSim_demo(R,refstairs,measured,KoopmanStates,PolyLiftingFunction,controller, ...
        ControlSetStr,Vinf,vinfStr,kTimestep,sol_update,OptYaw,KoopAIC,Rvec,RU,KoopFeat,AICnoV1);
    resCL = [resCL; table(sol_update, JQ*1e3, JR, 'VariableNames', {'Adaptive','TE_kW','AA'})]; %#ok<AGROW>
end
disp(VAFol)
disp(resCL)

%% Copy png and eps in one folder (short file names)
thesisDir = fullfile(mainDir,'figThesisAdaptiveAIC');
if ~isfolder(thesisDir)
    mkdir(thesisDir);
end
dirOL = fullfile(mainDir,'DataInOutWfSim','AdaptiveAIC');
dirCL = fullfile(mainDir,'DataT2OLWFSim','Vinf7dot0To9dot0_CL_K06_P0');
% {source, target name}
figs = {fullfile(dirOL,'ValPVal_Non_adaptive'), 'AdaptiveAIC_OL_NonAdaptive'
        fullfile(dirOL,'ValPVal_Adaptive'), 'AdaptiveAIC_OL_Adaptive'
        fullfile(dirOL,'AdaptiveAIC_OL_Sigma'), 'AdaptiveAIC_OL_Sigma'
        fullfile(dirOL,'AdaptiveAIC_Test9'), 'AdaptiveAIC_OL_Test9'
        fullfile(dirCL,'Stairs20R1dot0e_neg06_AICnoV1_Update0_Tex'), 'AdaptiveAIC_CL_NonAdaptive'
        fullfile(dirCL,'Stairs20R1dot0e_neg06_AICnoV1_Update1_Tex'), 'AdaptiveAIC_CL_Adaptive'};
filetype = {'.png','.eps'};
for idx = 1:size(figs,1)
    for idxT = 1:length(filetype)
        copyfile([figs{idx,1},filetype{idxT}], fullfile(thesisDir,[figs{idx,2},filetype{idxT}]));
    end
end
