% mainGeneratePlotThesis creates the thesis plots of the Koopman MPC
% without the wind speed V1 (input-only Koopman power model 'Ur n=6',
% KoopFeatModel3.mat):
%   - Open-loop identification: diagnostic plots of Ur1, Ur2 and the farm
%     power on the three validation sets (plotPowerStatesDiag)
%   - Closed loop: Baseline AIC, AIC with static yaw and WRC-AIC
%     (plotWFSimTimeseriesTex)
% Steps whose results exist are not repeated (open-loop data, Koopman
% models, closed-loop plot data <figure name>_plotData.mat), so a second
% call only replots. All figures are copied into the folder figThesis.
clc; clear; close all;
addpath('WFSimCode');
addpath('KoopmanIODMD');
mainDir = pwd;
dirModels = fullfile(mainDir,'DataInOutWfSim');
dirInputsOnly = fullfile(dirModels,'PowerInputsOnly');

%% 1. Create or load open loop data with WFSim (as mainGeneratePlots)
% Identification (noise_step, ct_levels) and validation sets (ct_levels_val,
% steps, steps_Ct_comb) of the input-only model, AllComb for the model in
% step 2; the run is skipped if the data exist
R = 1e-6; refstairs = 0; measured = 1; KoopmanStates = 6; PolyLiftingFunction = 0;
controller = 0; Vinf = 8; vinfstr = '';
ControlSetStrs = {'sowfa_2turb_yaw_noise_step','sowfa_2turb_yaw_ct_levels', ...
    'sowfa_2turb_yaw_ct_levels_val','sowfa_2turb_yaw_steps','sowfa_2turb_yaw_steps_Ct_comb', ...
    'sowfa_2turb_alm_turbl_AllComb'};
for idx = 1:numel(ControlSetStrs)
    WFSim_demo(R,refstairs,measured,KoopmanStates,PolyLiftingFunction,controller,ControlSetStrs{idx},Vinf,vinfstr);
    close all
end

%% 2. Koopman models
% MPCinitv loads the published Ur model (eDMDresults_UasOutput) and power
% model (eDMDresults_UasOutput_MIMO_PT) with V1 also for KoopFeat > 0; they
% are evaluated, but the applied input of the runs without wind (KoopAIC = 0)
% only depends on the input-only power model and the measured Ur
koopmanDirs = fullfile(dirModels,{'eDMDresults_UasOutput','eDMDresults_UasOutput_MIMO_PT'}, ...
    'Vin_Vinf8dot0_states06');
noStates = 6; PolyVec = 0; useVal = 1; percentTrain = .6; NTMstr = ''; t0 = 240;
if ~isfolder(koopmanDirs{1})
    filenameId = 'Vinf8dot0_sowfa_2turb_alm_turbl_AllComb';
    MainWfSimVin(0,filenameId,filenameId,noStates,PolyVec,useVal,percentTrain,NTMstr,t0);
    close all
end
if ~isfolder(koopmanDirs{2})
    MainWfSimVinYawPowerOut(1,'Vinf8dot0_sowfa_2turb_yaw_noise_step','Vinf8dot0_sowfa_2turb_yaw_steps', ...
        noStates,PolyVec,useVal,percentTrain,NTMstr,t0);
    close all
end

% Input-only power models (as runKoopStatesCL): states Ur n=4, 6, 8, 4+d,
% 6+d saved as KoopFeatModel2..6.mat; the thesis uses 'Ur n=6' (KoopFeat = 3)
idN = {'Vinf8dot0_sowfa_2turb_yaw_noise_step','Vinf8dot0_sowfa_2turb_yaw_ct_levels'};
valN = {'Vinf8dot0_sowfa_2turb_yaw_ct_levels_val','Vinf8dot0_sowfa_2turb_yaw_steps', ...
    'Vinf8dot0_sowfa_2turb_yaw_steps_Ct_comb'};
stateSets = {'Ur n=4','Ur n=6','Ur n=8','Ur n=4+d','Ur n=6+d'};
KoopFeat = 3;
if ~isfile(fullfile(dirInputsOnly,sprintf('KoopFeatModel%d.mat',KoopFeat)))
    addpath(fullfile('WFSimCode','bin','Kmpc'));
    [~, Mdl] = testPowerStates(idN, valN, {'NL+1+M+x'}, 240, 100, stateSets);
    close all
    for idxS = 1:numel(Mdl)
        model = Mdl(idxS);
        model.t0 = 240;
        model.filenameId = idN;
        save(fullfile(dirInputsOnly,sprintf('KoopFeatModel%d.mat',idxS+1)), 'model');
    end
end

%% 3. Diagnostic plots of the open-loop identification (input-only model)
% Per validation set Ur1, Ur2 and the farm power: WFSim vs. model
VAFdiag = plotPowerStatesDiag(sprintf('KoopFeatModel%d.mat',KoopFeat), valN, 240);

%% 4. Closed loop without wind
R = 1e-6; refstairs = -1; measured = 0; controller = 2; Vinf = 8;
ControlSetStr = 'sowfa_2turb_alm_turbl'; PolyLiftingFunction = 0; KoopmanStates = 6;
vinfStr = ''; kTimestep = 2; sol_update = 0; Rvec = [1,1,0.001]; RU = 100;
UrFromPower = 1; % 1 (default): Ur estimated from the power (urFromPower), 0: measured Ur
urStr = repmat('_UrP',1,UrFromPower); % file name suffix of WFSim_demo

% {name, OptYaw, KoopAIC, KoopFeat}
cases = {'Baseline AIC',          0, 0, KoopFeat
         'AIC static yaw',   2, 0, KoopFeat
         'WRC-AIC',               1, 0, KoopFeat};
resCL = table();
for idxC = 1:size(cases,1)
    [~,JR,JQ] = WFSim_demo(R,refstairs,measured,KoopmanStates,PolyLiftingFunction,...
        controller,ControlSetStr,Vinf,vinfStr,kTimestep,sol_update,cases{idxC,2},cases{idxC,3},Rvec,RU,cases{idxC,4},0,[],UrFromPower);
    fig = gcf; fig.Name = cases{idxC,1};
    resCL = [resCL; table(string(cases{idxC,1}), JQ*1e3, JR, 'VariableNames', {'Case','TE_kW','AA'})]; %#ok<AGROW>
end
disp(VAFdiag)
disp(resCL)

%% Copy png and eps in one folder
thesisDir = fullfile(mainDir,'figThesis');
if ~isfolder(thesisDir)
    mkdir(thesisDir);
end
origDirCL = fullfile(mainDir,'DataT2OLWFSim','Vinf8dot0_CL_K06_P0');
filenamesCL = {sprintf('R1dot0e_neg06_KoopAIC0_OptYaw0_KoopFeat%d%s_Update0_Tex',KoopFeat,urStr), ...
    sprintf('R1dot0e_neg06_KoopAIC0_OptYaw2_KoopFeat%d%s_Update0_Tex',KoopFeat,urStr), ...
    sprintf('R1dot0e_neg06_KoopAIC0_KoopFeat%d%s_Update0_Tex',KoopFeat,urStr)};
filetype = {'.png','.eps'};
for idx = 1:length(filenamesCL)
    for idxT = 1:length(filetype)
        copyfile(fullfile(origDirCL,[filenamesCL{idx},filetype{idxT}]), thesisDir);
    end
end
for diagDir = {sprintf('Diag_KoopFeatModel%d_Ur',KoopFeat), sprintf('Diag_KoopFeatModel%d_P',KoopFeat)}
    for idxT = 1:length(filetype)
        files = dir(fullfile(dirInputsOnly,diagDir{1},['*',filetype{idxT}]));
        for idx = 1:numel(files) % prefix Ur/P: same file names in both folders
            copyfile(fullfile(files(idx).folder,files(idx).name), ...
                fullfile(thesisDir,[diagDir{1},'_',files(idx).name]));
        end
    end
end

%% Copy png and eps with the file names of the thesis folder figures/figures_IFAC23
% (the folder can be copied to the thesis as it is)
ifacDir = fullfile(mainDir,'figures_IFAC23');
if ~isfolder(ifacDir)
    mkdir(ifacDir);
end
% Identification and validation set of the WRC/AIC diagnostic plots
idSetWRC = 'yaw_noise_step'; % identification data described in the thesis
valSetWRC = 'yaw_ct_levels_val';
diagFile = @(type, setName, var) fullfile(dirInputsOnly, sprintf('Diag_KoopFeatModel%d_%s',KoopFeat,var), ...
    ['Val', type, '_Yaw_', strrep(setName,'yaw_',''), '__n_K', ...
    sprintf('__%d__U_r___U_r_2___U_r_3__', KoopmanStates)]); % file name of plotEDMDinputsEffWind
% {source, thesis file name}
figsIFAC = {fullfile(origDirCL,filenamesCL{1}), '04FigTimeSeries0_tex'
            fullfile(origDirCL,filenamesCL{2}), '05FigTimeSeriesQLMPC_tex'
            fullfile(origDirCL,filenamesCL{3}), '06FigTimeSeriesLinMPC_tex'
            diagFile('Id', idSetWRC, 'Ur'), '03aFig_WRC_Id_Vr'
            diagFile('Id', idSetWRC, 'P'), '03aFig_WRC_Id_P'
            diagFile('Val', valSetWRC, 'Ur'), '03aFig_WRC_Val_Vr'
            diagFile('Val', valSetWRC, 'P'), '03aFig_WRC_Val_P'};
for idx = 1:size(figsIFAC,1)
    for idxT = 1:length(filetype)
        copyfile([figsIFAC{idx,1},filetype{idxT}], fullfile(ifacDir,[figsIFAC{idx,2},filetype{idxT}]));
    end
end
