% mainGeneratePlotThesis creates the thesis plots of the Koopman MPC
% without the wind speed V1 for three input-only Koopman power models
% (modelCases): the control-affine model 'Ur n=6' with the nonlinear
% features of the filtered inputs as states and the three controls ct1,
% ct2, phi1 as inputs (KoopFeatModel9.mat, figure names as before), 'Ur n=6'
% with the three controls only (KoopFeatModel7.mat, figure names + '_n6')
% and 'Ur n=6' with 13 input features incl. moving means (KoopFeatModel3.mat,
% figure names + '_nl'):
%   - Open-loop identification: VAF of Ur1, Ur2 and the farm power and
%     diagnostic plots of the control-affine model with the estimates of
%     the 3-input model as grey dashed lines (plotPowerStatesDiag)
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

% Input-only Koopman power model of the three closed loops (Baseline AIC,
% AIC static yaw, WRC-AIC): x(k+1) = A x(k) + B w(k), P_WF(k) = C x(k) + D w(k)
% States (state set, see testPowerStates), e.g.
%   'Ur n=6': the six states Ur1, Ur2, Ur1^2, Ur2^2, Ur1^3, Ur2^3
%   'Ur n=4', 'Ur n=8', 'Ur n=4+d', 'Ur n=6+d'
% Inputs w (feature set, see powerFeatures), e.g.
%   'NL+1+M+x': ct1, ct2, phi1, 1, two nonlinear power features, four
%     moving means of the inputs (tau = 100 s, i.e. four further filter
%     states) and three products (13 inputs)
%   'NL+1': ct1, ct2, phi1, 1 and the two nonlinear power features
%   'u': the three controls ct1, ct2, phi1 (one output P_WF)
%   'NL+1+M+x,aff': control-affine model with the features of 'NL+1+M+x'
%     as functions of the inputs filtered as in the turbine model (three
%     further filter states) and the three controls as the only inputs
% Each combination is saved as KoopFeatModel<KoopFeat>.mat ({state set,
% feature set, KoopFeat}); a missing model is identified
wrcModels = {'Ur n=4',   'NL+1+M+x', 2
             'Ur n=6',   'NL+1+M+x', 3 % 6 + 4 filter states, 13 inputs
             'Ur n=8',   'NL+1+M+x', 4
             'Ur n=4+d', 'NL+1+M+x', 5
             'Ur n=6+d', 'NL+1+M+x', 6
             'Ur n=6',   'u',        7 % 6 states, 3 inputs
             'Ur n=6',   'NL+1',     8
             'Ur n=6',   'NL+1+M+x,aff', 9}; % 6 + 3 + 4 filter states, 3 inputs (control affine)
% Cases of the loop: {state set, feature set, suffix of the figure names,
% case whose estimates are shown as grey dashed lines in the diagnostic
% plots (0: none), diagnostic plots (1) or VAF only (0)}
modelCases = {'Ur n=6', 'NL+1+M+x,aff', '',    2, 1
              'Ur n=6', 'u',            '_n6', 0, 0
              'Ur n=6', 'NL+1+M+x',     '_nl', 0, 0};
% names of the cases in the diagnostic plots (TeX, char(966): varphi, not supported by the TeX interpreter)
refNames = {['n_', char(966), ' = 13'], ['n_', char(966), ' = 6'], ['n_', char(966), ' = 10']};
tableNames = {'13 states, 3 inputs (control affine)', '6 states, 3 inputs', '10 states, 13 inputs'}; % names of the cases in the table
tableOrder = [2 3 1]; % rows of the table (cases)
showCase = 1; % only the figures of this case stay open
idN = {'Vinf8dot0_sowfa_2turb_yaw_noise_step','Vinf8dot0_sowfa_2turb_yaw_ct_levels'};
valN = {'Vinf8dot0_sowfa_2turb_yaw_ct_levels_val','Vinf8dot0_sowfa_2turb_yaw_steps', ...
    'Vinf8dot0_sowfa_2turb_yaw_steps_Ct_comb'};
% Data sets of the comparison table (VAF: ct_levels, ct_levels_val,
% steps_Ct_comb) and of the diagnostic plots (by default ct_levels and
% ct_levels_val); allDiagSets = 1 plots all identification and validation sets
allDiagSets = 0;
idPlot = {'Vinf8dot0_sowfa_2turb_yaw_ct_levels'};
valPlot = {'Vinf8dot0_sowfa_2turb_yaw_ct_levels_val','Vinf8dot0_sowfa_2turb_yaw_steps_Ct_comb'};
setsPlot = {'Vinf8dot0_sowfa_2turb_yaw_ct_levels','Vinf8dot0_sowfa_2turb_yaw_ct_levels_val'};
if allDiagSets
    idPlot = idN;
    valPlot = valN;
    setsPlot = [idN, valN];
end
addpath(fullfile('WFSimCode','bin','Kmpc'));

% Closed loop settings (step 4)
R = 1e-6; refstairs = -1; measured = 0; controller = 2; Vinf = 8;
ControlSetStr = 'sowfa_2turb_alm_turbl'; PolyLiftingFunction = 0; KoopmanStates = 6;
vinfStr = ''; kTimestep = 2; sol_update = 0; Rvec = [1,1,0.001]; RU = 100;
UrFromPower = 1; % 1 (default): Ur estimated from the power (urFromPower), 0: measured Ur
urStr = repmat('_UrP',1,UrFromPower); % file name suffix of WFSim_demo

% Folders of the figures: all in figThesis; thesis file names in figures_IFAC23
% (the folder can be copied to the thesis as it is)
thesisDir = fullfile(mainDir,'figThesis');
ifacDir = fullfile(mainDir,'figures_IFAC23');
for figDir = {thesisDir, ifacDir}
    if ~isfolder(figDir{1})
        mkdir(figDir{1});
    end
end
origDirCL = fullfile(mainDir,'DataT2OLWFSim','Vinf8dot0_CL_K06_P0');
filetype = {'.png','.eps'};
% Identification and validation set of the WRC/AIC diagnostic plots
idSetWRC = 'yaw_ct_levels'; % identification set with random CT and yaw levels (runCtLevelsOL)
valSetWRC = 'yaw_ct_levels_val'; % validation set of the same type (runCtLevelsOL)

VAFall = table();
resCL = table();
% Koopman models of the cases (identified if missing)
KoopFeatCases = zeros(size(modelCases,1),1);
for idxMC = 1:size(modelCases,1)
    [stateSetWRC, featureSetWRC] = modelCases{idxMC,1:2};
    idxM = find(strcmp(wrcModels(:,1),stateSetWRC) & strcmp(wrcModels(:,2),featureSetWRC));
    if isempty(idxM) % new combination: next free number
        wrcModels(end+1,:) = {stateSetWRC, featureSetWRC, max([wrcModels{:,3}]) + 1}; %#ok<SAGROW>
        idxM = size(wrcModels,1);
    end
    KoopFeat = wrcModels{idxM,3};
    KoopFeatCases(idxMC) = KoopFeat;
    if ~isfile(fullfile(dirInputsOnly,sprintf('KoopFeatModel%d.mat',KoopFeat)))
        [~, Mdl] = testPowerStates(idN, valN, {featureSetWRC}, 240, 100, {stateSetWRC});
        close all
        model = Mdl(1);
        model.t0 = 240;
        model.filenameId = idN;
        save(fullfile(dirInputsOnly,sprintf('KoopFeatModel%d.mat',KoopFeat)), 'model');
    end
    fprintf('Power model KoopFeatModel%d.mat: states %s, inputs %s\n', KoopFeat, stateSetWRC, featureSetWRC);
end
for idxMC = 1:size(modelCases,1)
    [stateSetWRC, featureSetWRC, figSuffix, refCase, diagPlots] = modelCases{idxMC,:};
    KoopFeat = KoopFeatCases(idxMC);

    % 3. Diagnostic plots of the open-loop identification (input-only model)
    % Per validation set Ur1, Ur2 and the farm power: WFSim vs. model
    % (estimates of the reference case as grey lines; old plots removed)
    refFile = ''; refName = '';
    if refCase > 0
        refFile = sprintf('KoopFeatModel%d.mat',KoopFeatCases(refCase));
        refName = refNames{refCase};
    end
    for var = {'Ur','P'}
        delete(fullfile(dirInputsOnly,sprintf('Diag_KoopFeatModel%d_%s',KoopFeat,var{1}),'*.*'));
    end
    setsPlotMC = setsPlot;
    if ~diagPlots, setsPlotMC = {}; end
    VAFdiag = plotPowerStatesDiag(sprintf('KoopFeatModel%d.mat',KoopFeat), valPlot, 240, refFile, refName, idPlot, ...
        refNames{idxMC}, setsPlotMC);
    VAFdiag.Model = repmat(string(featureSetWRC), height(VAFdiag), 1);
    VAFall = [VAFall; VAFdiag]; %#ok<AGROW>

    % 4. Closed loop without wind: {name, OptYaw, KoopAIC, KoopFeat}
    cases = {'Baseline AIC',          0, 0, KoopFeat
             'AIC static yaw',   2, 0, KoopFeat
             'WRC-AIC',               1, 0, KoopFeat};
    for idxC = 1:size(cases,1)
        [~,JR,JQ] = WFSim_demo(R,refstairs,measured,KoopmanStates,PolyLiftingFunction,...
            controller,ControlSetStr,Vinf,vinfStr,kTimestep,sol_update,cases{idxC,2},cases{idxC,3},Rvec,RU,cases{idxC,4},0,[],UrFromPower);
        fig = gcf; fig.Name = [cases{idxC,1},figSuffix];
        resCL = [resCL; table(string(featureSetWRC), string(cases{idxC,1}), JQ*1e3, JR, NaN, NaN, ...
            'VariableNames', {'Model','Case','TE_kW','AA','TE240_kW','AA240'})]; %#ok<AGROW>
    end

    % Copy png and eps in one folder (figure names + figSuffix)
    filenamesCL = {sprintf('R1dot0e_neg06_KoopAIC0_OptYaw0_KoopFeat%d%s_Update0_Tex',KoopFeat,urStr), ...
        sprintf('R1dot0e_neg06_KoopAIC0_OptYaw2_KoopFeat%d%s_Update0_Tex',KoopFeat,urStr), ...
        sprintf('R1dot0e_neg06_KoopAIC0_KoopFeat%d%s_Update0_Tex',KoopFeat,urStr)};
    % TE and AA from sample 240 on (after the start-up of the wake); the
    % plots and their titles: from 500 s (WFSim_demo, as JQ, JR above)
    for idx = 1:length(filenamesCL)
        tmp = load(fullfile(origDirCL,[strrep(filenamesCL{idx},'_Tex',''),'_plotData.mat']),'pd');
        [TE240, AA240] = evalClosedLoop(tmp.pd, 240);
        iRow = height(resCL) - length(filenamesCL) + idx; % rows of this model, order of cases
        resCL.TE240_kW(iRow) = TE240*1e3;
        resCL.AA240(iRow) = AA240;
    end
    for idx = 1:length(filenamesCL)
        for idxT = 1:length(filetype)
            copyfile(fullfile(origDirCL,[filenamesCL{idx},filetype{idxT}]), ...
                fullfile(thesisDir,[filenamesCL{idx},figSuffix,filetype{idxT}]));
        end
    end
    for diagDir = {sprintf('Diag_KoopFeatModel%d_Ur',KoopFeat), sprintf('Diag_KoopFeatModel%d_P',KoopFeat)}
        for idxT = 1:length(filetype)
            files = dir(fullfile(dirInputsOnly,diagDir{1},['*',filetype{idxT}]));
            for idx = 1:numel(files) % prefix Ur/P: same file names in both folders
                [~, nameF] = fileparts(files(idx).name);
                copyfile(fullfile(files(idx).folder,files(idx).name), ...
                    fullfile(thesisDir,[diagDir{1},'_',nameF,figSuffix,filetype{idxT}]));
            end
        end
    end

    % Copy png and eps with the file names of the thesis folder figures/figures_IFAC23
    % (file name of plotEDMDinputsEffWind from the title: without the states
    % if a reference model is shown)
    diagSuffix = sprintf('__n_K__%d__U_r___U_r_2___U_r_3__', KoopmanStates);
    if refCase > 0, diagSuffix = ''; end
    diagFile = @(type, setName, var) fullfile(dirInputsOnly, sprintf('Diag_KoopFeatModel%d_%s',KoopFeat,var), ...
        ['Val', type, '_Yaw_', strrep(setName,'yaw_',''), diagSuffix]);
    % {source, thesis file name}
    figsIFAC = {fullfile(origDirCL,filenamesCL{1}), '04FigTimeSeries0_tex'
                fullfile(origDirCL,filenamesCL{2}), '05FigTimeSeriesQLMPC_tex'
                fullfile(origDirCL,filenamesCL{3}), '06FigTimeSeriesLinMPC_tex'
                diagFile('Id', idSetWRC, 'Ur'), '03aFig_WRC_Id_Vr'
                diagFile('Id', idSetWRC, 'P'), '03aFig_WRC_Id_P'
                diagFile('Val', valSetWRC, 'Ur'), '03aFig_WRC_Val_Vr'
                diagFile('Val', valSetWRC, 'P'), '03aFig_WRC_Val_P'};
    for idx = 1:size(figsIFAC,1)
        if ~diagPlots && ~isfile([figsIFAC{idx,1},filetype{1}]) % no diagnostic plots of this case
            continue
        end
        for idxT = 1:length(filetype)
            copyfile([figsIFAC{idx,1},filetype{idxT}], fullfile(ifacDir,[figsIFAC{idx,2},figSuffix,filetype{idxT}]));
        end
    end
    if idxMC ~= showCase % keep the diagnostic plots (names 'Ur ..', 'P ..') open
        figs = findobj('Type','figure');
        close(figs(~startsWith(get(figs,'Name'),{'Ur ','P '})));
    end
end
disp(VAFall)
disp(resCL)

% LaTeX table of the comparison (figures_IFAC23/tab_WRC_ModelComparison.tex)
maxEig = zeros(numel(tableOrder),1);
for idx = 1:numel(tableOrder)
    tmp = load(fullfile(dirInputsOnly,sprintf('KoopFeatModel%d.mat',KoopFeatCases(tableOrder(idx)))),'model');
    maxEig(idx) = max(abs(eig(tmp.model.A)));
end
writeWRCModelTable(VAFall, resCL, modelCases(tableOrder,2), tableNames(tableOrder), maxEig, ...
    fullfile(ifacDir,'tab_WRC_ModelComparison.tex'));
