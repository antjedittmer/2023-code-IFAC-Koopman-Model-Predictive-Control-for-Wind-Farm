function mainGeneratePlotThesisAdaptiveAIC(opts)
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
% call only replots. The figures are copied into figThesisAdaptiveAIC and
% figures_CDC22.
%
% Name-value inputs: settings of the online update (koopmanAdaptStep)
% - 'h': window/batch length (samples), default 5
% - 'nu': trigger threshold (m/s), default 0.15
% With other settings than the default, only the figures that depend on
% them are created: open-loop adaptive validation, singular values, test
% at 9 m/s and adaptive closed loop. Their titles also show nu and h. In
% figThesisAdaptiveAIC and figures_CDC22 they get the suffix
% _nu<%2.2f>_h<%d>, e.g. Fig03_OL_Bode_tex_nu0.10_h5.png; in the data
% folders the suffix _nu<100*nu>_h<h> (_nu010_h5, as WFSim_demo).
% - 'showTitle': true shows, false hides the titles of all figures (title
%   switch ShowThesisTitle of runGenerateAllThesisFigures, restored at the
%   end). Default []: unchanged, i.e. as set by runGenerateAllThesisFigures,
%   titles shown if the function is called on its own.
% Examples:
%   mainGeneratePlotThesisAdaptiveAIC                    % thesis figures
%   mainGeneratePlotThesisAdaptiveAIC('nu',0.1)          % nu = 0.10 m/s, h = 5
%   mainGeneratePlotThesisAdaptiveAIC('h',10,'nu',0.05)
%   mainGeneratePlotThesisAdaptiveAIC('nu',0.1,'showTitle',false) % no titles
arguments
    opts.h (1,1) double {mustBeInteger, mustBePositive} = 5
    opts.nu (1,1) double {mustBePositive} = 0.15
    opts.showTitle {mustBeScalarOrEmpty, mustBeNumericOrLogical} = []
end
if ~isempty(opts.showTitle)
    hadSwitch = isappdata(groot, 'ShowThesisTitle');
    oldSwitch = [];
    if hadSwitch, oldSwitch = getappdata(groot, 'ShowThesisTitle'); end
    setappdata(groot, 'ShowThesisTitle', logical(opts.showTitle));
    restoreSwitch = onCleanup(@() restoreTitleSwitch(hadSwitch, oldSwitch)); % also after an error
end
clc; close all;
mainDir = fileparts(mfilename('fullpath'));
cd(mainDir);
addpath('WFSimCode');
addpath('KoopmanIODMD');
addpath(fullfile('WFSimCode','bin','Kmpc'));
warning('off','MATLAB:nearlySingularMatrix'); % online updates with little excitation

% Settings of the online update and file name suffixes
isDefault = opts.h == 5 && abs(opts.nu - 0.15) < 1e-12;
adOL = struct('h',opts.h,'nu',opts.nu,'resetState',true,'trigChannels',1);
if isDefault
    adCL = []; % WFSim_demo: settings of MPCinitv, file names without suffix
    srcSfx = '';
    thesisSfx = '';
else
    adCL = adOL;
    srcSfx = sprintf('_nu%03d_h%d', round(100*opts.nu), opts.h); % as WFSim_demo
    thesisSfx = sprintf('_nu%2.2f_h%d', opts.nu, opts.h);
end

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
loadData_Id_Vinf7to9 = 1;
if loadData_Id_Vinf7to9 && ...
        ~isfile(fullfile(mainDir,'DataInOutWfSim','AdaptiveAICdbs','KoopAICnoV1.mat'))
    identifyAICnoV1();
end
if ~isfile(fullfile(mainDir,'DataInOutWfSim','AdaptiveAIC','KoopAICnoV1_Vinf9.mat'))
    identifyAICnoV1('Vinf9dot0_sowfa_2turb_ct_levels', [], 'KoopAICnoV1_Vinf9.mat');
end

%% 3. Open-loop validation: non-adaptive and online updated model,
% singular values of the updated models vs. the offline model at 9 m/s
% models at 240, 5100, 6000, 8000 s (aicNoV1ModelTimes), marked in color in
% the adaptive validation plot and shown with their singular values
plotAICnoV1Sigma([], [], [], adOL, [], srcSfx);
VAFol = validateAICnoV1([], [], adOL, aicNoV1ModelTimes(), srcSfx);
testAICnoV1Vinf9([], [], adOL, srcSfx);
if isDefault
    % 7 m/s model on its identification set and the validation set up to
    % the ramp (does not depend on the update settings)
    [VAFidVal, filesIdVal] = plotAICnoV1IdVal();
end
cd(mainDir)

%% 4. Closed loop: non-adaptive and adaptive AIC without V1
R = 1e-6; refstairs = 20; measured = 0; KoopmanStates = 6; PolyLiftingFunction = 0;
controller = 2; ControlSetStr = 'sowfa_2turb_aic_5000'; vinfStr = ''; kTimestep = 2;
OptYaw = 0; KoopAIC = 1; Rvec = [1,1,0.001]; RU = 100; KoopFeat = 0; AICnoV1 = 1;
Vinf = min(7 + 0.1*max((1:5001) - 2000, 0), 9);
UrFromPower = 1; % 1 (default): Ur estimated from the power (urFromPower), 0: measured Ur
urStr = repmat('_UrP',1,UrFromPower); % file name suffix of WFSim_demo
resCL = table();
% the non-adaptive run does not depend on the update settings
updates = [0 1];
if ~isDefault, updates = 1; end
for sol_update = updates
    clear qLPVMPCcontroller_dU_quadprog % reset the persistent buffers of the controller
    rng(0); % excitation signal of the adaptive scheme
    [~,JR,JQ] = WFSim_demo(R,refstairs,measured,KoopmanStates,PolyLiftingFunction,controller, ...
        ControlSetStr,Vinf,vinfStr,kTimestep,sol_update,OptYaw,KoopAIC,Rvec,RU,KoopFeat,AICnoV1, ...
        adCL,UrFromPower);
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
figCLAd = fullfile(dirCL,['Stairs20R1dot0e_neg06_AICnoV1',srcSfx,urStr,'_Update1_Tex']);
% {source, target name}; without the suffix only the figures that depend on
% the update settings
figs = {fullfile(dirOL,['ValPVal_Adaptive',srcSfx]), 'AdaptiveAIC_OL_Adaptive'
        fullfile(dirOL,['AdaptiveAIC_OL_Sigma',srcSfx]), 'AdaptiveAIC_OL_Sigma'
        fullfile(dirOL,['AdaptiveAIC_Test9',srcSfx]), 'AdaptiveAIC_OL_Test9'
        figCLAd, 'AdaptiveAIC_CL_Adaptive'};
if isDefault
    figs = [figs
            {fullfile(dirOL,'ValPVal_Non_adaptive'), 'AdaptiveAIC_OL_NonAdaptive'
             fullfile(dirCL,['Stairs20R1dot0e_neg06_AICnoV1',urStr,'_Update0_Tex']), 'AdaptiveAIC_CL_NonAdaptive'}];
end
copyFigures(figs, thesisDir, thesisSfx);

%% Copy png and eps with the file names of the thesis folder figures/figures_CDC22
% (the folder can be copied to the thesis as it is)
cdcDir = fullfile(mainDir,'figures_CDC22');
if ~isfolder(cdcDir)
    mkdir(cdcDir);
end
% {source, thesis file name}
figsCDC = {fullfile(dirOL,['ValPVal_Adaptive',srcSfx]), 'Fig02_OL_Update1_tex'
           fullfile(dirOL,['AdaptiveAIC_OL_Sigma',srcSfx]), 'Fig03_OL_Bode_tex'
           figCLAd, 'Fig05_CL_Update1_tex'};
if isDefault
    figsCDC = [figsCDC
               {filesIdVal{1}, 'Fig00_AIC7mpers_OL_id'
                filesIdVal{2}, 'Fig00_AIC7mpers_OL_val'
                fullfile(dirOL,'ValPVal_Non_adaptive'), 'Fig01_OL_Update0_tex'
                fullfile(dirCL,['Stairs20R1dot0e_neg06_AICnoV1',urStr,'_Update0_Tex']), 'Fig04_CL_Update0_tex'}];
end
copyFigures(figsCDC, cdcDir, thesisSfx);
if isDefault
    disp(VAFidVal)
end
end

function restoreTitleSwitch(hadSwitch, oldSwitch)
% restoreTitleSwitch title switch ShowThesisTitle as before the call
if hadSwitch
    setappdata(groot, 'ShowThesisTitle', oldSwitch);
elseif isappdata(groot, 'ShowThesisTitle')
    rmappdata(groot, 'ShowThesisTitle');
end
end

function copyFigures(figs, targetDir, sfx)
% copyFigures copies the png and eps of the figures {source, target name}
% into targetDir, with the suffix sfx after the target name
filetype = {'.png','.eps'};
for idx = 1:size(figs,1)
    for idxT = 1:length(filetype)
        copyfile([figs{idx,1},filetype{idxT}], fullfile(targetDir,[figs{idx,2},sfx,filetype{idxT}]));
    end
end
end
