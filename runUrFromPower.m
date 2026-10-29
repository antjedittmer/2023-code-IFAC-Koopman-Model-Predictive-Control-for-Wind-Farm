% runUrFromPower closed loop with the rotor effective wind speeds estimated
% from the power (urFromPower) instead of measured, for
% - the adaptive AIC without V1 (nu = 0.15 m/s, h = 5, ramp 7 -> 9 m/s, as
%   mainGeneratePlotThesisAdaptiveAIC, step 4)
% - the WRC-AIC (8 m/s, as mainGeneratePlotThesis, step 4)
% Requires the models of mainGeneratePlotThesis and
% mainGeneratePlotThesisAdaptiveAIC. Figures (png, eps) and the table
% (csv) are saved in figUrFromPower
mainDir = fileparts(mfilename('fullpath'));
cd(mainDir);
addpath('WFSimCode'); addpath('KoopmanIODMD'); addpath(fullfile('WFSimCode','bin','Kmpc'));
dirOut = fullfile(mainDir,'figUrFromPower');
if ~isfolder(dirOut), mkdir(dirOut); end

R = 1e-6; measured = 0; KoopmanStates = 6; PolyLiftingFunction = 0; controller = 2;
vinfStr = ''; kTimestep = 2; Rvec = [1,1,0.001]; RU = 100;
res = table();
for UrFromPower = [0 1]
    UrStr = {'measured','fromPower'};

    % Adaptive AIC without V1 (default settings nu = 0.15 m/s, h = 5)
    clear qLPVMPCcontroller_dU_quadprog % reset the persistent buffers of the controller
    rng(0); % excitation signal of the adaptive scheme
    Vinf = min(7 + 0.1*max((1:5001) - 2000, 0), 9);
    [~,JR,JQ] = WFSim_demo(R,20,measured,KoopmanStates,PolyLiftingFunction,controller, ...
        'sowfa_2turb_aic_5000',Vinf,vinfStr,kTimestep,1,0,1,Rvec,RU,0,1,[],UrFromPower);
    close all;
    res = [res; table("Adaptive AIC", string(UrStr{UrFromPower+1}), JQ*1e3, JR, ...
        'VariableNames', {'Case','Ur','TE_kW','AA'})]; %#ok<AGROW>
    f = sprintf('Stairs20R1dot0e_neg06_AICnoV1%s_Update1_Tex', repmat('_UrP',1,UrFromPower));
    copyFig(fullfile(mainDir,'DataT2OLWFSim','Vinf7dot0To9dot0_CL_K06_P0',f), ...
        fullfile(dirOut,sprintf('CL_adaptiveAIC_Ur_%s',UrStr{UrFromPower+1})));

    % WRC-AIC (KoopFeat = 3: input-only power model 'Ur n=6')
    clear qLPVMPCcontroller_dU_quadprog
    [~,JR,JQ] = WFSim_demo(R,-1,measured,KoopmanStates,PolyLiftingFunction,controller, ...
        'sowfa_2turb_alm_turbl',8,vinfStr,kTimestep,0,1,0,Rvec,RU,3,0,[],UrFromPower);
    close all;
    res = [res; table("WRC-AIC", string(UrStr{UrFromPower+1}), JQ*1e3, JR, ...
        'VariableNames', {'Case','Ur','TE_kW','AA'})]; %#ok<AGROW>
    f = sprintf('R1dot0e_neg06_KoopAIC0_KoopFeat3%s_Update0_Tex', repmat('_UrP',1,UrFromPower));
    copyFig(fullfile(mainDir,'DataT2OLWFSim','Vinf8dot0_CL_K06_P0',f), ...
        fullfile(dirOut,sprintf('CL_WRCAIC_Ur_%s',UrStr{UrFromPower+1})));
end
disp(res)
writetable(res, fullfile(dirOut,'UrFromPower.csv'));

function copyFig(src, dst)
% copyFig copies the png and eps of a figure, if they exist
for ext = {'.png','.eps'}
    if isfile([src,ext{1}])
        copyfile([src,ext{1}], [dst,ext{1}]);
    else
        warning('Figure not found: %s', [src,ext{1}]);
    end
end
end
