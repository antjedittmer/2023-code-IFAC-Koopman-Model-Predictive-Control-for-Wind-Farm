% Open loop WFSim runs for the identification and validation sets with
% random CT_prime and yaw levels (settings as step 1 of mainGeneratePlots)
addpath('WFSimCode'); addpath('KoopmanIODMD');
R = 1e-6; refstairs = 0; measured = 1; KoopmanStates = 6; PolyLiftingFunction = 0;
controller = 0; Vinf = 8; vinfstr = '';
for ControlSetStr = {'sowfa_2turb_yaw_ct_levels', 'sowfa_2turb_yaw_ct_levels_val'}
    WFSim_demo(R,refstairs,measured,KoopmanStates,PolyLiftingFunction,controller,ControlSetStr{1},Vinf,vinfstr);
    close all
end
