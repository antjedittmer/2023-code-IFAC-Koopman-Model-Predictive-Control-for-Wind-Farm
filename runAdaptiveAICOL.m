% Open-loop WFSim data for the adaptive AIC without the wind speed V1
% (thrust control only, yaw 0, settings as step 1 of mainGeneratePlots):
% - identification: random CT_prime levels at constant V_inf = 7 m/s
% - validation: random CT_prime levels, V_inf ramp from 7 to 9 m/s with
%   0.1 m/s per sample at t = 5000 s (as the closed-loop test case)
addpath('WFSimCode'); addpath('KoopmanIODMD');
R = 1e-6; refstairs = 0; measured = 1; KoopmanStates = 6; PolyLiftingFunction = 0;
controller = 0; vinfstr = '';

WFSim_demo(R,refstairs,measured,KoopmanStates,PolyLiftingFunction,controller, ...
    'sowfa_2turb_ct_levels',7,vinfstr);
close all

VinfRamp = windRamp(7, 9, 0.1, 5000, 10001);
WFSim_demo(R,refstairs,measured,KoopmanStates,PolyLiftingFunction,controller, ...
    'sowfa_2turb_ct_levels_val',VinfRamp,vinfstr);
close all

function V = windRamp(V0, V1, rate, kStart, N)
% Free stream wind V0 up to sample kStart, then a ramp to V1 with rate per sample
V = min(V0 + rate*max((1:N) - kStart, 0), V1);
end
