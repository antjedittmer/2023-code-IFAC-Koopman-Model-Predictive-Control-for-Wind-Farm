function [TE, AA] = evalClosedLoop(pd, t0)
% evalClosedLoop tracking error and actuator activity of a closed-loop run
% from sample t0 on (to the end of the simulation), with the formulas of
% plotWFSimTimeseriesLatex2 (TE, from sample 60 at the earliest) and
% saveDataEDMD (AA)
%
% Inputs
% - pd: plot data of WFSim_demo (Pref, Power, CT_prime, Phi)
% - t0: first sample of the evaluation
%
% Outputs
% - TE: mean absolute power tracking error (MW, as JQ)
% - AA: actuator activity (as JR)

e = pd.Pref(:) - sum(pd.Power)';
TE = mean(abs(e(max(60,t0):end)))/10^6;
AA = mean(mean(abs(diff(pd.CT_prime(:,t0:end)'))))*2 + mean(mean(abs(diff(pd.Phi(:,t0:end)'))))*0.1 + ...
    mean(mean(abs(pd.Phi(:,t0:end)')))*10^(-3);
end
