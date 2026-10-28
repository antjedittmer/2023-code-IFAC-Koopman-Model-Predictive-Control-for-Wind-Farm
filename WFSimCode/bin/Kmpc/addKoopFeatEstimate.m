function pd = addKoopFeatEstimate(pd, dirModels, nAhead)
% addKoopFeatEstimate adds the n-step-ahead estimates of the input-only
% Koopman power model (koopFeatEstimate) to the closed-loop plot data pd of
% a run with KoopFeat > 0: V_est (Vr1, Vr2) and P_est (farm power, W). The
% logged V_est of these runs comes from the Ur model with V1, which the
% controller does not use, so it is replaced.
%
% Inputs
% - pd: plot data of WFSim_demo (V, CT_prime, Phi, KoopFeat)
% - dirModels: folder DataInOutWfSim with PowerInputsOnly/KoopFeatModel<k>.mat
% - nAhead: steps ahead (default 10 = MPC horizon mpc.Nh, the predictions
%   the control decisions are based on)

if nargin < 3
    nAhead = 10;
end
if ~isfield(pd,'KoopFeat') || pd.KoopFeat == 0 || ...
        (isfield(pd,'nAhead') && pd.nAhead == nAhead)
    return;
end
featFile = 'KoopFeatModel.mat'; % KoopFeat = 1 (as MPCinitv)
if pd.KoopFeat > 1
    featFile = sprintf('KoopFeatModel%d.mat', pd.KoopFeat);
end
tmp = load(fullfile(dirModels,'PowerInputsOnly',featFile),'model');
[pd.V_est, pd.P_est] = koopFeatEstimate(tmp.model, pd.V, [pd.CT_prime; pd.Phi(1,:)], nAhead);
pd.nAhead = nAhead;
end
