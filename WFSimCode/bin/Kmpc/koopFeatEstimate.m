function [VrEst, PEst] = koopFeatEstimate(model, Vr, U, nAhead)
% koopFeatEstimate n-step-ahead estimates of the input-only Koopman power
% model (predictPowerFeat) in closed loop: as in the MPC, the state is
% initialized from the measured rotor effective wind speeds (koopFeatState)
% and the moving means of the features are continued with the applied
% inputs (powerFeatures); the state is then propagated n steps
%   x(k|k-n) = A^n x(k-n) + sum A^(n-1-j) B w(k-n+j),  x(k-n) = koopFeatState(Vr(k-n))
%   VrEst(k) = x(k|k-n)(1:2),  PEst(k) = C x(k|k-n) + D w(k)
% With nAhead = mpc.Nh this is the prediction at the end of the MPC horizon
% (with the applied instead of the planned inputs).
%
% Inputs
% - model: struct of the power model (KoopFeatModel<k>.mat)
% - Vr: measured rotor effective wind speeds [Ur1; Ur2] (2 x K)
% - U: applied inputs [ct1; ct2; phi1] (3 x K)
% - nAhead: number of steps ahead (default 1)
%
% Outputs
% - VrEst: estimated Vr (2 x K, NaN in the first nAhead steps)
% - PEst: estimated farm power (W, 1 x K, NaN in the first nAhead steps)

if nargin < 4
    nAhead = 1;
end
K = size(Vr,2);
W = powerFeatures(U, model.featureSet, model.tau); % steady state before k = 1
urS = struct('M', Vr(:,1), 'UrPrev', Vr(:,1), 'k', 1);
VrEst = NaN(2,K);
PEst = NaN(1,K);
for k = 1:K-nAhead
    [x, urS] = koopFeatState(model, Vr(:,k), urS);
    for j = 0:nAhead-1
        x = model.A*x + model.B*W(:,k+j);
    end
    VrEst(:,k+nAhead) = x(1:2);
    PEst(k+nAhead) = (model.C*x + model.D*W(:,k+nAhead))*1e6;
end
end
