function [UrEst, upd, err, modelEnd, modelsAt] = simAdaptiveKoopman(model, Ur, W, adaptive, ad, saveAt)
% simAdaptiveKoopman open-loop estimation of the effective wind speeds with
% the Koopman model as in the MPC (wfmodelNL): the model is simulated from
% the measured initial state, the predicted Ur are lifted again in each
% step (K_psi). With adaptive = true, the model is updated online with the
% measured Ur (koopmanAdaptStep), as in the controller.
%
% Inputs
% - model: struct with A, B, L, Lu (identifyAICnoV1)
% - Ur: measured effective wind speeds (2 x N)
% - W: inputs (nw x N), here [CT1; CT2]
% - adaptive: online update on/off
% - ad: settings of the update (h, nu, maxEig), default h = 5, nu = 0.15 m/s;
%   ad.resetState = true: after an update the estimator continues from the
%   measured Ur (default false: free run as in wfmodelNL)
% - saveAt: samples after which the model is returned in modelsAt (optional)
%
% Outputs
% - UrEst: estimated effective wind speeds (2 x N)
% - upd: update flags (1 x N)
% - err: estimation error of the last h samples (1 x N)
% - modelEnd: model after the last sample
% - modelsAt: models (fields A, B, k) after the samples saveAt

if nargin < 5, ad = []; end
if nargin < 6, saveAt = []; end
modelsAt = struct('A',{},'B',{},'k',{});
if ismember(1, saveAt) % initial model
    modelsAt(1) = struct('A',model.A,'B',model.B,'k',1);
end
resetState = isstruct(ad) && isfield(ad,'resetState') && ad.resetState;
N = size(Ur,2);
UrEst = nan(2,N); UrEst(:,1) = Ur(:,1);
upd = zeros(1,N); err = zeros(1,N);
x = lift(Ur(:,1));
for k = 2:N
    xn = model.A*x + model.B*W(:,k-1);
    UrEst(:,k) = xn(1:2);
    x = lift(xn(1:2));
    if adaptive
        [model, ad, upd(k)] = koopmanAdaptStep(model, ad, Ur(:,k), W(:,k), UrEst(:,k));
        err(k) = ad.err;
        if upd(k) && resetState
            x = lift(Ur(:,k));
        end
    end
    if ismember(k, saveAt)
        modelsAt(end+1) = struct('A',model.A,'B',model.B,'k',k); %#ok<AGROW>
    end
end
modelEnd = model;
end

function phi = lift(Ur)
phi = [Ur; Ur.^2; Ur.^3];
end
