function [x0, urS] = koopFeatState(model, Ur, urS)
% koopFeatState initial state of the input-only Koopman power model
% (predictPowerFeat) from the measured rotor effective wind speeds
%
% States (model.stateSet, KoopmanIODMD/testPowerStates):
%   'Ur n=4/6/8'   first n of Ur1, Ur2, Ur1^2, Ur2^2, Ur1^3, Ur2^3, DUr1, DUr2
%                  (K_psi; DUr = Ur - moving mean of the past Ur over 25 s)
%   'Ur n=..+d'    plus the previous Ur1(k-1), Ur2(k-1)
% Models without stateSet (exportPowerFeatureModel): 'Ur n=<model.n>'
%
% Inputs
% - model: struct of the power model
% - Ur: current rotor effective wind speeds [Ur1; Ur2]
% - urS: past Ur and moving mean (struct M, UrPrev, k as structPrev of K_psi)
%
% Outputs
% - x0: initial state
% - urS: updated with Ur (call once per time step for the next step)

stateSet = sprintf('Ur n=%d', model.n);
if isfield(model,'stateSet')
    stateSet = model.stateSet;
end
nL = sscanf(stateSet, 'Ur n=%d');

% Moving mean of the past Ur as in K_psi
l = min(25, urS.k);
M = (l-1)/l*urS.M + urS.UrPrev/l;
DUr = Ur - M;
xAll = [Ur(1); Ur(2); Ur(1)^2; Ur(2)^2; Ur(1)^3; Ur(2)^3; DUr(1); DUr(2)];
x0 = xAll(1:nL);
if endsWith(stateSet, '+d')
    x0 = [x0; urS.UrPrev];
end
urS = struct('M', M, 'UrPrev', Ur, 'k', urS.k + 1);
end
