function z0 = koopFeatLinState(model, x0, zFilt)
% koopFeatLinState initial state of the linear power model (feature set
% 'NL+1+M+x,lin'): lifted rotor effective wind speeds and the feature
% states, calculated exactly from the filter states
%
% Inputs
% - model: struct of the power model (tau: time constant of the moving means)
% - x0: lifted Ur (koopFeatState, 6 x 1)
% - zFilt: filter states [moving means (4); filtered inputs s (3)], as the
%   second output of powerFeatures (mpc.featM)
%
% Output
% - z0: [x0; feature states (featStatesLin)] (21 x 1)

s0 = zFilt(5:7);
W0 = powerFeatures(s0, model.featureSet, model.tau, zFilt); % first column: features of s0, M0
z0 = [x0(:); featStatesLin(W0)];
end
