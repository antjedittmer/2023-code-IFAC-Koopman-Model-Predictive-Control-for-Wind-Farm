function P = predictPowerFeat(model, Ur, Mprev, U)
% predictPowerFeat farm power (W) of the input-only Koopman power model
% (exportPowerFeatureModel) over the horizon of the inputs U
%
% x(k+1) = A x(k) + B w(k), P(k) = C x(k) + D w(k) (P in MW in the model)
% with the lifted rotor effective wind speeds x and the features w
% (powerFeatures) of the inputs u(k) = [ct1; ct2; phi1].
%
% Inputs
% - model: struct with A, B, C, D, n, featureSet, tau
% - Ur: initial state (koopFeatState) or the current rotor effective wind
%   speeds [Ur1; Ur2] (lifted with the first n observables of K_psi)
% - Mprev: moving means of the features before the first step
% - U: inputs [ct1; ct2; phi1] (3 x Nh)

x = Ur(:);
if numel(x) ~= model.n
    x = liftUr(Ur, model.n);
end
W = powerFeatures(U, model.featureSet, model.tau, Mprev);
Nh = size(U,2);
P = zeros(Nh,1);
for k = 1:Nh
    P(k) = model.C*x + model.D*W(:,k);
    x = model.A*x + model.B*W(:,k);
end
P = P*1e6;
end

function x = liftUr(Ur, n)
% First n observables of K_psi: Ur1, Ur2, Ur1^2, Ur2^2, Ur1^3, Ur2^3
xAll = [Ur(1); Ur(2); Ur(1)^2; Ur(2)^2; Ur(1)^3; Ur(2)^3];
x = xAll(1:n);
end
