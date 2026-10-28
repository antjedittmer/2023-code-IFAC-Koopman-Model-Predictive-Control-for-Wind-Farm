function Kp = MatsysPfeat(model, Ur, Mprev, U0)
% MatsysPfeat linearizes the input-only Koopman power model
% (predictPowerFeat) around the input trajectory U0 over the horizon and
% returns it in the form of MatsysP: P = L*xinit + S*U with xinit = 1,
% S = dP/dU at U0 and L = P(U0) - S*U0. U is stacked per time step as
% [ct1; ct2; phi1] (as mpc.utemp in H_G_Kmatrix_dU).
%
% Inputs
% - model, Ur, Mprev: see predictPowerFeat (Ur: initial state of koopFeatState)
% - U0: linearization point [ct1; ct2; phi1] (3 x Nh)

[nu, Nh] = size(U0);
P0 = predictPowerFeat(model, Ur, Mprev, U0);
dU = [1e-4; 1e-4; 1e-2]; % finite difference steps (ct1, ct2, phi1 in deg)
S = zeros(Nh, nu*Nh);
for idx = 1:nu*Nh
    U = U0;
    U(idx) = U(idx) + dU(mod(idx-1, nu) + 1);
    S(:,idx) = (predictPowerFeat(model, Ur, Mprev, U) - P0)/dU(mod(idx-1, nu) + 1);
end
Kp.S = S;
Kp.L = P0 - S*U0(:);
Kp.xinit = 1;
Kp.P0 = P0;
end
