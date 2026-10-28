function T = testPowerStableVariants(filenameId,filenameVal,t0)
% testPowerStableVariants compares stabilized variants of the linear
% Koopman power model (state space form used in the MPC) for constant free
% stream wind. Identification and validation data and the lifting K_psi as
% in testPowerInputsOnly (lifted states of the measured Ur1, Ur2; inputs
% ct1, ct2, phi1).
%
% Variants (each with the input sets 'V1' = ct1, ct2, phi1, V1 and
% 'u+1' = ct1, ct2, phi1, constant one):
%   - least squares as in eDMD_RT_Vin_P, n = 2, 4, 6 lifted states
%   - ridge regression (scaled regressors), n = 6
%   - eigenvalues of A limited to |lambda| <= 0.99, B refitted, n = 6
% and, for 'u+1' only, an FIR model P(k) = c + sum_i theta_i'*u(k-i) with
% ridge regression (no state, always stable).
%
% All models are simulated over the validation set from the measured
% initial state (lifted Ur1, Ur2 at the first sample); output: farm power.

if ~nargin
    filenameId = 'Vinf8dot0_sowfa_2turb_yaw_noise_step';
    filenameVal = {'Vinf8dot0_sowfa_2turb_yaw_steps'};
    t0 = 240;
end
if ischar(filenameVal)
    filenameVal = {filenameVal};
end
nMax = 6;

% Data sets (lifted states with nMax states, the first n rows are the
% lifting with n states as K_psi stacks them in a fixed order)
Dall = cell(numel(filenameVal),1);
for idxV = 1:numel(filenameVal)
    [~, Dall{idxV}] = testPowerInputsOnly(filenameId, filenameVal{idxV}, nMax, t0);
    close all
end
Did = Dall{1}.id;

variants = {'LS n=2', 2, 0, Inf; 'LS n=4', 4, 0, Inf; 'LS n=6', 6, 0, Inf; ...
    'ridge 1e-6 n=6', 6, 1e-6, Inf; 'ridge 1e-5 n=6', 6, 1e-5, Inf; 'ridge 1e-4 n=6', 6, 1e-4, Inf; ...
    'eig<=0.99 n=6', 6, 0, 0.99; 'eig<=0.95 n=6', 6, 0, 0.95};
inputSets = {'V1','u+1'};

T = table();
for idxI = 1:numel(inputSets)
    for idxVar = 1:size(variants,1)
        [name, n, lambda, rhoMax] = variants{idxVar,:};
        sys = fitLinear(Did, n, inputSets{idxI}, lambda, rhoMax);
        for idxV = 1:numel(filenameVal)
            Dv = Dall{idxV}.val;
            [w, P] = regressors(Dv, inputSets{idxI});
            y = lsim(sys, w', [], Dv.lifted(1:n,1))';
            T = [T; evalRow(inputSets{idxI}, name, filenameVal{idxV}, P, y, max(abs(eig(sys.A))))]; %#ok<AGROW>
        end
    end
end

% FIR models (inputs only + offset)
for L = [60 120 200]
    for lambda = [1e-6 1e-4]
        theta = fitFIR(Did, L, lambda);
        for idxV = 1:numel(filenameVal)
            Dv = Dall{idxV}.val;
            [y, P] = simFIR(Dv, L, theta);
            T = [T; evalRow('u+1', sprintf('FIR L=%d ridge %g', L, lambda), filenameVal{idxV}, P, y, 0)]; %#ok<AGROW>
        end
    end
end
disp(T)
dirOut = fullfile(fileparts(fileparts(mfilename('fullpath'))),'DataInOutWfSim','PowerInputsOnly');
save(fullfile(dirOut,'PowerStableVariants.mat'), 'T');
writetable(T, fullfile(dirOut,'PowerStableVariants.csv'));
end

function [w, P] = regressors(S, inputSet)
% Inputs [ct1;ct2;phi1;extra] and power of the data set S
nIn = 3;
inp = S.lifted(end-nIn+1:end,:);
switch inputSet
    case 'V1',  w = [inp; S.V1];
    case 'u+1', w = [inp; ones(1,numel(S.P))];
end
P = S.P/1e6; % MW
end

function sys = fitLinear(S, n, inputSet, lambda, rhoMax)
% fitLinear x(k+1) = A x(k) + B w(k), P(k) = C x(k) + D w(k) (as
% eDMD_RT_Vin_P) by (ridge) least squares on scaled regressors
[w, P] = regressors(S, inputSet);
X = S.lifted(1:n,:);
Psi = [X(:,1:end-1); w(:,1:end-1)];
Y = [X(:,2:end); P(1:end-1)];
s = sqrt(mean(Psi.^2,2)); % scaling of the regressors
Psis = Psi./s;
N = size(Psis,2);
K = (Y*Psis')/(Psis*Psis' + lambda*N*eye(size(Psis,1)));
K = K./s';
A = K(1:n,1:n); B = K(1:n,n+1:end); C = K(n+1,1:n); D = K(n+1,n+1:end);
if isfinite(rhoMax)
    [V, E] = eig(A);
    e = diag(E);
    e(abs(e) > rhoMax) = rhoMax*e(abs(e) > rhoMax)./abs(e(abs(e) > rhoMax));
    A = real(V*diag(e)/V);
    % refit B with the limited A, and C, D
    B = ((X(:,2:end) - A*X(:,1:end-1))*w(:,1:end-1)')/(w(:,1:end-1)*w(:,1:end-1)');
end
sys = ss(A, B, C, D, 1);
end

function theta = fitFIR(S, L, lambda)
% fitFIR P(k) = c + sum_{i=0}^{L-1} theta_i' u(k-i), ridge on scaled regressors
[Phi, P] = firRegressors(S, L);
s = sqrt(mean(Phi.^2,1));
Phis = Phi./s;
N = size(Phis,1);
theta = ((Phis'*Phis + lambda*N*eye(size(Phis,2)))\(Phis'*P))./s';
end

function [y, P] = simFIR(S, L, theta)
[Phi, P] = firRegressors(S, L);
y = (Phi*theta)';
P = P';
end

function [Phi, P] = firRegressors(S, L)
% Rows: [1, u(k)', u(k-1)', ..., u(k-L+1)'] for k = L..end (inputs ct1,
% ct2, phi1 in the lifted data); before the first sample u is held constant
inp = S.lifted(end-2:end,:);
inp = [repmat(inp(:,1),1,L-1), inp];
N = numel(S.P);
Phi = ones(N, 1 + 3*L);
for i = 0:L-1
    Phi(:, 2 + 3*i : 4 + 3*i) = inp(:, (L:L+N-1) - i)';
end
P = S.P(:)/1e6;
end

function row = evalRow(inputSet, name, valName, P, y, rho)
P = P(1:numel(y)); P = P(:)'; y = y(:)';
row = table(string(inputSet), string(name), string(strrep(valName,'Vinf8dot0_sowfa_2turb_','')), ...
    vaf(P, y), sqrt(mean((y - P).^2)), mean(y - P), max(abs(y - P)), rho, ...
    'VariableNames', {'Inputs','Variant','ValSet','VAF','RMSE_MW','MeanErr_MW','MaxErr_MW','maxAbsEig'});
end
