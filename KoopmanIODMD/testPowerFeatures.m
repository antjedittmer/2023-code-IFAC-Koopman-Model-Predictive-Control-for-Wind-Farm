function T = testPowerFeatures(filenameId, filenameVal, t0, tau)
% testPowerFeatures compares input features for the input-only farm power
% model (constant free stream wind, no wind speed measurement). Data and
% lifting of the measured Ur1, Ur2 as in testPowerInputsOnly.
%
% Feature sets w(k) (inputs ct1, ct2, phi1 in deg; c1 = cos(phi1)):
%   'u+1'        ct1, ct2, phi1, 1
%   'u+1+M'      as 'u+1' plus the moving means M(ct1), M(phi1)
%   'NL+1'       as 'u+1' plus g1 = ct1/(1+ct1/4)^3*c1^3, g2 = ct2/(1+ct2/4)^3
%                (power of a turbine with CT_prime and yaw for a fixed
%                upstream wind speed: P ~ CT'*(cos(phi)*U/(1+CT'/4))^3)
%   'NL+1+M'     as 'NL+1' plus M(ct1), M(phi1), M(a1), M(a1*c1^2) with the
%                induction a1 = (ct1/4)/(1+ct1/4) of turbine 1
%   'NL+1+M+x'   as 'NL+1+M' plus g2*M(a1), g2*M(phi1), g2*M(a1*c1^2)
%                (wake of turbine 1 at turbine 2 times its power feature)
% The moving mean M is an exponential moving average with time constant
% tau (default 100 s): M(k+1) = (1-1/tau)*M(k) + 1/tau*x(k), i.e. a linear
% state that can be added to the MPC prediction model.
%
% Model structures (output PT1 + PT2):
%   'SS n=4/6'  Koopman state space model: lifted Ur states (K_psi) and the
%               features as inputs (as eDMD_RT_Vin_P, linear simulation as
%               in the MPC, initial state from the data)
%   'FIR'       P(k) = sum_j theta_j'*w(k-5j), j = 0..29 (150 s), no state
%   'static'    P(k) = theta'*w(k)

if nargin < 3
    t0 = 240;
end
if nargin < 4
    tau = 100;
end
if ischar(filenameVal)
    filenameVal = {filenameVal};
end

addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))),'WFSimCode','bin','Kmpc')); % powerFeatures
Dall = cell(numel(filenameVal),1);
for idxV = 1:numel(filenameVal)
    [~, Dall{idxV}] = testPowerInputsOnly(filenameId, filenameVal{idxV}, 6, t0);
    close all
end
Did = Dall{1}.id;

featureSets = {'u+1','u+1+M','NL+1','NL+1+M','NL+1+M+x'};
structures = {'SS n=4','SS n=6','FIR','static'};

T = table();
for idxF = 1:numel(featureSets)
    Wid = features(Did, featureSets{idxF}, tau);
    for idxS = 1:numel(structures)
        mdl = fitModel(Did, Wid, structures{idxS});
        for idxV = 1:numel(filenameVal)
            Dv = Dall{idxV}.val;
            y = simModel(mdl, Dv, features(Dv, featureSets{idxF}, tau));
            P = Dv.P/1e6;
            k = 1:min(numel(y), numel(P));
            T = [T; table(string(featureSets{idxF}), string(structures{idxS}), ...
                string(strrep(filenameVal{idxV},'Vinf8dot0_sowfa_2turb_','')), ...
                vaf(P(k), y(k)), sqrt(mean((y(k) - P(k)).^2)), mean(y(k) - P(k)), mdl.rho, ...
                'VariableNames', {'Features','Model','ValSet','VAF','RMSE_MW','MeanErr_MW','maxAbsEig'})]; %#ok<AGROW>
        end
    end
end
disp(T)
dirOut = fullfile(fileparts(fileparts(mfilename('fullpath'))),'DataInOutWfSim','PowerInputsOnly');
save(fullfile(dirOut,'PowerFeatures.mat'), 'T', 'filenameId', 'filenameVal', 'tau');
writetable(T, fullfile(dirOut,'PowerFeatures.csv'));
end

function W = features(S, name, tau)
% Features of the data set S (powerFeatures, also used in the MPC)
W = powerFeatures(S.lifted(end-2:end,:), name, tau);
end

function mdl = fitModel(S, W, structure)
P = S.P/1e6;
lambda = 1e-8; % small ridge on scaled regressors (collinear features)
switch structure
    case {'SS n=4','SS n=6'}
        n = sscanf(structure, 'SS n=%d');
        X = S.lifted(1:n,:);
        Psi = [X(:,1:end-1); W(:,1:end-1)];
        Y = [X(:,2:end); P(1:end-1)];
        K = ridgeLS(Psi', Y', lambda)';
        mdl.sys = ss(K(1:n,1:n), K(1:n,n+1:end), K(n+1,1:n), K(n+1,n+1:end), 1);
        mdl.n = n;
        mdl.rho = max(abs(eig(mdl.sys.A)));
    case 'FIR'
        mdl.lags = 0:5:145;
        mdl.theta = ridgeLS(lagged(W, mdl.lags), P(:), lambda);
        mdl.rho = 0;
    case 'static'
        mdl.lags = 0;
        mdl.theta = ridgeLS(W', P(:), lambda);
        mdl.rho = 0;
end
mdl.structure = structure;
end

function y = simModel(mdl, S, W)
if isfield(mdl, 'sys')
    y = lsim(mdl.sys, W', [], S.lifted(1:mdl.n,1))';
else
    y = (lagged(W, mdl.lags)*mdl.theta)';
end
end

function Phi = lagged(W, lags)
% Rows [w(k)', w(k-l2)', ...]; before the first sample w is held constant
[nw, N] = size(W);
Wp = [repmat(W(:,1), 1, max(lags)), W];
Phi = zeros(N, nw*numel(lags));
for j = 1:numel(lags)
    Phi(:, (j-1)*nw + (1:nw)) = Wp(:, (1:N) + max(lags) - lags(j))';
end
end

function theta = ridgeLS(Phi, Y, lambda)
% (ridge) least squares Phi*theta = Y on regressors scaled to unit rms
s = sqrt(mean(Phi.^2, 1));
s(s == 0) = 1;
Phis = Phi./s;
theta = ((Phis'*Phis + lambda*size(Phis,1)*eye(size(Phis,2)))\(Phis'*Y))./s';
end
