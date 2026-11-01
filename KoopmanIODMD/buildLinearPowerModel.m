function model = buildLinearPowerModel(modelAff, filenameId, filenameVal, t0)
% buildLinearPowerModel linear Koopman power model of the WRC-AIC from the
% control-affine model (feature set 'NL+1+M+x,aff'): the features become
% states, so that the MPC needs no linearization (MatsysPfeat)
%
%   z(k+1) = A z(k) + B u(k),  P_WF(k) = C z(k) + D u(k)  (MW)
%
% States z (21): lifted Ur (Ur1, Ur2, Ur1^2, Ur2^2, Ur1^3, Ur2^3) and the
% feature states of featStatesLin: filtered inputs s = [ct1; ct2; phi1],
% offset 1, g1, g2, moving means Mc, Mgamma, Ma, Mac, products g2*Ma,
% g2*Mgamma, g2*Mac, induction a1 and ac = a1*cos^2(phi1) (all of s);
% inputs u (3): the controls before the filter of the turbine model.
% Rows of A, B, C, D:
% - lifted Ur and output: A_w, B_w, B_u, C_w, D_w, D_u of the control-affine
%   model (linear in its features and controls, i.e. unchanged)
% - filtered inputs: s(k+1) = (1-tauU) s(k) + tauU u(k) (inputFilterCoeff)
% - offset: 1(k+1) = 1(k)
% - moving means: M(k+1) = (1-1/tau) M(k) + 1/tau [s1; s3; a1; ac](k)
% - g1, g2, products, a1, ac: nonlinear functions of s, i.e. of the inputs;
%   identified by least squares from the identification data (EDMD) with
%   all states and inputs as regressors (the only approximated rows)
%
% Inputs
% - modelAff: control-affine model (KoopFeatModel<k>.mat, 'NL+1+M+x,aff')
% - filenameId, filenameVal: data sets as in testPowerStates (identification
%   from t0 on, the validation set is only needed by testPowerInputsOnly)
% - t0: first sample

assert(strcmp(modelAff.featureSet, 'NL+1+M+x,aff'), 'buildLinearPowerModel: control-affine model expected');
mainDir = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(mainDir,'WFSimCode','bin','Kmpc'));
if iscell(filenameVal), filenameVal = filenameVal{1}; end
evalc('[~, D] = testPowerInputsOnly(filenameId, filenameVal, 8, t0);');
close(findobj('Type','figure','Name','PowerInputsOnly'));
S = D.id;

% State data as in testPowerStates (from the second sample on)
featureSet = 'NL+1+M+x,lin';
tau = modelAff.tau;
W = powerFeatures(S.lifted(end-2:end,:), featureSet, tau);
[F, U] = featStatesLin(W);
nL = modelAff.n; % lifted Ur
Z = [S.lifted(1:nL,:); F];
Z = Z(:,2:end); U = U(:,2:end);
nz = size(Z,1); nu = size(U,1);

% Indices of the states
iS = nL + (1:3); iOne = nL + 4; iG = nL + [5 6]; iM = nL + (7:10);
iX = nL + (11:13); iA = nL + [14 15];
iNL = [iG, iX, iA]; % rows identified by least squares

A = zeros(nz); B = zeros(nz, nu);
A(1:nL, 1:nL + 13) = [modelAff.A, modelAff.B(:,1:13)]; % lifted Ur: A_w, B_w
B(1:nL,:) = modelAff.B(:,14:16); % B_u
C = [modelAff.C, modelAff.D(1:13), 0, 0]; % C_w, D_w
Dm = modelAff.D(14:16); % D_u
tauU = inputFilterCoeff;
A(iS,iS) = (1 - tauU)*eye(3);
B(iS,:) = tauU*eye(3);
A(iOne,iOne) = 1;
A(iM,iM) = (1 - 1/tau)*eye(4);
A(iM,[iS(1), iS(3), iA]) = 1/tau*eye(4); % h = [s1; s3; a1; ac]

% Nonlinear rows: least squares (scaled regressors, small ridge term as
% testPowerStates); regressors: the feature states and the inputs, not the
% lifted Ur, as the features are functions of the filtered inputs and the
% moving means only (the feature states are then driven by u only)
iF = nL+1:nz;
Psi = [Z(iF,1:end-1); U(:,1:end-1)];
Y = Z(iNL,2:end);
sc = sqrt(mean(Psi.^2, 2));
K = (Y*(Psi./sc)')/((Psi./sc)*(Psi./sc)' + 1e-8*size(Psi,2)*eye(size(Psi,1)));
K = K./sc';
A(iNL,iF) = K(:,1:numel(iF));
B(iNL,:) = K(:,numel(iF)+1:end);

% Check: one-step prediction of the states and the output on the data
Zp = A*Z(:,1:end-1) + B*U(:,1:end-1);
errNL = sqrt(mean((Zp(iNL,:) - Z(iNL,2:end)).^2, 2))';
errExact = max(max(abs(Zp(setdiff(nL+1:nz, iNL),:) - Z(setdiff(nL+1:nz, iNL),2:end))));
eigA = eig(A(setdiff(1:nz, iOne), setdiff(1:nz, iOne))); % without the offset (eigenvalue 1)
fprintf('Linear power model: %d states, %d inputs, max|eig(A)| without offset %.4f\n', nz, nu, max(abs(eigA)));
fprintf('  one-step RMSE of the identified rows (g1 g2 g2Ma g2Mg g2Mac a1 ac): %s\n', sprintf('%.2e ', errNL));
fprintf('  max one-step error of the known filter rows: %.2e\n', errExact);

model = struct('A',A,'B',B,'C',C,'D',Dm,'n',nL,'nz',nz,'stateSet',modelAff.stateSet, ...
    'featureSet',featureSet,'tau',tau,'maxAbsEig',max(abs(eigA)), ...
    't0',t0,'filenameId',{filenameId},'modelAff',modelAff);
end
