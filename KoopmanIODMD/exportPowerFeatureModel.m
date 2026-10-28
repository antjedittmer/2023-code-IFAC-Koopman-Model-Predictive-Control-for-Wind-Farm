function model = exportPowerFeatureModel(filenameId, n, featureSet, tau, t0)
% exportPowerFeatureModel identifies the input-only Koopman power model
% (lifted Ur states, input features of powerFeatures, output PT1 + PT2 in
% MW) by least squares as testPowerFeatures and saves it for the MPC
% (WFSim_demo with KoopFeat = 1) in DataInOutWfSim/PowerInputsOnly/KoopFeatModel.mat
%
% x(k+1) = A x(k) + B w(k), P(k) = C x(k) + D w(k)

if nargin < 1 || isempty(filenameId)
    filenameId = 'Vinf8dot0_sowfa_2turb_yaw_noise_step';
end
if nargin < 2, n = 4; end
if nargin < 3, featureSet = 'NL+1+M+x'; end
if nargin < 4, tau = 100; end
if nargin < 5, t0 = 240; end

mainDir = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(mainDir,'WFSimCode','bin','Kmpc')); % powerFeatures

% Identification data as testPowerInputsOnly (validation set not used here)
idName = filenameId;
if iscell(idName), idName = idName{1}; end
[~, D] = testPowerInputsOnly(filenameId, idName, 6, t0);
close all
S = D.id;

W = powerFeatures(S.lifted(end-2:end,:), featureSet, tau);
X = S.lifted(1:n,:);
P = S.P/1e6;
Psi = [X(:,1:end-1); W(:,1:end-1)];
Y = [X(:,2:end); P(1:end-1)];
sc = sqrt(mean(Psi.^2, 2)); % scaled regressors, small ridge (collinear features)
Psis = Psi./sc;
K = (Y*Psis')/(Psis*Psis' + 1e-8*size(Psis,2)*eye(size(Psis,1)));
K = K./sc';

model = struct('A', K(1:n,1:n), 'B', K(1:n,n+1:end), 'C', K(n+1,1:n), 'D', K(n+1,n+1:end), ...
    'n', n, 'featureSet', featureSet, 'tau', tau, 't0', t0, 'filenameId', {filenameId}, ...
    'maxAbsEig', max(abs(eig(K(1:n,1:n)))));
dirOut = fullfile(mainDir,'DataInOutWfSim','PowerInputsOnly');
if ~isfolder(dirOut)
    mkdir(dirOut);
end
save(fullfile(dirOut,'KoopFeatModel.mat'), 'model');
fprintf('Saved KoopFeatModel.mat: n = %d, %s, tau = %g s, max|eig(A)| = %.4f\n', n, featureSet, tau, model.maxAbsEig);
end
