function model = identifyAICnoV1(filenameId, t0)
% identifyAICnoV1 Koopman model of the effective wind speeds for the
% adaptive AIC without the measured wind speed V1 (thesis Section
% "Koopman Matrix Creation" with the thrust control signals only,
% w = [CT1; CT2] instead of [V1; CT1; CT2], no yaw):
% phi(k+1) = A phi(k) + B w(k), phi = [Ur1 Ur2 Ur1^2 Ur2^2 Ur1^3 Ur2^3].
% Saved as DataInOutWfSim/AdaptiveAIC/KoopAICnoV1.mat with the correlation
% matrices L = G+ Gw', Lu = Gw Gw' (sums over the samples, as in the thesis)
% as initial values of the online update.
%
% Inputs
% - filenameId: open-loop data (default: random CT levels at V_inf = 7 m/s)
% - t0: first sample (s), default 2: the start-up of the wake (Vr1 = Vr2)
%   is needed, the closed loop starts from the undisturbed flow

if nargin < 1 || isempty(filenameId)
    filenameId = 'Vinf7dot0_sowfa_2turb_ct_levels';
end
if nargin < 2, t0 = 2; end
mainDir = fileparts(fileparts(mfilename('fullpath')));
S = load(fullfile(mainDir,'DataT2OLWFSim','Vinf7dot0_OL_Ct',filenameId),'Ur1','Ur2','Ct1','Ct2');
k = t0:numel(S.Ur1);
Ur = [S.Ur1(k); S.Ur2(k)];
W = [S.Ct1(k); S.Ct2(k)];

G = [Ur; Ur.^2; Ur.^3];
Gw = [G(:,1:end-1); W(:,1:end-1)];
Gp = G(:,2:end);
N = size(Gw,2);
L = Gp*Gw';
Lu = Gw*Gw';
Kmat = L/Lu;
n = size(G,1);
model = struct('A',Kmat(:,1:n),'B',Kmat(:,n+1:end),'L',L,'Lu',Lu,'n',n, ...
    'inputs',{{'CT1','CT2'}},'filenameId',filenameId,'t0',t0);
model.maxAbsEig = max(abs(eig(model.A)));

dirOut = fullfile(mainDir,'DataInOutWfSim','AdaptiveAIC');
if ~isfolder(dirOut), mkdir(dirOut); end
save(fullfile(dirOut,'KoopAICnoV1.mat'),'model');
fprintf('KoopAICnoV1: %d samples, max|eig(A)| = %.4f\n', N, model.maxAbsEig);
end
