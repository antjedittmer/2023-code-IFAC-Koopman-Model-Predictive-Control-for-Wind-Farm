function [T, Mdl] = testPowerStates(filenameId, filenameVal, featureSets, t0, tau, stateSets)
% testPowerStates compares Koopman states of the input-only farm power
% model: the lifted rotor effective wind speeds Ur (as K_psi) and liftings
% of the turbine powers P1, P2 of the previous time step (MW), which are
% measured directly. Inputs: features of powerFeatures (ct1, ct2, phi1).
%
% Model: x(k+1) = A x(k) + B w(k), P(k) = C x(k) + D w(k), identified by
% least squares as testPowerFeatures, simulated linearly over the
% validation set from the measured initial state (as in the MPC).
%
% Feature sets: see powerFeatures; '<set>+V1' appends the measured wind
% speed V1 (u in front of turbine 1) as input, 'u+V1' = ct1, ct2, phi1, V1
% (inputs of MainWfSimVinYawPowerOut, no constant).
%
% States x(k) (Pi = Pi(k-1) in MW, cti = cti(k-1)):
%   'Ur n=4'        Ur1, Ur2, Ur1^2, Ur2^2 (reference, as exportPowerFeatureModel)
%   'Ur n=6'        'Ur n=4', Ur1^3, Ur2^3 (K_psi)
%   'Ur n=8'        'Ur n=6', DUr1, DUr2 (K_psi: deviation from the moving
%                   mean of the past Ur over 25 s)
%   'Ur n=4+d', 'Ur n=6+d', 'Ur n=8+d'   plus the delayed Ur1(k-1), Ur2(k-1)
%   'Ur n=6+D<N>'   'Ur n=6' plus DUr1, DUr2 with the moving mean over N
%                   samples (recursion of K_psi, which uses N = 25)
%   'Ur n=6+D<N>,M' as 'Ur n=6+D<N>' plus the moving means M1, M2
%   'Ur n=6+M<N>fix' 'Ur n=6' plus M1, M2; the rows of M in A are not fitted
%                   but set to the known recursion M(k+1) = (1-1/N)*M(k) + Ur(k)/N
%   'Ur n=4+Ur1(25,50,100)'  'Ur n=4' plus the shift register Ur1(k-1..k-dmax)
%                   (fixed shift rows); only the listed delays enter the Ur rows
%                   and the output; Ur1 before the first sample held at Ur1(1)
%   'P'             P1, P2
%   'P,sqrtP'       P1, P2, sqrt(P1), sqrt(P2)
%   'P,cbrtP'       P1, P2, P1^(1/3), P2^(1/3)
%   'P,Ue'          P1, P2, (P1/ct1)^(1/3), (P2/ct2)^(1/3)  (~ Ur, P = c*ct'*Ue^3)
%   'P,sqrtP,cbrtP' P1, P2, sqrt(P1), sqrt(P2), P1^(1/3), P2^(1/3)

if nargin < 1 || isempty(filenameId)
    filenameId = 'Vinf8dot0_sowfa_2turb_yaw_noise_step';
end
if nargin < 2 || isempty(filenameVal)
    filenameVal = {'Vinf8dot0_sowfa_2turb_yaw_steps','Vinf8dot0_sowfa_2turb_yaw_steps_Ct_comb'};
end
if ischar(filenameVal), filenameVal = {filenameVal}; end
if nargin < 3 || isempty(featureSets), featureSets = {'u+1','NL+1+M+x'}; end
if nargin < 4, t0 = 240; end
if nargin < 5, tau = 100; end

mainDir = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(mainDir,'WFSimCode','bin','Kmpc')); % powerFeatures
Dall = cell(numel(filenameVal),1);
for idxV = 1:numel(filenameVal)
    evalc('[~, Dall{idxV}] = testPowerInputsOnly(filenameId, filenameVal{idxV}, 8, t0);');
    close all
end

if nargin < 6 || isempty(stateSets)
    stateSets = {'Ur n=4','P','P,sqrtP','P,cbrtP','P,Ue','P,sqrtP,cbrtP'};
end
T = table();
Mdl = struct('A',{},'B',{},'C',{},'D',{},'n',{},'stateSet',{},'featureSet',{},'tau',{},'maxAbsEig',{});
for idxF = 1:numel(featureSets)
    for idxS = 1:numel(stateSets)
        [X, W, P] = dataMatrices(Dall{1}.id, stateSets{idxS}, featureSets{idxF}, tau);
        n = size(X,1);
        Psi = [X(:,1:end-1); W(:,1:end-1)];
        Y = [X(:,2:end); P(1:end-1)];
        sc = sqrt(mean(Psi.^2, 2));
        K = (Y*(Psi./sc)')/((Psi./sc)*(Psi./sc)' + 1e-8*size(Psi,2)*eye(size(Psi,1)));
        K = K./sc';
        tokFix = regexp(stateSets{idxS},'\+M(\d+)fix$','tokens','once');
        if ~isempty(tokFix) % known moving mean dynamics (last two states)
            N = str2double(tokFix{1});
            K(n-1:n,:) = 0;
            K(n-1,[1,n-1]) = [1/N, 1-1/N];
            K(n,[2,n]) = [1/N, 1-1/N];
        end
        tokDel = regexp(stateSets{idxS},'^Ur n=(\d+)\+Ur1\(([\d,]+)\)$','tokens','once');
        if ~isempty(tokDel) % shift register Ur1(k-1..k-dmax): fixed shift rows,
            % only the delays dl enter the Ur rows and the output (refit)
            nL = str2double(tokDel{1});
            dl = str2double(split(tokDel{2},','))';
            iSel = [1:nL, nL + dl, n+1:size(Psi,1)];
            iFree = [1:nL, n+1];
            sc = sqrt(mean(Psi(iSel,:).^2, 2));
            Ps = Psi(iSel,:)./sc;
            K(iFree,:) = 0;
            K(iFree,iSel) = ((Y(iFree,:)*Ps')/(Ps*Ps' + 1e-8*size(Ps,2)*eye(numel(iSel))))./sc';
            K(nL+1:n,:) = 0;
            K(nL+1,1) = 1;
            K(sub2ind(size(K), nL+2:n, nL+1:n-1)) = 1;
        end
        sys = ss(K(1:n,1:n), K(1:n,n+1:end), K(n+1,1:n), K(n+1,n+1:end), 1);
        Mdl(end+1) = struct('A',sys.A,'B',sys.B,'C',sys.C,'D',sys.D,'n',n, ... %#ok<AGROW>
            'stateSet',stateSets{idxS},'featureSet',featureSets{idxF},'tau',tau,'maxAbsEig',max(abs(eig(sys.A))));
        for idxV = 1:numel(filenameVal)
            [Xv, Wv, Pv] = dataMatrices(Dall{idxV}.val, stateSets{idxS}, featureSets{idxF}, tau);
            y = lsim(sys, Wv', [], Xv(:,1))';
            T = [T; table(string(featureSets{idxF}), string(stateSets{idxS}), ...
                string(strrep(filenameVal{idxV},'Vinf8dot0_sowfa_2turb_','')), ...
                vaf(Pv, y), sqrt(mean((y - Pv).^2)), mean(y - Pv), max(abs(eig(sys.A))), ...
                'VariableNames', {'Features','States','ValSet','VAF','RMSE_MW','MeanErr_MW','maxAbsEig'})]; %#ok<AGROW>
        end
    end
end
disp(T)
dirOut = fullfile(mainDir,'DataInOutWfSim','PowerInputsOnly');
writetable(T, fullfile(dirOut,'PowerStates.csv'));
end

function [X, W, P] = dataMatrices(S, stateSet, featureSet, tau)
% States from k = 2 on (power states use the previous time step)
inp = S.lifted(end-2:end,:);
if endsWith(featureSet,'+V1') % measured wind speed V1 as additional input
    baseSet = extractBefore(featureSet, strlength(featureSet) - 2);
    if strcmp(baseSet,'u') % as the paper: ct1, ct2, phi1, V1 (no constant)
        W = [inp; S.V1];
    else
        W = [powerFeatures(inp, baseSet, tau); S.V1];
    end
else
    W = powerFeatures(inp, featureSet, tau);
end
P12 = S.P12/1e6;
Pp = P12(:,1:end-1);      % P(k-1)
ctp = inp(1:2,1:end-1);   % ct(k-1)
switch stateSet
    case {'Ur n=4','Ur n=6','Ur n=8'}
        X = S.lifted(1:sscanf(stateSet,'Ur n=%d'),2:end);
    case {'Ur n=4+d','Ur n=6+d','Ur n=8+d'}
        X = [S.lifted(1:sscanf(stateSet,'Ur n=%d'),2:end); S.lifted(1:2,1:end-1)];
    case 'P'
        X = Pp;
    case regexp(stateSet,'^Ur n=\d+\+Ur1\([\d,]+\)$','match','once')
        tok = regexp(stateSet,'^Ur n=(\d+)\+Ur1\(([\d,]+)\)$','tokens','once');
        dmax = max(str2double(split(tok{2},',')));
        Ur1 = [repmat(S.lifted(1,1), 1, dmax), S.lifted(1,:)]; % held before the start
        N = size(S.lifted,2);
        X = [S.lifted(1:str2double(tok{1}),:); zeros(dmax, N)];
        for j = 1:dmax
            X(str2double(tok{1}) + j,:) = Ur1((1:N) + dmax - j); % Ur1(k-j)
        end
        X = X(:,2:end);
    case regexp(stateSet,'^Ur n=\d+\+M\d+fix$','match','once')
        tok = regexp(stateSet,'^Ur n=(\d+)\+M(\d+)','tokens','once');
        [~, M] = movingMeanDev(S.lifted(1:2,:), str2double(tok{2}));
        X = [S.lifted(1:str2double(tok{1}),:); M];
        X = X(:,2:end);
    case regexp(stateSet,'^Ur n=\d+\+D\d+(,M)?$','match','once')
        tok = regexp(stateSet,'^Ur n=(\d+)\+D(\d+)','tokens','once');
        [DUr, M] = movingMeanDev(S.lifted(1:2,:), str2double(tok{2}));
        X = S.lifted(1:str2double(tok{1}),:);
        X = [X; DUr];
        if endsWith(stateSet,',M'), X = [X; M]; end
        X = X(:,2:end);
    case 'P,sqrtP'
        X = [Pp; sqrt(Pp)];
    case 'P,cbrtP'
        X = [Pp; nthroot(Pp,3)];
    case 'P,Ue'
        X = [Pp; nthroot(Pp./ctp,3)];
    case 'P,sqrtP,cbrtP'
        X = [Pp; sqrt(Pp); nthroot(Pp,3)];
end
W = W(:,2:end);
P = S.P(2:end)/1e6;
end

function [DUr, M] = movingMeanDev(Ur, N)
% Deviation of Ur from the moving mean of the past Ur over N samples, with
% the recursion of K_psi (initial mean Ur(:,1), k = 1 at the first sample)
M = zeros(size(Ur));
Mk = Ur(:,1);
UrPrev = Ur(:,1);
for k = 1:size(Ur,2)
    l = min(N,k);
    Mk = (l-1)/l*Mk + UrPrev/l;
    M(:,k) = Mk;
    UrPrev = Ur(:,k);
end
DUr = Ur - M;
end
