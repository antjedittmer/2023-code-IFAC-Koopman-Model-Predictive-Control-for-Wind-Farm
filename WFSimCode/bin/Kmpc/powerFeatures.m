function [W, M] = powerFeatures(inp, featureSet, tau, Mprev)
% powerFeatures input features of the input-only Koopman power model
% (identification in KoopmanIODMD/testPowerFeatures, prediction in the MPC)
%
% Inputs
% - inp: [ct1; ct2; phi1] (phi1 in deg), one column per time step
% - featureSet: 'u', 'u+1', 'u+1+M', 'NL+1', 'NL+1+M' or 'NL+1+M+x', or
%   one of them with ',aff' (control-affine model, see below)
% - tau: time constant of the moving means (s, time step 1 s)
% - Mprev: moving means [M(ct1); M(phi1); M(a1); M(a1*c1^2)] before the
%   first column (default: the values of the first column, i.e. steady state);
%   ',aff': [moving means; filtered inputs s] (7 x 1)
%
% Outputs
% - W: features, one column per time step
% - M: moving means after the last column (state for the next call);
%   ',aff': [moving means; filtered inputs s] after the last column
%
% Features (c1 = cos(phi1), induction a1 = (ct1/4)/(1+ct1/4)):
%   'u'          ct1, ct2, phi1 (the three controls only)
%   'u+1'        ct1, ct2, phi1, 1
%   'u+1+M'      'u+1', M(ct1), M(phi1)
%   'NL+1'       'u+1', g1 = ct1/(1+ct1/4)^3*c1^3, g2 = ct2/(1+ct2/4)^3
%   'NL+1+M'     'NL+1', M(ct1), M(phi1), M(a1), M(a1*c1^2)
%   'NL+1+M+x'   'NL+1+M', g2*M(a1), g2*M(phi1), g2*M(a1*c1^2)
% The moving mean is an exponential moving average including the current
% step: M(k) = (1-1/tau)*M(k-1) + 1/tau*h(k).
%
% Control-affine model ('<set>,aff'): the features of <set> are functions
% of the filtered inputs s, i.e. states, and the three controls u enter
% linearly (W = [features(s, M); u]):
%   s(k+1) = (1-tauU)*s(k) + tauU*u(k)       (filter of the turbine model,
%                                             tauU = inputFilterCoeff)
%   M(k+1) = (1-1/tau)*M(k) + 1/tau*h(s(k))
% inp are the applied inputs, i.e. the filtered controls: s(k+1) = inp(k),
% and u(k) = (inp(k) - (1-tauU)*s(k))/tauU are the controls before the
% filter (in the MPC the decision variables, predictPowerFeat).

if nargin < 4
    Mprev = [];
end
if endsWith(featureSet, ',aff')
    [W, M] = affineFeatures(inp, extractBefore(featureSet, ',aff'), tau, Mprev);
    return
end
h = meanSignals(inp); % signals of the moving means
if isempty(Mprev)
    Mprev = h(:,1);
end
beta = 1 - 1/tau;
Mall = zeros(size(h));
for idx = 1:size(h,1)
    Mall(idx,:) = filter(1/tau, [1, -beta], h(idx,:), beta*Mprev(idx));
end
M = Mall(:,end);
W = features(inp, Mall, featureSet);
end

function [W, z] = affineFeatures(inp, featureSet, tau, z0)
% Features of the filtered inputs s and the controls u (control-affine model)
tauU = inputFilterCoeff;
if isempty(z0) % steady state of the first column
    z0 = [meanSignals(inp(:,1)); inp(:,1)];
end
N = size(inp,2);
s = z0(5:7); Mk = z0(1:4);
S = zeros(3,N); Mall = zeros(4,N); U = zeros(3,N);
for k = 1:N
    S(:,k) = s;
    Mall(:,k) = Mk;
    U(:,k) = (inp(:,k) - (1 - tauU)*s)/tauU;
    Mk = (1 - 1/tau)*Mk + 1/tau*meanSignals(s);
    s = inp(:,k);
end
W = [features(S, Mall, featureSet); U];
z = [Mk; s];
end

function h = meanSignals(inp)
% Signals of the moving means: ct1, phi1, a1, a1*c1^2
ct1 = inp(1,:); phi1 = inp(3,:);
c1 = cosd(phi1);
a1 = (ct1/4)./(1 + ct1/4);
h = [ct1; phi1; a1; a1.*c1.^2];
end

function W = features(inp, Mall, featureSet)
% Features of the feature set from the inputs and their moving means
ct1 = inp(1,:); ct2 = inp(2,:); phi1 = inp(3,:);
c1 = cosd(phi1);
g1 = ct1./(1 + ct1/4).^3.*c1.^3;
g2 = ct2./(1 + ct2/4).^3;
W = [ct1; ct2; phi1; ones(size(ct1))];
switch featureSet
    case 'u'
        W = W(1:3,:);
    case 'u+1'
    case 'u+1+M'
        W = [W; Mall(1:2,:)];
    case 'NL+1'
        W = [W; g1; g2];
    case 'NL+1+M'
        W = [W; g1; g2; Mall];
    case 'NL+1+M+x'
        W = [W; g1; g2; Mall; g2.*Mall(3,:); g2.*Mall(2,:); g2.*Mall(4,:)];
    otherwise
        error('Unknown feature set %s', featureSet);
end
end
