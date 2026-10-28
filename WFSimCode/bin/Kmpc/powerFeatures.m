function [W, M] = powerFeatures(inp, featureSet, tau, Mprev)
% powerFeatures input features of the input-only Koopman power model
% (identification in KoopmanIODMD/testPowerFeatures, prediction in the MPC)
%
% Inputs
% - inp: [ct1; ct2; phi1] (phi1 in deg), one column per time step
% - featureSet: 'u+1', 'u+1+M', 'NL+1', 'NL+1+M' or 'NL+1+M+x'
% - tau: time constant of the moving means (s, time step 1 s)
% - Mprev: moving means [M(ct1); M(phi1); M(a1); M(a1*c1^2)] before the
%   first column (default: the values of the first column, i.e. steady state)
%
% Outputs
% - W: features, one column per time step
% - M: moving means after the last column (state for the next call)
%
% Features (c1 = cos(phi1), induction a1 = (ct1/4)/(1+ct1/4)):
%   'u+1'        ct1, ct2, phi1, 1
%   'u+1+M'      'u+1', M(ct1), M(phi1)
%   'NL+1'       'u+1', g1 = ct1/(1+ct1/4)^3*c1^3, g2 = ct2/(1+ct2/4)^3
%   'NL+1+M'     'NL+1', M(ct1), M(phi1), M(a1), M(a1*c1^2)
%   'NL+1+M+x'   'NL+1+M', g2*M(a1), g2*M(phi1), g2*M(a1*c1^2)
% The moving mean is an exponential moving average including the current
% step: M(k) = (1-1/tau)*M(k-1) + 1/tau*h(k).

ct1 = inp(1,:); ct2 = inp(2,:); phi1 = inp(3,:);
c1 = cosd(phi1);
g1 = ct1./(1 + ct1/4).^3.*c1.^3;
g2 = ct2./(1 + ct2/4).^3;
a1 = (ct1/4)./(1 + ct1/4);
h = [ct1; phi1; a1; a1.*c1.^2]; % signals of the moving means
if nargin < 4 || isempty(Mprev)
    Mprev = h(:,1);
end
beta = 1 - 1/tau;
Mall = zeros(size(h));
for idx = 1:size(h,1)
    Mall(idx,:) = filter(1/tau, [1, -beta], h(idx,:), beta*Mprev(idx));
end
M = Mall(:,end);

W = [ct1; ct2; phi1; ones(size(ct1))];
switch featureSet
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
