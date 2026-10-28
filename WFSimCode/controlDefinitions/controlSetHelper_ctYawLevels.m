function turbInputSet = controlSetHelper_ctYawLevels(seed, tEnd)
% controlSetHelper_ctYawLevels open loop excitation for the Koopman sys id:
% the thrust coefficients CT_prime of both turbines and the yaw angle of
% turbine 1 jump independently between random levels, held for 150 s to
% 300 s each (longer than the wake travel time from WT1 to WT2 of roughly
% 100 s), with small noise on top and the rate limits of the MPC.
%
% Inputs
% - seed: random seed (identification and validation sets differ)
% - tEnd: length of the input sequence (s)

rng(seed);
holdRange = [150 300]; % (s)
ctRange = [0.5 2];
phiRange = [0 30];      % (deg)
ctRLim = 0.05;          % (1/s)
phiRLim = 1;            % (deg/s)
nNoise = 10;            % noise held for nNoise s
ctNoise = 0.05;
phiNoise = 0.5;         % (deg)

t = 0:tEnd-1;
CT_prime = [randomLevels(t, holdRange, ctRange); randomLevels(t, holdRange, ctRange)];
phi1 = randomLevels(t, holdRange, phiRange);

CT_prime = CT_prime + ctNoise*heldNoise(2, numel(t), nNoise);
phi1 = phi1 + phiNoise*heldNoise(1, numel(t), nNoise);
CT_prime = min(max(CT_prime, 0.2), 2);
phi1 = min(max(phi1, 0), 30);

CT_prime = rateLimit(CT_prime, ctRLim);
phi1 = rateLimit(phi1, phiRLim);

turbInputSet.t = t;
turbInputSet.CT_prime = CT_prime;
turbInputSet.phi = [phi1; zeros(size(phi1))];
turbInputSet.interpMethod = 'lin'; % Linear interpolation over time
end

function x = randomLevels(t, holdRange, valRange)
% Piecewise constant signal: random levels in valRange, random hold times
x = nan(size(t));
k = 1;
while k <= numel(t)
    tHold = randi(holdRange);
    x(k:min(k+tHold-1, numel(t))) = valRange(1) + diff(valRange)*rand;
    k = k + tHold;
end
end

function w = heldNoise(nRows, nCols, nHold)
% Normally distributed noise, each value held for nHold samples
w = repelem(randn(nRows, ceil(nCols/nHold)), 1, nHold);
w = w(:, 1:nCols);
end

function x = rateLimit(x, rLim)
for idx = 2:size(x,2)
    x(:,idx) = x(:,idx-1) + min(max(x(:,idx) - x(:,idx-1), -rLim), rLim);
end
end
