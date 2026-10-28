function turbInputSet = controlSet_sowfa_2turb_aic_5000(Wp) %#ok<INUSD>
% Closed-loop test case of the adaptive AIC (5000 s): the inputs are only
% applied before the controller starts (CT_prime = 2, yaw 0)
t = 0:5000;
turbInputSet.t = t;
turbInputSet.CT_prime = 2*ones(2,numel(t));
turbInputSet.phi = zeros(2,numel(t));
turbInputSet.interpMethod = 'lin';
end
