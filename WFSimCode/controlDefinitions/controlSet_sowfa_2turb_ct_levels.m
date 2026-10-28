function turbInputSet = controlSet_sowfa_2turb_ct_levels(Wp) %#ok<INUSD>
% Identification set for the AIC (thrust control only): random CT_prime
% levels of both turbines, yaw 0 (controlSetHelper_ctYawLevels, seed 3)
turbInputSet = controlSetHelper_ctYawLevels(3, 30000, [0 0]);
end
