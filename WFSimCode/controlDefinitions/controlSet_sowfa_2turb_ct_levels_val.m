function turbInputSet = controlSet_sowfa_2turb_ct_levels_val(Wp) %#ok<INUSD>
% Validation set for the AIC (thrust control only): random CT_prime
% levels of both turbines, yaw 0 (controlSetHelper_ctYawLevels, seed 4)
turbInputSet = controlSetHelper_ctYawLevels(4, 10000, [0 0]);
end
