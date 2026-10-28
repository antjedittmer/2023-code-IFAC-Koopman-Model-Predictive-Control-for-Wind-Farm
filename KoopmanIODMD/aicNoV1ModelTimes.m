function sel = aicNoV1ModelTimes()
% aicNoV1ModelTimes times of the models of the adaptive AIC without V1 that
% are shown with their singular values (plotAICnoV1Sigma) and marked in the
% update subplot of the adaptive open-loop validation (validateAICnoV1):
% the initial model at t0 = 240 s and the updated models at 5100, 6000 and
% 8000 s (ramp from 7 to 9 m/s at 5000 s, last update at 6235 s).
% Colors: initial model gray, the others in colors not used for WT1, WT2
% (lines 1, 2).
%
% Output
% - sel: table with Time_s and Color

cl = lines(5);
sel = table([240; 5100; 6000; 8000], [0.5 0.5 0.5; cl(3:5,:)], ...
    'VariableNames', {'Time_s','Color'});
end
