function ad = adaptSettings(ad)
% adaptSettings settings of the online update of the adaptive AIC without
% V1 (koopmanAdaptStep): window/batch length h = 5, trigger threshold
% nu = 0.15 m/s, the estimator continues from the measured Ur after an
% update (resetState), trigger on Vr1 only (trigChannels = 1). Fields of
% the input ad replace the defaults.
%
% Input
% - ad: struct with some of the fields h, nu, resetState, trigChannels,
%   maxEig (optional)
%
% Output
% - ad: complete settings

adDefault = struct('h',5,'nu',0.15,'resetState',true,'trigChannels',1);
if nargin < 1 || isempty(ad)
    ad = adDefault;
    return
end
for fn = fieldnames(adDefault)'
    if ~isfield(ad, fn{1})
        ad.(fn{1}) = adDefault.(fn{1});
    end
end
end
