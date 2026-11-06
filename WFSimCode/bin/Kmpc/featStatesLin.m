function [F, U] = featStatesLin(W)
% featStatesLin feature states and controls of the linear power model
% (feature set 'NL+1+M+x,lin', buildLinearPowerModel) from the features of
% the control-affine model W = [w; u] (powerFeatures, 'NL+1+M+x,aff')
%
% Outputs
% - F: feature states [s1; s2; s3; 1; g1; g2; Mc; Mgamma; Ma; Mac;
%   g2*Ma; g2*Mgamma; g2*Mac; a1; ac] (15 x N) with the filtered inputs
%   s = [ct1; ct2; phi1], the induction a1 and ac = a1*cos^2(phi1) of s
% - U: controls before the filter (3 x N), the inputs of the model

assert(size(W,1) == 16, 'featStatesLin: features of ''NL+1+M+x,aff'' expected');
s = W(1:3,:);
a1 = (s(1,:)/4)./(1 + s(1,:)/4);
F = [W(1:13,:); a1; a1.*cosd(s(3,:)).^2];
U = W(14:16,:);
end
