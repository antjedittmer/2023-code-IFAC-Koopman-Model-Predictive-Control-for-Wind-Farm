function [model, ad, updated] = koopmanAdaptStep(model, ad, Ur, w, UrEst)
% koopmanAdaptStep online update of the Koopman effective wind speed model
% (adaptive AIC, thesis Section "Koopman Matrix Update", cisneros2020data):
% the model is updated by recursive least squares with a batch of the last
% h samples if the estimation error of the last h samples exceeds nu, and
% the update counter is reset to 1 otherwise.
%
% Model: phi(k+1) = A phi(k) + B w(k), phi = [Ur1 Ur2 Ur1^2 Ur2^2 Ur1^3 Ur2^3]
% K = [A B] = L/Lu, L = sum phi(k+1) phi_w(k)', Lu = sum phi_w(k) phi_w(k)'
%
% Inputs
% - model: struct with A, B, L, Lu (identifyAICnoV1)
% - ad: state of the update (struct with h, nu, nUpd, buffers); [] initializes;
%   optional maxEig: updates resulting in max|eig(A)| >= maxEig are rejected
%   (default inf: all updates are accepted, as in the thesis)
% - Ur: measured effective wind speeds [Ur1; Ur2] at sample k
% - w: inputs at sample k ([CT1; CT2])
% - UrEst: estimated effective wind speeds at sample k
%
% Outputs
% - model: updated model (A, B, L, Lu)
% - ad: updated state (err: estimation error of the last h samples)
% - updated: true if the model was updated in this step

if isempty(ad) || ~isfield(ad,'UrBuf')
    if isempty(ad), ad = struct(); end
    if ~isfield(ad,'h'), ad.h = 5; end          % window/batch length
    if ~isfield(ad,'nu'), ad.nu = 0.15; end     % trigger threshold (m/s)
    if ~isfield(ad,'maxEig'), ad.maxEig = inf; end
    ad.nUpd = 1;
    ad.UrBuf = zeros(2,0); ad.wBuf = zeros(numel(w),0); ad.estBuf = zeros(2,0);
    ad.err = 0;
end
h = ad.h;
ad.UrBuf = [ad.UrBuf(:,max(1,end-h+1):end), Ur(:)];
ad.wBuf = [ad.wBuf(:,max(1,end-h+1):end), w(:)];
ad.estBuf = [ad.estBuf(:,max(1,end-h+1):end), UrEst(:)];

updated = false;
if size(ad.UrBuf,2) < h + 1
    return
end
ad.err = max(abs(ad.UrBuf(:,2:end) - ad.estBuf(:,2:end)),[],'all');
if ad.err >= ad.nu
    % Batch of h samples z(k-h..k-1) and successors x(k-h+1..k)
    Z = [liftUr6(ad.UrBuf(:,1:h)); ad.wBuf(:,1:h)];
    Xp = liftUr6(ad.UrBuf(:,2:h+1));
    nUpd = ad.nUpd + 1;
    Lu = ((nUpd-1)*model.Lu + Z*Z')/nUpd;
    L = ((nUpd-1)*model.L + Xp*Z')/nUpd;
    Kmat = L/Lu;
    n = size(model.A,1);
    if max(abs(eig(Kmat(:,1:n)))) < ad.maxEig
        ad.nUpd = nUpd;
        model.Lu = Lu;
        model.L = L;
        model.A = Kmat(:,1:n);
        model.B = Kmat(:,n+1:end);
        updated = true;
    end
else
    ad.nUpd = 1;
end
end

function phi = liftUr6(Ur)
% Lifted states of K_psi (n = 6): Ur1, Ur2, Ur1^2, Ur2^2, Ur1^3, Ur2^3
phi = [Ur; Ur.^2; Ur.^3];
end
