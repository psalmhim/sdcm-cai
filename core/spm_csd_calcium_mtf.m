function [y,w,S,Gu,Gn] = spm_csd_calcium_mtf(P,M,U)
%==========================================================================
% spm_csd_calcium_mtf
% -------------------------------------------------------------------------
% Spectral response of a Calcium-DCM (transfer function × noise spectrum)
%
% FORMAT [y,w,S,Gu,Gn] = spm_csd_calcium_mtf(P,M,U)
%
% P - model parameters (includes A,B,C, tau_Ca, kappa_Ca, alpha, beta)
% M - model structure (fields: .Hz, .u, .autoscale)
% U - model inputs (expects U.csd as complex cross spectra)
%
% y  - cross-spectral density (nw × nn × nn)
% w  - frequencies (Hz)
% S  - directed transfer functions (complex)
% Gu - neuronal (endogenous) fluctuation spectra
% Gn - observation noise spectra
%
% This function parallels spm_csd_fmri_mtf but implements calcium kinetics
% via spm_fx_calcium, with neuronal-to-calcium transfer rather than BOLD.
%
%--------------------------------------------------------------------------
% Author: HJ Park & ChatGPT, 2025
%==========================================================================

% Multiple sessions
[m,n,s] = size(P.B);
if s > 1
    % centre session-specific parameters
    for i = 1:n
        for j = 1:n
            P.B(i,j,:) = P.B(i,j,:) - mean(P.B(i,j,:));
        end
    end
    % session-wise predictions
    for i = 1:s
        Q    = P;
        Q.A  = Q.A + Q.B(:,:,i);
        Q.B  = zeros(n,n,0);
        y{i} = spm_csd_calcium_mtf(Q,M,U);
    end
    return
end

%==========================================================================
% Frequency grid and dimensional setup
%==========================================================================
w    = M.Hz(:);
nw   = length(w);
nn   = size(P.A,1);
nu   = length(M.u);
form = '1/f';

if ~isfield(M,'autoscale'), M.autoscale = false; end

%==========================================================================
% Spectra of neuronal fluctuations (Gu) and observation noise (Gn)
%==========================================================================
Gu = zeros(nw,nu,nu);
Gn = zeros(nw,nn,nn);

% experimental (driven) components
if any(any(P.C))
    for i = 1:nu
        for j = 1:nu
            for k = 1:nw
                Gu(k,i,j) = P.C(i,:)*squeeze(U.csd(k,:,:))*P.C(j,:)';
            end
        end
    end
end

% endogenous neuronal fluctuations
for i = 1:nu
    if strcmp(form,'1/f')
        G = w.^(-exp(P.a(2,1)));
    else
        G = spm_mar2csd(exp(P.a(2,1)),w);
    end
    Gu(:,i,i) = Gu(:,i,i) + exp(P.a(1,1))*G/sum(G);
end

% region-specific observation noise
for i = 1:nn
    if strcmp(form,'1/f')
        G = w.^(-exp(P.b(2,1))/2);
    else
        G = spm_mar2csd(exp(P.b(2,1))/2,w);
    end
    Gn(:,i,i) = Gn(:,i,i) + exp(P.c(1,i))*G/sum(G);
end

% global (shared) observation noise
if strcmp(form,'1/f')
    G = w.^(-exp(P.b(2,1))/2);
else
    G = spm_mar2csd(exp(P.b(2,1))/2,w);
end
for i = 1:nn
    for j = i:nn
        Gn(:,i,j) = Gn(:,i,j) + exp(P.b(1,1))*G/sum(G);
        Gn(:,j,i) = Gn(:,i,j);
    end
end

%==========================================================================
% Directed transfer functions (Calcium-specific)
%==========================================================================
% Only set an identity `P.C` if caller did not provide one. Overwriting a
% supplied `P.C` would remove any driven/input mapping used above.
if ~isfield(P,'C') || isempty(P.C)
    P.C = speye(nn,nu);
end

% Use custom transfer function that handles calcium model properly
[S, Hz_tf] = spm_dcm_mtf_calcium(P, M);

% Interpolate to match requested frequency vector if needed
if length(Hz_tf) ~= nw || any(abs(Hz_tf - w) > 1e-6)
    S_interp = zeros(nw, nn, nn);
    for i = 1:nn
        for j = 1:nn
            S_interp(:,i,j) = interp1(Hz_tf, squeeze(S(:,i,j)), w, 'linear', 'extrap');
        end
    end
    S = S_interp;
end

%==========================================================================
% Predicted cross-spectral densities
%==========================================================================
G = zeros(nw,nn,nn);
for i = 1:nw
    S_i = reshape(S(i,:,:),nn,nn);
    Gu_i = reshape(Gu(i,:,:),nu,nu);
    G(i,:,:) = S_i * Gu_i * S_i';
end

% Add observation noise
if isfield(M,'g')
    y = G + Gn;
else
    y = G;
end

%==========================================================================
% Optional normalization (for cross-modality comparability)
%==========================================================================
if M.autoscale
    for k = 1:nw
        y(k,:,:) = y(k,:,:)/max(abs(y(k,:,:)),[],'all');
    end
end

%==========================================================================
% Autoregressive parameterisation (same as fMRI)
%==========================================================================
y = spm_mar2csd(spm_csd2mar(y,M.Hz,M.p - 1),M.Hz);

end
