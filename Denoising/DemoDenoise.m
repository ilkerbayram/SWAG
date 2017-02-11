% Demo file for the denoising experiment from the manuscript
% 'A Penalty Promoting Sparsity Within and Across Groups', by I. Bayram and
% S. Bulek, 2016. 
%
% Matlab code by I. Bayram and S. Bulek, 2016.

clear all;
close all;

[x,fs] = audioread('sp.wav'); % the clean signal

[n,fs2] = audioread('ambience_casino.wav'); % the noise signal
n = n((1:length(x)),1);
n = sqrt(length(n)) * n / sqrt(sum(n.^2));
SNR = 5;
sig = sqrt( sum(x.^2)/(10^(SNR/10))/length(x) );

y = x + sig * n; % observed signal

%%% STFT parameters
N = round(60*fs/1000);
Hop = round(N/4);
win = hamming(N);
param.win = NormalizeW(win,Hop);
param.hop = Hop;



% penalty function parameters
% make sure that  (K*lam*gam < 1) to ensure convexity
param.K = 16; % group size (i.e., we are using a K-variate penalty)
param.Gs = 8; % group size for the mixed norm


% Maximum number of iterations for the Douglas-Rachford algorithm
param.MAX_ITER = 100;


%% l_{2,1} norm regularization
lamlist3 = (3) *  sig * 1.5.^(-20:5); % sweep search for the best lambda around sigma
SNRlistMixed = [];
bestSNR = - inf;
bestxMixed = zeros(size(x));
for lam = lamlist3,
    param.lam = lam;
    
    [rec] = DenoiseMixed(y,param);
    
    %oracle normalization of energy
    enr = sqrt(sum(rec.^2));
    rec = enx * rec / enr; % normalize for a better comparison of SNR
    
    SNRlistMixed = [SNRlistMixed snr(x,rec-x)]
    if SNRlistMixed(end) > bestSNR,
        bestSNR = SNRlistMixed(end);
        bestxMixed = rec;
    end
end

%% Proposed hybrid penalty
[m,k] = max(SNRlistMixed);
param.lam = lamlist3(k) * 0.5; % use half of the lambda value found in l_{2,1} regularization
gamlist = 0.1 * 1.2.^(-15:0) / param.lam; % sweep search for the best gamma
SNRlistHybrid = [];
bestSNR = - inf;
bestxHybrid = zeros(size(y)); % will hold the best reconstruction for the hybrid 
for gam = gamlist,
    param.gam = gam;
    
    [rec] = DenoiseHybrid(y,param);
    
    % oracle normalization of energy
    enr = sqrt(sum(rec.^2));
    rec = enx * rec / enr; % normalize for a better comparison of SNR
    
    SNRlistHybrid = [SNRlistHybrid snr(x,rec-x)]
    if SNRlistHybrid(end) > bestSNR,
        bestSNR = SNRlistHybrid(end);
        bestxHybrid = rec;
    end
end

%% l_1 regularization

lamlist = 3 *  sig * 1.2.^(-20:0); % search for the best lambda value...
SNRlistSoft = [];
bestSNR = - inf;
bestxSoft = zeros(size(x)); % will hold the best reconstruction for l1 regularization

for lam = lamlist,
    param.lam = lam;
    
    [rec] = DenoiseL1(y,param);
    
    %oracle normalization of energy
    enr = sqrt(sum(rec.^2));
    rec = enx * rec / enr; % normalize for a better comparison of SNR
    
    SNRlistSoft = [SNRlistSoft snr(x,rec-x)]
    if SNRlistSoft(end) > bestSNR,
        bestSNR = SNRlistSoft(end);
        bestxSoft = rec;
    end
end



%% Proposed penalty P_{\gamma}

[m,k] = max(SNRlistSoft);
param.lam = lamlist(k) * 0.5; % use half of the lambda value found in l_1 regularization

gamlist = 0.1 * 1.2.^(-15:0) / param.lam; % sweep search for the best gamma
SNRlist = [];
bestSNR = - inf;
bestx = zeros(size(y)); % will hold the best reconstruction for the proposed method

for gam = gamlist,
    param.gam = gam;
    
    [rec] = DenoiseProposed(y,param);
    
    % oracle normalization of energy
    enr = sqrt(sum(rec.^2));
    rec = enx * rec / enr; % normalize for a better comparison of SNR
    
    SNRlist = [SNRlist snr(x,rec-x)]
    if SNRlist(end) > bestSNR,
        bestSNR = SNRlist(end);
        bestx = rec;
    end
end



%% E-Lasso penalty
lamlist2 = 2 *  sig * 1.2.^(-15:0); % sweep search for the best lambda value...
SNRlistE = [];
bestSNR = - inf;
bestxE = zeros(size(x)); % will hold the best reconstruction for E-Lasso

for lam = lamlist2,
    param.lam = lam;
    [rec] = DenoiseELasso(y,param);
    
    %oracle normalization of energy
    enr = sqrt(sum(rec.^2));
    rec = enx * rec / enr; % normalize for a better comparison of SNR
    
    SNRlistE = [SNRlistE snr(x,rec-x)]
    if SNRlistE(end) > bestSNR,
        bestSNR = SNRlistE(end);
        bestxE = rec;
    end
end

%% list the highest SNR for each method ... 
SNR_nvar = max(SNRlist)
SNR_soft = max(SNRlistSoft)
SNR_E = max(SNRlistE)
SNR_Mixed = max(SNRlistMixed)
SNR_MixedNvar = max(SNRlistHybrid)


%% display the spectrograms
% set the parameters for visualization
param.N = round(50*fs/1000);
param.hop = round(param.N/4);
param.win = hamming(param.N);
param.win = NormalizeW(param.win,param.hop);

param.Fr = [0 2500]; % frequency range for STFT display
param.clim = [-50 0]; % decibel range


Y = STFT(n,win,Hop);
Norm = max(abs(Y(:)));
figure;
subplot(3,1,[1 2]);
DispSTFT(Y/Norm, param);
title('Spectrogram of Noise');

Y = STFT(x,win,Hop);
Norm = max(abs(Y(:)));
figure;
subplot(3,1,[1 2]);
DispSTFT(Y/Norm, param);
title('Clean Signal');

Y = STFT(y,win,Hop);
Norm = max(abs(Y(:)));
figure;
subplot(3,1,[1 2]);
DispSTFT(Y/Norm, param);
title('Observed Noisy Signal');

Y = STFT(bestx,win,Hop);
figure;
subplot(3,1,[1 2]);
DispSTFT(Y/Norm, param);
title('Proposed');

Y = STFT(bestxSoft,win,Hop);
figure;
subplot(3,1,[1 2]);
DispSTFT(Y/Norm, param);
title('l_1 norm');

Y = STFT(bestxE,win,Hop);
figure;
subplot(3,1,[1 2]);
DispSTFT(Y/Norm, param);
title('E-Lasso');

Y = STFT(bestxMixed,win,Hop);
figure;
subplot(3,1,[1 2]);
DispSTFT(Y/Norm, param);
title('Mixed Norm');

Y = STFT(bestxMixedNvar,win,Hop);
figure;
subplot(3,1,[1 2]);
DispSTFT(Y/Norm, param);
title('Proposed Hybrid Penalty');