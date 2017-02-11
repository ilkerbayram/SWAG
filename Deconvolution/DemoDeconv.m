% Main program for the deconvolution simulations in Sec.IV of the manuscript
% 'A Penalty Function Promoting  Sparsity Within and Across Groups',
% by I. Bayram and S. Bulek, 2016.
%
% compares four algorithms:
%  - iterative p-shrinkage algorithm (IPS)
%  - mixed norms
%  - sparse Group Lasso (SGL)
%  - proposed  penalty

% %% sweep search
% % first determine the parameters for different SNRs with a sweep search
% % the script files below perform a search and save the search results in
% % .mat files for later access
% SweepSearch_ips;
% SweepSearch_mixed;
% SweepSearch_sgl;
% SweepSearch_Proposed;

%% deconvolution with the optimal parameters found by sweep search

clear all;
close all;

clc;

% Change to choose another realization of gaussian noise
seed = 1 ;
randn('state',seed);
rand('state',seed);

%% the wavelet
load('ricker.mat') ;
wL = length(h); % length of the blur function

%% the sparse signal
load('sparse_sig.mat')
xL = length(x);

figure
subplot(212)
plot(x)
title(['Sparse signal with nonzero samples = ' num2str(sum(x~=0))])

subplot(211)
plot(h)
title('Wavelet')

%% produce the noise-free observation
y_noisefree = convH(x,h); % noise-free observations

MAX_ITER = 2000;

NTrial = 100; % number of trials

Gs = 8;  % group size for SGL, mixed norm and the propose


SNRlist = [5 10 15 20]; % SNR can be obe of these values (sweep search is performed for these values)
SNRind = 1;

SNR = SNRlist(SNRind);


% parameters

sig = sqrt( sum(abs(y_noisefree).^2) / (xL * 10^(SNR/10))); % noise std

alp = 2 * sum(abs(xcorr(h))); % this is an upper bound for || H' * H ||


%% Deconvolution with ISTA using the proposed penalty

t = strcat('SweepPropSNR',num2str(SNR),'.mat');
load(t);
[a,b] = max(mean_srer);
lam_Proposed =  beta_set(b) * sig; % find the best lambda...

gam = 0.9 * min(alp,1) / lam_Proposed ; % this ensures that the denoising step is convex...

%% SGL parameters
alpSGL = 0.95; % this should be in the range (0,1)
t = strcat('SweepSGLSNR',num2str(SNR),'.mat');
load(t);
[a,b] = max(mean_srer);
lamSGL = beta_set(b) * sig;

%% mixed norm parameters
t = strcat('SweepMixedSNR',num2str(SNR),'.mat');
load(t);
[a,b] = max(mean_srer);
lamM =  beta_set(b) * sig;


%% IPS parameters
p_pShrink = -0.5;
t = strcat('SweepIPSSNR',num2str(SNR),'.mat');
load(t);
[a,b] = max(mean_srer);
lam_pShrink =  beta_set(b) * sig; % 0.21 for SNR=5, 0.25 for SNR=10 & SNR=15, 0.32 for SNR=20


%%
disp(['SNR = ' num2str(SNR)])
disp(['No of iterations = ' num2str(MAX_ITER)])
disp(['Group size = ' num2str(Gs)])


snr_store = zeros(NTrial, 8); % store steady-state performance value

snr_flag = 1; % store performance evolution

snr_iter_sgl_store = zeros(MAX_ITER, NTrial); % store transient performance values
snr_iter_mixed_store = zeros(MAX_ITER, NTrial);
snr_iter_pShrink_store = zeros(MAX_ITER, NTrial);
snr_iter_prop_store = zeros(MAX_ITER, NTrial);


for trial = 1:NTrial,
    
    noise = randn(1, xL);
    y = y_noisefree + noise * sig; % the noisy observations
    
    [x_sgl, cost_sgl, snr_iter_sgl] = deconv_SGL(y, h, alp, lamSGL, alpSGL, Gs, MAX_ITER, x, snr_flag); % sparse group lasso
    
    [x_mixed, cost_mixed, snr_iter_mixed] = deconv_SGL(y, h, alp, lamM, 0, Gs, MAX_ITER, x, snr_flag); % l_{2,1} norm
    
    [x_pShrink, snr_iter_pShrink] = deconv_pShrink(y, h, alp, lam_pShrink, p_pShrink, MAX_ITER, x, snr_flag); % iterative p-shrinkage
    
    [x_Proposed, cost_nVariate, snr_iter_prop] = deconv_Proposed(y, h, alp, lam_Proposed, gam, Gs, MAX_ITER, x, snr_flag); % proposed
    
    % debiasing
    x_sglDb = DeBias(y,h,x_sgl);
    
    x_mixedDb = DeBias(y,h,x_mixed);
    
    x_pShrinkDb = DeBias(y,h,x_pShrink);
    
    x_ProposedDb = DeBias(y,h,x_Proposed);
    
    % performance measurement
    snr_sgl = snr(x, x - x_sgl);
    snr_mixed = snr(x, x - x_mixed);
    snr_pShrink = snr(x, x - x_pShrink);
    snr_Prop = snr(x, x - x_Proposed);
    
    snr_sglDb = snr(x, x - x_sglDb);
    snr_mixedDb = snr(x, x - x_mixedDb);
    snr_pShrinkDb = snr(x, x - x_pShrinkDb);
    snr_ProposedDb = snr(x, x - x_Proposed);
    
    snr_store(trial,:) = [snr_sgl, snr_mixed, snr_Prop, snr_pShrink, snr_sglDb, snr_mixedDb, snr_ProposedDb, snr_pShrinkDb];  % store steady-state performance values
    
    snr_iter_sgl_store(:,trial) = snr_iter_sgl;  % store transient performance values
    snr_iter_mixed_store(:,trial) = snr_iter_mixed;
    snr_iter_pShrink_store(:,trial) = snr_iter_pShrink;
    snr_iter_prop_store(:,trial) = snr_iter_prop;
    
end

% performance statistics
mean_snr = mean(snr_store,1);
std_snr = std(snr_store,0,1);


mean_snr_iter_sgl = mean(snr_iter_sgl_store,2);
mean_snr_iter_mixed = mean(snr_iter_mixed_store,2);
mean_snr_iter_ips = mean(snr_iter_pShrink_store,2);
mean_snr_iter_nvt = mean(snr_iter_prop_store,2);

std_snr_iter_sgl = std(snr_iter_sgl_store,0,2);
std_snr_iter_mixed = std(snr_iter_mixed_store,0,2);
std_snr_iter_ips = std(snr_iter_pShrink_store,0,2);
std_snr_iter_nvt = std(snr_iter_prop_store,0,2);


figure
hold
p1 = plot(mean_snr_iter_sgl, 'r');
p1b = plot(mean_snr_iter_sgl + 3*std_snr_iter_sgl,'r--');
p1c = plot(mean_snr_iter_sgl - 3*std_snr_iter_sgl,'r--');

pm = plot(mean_snr_iter_mixed, 'y');
pmb = plot(mean_snr_iter_mixed + 3*std_snr_iter_mixed,'y--');
pmc = plot(mean_snr_iter_mixed - 3*std_snr_iter_mixed,'y--');

p2 = plot(mean_snr_iter_nvt, 'b');
p2b = plot(mean_snr_iter_nvt + 3*std_snr_iter_nvt,'b--');
p2c = plot(mean_snr_iter_nvt - 3*std_snr_iter_nvt,'b--');

p3 = plot(mean_snr_iter_ips, 'g');
p3b = plot(mean_snr_iter_ips + 3*std_snr_iter_ips,'g--');
p3c = plot(mean_snr_iter_ips - 3*std_snr_iter_ips,'g--');


title(['Reliability plots at SNR = ' num2str(SNR) 'dB' ])
xlabel('Iteration index')
ylabel('SRER (dB)')
legend([p1,pm,p2,p3],'SGL','Mixed','Proposed','IPS')

