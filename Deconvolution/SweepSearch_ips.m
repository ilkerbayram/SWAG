%% Sparse deconvolution with IPS
%% Find regularization parameter lambda with a sweep search

clear all;
close all;
clc;

savedir = cd; % save into the current folder...
%% the wavelet
load('ricker.mat') ;
wL = length(h); % length of the blur function

%% the sparse signal
load('sparse_sig.mat')
xL = length(x);

%% produce the noise-free observation
y_noisefree = convH(x,h); % noise-free observations

MAX_ITER = 1e3; % max num of iterations.

NTrial = 50; % number of trials, each having a different noise realization
%% IPS parameters
alp = 2 * sum(abs(xcorr(h))); % this is an upper bound for || H' * H ||

p_pShrink = -0.5; % the p-parameter

%
beta_min = 0.05;
beta_max = 1.5;
beta_set = exp(linspace(log(beta_min), log(beta_max), 20));

SNRlist = 5:5:20;

for SNR = SNRlist,
    sig = sqrt( sum(abs(y_noisefree).^2) / (xL * 10^(SNR/10))); % noise std
    
    disp(['input SNR = ' num2str(SNR)])
    disp(['No of iterations = ' num2str(MAX_ITER)])
    
    mean_srer = zeros(length(beta_set), 1);
    std_srer = zeros(length(beta_set), 1);
    
    tic
    for beta_index = 1:length(beta_set)
        
        beta = beta_set(beta_index);
        
        lam = beta * sig;
        
        disp(['beta = ' num2str(beta)])
        
        % Change to choose another realization of gaussian noise
        seed = 1 ;
        randn('state',seed);
        
        srer_store = zeros(NTrial,1);
        
        for trial = 1:NTrial,
            
            noise = randn(1, xL);
            y = y_noisefree + noise * sig; % the noisy observations
            
            [x_pShrink, snr_iter_pShrink] = deconv_pShrink(y, h, alp, lam, p_pShrink, MAX_ITER, x, 0);
            
            srer_store(trial) = snr(x, x - x_pShrink);
            
        end
        
        mean_srer(beta_index) = mean(srer_store);
        std_srer(beta_index) = std(srer_store);
        
        toc
    end
    
    sd = strcat(savedir,'\SweepIPSSNR',num2str(SNR),'.mat');
    save(sd, 'srer_store', 'beta_set', 'mean_srer', 'std_srer');
end