% Demo file for the experiment showing the dependence on the product
% gamma * lambda as the number of non-zeros in the signal changes (fig 4. in the manuscript)

clear all;
close all;

L = 10; % length of the signal
N = 10000; % number of trials

Klist = 1:5; % list of K values to consider

indgam = 1.05.^(-100:2:-1);

SNRin = 5; % input SNR

SNRlist = zeros(length(indgam),50); % output SNR wrt varying (lambda * gamma)


n = 0;
for K = Klist, % repeat for each K value
    n = n+1;
    
    % produce the clean signal
    x = [randn(K,N); zeros(L-K,N) ];
    
    % determine the signal std to achieve observations with input SNR = SNRin
    ex = sum(abs(x(:)).^2);
    en = ex / 10^(SNRin/10);
    sig = sqrt(en / length(x(:)));
    
    % the noisy observations
    y = x + sig * randn(size(x));
    
    % lambda is fixed to sigma / 2
    lam = sig/2;
    
    gamlist = (1/lam) * indgam;

    k = 0;
    for gam = gamlist,
        k = k+1;
        
        z = TholdProposedLinear(y,lam,gam);
        
        SNRlist(n,k) = snr(x,z - x) - SNRin;
    end    
end


%% display

figure;
for n = 1:5,
    plot(indgam(:),SNRlist(n,:),'color',rand(1,3)); hold on;
end
axis([0,0.9, 0, 8]);
legend('K=1','K=2','K=3','K=4','K=5','Location','East');
xlabel('\lambda \gamma');
ylabel('SNR Gain (dB)');
