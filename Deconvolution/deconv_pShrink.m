function [z, snr_iter] = deconv_pShrink(y, h, alp, lam, p, Nit, x, snr_flag)
% [z] = deconv_pShrink(y, h, alp, lam, p, Nit)
% Sparse deconvolution using the pShrink threshold function
%
% INPUT
%   y - noisy data
%   h - filter impulse response
%   alp - stepsize for the algorithm
%   lam - weight of the regularization parameter
%   p - the parameter for the threshold
%   Nit - number of iterations
%   x - clean sparse signal (reference for performance measurement)
%   snr_flag - flag for storing performance evolution with iterations
%
% OUTPUT
%   x - deconvolved sparse signal
%   snr_iter - performance evolution with iterations
% 
% I. Bayram and S. Bulek, 2016


Thold = @(x,tau,p) sign(x) .* max( abs(x) - tau^(2-p) * abs(x).^(p-1), 0);


z = zeros(size(y)); % z will hold the reconstruction
snr_iter = zeros(Nit,1);


for iter = 1 : Nit,
    
    z = z + (1/alp) * convHT(y - convH(z,h), h);
    z = Thold(z,lam,p);
    
    if snr_flag
        snr_iter(iter) = snr(x, x - z);
    end
    
    
end
