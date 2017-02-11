function [z, cost, snr_iter] = deconv_Proposed(y, h, alp, lam, gam, Gs, Nit, x, snr_flag)
% [x, cost] = deconv_L1(y, lam, h, Nit)
% Sparse deconvolution with the nvariate penalty
% Cost function : 0.5 * sum(abs(y-conv(h,x)).^2) + lam * sum(abs(x));
%
% INPUT
%   y - noisy data
%   h - filter impulse response
%   alp - the step-size parameter
%   lam, gam - weights for the regularization function
%   Gs - group size
%   Nit - number of iterations
%   x - clean sparse signal (reference for performance measurement)
%   snr_flag - flag for storing performance evolution with iterations
%
% OUTPUT
%   x - deconvolved sparse signal
%   cost - cost function history
%   snr_iter - performance evolution with iterations
% 
% I. Bayram and S. Bulek, 2016

xL = length(y);

z = zeros(size(y)); % z will hold the reconstruction
cost = zeros(1,Nit);
snr_iter = zeros(Nit,1);


for iter = 1:Nit,
    
    z = z + (1/alp) * convHT(y - convH(z,h), h);
    z = reshape(z,Gs,xL/Gs);
    z = TholdProposedLinear(z,lam/alp,gam);
    cost(iter) = ComputeCost(y,z,h,lam,gam);
    z = reshape(z,1,xL);
    
    if snr_flag
        snr_iter(iter) = snr(x, x - z);
    end
      
    if (iter > 1),
        if (cost(iter) > (cost(iter-1) + 10^(-5))),
            disp('Unexpected increase in cost');
        end
    end
    
end
