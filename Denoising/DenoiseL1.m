function [x] = DenoiseL1(y,param)
% input parameters : 
%
% y : noisy signal in the time domain
% param.win : window for the STFT
% param.hop : hop-size for the STFT
% param.MAX_ITER : maximum number of iterations
% param.K : group-size for E-Lasso, and the proposed methods
% param.lam : lambda parameter for the regularizer
% param.Gs : group-size for the mixed norm
%
% output parameters :
%
% x : denoised output

SoftThold = @(x,tau) max( abs(x) - tau, 0) .* x ./ (abs(x) + 10^(-10)); % soft threshold function

bet = 0.5; % beta parameter used in the Douglas-Rachford algorithm

Y = STFT(y,param.win,param.hop);
T = Y;
U = zeros(size(Y));
Z = zeros(size(Y));

wb = waitbar(0,'Please wait');

for iter = 1:param.MAX_ITER,
    waitbar(iter/param.MAX_ITER,wb)
    
    % project T onto the range of S
    u = ISTFT(T,param.win,param.hop);
    u = real(u(1:length(y)));
    U = STFT(u,param.win,param.hop);
    
    % proximal step for the penalty
    Z = ( 2*U - T + bet*Y )/(1+bet);
    
    % threshold step
    Z = SoftThold(Z,param.lam * bet / (1+bet));
    
    % update T
    T = Z + T - U;
end

close(wb)

x = ISTFT(T,param.win,param.hop);
x = real(x(1:length(y)));