function [g] = STFT(x,win,Hop)
% STFT
% x : input signal
% win : Window for STFT
% Hop : Hop-Size

% Ilker Bayram
% ibayram@itu.edu.tr
% Istanbul Teknik Universitesi, 2011


% make the inputs column-vectors
win = win(:);
x = x(:);

N = length(x);
W = length(win); 
K = floor((N - W)/Hop);
len2 = W + K*Hop;
K = K+1;
if len2 < N,
    x = [x; zeros(W - (N - len2),1)];
    K = K + 1;
end

g = zeros(W,K);
for k = 1:W,
    g(k,:) = win(k).*x(k:Hop:k + Hop*(K-1));
end
g = fft(g)/sqrt(length(win));