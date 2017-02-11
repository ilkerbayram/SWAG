function [y] = ISTFT(g,win,Hop)
% Inverse STFT
% g : input STFT (output of the function STFT.m)
% win : Window for STFT
% Hop : Hop-Size

% Ilker Bayram
% ibayram@itu.edu.tr
% Istanbul Teknik Universitesi, 2011

win = win(:).';

c = ifft(g)*sqrt(length(win));
K = size(g,2) - 1;
W = size(g,1);
y = zeros(K*Hop + W,1);

c = c.';

for k = 1:size(g,1),
    y(k:Hop:K*Hop+k) = y(k:Hop:K*Hop+k) + win(k)*c(:,k);
end