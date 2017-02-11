function [win2] = NormalizeW(win,Hop)
% Normalizes the window so that the STFT is a Parseval frame
% win : input window
% Hop : Hop-size
% 
% Matlab code by Ilker Bayram, 2011

N = length(win);
K = floor(N/Hop);
win2 = win.^2;
z = win2;
for n = 1:K,%shift to the left
    z(1:end-n*Hop) = z(1:end-n*Hop) + win2(n*Hop+1:end);
end
for n = 1:K,%shift to the right
    z(n*Hop+1:end) = z(n*Hop+1:end) + win2(1:end-n*Hop);
end
win2 = win./sqrt(z);