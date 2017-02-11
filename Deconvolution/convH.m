function [x] = convH(y,h)
% applies H to y, produces a signal of the same size

x = conv(y,h);
x = x(1:length(y));