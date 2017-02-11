function [x] = convHT(y,h)
% this is the transpose of convH

x = conv(y,h(end:-1:1));
x = x(length(h):end);