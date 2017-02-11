function [x] = DeBias(y,h,x)
% Debiasing

th = 10^(-5);
ind = abs(x) > th;

y = y(:);
x = x(:);
h = h(:);

H = convmtx(h,length(x));
H = H(1:length(x),:);
H = H(:,ind);

z = (H' * H + 0.001 * eye(sum(ind))) \ (H' * y) ;

x(ind) = z;
x = x';