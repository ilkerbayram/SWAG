function [x] = TholdMixed(z,lam)
% applies the l2 norm threshold to each column of z

eps = 10^(-10);

mz = sqrt(sum(abs(z).^2));
mz2 = max(0, mz - lam);

mz2 = mz2 ./ (mz + eps);

x = bsxfun(@times, z, mz2);