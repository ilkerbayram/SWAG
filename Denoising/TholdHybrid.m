function [x] = TholdHybrid(z,lam,gam,K)
% hybrid threshold function obtained by applying the proposed threshold to
% the norms of the groups in each supergroup of size K

eps = 10^(-10);

mz = sqrt(sum(abs(z).^2));
mz2 = reshape(mz,K,length(mz(:)) / K);

mz2 = TholdProposedLinear(mz2,lam,gam);
mz2 = reshape(mz2,1,length(mz(:)));
mz2 = mz2 ./ (mz + eps);

x = bsxfun(@times, z, mz2);