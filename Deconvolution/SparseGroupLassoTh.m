function [x] = SparseGroupLassoTh(y,alp,lam)
% applies the Sparse-Group Lasso threshold function to each column of y

x = sign(y) .* max(abs(y) - alp * lam, 0);
mag = sqrt(sum(abs(x).^2));
nx = max( mag - (1-alp)*lam, 0);
nx = nx ./ (mag + 10^(-10));
x = bsxfun(@times,x,nx);