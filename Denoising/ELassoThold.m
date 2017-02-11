function [x] = ELassoThold(z,lam)
% applies the E-lasso threshold to each column of z

eps = 10^(-10);

N = size(z,1); % N-variate threshold

mag = sort(abs(z),1,'descend'); % sort the magnitudes in descending order
S = cumsum(mag,1);

x = z;
for c = 1:size(z,2),
    h = 2 * lam * S(1,c);
    k = 1;
    while k < N,
        if mag(k+1,c) > h,
            k = k + 1;
            h = 2 * lam * S(k,c);
        else            
            break;
        end
    end
    col = z(:,c);
    col = max(1 - h ./ (abs(col) + eps), 0) .* col;
    x(:,c) = col;
end