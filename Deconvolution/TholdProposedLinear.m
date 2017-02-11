function [x] = TholdProposedLinear(z,lam,gam)
% applies the proposed threshold to each column of z
% in this implementation, a linear search for 'k' is employed

eps = 10^(-10);

N = size(z,1); % N-variate threshold

mag = sort(abs(z),1,'descend'); % sort the magnitudes in descending order
S = cumsum(mag,1);
lg = lam*gam;
cnst = lam*(1 - lg);
x = z;
for c = 1:size(z,2),
    h = lam;
    k = 0;
    while k < N,
        if mag(k+1,c) > h,
            k = k + 1;
            h = (cnst + lg * S(k,c)) / (1 + (k-1) * lg );
        else            
            break;
        end
    end
    col = z(:,c);
    col = max(1 - h ./ (abs(col) + eps), 0) .* col / (1-lg);
    x(:,c) = col;
end