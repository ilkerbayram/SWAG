function [c] = ComputeCostSGL(y,z,h,alp,lam)
% computes cost for DemoDeconv
c = 0;
% the penalty term
c = sum( lam * (1-alp) * sqrt(sum(abs(z).^2,1)) + lam * alp * sum(abs(z),1));
% data fidelity term
z = reshape(z,size(y,1),size(y,2));
c = c + 0.5 * sum( abs(y - convH(z,h)).^2);