function [c] = ComputeCost(y,z,h,lam,gam)
% computes cost for the proposed formulation of deconvolution
c = 0;
% the penalty term
c = sum( sum(abs(z),1).^2 - sum(abs(z).^2,1) );
c = lam * ( gam * c / 2 + sum(abs(z(:))) );
% data fidelity term
z = reshape(z,size(y,1),size(y,2));
c = c + 0.5 * sum( abs(y - convH(z,h)).^2);