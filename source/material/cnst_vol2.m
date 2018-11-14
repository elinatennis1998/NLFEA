%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,P,var]=cnst_vol(F,var,mat,ndim)
% NOTE: THIS ONE DOES NOT WORK YET; IT WAS AN ATTEMPT TO CONVERT A TO C
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% determine tangent, stress and internal variable 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
emod = mat(1);   nue = mat(2);   kt = mat(3);  kc = mat(4);
mt   = mat(5);   mc  = mat(6);   tt = mat(7);  tc = mat(8);
dt   = mat(9);

xmu = emod / 2.0 / (1.0+nue);
xlm = emod * nue / (1.0+nue) / ( 1.0-2.0*nue );

%%% update internal variable %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[var,fac1,fac2]=updt_vol2(F,var,mat,ndim);

theta  = var(1) + 1;
delta  = eye(ndim);
Fe = F / theta;
Fe_inv = inv(Fe);
Je     = det(Fe);
G_inv = delta;
g = delta;
% C = F'*F;
Ce_inv = Fe_inv*Fe_inv';
  
%%% first piola kirchhoff stress %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S = (xmu * G_inv + (xlm * log(Je) - xmu) * Ce_inv);
P =  Fe * S; 
Pe   = xmu * Fe +            (xlm*log(Je) - xmu) * Fe_inv';

%%% tangent %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = zeros(ndim,ndim,ndim,ndim);
Cmat = A;
CeCe = xmu * Fe + (ndim*xlm - xlm*log(Je) + xmu) * Fe_inv';
for I=1:ndim
for J=1:ndim
for K=1:ndim
for L=1:ndim
  Cmat(I,J,K,L) =  (xlm                 * Ce_inv(I,J)*Ce_inv(K,L) ...
			    - (xlm * log(JF) - xmu) * (Ce_inv(I,K)*Ce_inv(J,L) + Ce_inv(I,L)*Ce_inv(J,K))) ...
                + fact * S(I,J)*S(K,L);		 
end, end, end, end

for i=1:ndim
for I=1:ndim
for k=1:ndim
for K=1:ndim
  A(i,I,k,K) =  S(I,K)*g(i,k) + F(i,1)*Cmat(1,I,1,K)*F(k,1) ...
             + F(i,2)*Cmat(2,I,1,K)*F(k,1) + F(i,3)*Cmat(3,I,1,K)*F(k,1) ...
             + F(i,1)*Cmat(1,I,2,K)*F(k,2) + F(i,2)*Cmat(2,I,2,K)*F(k,2) ...
             + F(i,3)*Cmat(3,I,2,K)*F(k,2) + F(i,1)*Cmat(1,I,3,K)*F(k,3) ...
             + F(i,2)*Cmat(2,I,3,K)*F(k,3) + F(i,3)*Cmat(3,I,3,K)*F(k,3);	 
end, end, end, end
%%% tangent %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pe   = xmu * Fe +            (xlm*log(Je) - xmu) * Fe_inv';
% AeFe = xmu * Fe + (ndim*xlm - xlm*log(Je) + xmu) * Fe_inv';
% 
% ten1 = fac1 *  AeFe;
% ten2 = fac2 * (Pe+AeFe);
for i=1:ndim
for j=1:ndim
for k=1:ndim
for l=1:ndim 
  A(i,j,k,l) =  xlm                  * Fe_inv(j,i)*Fe_inv(l,k) ...
			 - (xlm * log(Je) - xmu) * Fe_inv(l,i)*Fe_inv(j,k) ...
			 +                  xmu  *  delta(i,k)* delta(j,l) ...
             +                           ten1(i,j)*  ten2(k,l);			 
end, end, end, end
A = A / theta; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

