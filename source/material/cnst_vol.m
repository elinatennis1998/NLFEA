%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,P,var]=cnst_vol(F,var,mat,ndim)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% determine tangent, stress and internal variable 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
emod = mat(1);   nue = mat(2);   kt = mat(3);  kc = mat(4);
mt   = mat(5);   mc  = mat(6);   tt = mat(7);  tc = mat(8);
dt   = mat(9);

xmu = emod / 2.0 / (1.0+nue);
xlm = emod * nue / (1.0+nue) / ( 1.0-2.0*nue ); 

%%% update internal variable %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[var,ten1,ten2]=updt_vol(F,var,mat,ndim);

theta  = var(1) + 1;
Fe     = F / theta;
Fe_inv = inv(Fe);
Je     = det(Fe);
delta  = eye(ndim); 

%%% first piola kirchhoff stress %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P = xmu * Fe + (xlm * log(Je) - xmu) * Fe_inv';

%%% tangent %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

