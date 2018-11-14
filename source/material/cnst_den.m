%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,P,var]=cnst_den(F,var,mat,ndim)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% determine tangent, stress and internal variable 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

emod = mat(1);   nue  = mat(2);   rho0 = mat(3);   psi0 = mat(4);   
expm = mat(5);   expn = mat(6);   dt   = mat(7);

xmu = emod / 2.0 / (1.0+nue);
xlm = emod * nue / (1.0+nue) / ( 1.0-2.0*nue ); 

%%% update internal variable %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[var,facs,fact]=updt_den(F,var,mat,ndim);

F_inv = inv(F);
J     = det(F);
delta = eye(ndim);
  
%%% first piola kirchhoff stress %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P =  facs * (xmu * F + (xlm * log(J) - xmu) * F_inv'); 
   
%%% tangent %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:ndim
for j=1:ndim
for k=1:ndim
for l=1:ndim
  A(i,j,k,l) =  xlm                 * F_inv(j,i)*F_inv(l,k) ...
			 - (xlm * log(J) - xmu) * F_inv(l,i)*F_inv(j,k) ...
			 +                 xmu  * delta(i,k)*delta(j,l);
  A(i,j,k,l) = facs * A(i,j,k,l)    +  fact * P(i,j)*P(k,l);			 
end, end, end, end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
