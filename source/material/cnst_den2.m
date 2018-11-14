%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,P,var]=cnst_den2(F,var,mat,ndim)
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
JF     = det(F);
G_inv = eye(ndim);
g = eye(ndim);
% C = F'*F;
C_inv = F_inv*F_inv';
  
%%% first piola kirchhoff stress %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S = facs * (xmu * G_inv + (xlm * log(JF) - xmu) * C_inv);
P =  F * S; 
   
%%% tangent %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = zeros(ndim,ndim,ndim,ndim);
Cmat = A;
for I=1:ndim
for J=1:ndim
for K=1:ndim
for L=1:ndim
  Cmat(I,J,K,L) =  facs * (xlm                 * C_inv(I,J)*C_inv(K,L) ...
			    - (xlm * log(JF) - xmu) * (C_inv(I,K)*C_inv(J,L) + C_inv(I,L)*C_inv(J,K))) ...
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
