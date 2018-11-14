%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,P,var]=cnst_den3(F,var,mat,ndim)
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
g = eye(ndim);
b = F*F';
  
%%% first piola kirchhoff stress %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = facs * (xmu * b + (xlm * log(JF) - xmu) * g);
P =  t * F_inv'; 
   
%%% tangent %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = zeros(ndim,ndim,ndim,ndim);
c = A;
for i=1:ndim
for j=1:ndim
for k=1:ndim
for l=1:ndim
  c(i,j,k,l) =  facs * (xlm                 * g(i,j)*g(k,l) ...
	         - (xlm * log(JF) - xmu) * (g(i,k)*g(j,l) + g(i,l)*g(j,k))) ...
             + fact * t(i,j)*t(k,l);		 
end, end, end, end

for i=1:ndim
for I=1:ndim
for k=1:ndim
for K=1:ndim
  A(i,I,k,K) = F_inv(I,2)*g(i,k)*t(2,1)*F_inv(K,1) + F_inv(I,3)*g(i,k)*t(3,1)*F_inv(K,1) ...
             + F_inv(I,1)*g(i,k)*t(1,2)*F_inv(K,2) + F_inv(I,2)*g(i,k)*t(2,2)*F_inv(K,2) ...
             + F_inv(I,3)*g(i,k)*t(3,2)*F_inv(K,2) + F_inv(I,1)*g(i,k)*t(1,3)*F_inv(K,3) ...
             + F_inv(I,2)*g(i,k)*t(2,3)*F_inv(K,3) + F_inv(I,3)*g(i,k)*t(3,3)*F_inv(K,3) ...
             + F_inv(I,1)*g(i,k)*t(1,1)*F_inv(K,1) + F_inv(I,1)*c(1,i,1,k)*F_inv(K,1) ...
             + F_inv(I,2)*c(2,i,1,k)*F_inv(K,1) + F_inv(I,3)*c(3,i,1,k)*F_inv(K,1) ...
             + F_inv(I,1)*c(1,i,2,k)*F_inv(K,2) + F_inv(I,2)*c(2,i,2,k)*F_inv(K,2) ...
             + F_inv(I,3)*c(3,i,2,k)*F_inv(K,2) + F_inv(I,1)*c(1,i,3,k)*F_inv(K,3) ...
             + F_inv(I,2)*c(2,i,3,k)*F_inv(K,3) + F_inv(I,3)*c(3,i,3,k)*F_inv(K,3);	 
end, end, end, end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
