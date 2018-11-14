%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [var,facs,fact]=updt_den(F,var,mat,ndim)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% density update for functional adaption
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input:  F        -> deformation gradient
%         var      -> var = rho / rho_0 - 1 internal variable density
%         mat      -> material parameters
% output: var      -> var = rho / rho_0 - 1 internal variable density
%         facs     -> factor for stresses    
%         fact     -> factor for tangent operator    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tol = 1e-8;1e-12;1e-17;      

emod = mat(1);   nue  = mat(2);   rho0 = mat(3);   psi0 = mat(4);   
expm = mat(5);   expn = mat(6);   dt   = mat(7);

xmu = emod / 2.0 / (1.0+nue);
xlm = emod * nue / (1.0+nue) / ( 1.0-2.0*nue ); 

J  = det(F);
% if abs(J-1) > 1e-14
% C  = F'*F;
% J
% F(1,1)*F(2,2)*F(3,3)
% J-1
% F(1,1)*F(2,2)*F(3,3)-1
% I1 = trace(C)
% t1 = xlm/2 * log(J)^2
% t2 = I1 - ndim - 2*log(J)
% t3 = xmu/2 * t2
% psi0_neo = t1 + t3
% else
C  = F'*F;
I1 = trace(C);
psi0_neo = xlm/2 * log(J)^2  +  xmu/2 *(I1 - ndim - 2*log(J));
% end
       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% euler backward - implicit time integration (quad.conv.@large str)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if abs(psi0_neo) > tol*1e3 
          
rho_k0 = (1 + var) * rho0;
rho_k1 = (1 + var) * rho0; 

iter = 0;   res  = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%% local newton-raphson iteration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    while abs(res) > tol
      iter=iter+1;
  	  
      res =rho_k1-rho_k0-((rho_k1/rho0)^(expn-expm)*psi0_neo-psi0)*dt; 
	  dres=1-(expn-expm)*(rho_k1/rho0)^(expn-expm)*psi0_neo*dt/rho_k1;
	  
      drho   =-res/dres;  
	  rho_k1 = rho_k1+drho; 
      if iter == 1
          tol = tol*abs(res);
      end
      
	  if(iter>20); disp(['*** NO LOCAL CONVERGENCE ***']); return; end;      
    end
%%% local newton-raphson iteration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
rho = rho_k1;
else
rho = rho0;
end

var = rho / rho0 - 1;
facs= (rho/rho0)^expn;

if (dt~=0)
  dres=  1    - (expn-expm) * (rho/rho0)^(expn-expm) / rho *psi0_neo*dt ;
  fact=          expn / rho * (rho/rho0)^(    -expm) * dt / dres;
else
  fact = 0.0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
