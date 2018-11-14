%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [var,fac1,fac2]=updt_vol2(F,var,mat,ndim)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% volume update for functional adaption
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input:  F        -> deformation gradient
%         var      -> var = theta    ... internal variable growth factor
%         mat      -> material parameters
% output: var      -> var = theta    ... internal variable growth factor
%         ten1     -> tensor 1 for tangent operator    
%         ten2     -> tensor 2 for tangent operator    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Equation numbers (#) refer to Computational Modelling of Isotropic Multiplicative Growth
% Himpel, Kuhl, Menzel, Steinmann, CMES, vol.8, no.2, pp.119-134, 2005
%
tol = 1e-10;      

emod = mat(1);   nue = mat(2);   kt = mat(3);  kc = mat(4);
mt   = mat(5);   mc  = mat(6);   tt = mat(7);  tc = mat(8);
dt   = mat(9);

xmu = emod / 2.0 / (1.0+nue);
xlm = emod * nue / (1.0+nue) / ( 1.0-2.0*nue ); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% euler backward - implicit time integration (quad.conv.@large str)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
the_k0 = var(1) + 1;
the_k1 = var(1) + 1;

iter = 0;   res  = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%% local newton-raphson iteration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    while abs(res) > tol
      iter=iter+1;
  	  
	  Fe = F / the_k1; % (6)
	  Ce = Fe' * Fe; % (7)
	  delta  = eye(ndim); 
	  Fe_inv = inv(Fe);
	  Ce_inv = inv(Ce);
	  Je = det(Fe);
      
      Se = xmu * delta + (xlm * log(Je) - xmu) * Ce_inv; % (61)	 
      Me = Ce*Se; % (53)
	  tr_Me = trace(Me);
	  
      CeLeCe = ndim * ndim * xlm - 2 * ndim * (xlm * log(Je) - xmu); % (54), after evaluation of the model given by (61)
	  dtrM_dthe = - 1/the_k1 * ( 2*tr_Me + CeLeCe ); % (54)

      if tr_Me > 0
	    k       = kt*((tt-the_k1)/(tt-1))^mt; % (35)
	    dk_dthe = k  /(the_k1-tt)        *mt; % (53)
      else 
	    k       = kc*((the_k1-tc)/(1-tc))^mc; % (35)
	    dk_dthe = k  /(the_k1-tc)        *mc; % (53)
      end
      
      res  = k * tr_Me * dt - the_k1 + the_k0; % (53)
	  dres =(dk_dthe * tr_Me + k * dtrM_dthe)*dt -1; % (52)
	  
	  dthe   = -res/dres;
	  the_k1 =  the_k1 + dthe;
      
	  if(iter>20); disp(['*** NO LOCAL CONVERGENCE ***']); return; end;      
    end
%%% local newton-raphson iteration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
if (abs(dres)<1e-12 | abs(the_k1-1)<1e-12)
  fac1 = 0.0;            
  fac2 = 0.0;
else
  fac1 = -1 / the_k1; % (46) 
  fac2 = -k / dres * dt; % (51,56)
end

% Pe   = xmu * Fe +            (xlm*log(Je) - xmu) * Fe_inv';
% AeFe = xmu * Fe + (ndim*xlm - xlm*log(Je) + xmu) * Fe_inv';
% 
% ten1 = fac1 *  AeFe;
% ten2 = fac2 * (Pe+AeFe);

var(1) = the_k1 - 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sca  = Ce * Le * Ce;
% ten1 = -1/the_k1  * Ae * Fe;
% ten2 = -k/dres*dt * dtrMdC * dCdF;
% ten2 = -k/dres*dt * Fe * [2 Se + Ce:Le] % with Fe*Se=Pe
% ten2 = -k/dres*dt *[Fe * Se + Fe:Ae]