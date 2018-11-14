%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [var,Pe,Aeg]=updt_cyl(F,var,mat,ndim)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cylindrical volume update for functional adaption
% Tim Truster 10/1/2013
% uses "IJSS39-Lubarda-Mechanics of solids, growing mass".pdf
% Equation numbers correspond to AF report dated 2012-10-14
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input:  F        -> deformation gradient
%         var      -> var = theta    ... internal variable growth factor
%         mat      -> material parameters
% output: var      -> var = theta    ... internal variable growth factor
%         ten1     -> tensor 1 for tangent operator    
%         ten2     -> tensor 2 for tangent operator    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tol = 1e-10;      

emod = mat(1);   nue = mat(2);   kt = mat(3);  kc = mat(4);
mt   = mat(5);   mc  = mat(6);   tt = mat(7);  tc = mat(8);
dt   = mat(9);   mo  = mat(10:12); %fiber orientation vector

xmu = emod / 2.0 / (1.0+nue);
xlm = emod * nue / (1.0+nue) / ( 1.0-2.0*nue ); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% euler backward - implicit time integration (quad.conv.@large str)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
the_k0 = var(1) + 1;
the_k1 = var(1) + 1;

[~, mo] = VecNormalize(mo); % ensure mo is a unit vector
mo = mo';
C = F' * F;
mCmpmCm = mo*(C*mo)' + (C*mo)*mo';
mCm = mo'*C*mo;
mCmmm = mCm*(mo*mo');
one = [1; 1; 1; 0; 0; 0; 0; 0; 0];

F_m = F*mo;

iter = 0;   res  = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%% local newton-raphson iteration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    while abs(res) > tol
      iter=iter+1;
  	  
	  delta  = eye(ndim); 
      Fg_inv = (1/the_k1)*delta + (1-1/the_k1)*(mo*mo'); % (2)
	  Fe = F * Fg_inv; % (1)
	  Ce = Fe' * Fe;
	  Fe_inv = inv(Fe);
	  Je = det(Fe);
      
      % compute Se, Cijkl_e
      [taue, cmat] = SigmaCmat3(Fe,Je,[1 2 0 emod nue],xlm); % taue is the Kirchhoff stress tensor (i.e. no Jacobian)
      P = tranr4(Fe_inv,Fe_inv); % transformation tensor F_iI*F_jJ
      P = [P' zeros(6,3); zeros(3,9)];
      Se = P*taue; % second PK stress
      Se_m = [Se(1) Se(4) Se(6); Se(4) Se(2) Se(5); Se(6) Se(5) Se(3)];
      Cmat_e = P*cmat*P';
      Pe = Fe*Se_m; % first PK stress % (5)
      Ae = CSFtoA(Cmat_e,Se_m,Fe,ndim);
      Me = Ce*Se_m; % (53)
	  tr_Me = trace(Me);
      mMem = mo'*Me*mo;
      
      dtrMe_dFe = zeros(ndim);
      dmMem_dFe = zeros(ndim);
      Fe_m = Fe*mo;
      Pe_m = Pe*mo;
      for k = 1:ndim
          for K = 1:ndim
              Ae_lLkK_Fe_lL = 0;
              Ae_lLkK_Fe_lM_m_M_m_L = 0;
              for l = 1:ndim
                  for L = 1:ndim
                      Ae_lLkK_Fe_lL = Ae_lLkK_Fe_lL + Ae(l,L,k,K)*Fe(l,L);
                      Ae_lLkK_Fe_lM_m_M_m_L = Ae_lLkK_Fe_lM_m_M_m_L + Ae(l,L,k,K)*Fe_m(l)*mo(L);
                  end
              end
              dtrMe_dFe(k,K) = Pe(k,K) + Ae_lLkK_Fe_lL; % (8)
              dmMem_dFe(k,K) = Pe_m(k)*mo(K) + Ae_lLkK_Fe_lM_m_M_m_L; % (9)
          end
      end
      dtrMe_dthe = 1/the_k1^2*(-trace(dtrMe_dFe'*F) + F_m'*dtrMe_dFe*mo); % (11)
      dmMem_dthe = 1/the_k1^2*(-trace(dmMem_dFe'*F) + F_m'*dmMem_dFe*mo); % (11)

      if (tr_Me - mMem) > 0 % (4)
	    kv      = kt*((tt-the_k1)/(tt-1))^mt;
	    dk_dthe = kv  /(the_k1-tt)        *mt;
      else 
	    kv      = kc*((the_k1-tc)/(1-tc))^mc;
	    dk_dthe = kv  /(the_k1-tc)        *mc;                    
      end
      
      res  = kv * (tr_Me - mMem) * dt - the_k1 + the_k0; % (3+6)
	  dres =(dk_dthe * (tr_Me - mMem) + kv * (dtrMe_dthe - dmMem_dthe))*dt -1; % (7)
	  
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
  fac1 = 1;    
  fac2 = -kv / dres * dt; 
  if iter>7
      iter
  end
end

if abs(the_k1-1)<1e-12
    
Aeg = Ae;
    
else
    
dFe_dthe = 1/the_k1^2*F*(-delta + mo'*mo); % (10)
AeFe = zeros(ndim);
for i = 1:ndim
    for I = 1:ndim
        AeFejJ = 0;
        for j = 1:ndim
            for J = 1:ndim
                AeFejJ = AeFejJ + Ae(i,I,j,J)*dFe_dthe(j,J);
            end
        end
        AeFe(i,I) = AeFejJ; % (14)
    end
end

ten1 = fac1 *  AeFe; % (16)
ten2 = fac2 * ((dtrMe_dFe - dmMem_dFe)*Fg_inv); % (17)

Aeg = zeros(ndim,ndim,ndim,ndim);
for i = 1:ndim
    for I = 1:ndim
        for j = 1:ndim
            for J = 1:ndim
                Ae_iIjK_Fg_inv_JK = 0;
                for K = 1:ndim
                    Ae_iIjK_Fg_inv_JK = Ae_iIjK_Fg_inv_JK + Ae(i,I,j,K)*Fg_inv(J,K);
                end
                Aeg(i,I,j,J) = Ae_iIjK_Fg_inv_JK + ten1(i,I)*ten2(j,J);
            end
        end
    end
end

end

var(1) = the_k1 - 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sca  = Ce * Le * Ce;
% ten1 = -1/the_k1  * Ae * Fe;
% ten2 = -k/dres*dt * dtrMdC * dCdF;
% ten2 = -k/dres*dt * Fe * [2 Se + Ce:Le] % with Fe*Se=Pe
% ten2 = -k/dres*dt *[Fe * Se + Fe:Ae]