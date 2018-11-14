%
% Driver to compute the reference and spatial material moduli by finite differencing
%
% From Sanjay's website on FEAP
% Modified 06/16/2014 to compute the numerical 6th order tensor
% It works, I think.
clear all;
format short e;

% Sample total strain tensor
plasversion = 2;

% if plasversion == 1
% strain = [
%    0.003481952464967
%    0.002811361626954
%                    0
%   -0.002486099946961
%                    0
%                    0]*2;
% mu = 28;
% bulk = 46.666666666666667;
% Khard = 2;
% Hhard = 1;
% sigy = 0.2;
% ep_n = zeros(6,1);
% beta_n = zeros(6,1);
% % beta_n(2) = -.1;
% a_n = 0.11;
%                 
% [sig,D] = J2RadialReturn0(strain,ep_n,beta_n,a_n,mu,bulk,Khard,Hhard,sigy);
% 
% else
    
deps_n1 = 1.0e-03 *[
   
  -0.657595020273064
   0.763239674371468
  -0.127234334914105
                   0];
ee_n = [
  -0.001315190040546
   0.001526479348743
  -0.000254468669828
                   0];
a_n = 0;
strain = ee_n + deps_n1;
ElemYM = 70;
Elemv = 0.25;
sigy = 0.2;
Khard = 0;1;
plasmodel = 1;

if plasmodel == 1
[DGAMA,LALGVA,RSTAVA,STRES] = SUVM(2,2,[ElemYM,Elemv],[0,100;sigy,Khard*100+sigy],[ee_n; a_n],strain);
C_n1 = CTVM(DGAMA,LALGVA(1),2,2,[ElemYM,Elemv],[0,100;sigy,Khard*100+sigy],RSTAVA,STRES);
elseif plasmodel == 2
[DGAM,LALGVA,RSTAVA,STRES] = SUTR(2,2,[ElemYM,Elemv],[0,100;sigy,Khard*100+sigy],[ee_n; a_n],strain);
C_n1 = CTTR(LALGVA(1),2,LALGVA,2,[ElemYM,Elemv],[0,100;sigy,Khard*100+sigy],RSTAVA,eps_n1,STRES);
end
sig = STRES;
D = C_n1;

% end
% [sig,d] = model(F);

D = D(1:3,1:3);

% pull back 6th order tensor
% Eref = computed value that I don't know


% Perturb right Cauchy-Green deformation tensor to generate columns of material stiffness
sqeps = (sqrt(eps)); % this adjusts the perturbation size to modulate the accuracy

Enumerical = zeros(9,3);
for i = 1:3

  de = zeros(4,1);
  if (i < 3 )
    if strain(i) == 0
     e = sqeps;
    else
     e = sqeps*strain(i);
    end
    de(i) = 1;
  else
    if strain(i) == 0
     e =sqeps;
    else
     e = sqeps*strain(i);
    end
    de(i) = 1;0.5;
  end

  % Perturbed strain tensor
  strainp = strain+e*de;

  % Get perturbed stresses from the model
% if plasversion == 1
%                 
% sigp = J2RadialReturn0(strainp,ep_n,beta_n,a_n,mu,bulk,Khard,Hhard,sigy);
% 
% else

if plasmodel == 1
[DGAMA,LALGVA,RSTAVA,STRES] = SUVM(2,2,[ElemYM,Elemv],[0,100;sigy,Khard*100+sigy],[ee_n; a_n],strainp);
Dp = CTVM(DGAMA,LALGVA(1),2,2,[ElemYM,Elemv],[0,100;sigy,Khard*100+sigy],RSTAVA,STRES);
elseif plasmodel == 2
[DGAM,LALGVA,RSTAVA,STRES] = SUTR(2,2,[ElemYM,Elemv],[0,100;sigy,Khard*100+sigy],[ee_n; a_n],strainp);
Dp = CTTR(LALGVA(1),2,LALGVA,2,[ElemYM,Elemv],[0,100;sigy,Khard*100+sigy],RSTAVA,eps_n1,STRES);
end

% end
%   [sigp, dp] = model(Fp);

  Dp = Dp(1:3,1:3);

  % Finite difference to get columns of 6x6 moduli matrix
  Enumerical((i-1)*3+1:3*i,:) = (1/e)*(Dp - D);

end
% 
% fprintf('Reference Frame Moduli from Model'); 
% Eref
fprintf('Reference Frame Moduli from Perturbation');
Enumerical
