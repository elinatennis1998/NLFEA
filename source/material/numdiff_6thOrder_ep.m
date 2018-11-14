%
% Driver to compute the reference and spatial material moduli by finite differencing
%
% From Sanjay's website on FEAP
% Modified 06/16/2014 to compute the numerical 6th order tensor
% It works, I think.
clear all;
format short e;

% Sample total strain tensor
strain = [
   0.003481952464967
   0.002811361626954
                   0
  -0.002486099946961
                   0
                   0]*2;
mu = 28;
bulk = 46.666666666666667;
Khard = 2;
Hhard = 1;
sigy = 0.2;
ep_n = zeros(6,1);
beta_n = zeros(6,1);
% beta_n(2) = -.1;
a_n = 0.11;

% Call model to get model cauchy stress and spatial moduli
[sig,D] = J2RadialReturn0(strain,ep_n,beta_n,a_n,mu,bulk,Khard,Hhard,sigy);
% [sig,d] = model(F);

% pull back 6th order tensor
% Eref = computed value that I don't know


% Perturb right Cauchy-Green deformation tensor to generate columns of material stiffness
sqeps = (sqrt(eps)); % this adjusts the perturbation size to modulate the accuracy

Enumerical = zeros(36,6);
for i = 1:6

  de = zeros(6,1);
  if (i < 4 )
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
  [sigp,Dp] = J2RadialReturn0(strainp,ep_n,beta_n,a_n,mu,bulk,Khard,Hhard,sigy);
%   [sigp, dp] = model(Fp);

  % Finite difference to get columns of 6x6 moduli matrix
  Enumerical((i-1)*6+1:6*i,:) = (1/e)*(Dp - D);

end
% 
% fprintf('Reference Frame Moduli from Model'); 
% Eref
fprintf('Reference Frame Moduli from Perturbation');
Enumerical
