%
% Driver to compute the elasto-plastic material moduli by finite differencing
%
% Based off of Sanjay's website on FEAP
% Only works currently for isotropic hardening; but Miehe says nothing
% special about history-dependent materials. so something must be wrong.
% Maybe my constitutive tangent has an error.
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
                
[sig,D] = J2RadialReturn0(strain,ep_n,beta_n,a_n,mu,bulk,Khard,Hhard,sigy);
% [sig,d] = model(F);

% Store stress in vector form
Sv = sig;


% Perturb right Cauchy-Green deformation tensor to generate columns of material stiffness
sqeps = (sqrt(eps));

Dnumerical = zeros(6,6);
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
  sigp = J2RadialReturn0(strainp,ep_n,beta_n,a_n,mu,bulk,Khard,Hhard,sigy);
%   [sigp, dp] = model(Fp);

  % Compute perturbed 2nd PK
  Spv = sigp;

  % Finite difference to get columns of 6x6 moduli matrix
  Dnumerical(:,i) = (1/e)*(Spv - Sv);
end
% Symmetrize
Dnumerical = 1/2*(Dnumerical + Dnumerical');

fprintf('Reference Frame Moduli from Model'); 
D
fprintf('Reference Frame Moduli from Perturbation');
Dnumerical
