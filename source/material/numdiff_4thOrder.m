%
% Driver to compute the reference and spatial material moduli by finite differencing
%
% From Sanjay's website on FEAP
clear all;
format short e;

% Input a deformation gradient
F = [ .4 1 0; 
      0 2 2;
      0 0 1];
J = det(F);
mu = 11;
K = 999;
E = mu*(3*K+2*mu)/(K+mu);
v = K/(2*(K+mu));
mateprop = [1 3 0 E v]; % all other combos of material flags work as well.
lam = getlam(mateprop);

% Call model to get model cauchy stress and spatial moduli
[sig,d] = SigmaCmat3i(F,J,mateprop,lam);
d = d(1:6,1:6);
% [sig,d] = model(F);

% Convert to matrix for for convenience
sigm = [sig(1) sig(4) sig(6);
        sig(4) sig(2) sig(5);
        sig(6) sig(5) sig(3)];

% Compute some useful quantitites
J = det(F);
C = F'*F;
[vec,val] = eig(C);
U = zeros(3,3);
for i = 1:3
 U = U + sqrt(val(i,i))*vec(:,i)*vec(:,i)';
end
R = F*inv(U);
Fi = inv(F);

% Compute 2nd PK stress
Sm = J*Fi*sigm*Fi';

% Store 2nd PK in vector form
Sv = [Sm(1,1);
      Sm(2,2);
      Sm(3,3);
      Sm(1,2);
      Sm(2,3);
      Sm(3,1)];

% Convenience arrays for mapping between Voigt and tensor notation
% Warning: my ordering is [11 22 33 12 23 31] !
io = [1 2 3 1 2 3];
it = [1 2 3 2 3 1];

% Compute reference moduli (by pull back and multiplication by J)
for i = 1:6
   ioi = io(i); iti = it(i);
   Mf = [ Fi(ioi, 1 )*Fi(iti, 1);
          Fi(ioi, 2 )*Fi(iti, 2);
          Fi(ioi, 3 )*Fi(iti, 3);
          Fi(ioi, 1 )*Fi(iti, 2) +  Fi(ioi, 2 )*Fi(iti, 1);
          Fi(ioi, 2 )*Fi(iti, 3) +  Fi(ioi, 3 )*Fi(iti, 2);
          Fi(ioi, 3 )*Fi(iti, 1) +  Fi(ioi, 1 )*Fi(iti, 3)];
          
          df = Mf'*d*J;

 for j = i:6

   ioj = io(j); itj = it(j);
   Mb = [ Fi(ioj, 1 )*Fi(itj, 1);
          Fi(ioj, 2 )*Fi(itj, 2);
          Fi(ioj, 3 )*Fi(itj, 3);
          Fi(ioj, 1 )*Fi(itj, 2) +  Fi(ioj, 2 )*Fi(itj, 1);
          Fi(ioj, 2 )*Fi(itj, 3) +  Fi(ioj, 3 )*Fi(itj, 2);
          Fi(ioj, 3 )*Fi(itj, 1) +  Fi(ioj, 1 )*Fi(itj, 3)];

   D(i,j) = df*Mb;

 end
end

% Fill by symmetry
for i = 1:6
 for j = (i+1):6
  D(j,i) = D(i,j);
 end
end


% Perturb right Cauchy-Green deformation tensor to generate columns of material stiffness
sqeps = (sqrt(eps));

for i = 1:6

  dC = zeros(3,3);
  if (i < 4 )
    if C(i,i) == 0
     e = sqeps;
    else
     e = sqeps*C(i,i);
    end
    dC(i,i) = 1;
  else
    if C(io(i),it(i)) == 0
     e =sqeps;
    else
     e = sqeps*C(io(i),it(i));
    end
    dC(io(i),it(i)) = 0.5;
    dC(it(i),io(i)) = 0.5;
  end

  % Perturbed right Cauchy-Green Deformation tensor
  Cp = C+e*dC;

  % Related perturbed kinematic measures
  [vecp,valp] = eig(Cp);
  Up = zeros(3,3);
  for k = 1:3
   Up = Up + sqrt(valp(k,k))*vecp(:,k)*vecp(:,k)';
  end
  Fp = R*Up;
  Jp = det(Fp);
  Fpi = inv(Fp);

  % Get perturbed stresses from the model
  sigp = SigmaCmat3i(Fp,Jp,mateprop,lam);
%   [sigp, dp] = model(Fp);
  sigpm = [sigp(1) sigp(4) sigp(6);
           sigp(4) sigp(2) sigp(5);
           sigp(6) sigp(5) sigp(3)];

  % Compute perturbed 2nd PK
  Spm = (Jp)*Fpi*sigpm*Fpi';
  Spv = [Spm(1,1);
         Spm(2,2);
         Spm(3,3);
         Spm(1,2);
         Spm(2,3);
         Spm(3,1)];

  % Finite difference to get columns of 6x6 moduli matrix
  Dnumerical(:,i) = (2/e)*(Spv - Sv);
end

fprintf('Reference Frame Moduli from Model'); 
D
fprintf('Reference Frame Moduli from Perturbation');
Dnumerical

% Report in pushed state also
for i = 1:6
   ioi = io(i); iti = it(i);
   Mf = [ F(ioi, 1 )*F(iti, 1);
          F(ioi, 2 )*F(iti, 2);
          F(ioi, 3 )*F(iti, 3);
          F(ioi, 1 )*F(iti, 2) +  F(ioi, 2 )*F(iti, 1);
          F(ioi, 2 )*F(iti, 3) +  F(ioi, 3 )*F(iti, 2);
          F(ioi, 3 )*F(iti, 1) +  F(ioi, 1 )*F(iti, 3)];

          df = Mf'*Dnumerical/J;
 for j = i:6

   ioj = io(j); itj = it(j);
   Mb = [ F(ioj, 1 )*F(itj, 1);
          F(ioj, 2 )*F(itj, 2);
          F(ioj, 3 )*F(itj, 3);
          F(ioj, 1 )*F(itj, 2) +  F(ioj, 2 )*F(itj, 1);
          F(ioj, 2 )*F(itj, 3) +  F(ioj, 3 )*F(itj, 2);
          F(ioj, 3 )*F(itj, 1) +  F(ioj, 1 )*F(itj, 3)];

   dnumerical(i,j) = df*Mb;

 end
end

% Fill by symmetry
for i = 1:6
 for j=(i+1):6
  dnumerical(j,i)=dnumerical(i,j);
 end
end


fprintf('Spatial Frame Moduli from Model'); 
d
fprintf('Spatial Frame Moduli from Perturbation'); 
dnumerical
