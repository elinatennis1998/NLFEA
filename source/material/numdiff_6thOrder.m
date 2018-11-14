%
% Driver to compute the reference and spatial material moduli by finite differencing
%
% From Sanjay's website on FEAP
% Modified 06/16/2014 to compute the numerical 6th order tensor
% It works!!!!!!!!!!!
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
[dmat]=dmat3(J,mateprop,lam);
dmat = dmat/J;
[sig,cmat] = SigmaCmat3i(F,J,mateprop,lam);
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
D = zeros(6,6);
for i = 1:6
   ioi = io(i); iti = it(i);
   Mf = [ Fi(ioi, 1 )*Fi(iti, 1);
          Fi(ioi, 2 )*Fi(iti, 2);
          Fi(ioi, 3 )*Fi(iti, 3);
          Fi(ioi, 1 )*Fi(iti, 2) +  Fi(ioi, 2 )*Fi(iti, 1);
          Fi(ioi, 2 )*Fi(iti, 3) +  Fi(ioi, 3 )*Fi(iti, 2);
          Fi(ioi, 3 )*Fi(iti, 1) +  Fi(ioi, 1 )*Fi(iti, 3)];
          
          df = Mf'*cmat(1:6,1:6)*J;

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

% pull back 6th order tensor
Eref = zeros(36,6);
Eref2 = zeros(36,6);
t = tranr4(Fi,Fi);
% get two of the contractions
for i = 1:6
    Eref2((i-1)*6+1:6*i,:) = J*t'*dmat((i-1)*6+1:6*i,1:6)*t;
end
% do the third contraction
for i = 1:6
    Erows = Eref2(i:6:30+i,:);
    Eref(i:6:30+i,:) = t'*Erows;
end


% Perturb right Cauchy-Green deformation tensor to generate columns of material stiffness
sqeps = (sqrt(eps)); % this adjusts the perturbation size to modulate the accuracy

Enumerical = zeros(36,6);
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

  % Get perturbed moduli from the model
%   [sigp, dp] = model(Fp);
  [sigma2R, cmatp] = SigmaCmat3i(Fp,Jp,mateprop,lam);
%   sigpm = [sigp(1) sigp(4) sigp(6);
%            sigp(4) sigp(2) sigp(5);
%            sigp(6) sigp(5) sigp(3)];

  % Compute perturbed 2nd PK
%   Spm = (Jp)*Fpi*sigpm*Fpi';
%   Spv = [Spm(1,1);
%          Spm(2,2);
%          Spm(3,3);
%          Spm(1,2);
%          Spm(2,3);
%          Spm(3,1)];
    % Compute perturbed reference moduli (by pull back and multiplication by J)
    Dp = zeros(6,6);
    for k = 1:6
       ioi = io(k); iti = it(k);
       Mf = [ Fpi(ioi, 1 )*Fpi(iti, 1);
              Fpi(ioi, 2 )*Fpi(iti, 2);
              Fpi(ioi, 3 )*Fpi(iti, 3);
              Fpi(ioi, 1 )*Fpi(iti, 2) +  Fpi(ioi, 2 )*Fpi(iti, 1);
              Fpi(ioi, 2 )*Fpi(iti, 3) +  Fpi(ioi, 3 )*Fpi(iti, 2);
              Fpi(ioi, 3 )*Fpi(iti, 1) +  Fpi(ioi, 1 )*Fpi(iti, 3)];

              df = Mf'*cmatp(1:6,1:6)*Jp;

     for j = k:6

       ioj = io(j); itj = it(j);
       Mb = [ Fpi(ioj, 1 )*Fpi(itj, 1);
              Fpi(ioj, 2 )*Fpi(itj, 2);
              Fpi(ioj, 3 )*Fpi(itj, 3);
              Fpi(ioj, 1 )*Fpi(itj, 2) +  Fpi(ioj, 2 )*Fpi(itj, 1);
              Fpi(ioj, 2 )*Fpi(itj, 3) +  Fpi(ioj, 3 )*Fpi(itj, 2);
              Fpi(ioj, 3 )*Fpi(itj, 1) +  Fpi(ioj, 1 )*Fpi(itj, 3)];

       Dp(k,j) = df*Mb;

     end
    end

    % Fill by symmetry
    for k = 1:6
     for j = (k+1):6
      Dp(j,k) = Dp(k,j);
     end
    end
    t = tranr4(Fpi,Fpi);
    Dp2 = Jp*t*cmatp(1:6,1:6)*t';

  % Finite difference to get columns of 6x6 moduli matrix
  Enumerical((i-1)*6+1:6*i,:) = (2/e)*(Dp - D);
end

fprintf('Reference Frame Moduli from Model'); 
Eref
fprintf('Reference Frame Moduli from Perturbation');
Enumerical

% Report in pushed state also
enumerical = zeros(36,6);
enumerical2 = zeros(36,6);
t = tranr4(F,F);
% get two of the contractions
for i = 1:6
    enumerical2((i-1)*6+1:6*i,:) = t'*Enumerical((i-1)*6+1:6*i,1:6)*t/J;
end
% do the third contraction
for i = 1:6
    Erows = enumerical2(i:6:30+i,:);
    enumerical(i:6:30+i,:) = t'*Erows;
end


fprintf('Spatial Frame Moduli from Model'); 
dmat
fprintf('Spatial Frame Moduli from Perturbation'); 
enumerical
