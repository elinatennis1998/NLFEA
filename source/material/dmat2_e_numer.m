function Enumerical = dmat2_e_numer(E,Patcha,Patchb,Patchc,cmat)
%
% Tim Truster
% 08/17/2014
% Routine for computing numerical version of 6th order constitutive tensor
% Test for small strain elasticity, goes with NL_Elem11_2dDG.m.
% Gives results similar to the analytical one.

ind3to2 = [1 2 4];
strain = E;

D = cmat;

% Perturb right Cauchy-Green deformation tensor to generate columns of material stiffness
sqeps = (sqrt(eps)); % this adjusts the perturbation size to modulate the accuracy
Enumerical = zeros(9,3);


%% evaluate perturbed moduli
for i = 1:3

        de = zeros(3,1);
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
        de(i) = 1;
        end

      % Perturbed strain tensor
      Ep = strain+e*de;
    Ep(3) = Ep(3)/2; %for extra factor of 2 on shear strain
    Ekk = Ep(1) + Ep(2);
    E2 = [Ep(1) Ep(3); Ep(3) Ep(2)];
    E2 = E2*E2;
    E2 = [E2(1,1); E2(2,2); E2(1,2)];

    Dp = [Patcha+Patchb+2*Patchc*Ep(1) Patcha Patchc*Ep(3) 
            Patcha Patcha+Patchb+2*Patchc*Ep(2) Patchc*Ep(3)
            Patchc*Ep(3) Patchc*Ep(3)  Patchb/2+Patchc/2*Ekk];


  % Finite difference to get columns of 6x6 moduli matrix
  Enumerical((i-1)*3+1:3*i,:) = (1/e)*(Dp - D);

end

% Symmetrize E
ind1 = [2 4 1]; ind2 = [1 1 2]; %xx-xx-yy
ind3 = sub2ind([9 3], ind1, ind2);
Enumerical(ind3) = sum(diag(Enumerical(ind1,ind2)))/3;
ind1 = [5 2 4]; ind2 = [1 2 2]; %xx-yy-yy
ind3 = sub2ind([9 3], ind1, ind2);
Enumerical(ind3) = sum(diag(Enumerical(ind1,ind2)))/3;
ind1 = [3 7 1]; ind2 = [1 1 3]; %xx-xx-xy
ind3 = sub2ind([9 3], ind1, ind2);
Enumerical(ind3) = sum(diag(Enumerical(ind1,ind2)))/3;
ind1 = [9 3 7]; ind2 = [1 3 3]; %xx-xy-xy
ind3 = sub2ind([9 3], ind1, ind2);
Enumerical(ind3) = sum(diag(Enumerical(ind1,ind2)))/3;
ind1 = [6 8 5]; ind2 = [2 2 3]; %yy-yy-xy
ind3 = sub2ind([9 3], ind1, ind2);
Enumerical(ind3) = sum(diag(Enumerical(ind1,ind2)))/3;
ind1 = [9 6 8]; ind2 = [2 3 3]; %yy-xy-xy
ind3 = sub2ind([9 3], ind1, ind2);
Enumerical(ind3) = sum(diag(Enumerical(ind1,ind2)))/3;
ind1 = [6 8 3 7 2 4]; ind2 = [1 1 2 2 3 3];
ind3 = sub2ind([9 3], ind1, ind2); %xx-yy-xy
Enumerical(ind3) = sum(diag(Enumerical(ind1,ind2)))/6;