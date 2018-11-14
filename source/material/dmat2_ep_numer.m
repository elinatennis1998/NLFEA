function Enumerical = dmat2_ep_numer(du,eep,beta_n,a_n,ElemYM,Elemv,Khard,Hhard,sigy,C_n1,plasversion,plasmodel)
%
% Tim Truster
% 07/29/2014
% Routine for computing numerical version of 6th order constitutive tensor
% Material routines copied from NL_Elem52_2d.m, isw=3.
% C_n1 is the algorithmic moduli for this point, pre-computed outside this
% routine.

ind3to2 = [1 2 4];

if plasversion == 1
                    
    bulk = ElemYM/(3*(1-2*Elemv));%Elemv*ElemYM/((1+Elemv)*(1-2*Elemv));
    mu = ElemYM/(2*(1+Elemv));

    eps2d = du(1:3); %2D enhanced strain
    eps3d = [eps2d(1); eps2d(2); 0; eps2d(3); 0; 0];
    ep_n = eep;
    strain = eps3d;

    D = C_n1;
    D = D(ind3to2,ind3to2);
                
else
    
    bulk = ElemYM/(3*(1-2*Elemv));%Elemv*ElemYM/((1+Elemv)*(1-2*Elemv));
    
    deps2d = du(1:3);
    deps_n1 = [deps2d(1); deps2d(2); deps2d(3); 0];
    ee_n = eep;
    strain = ee_n + deps_n1;

    D = C_n1;
    D = D(1:3,1:3);
    
end

% Perturb right Cauchy-Green deformation tensor to generate columns of material stiffness
sqeps = (sqrt(eps)); % this adjusts the perturbation size to modulate the accuracy
Enumerical = zeros(9,3);


%% evaluate perturbed moduli
for i = 1:3

    if plasversion == 1

      i2 = ind3to2(i);
      de = zeros(6,1);
      if (i < 3 )
        if strain(i2) == 0
         e = sqeps;
        else
         e = sqeps*strain(i2);
        end
        de(i2) = 1;
      else
        if strain(i2) == 0
         e =sqeps;
        else
         e = sqeps*strain(i2);
        end
        de(i2) = 1;
      end

      % Perturbed strain tensor
      strainp = strain+e*de;

        [~,~,~,~,~,~,Dp] = J2RadialReturn0(strainp,ep_n,beta_n,a_n,mu,bulk,Khard,Hhard,sigy);

        Dp = Dp(ind3to2,ind3to2);

    else

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
        de(i) = 1;
        end

        % Perturbed strain tensor
        strainp = strain+e*de;

        if plasmodel == 1
        [DGAMA,LALGVA,RSTAVA,STRES] = SUVM(2,2,[ElemYM,Elemv],[0,100;sigy,Khard*100+sigy],[ee_n; a_n],strainp);
        Dp = CTVM(DGAMA,LALGVA(1),2,2,[ElemYM,Elemv],[0,100;sigy,Khard*100+sigy],RSTAVA,STRES);
        elseif plasmodel == 2
        [DGAM,LALGVA,RSTAVA,STRES] = SUTR(2,2,[ElemYM,Elemv],[0,100;sigy,Khard*100+sigy],[ee_n; a_n],strainp);
        Dp = CTTR(LALGVA(1),2,LALGVA,2,[ElemYM,Elemv],[0,100;sigy,Khard*100+sigy],RSTAVA,strainp,STRES);
        end

        Dp = Dp(1:3,1:3) - bulk*[1 1 0]'*[1 1 0];

    end

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