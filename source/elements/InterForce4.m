% Tim Truster
% 07/14/2015
%
% Evaluate interactive force, from Hari and Rick's equation, using gradient
% of F

IntFor = zeros(ndf,1); % approximate value, for plotting


%% Branch on interactive force type
if InFoType == 2 % Drag force
%     InterForce3Both
    InterForce3Drag
elseif InFoType == 3 % Cohesive force
    InterForceCohesive
else % Incompatibility force
    
% Set densities
            
if iel1 == 33
    rho1 = mateprop(29)/JxX1;
elseif iel1 == 26 
    var = hr(nhb-1+ll);
    rhoM = (1 + var)*mateprop1(3);
    rho1 = rhoM/JxX1;
elseif iel1 == 32
%     var = hr(nhb-1+ll);
%     rhoM = (1 + var)^2*rhoS; % See eq (79) in my paper
%     rho1 = rhoM/JxX1;
end

if iel2 == 33
    rho2 = mateprop(30)/JxX2;
elseif iel2 == 26 
    var = hr(nhb-1+ll+8);
    rhoM = (1 + var)*mateprop2(3);
    rho2 = rhoM/JxX2;
elseif iel2 == 32
    rho2 = rho1;
%     var = hr(nhb-1+ll);
%     rhoM = (1 + var)^2*rhoS; % See eq (79) in my paper
%     rho2 = rhoM/JxX2;
end

% This is approximate
rhoMix = rho1 + rho2;

IntFor = (1-rho1)/rhoMix*fi1*dF1*Pe_vec1 - (1-rho2)/rhoMix*fi2*dF2*Pe_vec2;


%% Easy terms: 1-1 and 2-2

% helper matrix for variation of dF such that 
% dF1*Pe_vec1 = PN1*reshape(ul1,24,1) with dF1 evaluated using ul1
PN1 = Pe1*dN1;
% for i = 1:nel1
%     PN1(1:3,3*i-2:3*i) = PN1(1:3,3*i-2:3*i)';
% end
dN1b = reshape(PN1,9,8); % This is a little faster
dN1b = dN1b([1 4 7 2 5 8 3 6 9],:);
PN1 = reshape(dN1b,3,24);
ElemF(1:nel1*ndf) = ElemF(1:nel1*ndf) - c1*(Nmat1'*(1-rho1)/rhoMix*fi1*dF1*Pe_vec1); %(Pe_vec1'*dF1'*fi1')'

ElemK(1:nel1*ndf,1:nel1*ndf) = ElemK(1:nel1*ndf,1:nel1*ndf) + c1*(Nmat1'*(1-rho1)/rhoMix*fi1*dF1*Aeg_mat1*Bmat1 ...
             + Nmat1'*(1-rho1)/rhoMix*fi1*PN1);

% helper matrix for variation of dF such that 
% dF2*Pe_vec2 = PN2*reshape(ul2,24,1) with dF2 evaluated using ul2
PN2 = Pe2*dN2;
% for i = 1:nel2
%     PN2(1:3,3*i-2:3*i) = PN2(1:3,3*i-2:3*i)';
% end
dN1b = reshape(PN2,9,8);
dN1b = dN1b([1 4 7 2 5 8 3 6 9],:);
PN2 = reshape(dN1b,3,24);
ElemF(nel1*ndf+1:nel*ndf) = ElemF(nel1*ndf+1:nel*ndf) - c2*(Nmat2'*(1-rho2)/rhoMix*fi2*dF2*Pe_vec2);

ElemK(nel1*ndf+1:nel*ndf,nel1*ndf+1:nel*ndf) = ElemK(nel1*ndf+1:nel*ndf,nel1*ndf+1:nel*ndf) + c2*(Nmat2'*(1-rho2)/rhoMix*fi2*dF2*Aeg_mat2*Bmat2 ...
             + Nmat2'*(1-rho2)/rhoMix*fi2*PN2);


%% Projection terms
% 1-2
    
% Current physical location of int pt
xint = (xl2(1,1:nel2)+ul2(1,1:nel2))*shl;
yint = (xl2(2,1:nel2)+ul2(2,1:nel2))*shl;
zint = (xl2(3,1:nel2)+ul2(3,1:nel2))*shl;

% Find current point's location in deformed config of other
% constituent
xi = POU_Coord3(xint,yint,zint,xl1(:,1:nel1)+ul1(1:3,1:nel1),1,nel1);

% Evaluate  basis functions at integration point
shl1 = shlb(xi,nel1,nel1,0,0);
Nmat1p = zeros(3,nel*ndf/2);
for mm = 1:nel1
  Nmat1p(:,3*mm-2:3*mm) = [shl1(mm) 0        0         
                           0        shl1(mm) 0         
                           0        0        shl1(mm)];
end

ElemF(1:nel1*ndf) = ElemF(1:nel1*ndf) - c2*(-Nmat1p'*(1-rho2)/rhoMix*fi2*dF2*Pe_vec2); %(Pe_vec2'*dF2'*fi2')'

ElemK(1:nel1*ndf,nel1*ndf+1:nel*ndf) = ElemK(1:nel1*ndf,nel1*ndf+1:nel*ndf) - c2*(Nmat1p'*(1-rho2)/rhoMix*fi2*dF2*Aeg_mat2*Bmat2 ...
             + Nmat1p'*(1-rho2)/rhoMix*fi2*PN2);


% 2-1

% Current physical location of int pt
xint = (xl1(1,1:nel1)+ul1(1,1:nel1))*shl;
yint = (xl1(2,1:nel1)+ul1(2,1:nel1))*shl;
zint = (xl1(3,1:nel1)+ul1(3,1:nel1))*shl;

% Find current point's location in deformed config of other
% constituent
xi = POU_Coord3(xint,yint,zint,xl2(:,1:nel2)+ul2(1:3,1:nel2),1,nel2);

% Evaluate  basis functions at integration point
shl2 = shlb(xi,nel2,nel2,0,0);
Nmat2p = zeros(3,nel*ndf/2);
for mm = 1:nel1
  Nmat2p(:,3*mm-2:3*mm) = [shl2(mm) 0        0         
                           0        shl2(mm) 0         
                           0        0        shl2(mm)];
end

ElemF(nel1*ndf+1:nel*ndf) = ElemF(nel1*ndf+1:nel*ndf) - c1*(-Nmat2p'*(1-rho1)/rhoMix*fi1*dF1*Pe_vec1);

ElemK(nel1*ndf+1:nel*ndf,1:nel1*ndf) = ElemK(nel1*ndf+1:nel*ndf,1:nel1*ndf) - c1*(Nmat2'*(1-rho1)/rhoMix*fi1*dF1*Aeg_mat1*Bmat1 ...
             + Nmat2p'*(1-rho1)/rhoMix*fi1*PN1);

end
