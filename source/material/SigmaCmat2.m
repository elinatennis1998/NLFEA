function [sigma, cmat] = SigmaCmat2(F,J,mateprop,lam)
%
% Constitutive material model evaluation for NL elasticity
% 10/29/2011
% Output: sigma = J*sigma = tau = Kirchhoff stress
%         cmat = F_iI*F_jJ*F_kK*F_lL*C_IJKL

one = [1; 1; 0; 0];
mat1 = one*one';
matE = diag([2,2,1,0]);

if mateprop(1) > 0
    
    [sigmai, cmati] = SigmaCmatNSCST2i(F,J,mateprop);
    [theta1,theta2] = ThetaNS(J,mateprop);
    sigmap = lam*theta1*J*one;
    cpmat = lam*((theta2*J^2 + theta1*J)*mat1 - theta1*J*matE);
    sigma = sigmai + sigmap;
    cmat = cmati + cpmat;
    
else % Saint Venant-Kirchhoff Material
    
    PatchE = mateprop(4);
    Patchv = mateprop(5);
    mu = PatchE/(2*(1+Patchv));%80.19;
            
    C = F'*F;
    Ev = 1/2*[C(1,1)-1; C(2,2)-1; 2*C(1,2); 0];
    trE = one'*Ev;
    S = lam*trE*one + mu*matE*Ev;
    D = lam*(one*one') + mu*matE;
    P = [F(1,1)^2 F(1,2)^2 2*F(1,1)*F(1,2) 0
         F(2,1)^2 F(2,2)^2 2*F(2,1)*F(2,2) 0
         F(1,1)*F(2,1) F(1,2)*F(2,2) F(1,1)*F(2,2)+F(1,2)*F(2,1) 0
         zeros(1,4)];
    sigma = P*S;
    cmat = P*D*P';
    
end