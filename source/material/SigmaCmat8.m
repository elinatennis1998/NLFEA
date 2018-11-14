function [S, D] = SigmaCmat8(F,mu,lam)
%
% Constitutive material model evaluation for NL elasticity
% 10/29/2011
% Output: sigma = J*sigma = tau = Kirchhoff stress
%         cmat = F_iI*F_jJ*F_kK*F_lL*C_IJKL

one = [1; 1; 0; 0];
mat1 = one*one';
matE = diag([2,2,1,0]);

% if mateprop(1) > 0
%     
%     [sigmai, cmati] = SigmaCmatNSCST2i(F,J,mateprop);
%     [theta1,theta2] = ThetaNS(J,mateprop);
%     sigmap = lam*theta1*J*one;
%     cpmat = lam*((theta2*J^2 + theta1*J)*mat1 - theta1*J*matE);
%     sigma = sigmai + sigmap;
%     cmat = cmati + cpmat;
%     
% else % Saint Venant-Kirchhoff Material
            
    C = F'*F;
    Ev = 1/2*[C(1,1)-1; C(2,2)-1; 2*C(1,2); 0];
    trE = one'*Ev;
    S = lam*trE*one + mu*matE*Ev;
    D = lam*mat1 + mu*matE;
    
% end