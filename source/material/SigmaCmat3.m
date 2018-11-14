function [sigma, cmat] = SigmaCmat3(F,J,mateprop,lam)
%
% Constitutive material model evaluation for NL elasticity, 3D
% 01/04/2012
% Output: sigma = J*sigma = tau = Kirchhoff stress
%         cmat = F_iI*F_jJ*F_kK*F_lL*C_IJKL

one = [1; 1; 1; 0; 0; 0; 0; 0; 0];
mat1 = one*one';
matE = diag([2,2,2,1,1,1,0,0,0]);

[sigmai, cmati] = SigmaCmatNSCST3i(F,J,mateprop);
[theta1,theta2] = ThetaNS(J,mateprop);
sigmap = lam*theta1*J*one;
cpmat = lam*((theta2*J^2 + theta1*J)*mat1 - theta1*J*matE);
sigma = sigmai + sigmap;
cmat = cmati + cpmat;