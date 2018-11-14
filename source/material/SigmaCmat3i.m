function [sigma, cmat] = SigmaCmat3i(F,J,mateprop,lam)
% modified by Pinlei Chen for interface part
% Constitutive material model evaluation for NL elasticity, 3D
% 01/04/2012 05/20/2013
% Output: sigma = sigma = cauchy stress
%         cmat = F_iI*F_jJ*F_kK*F_lL*C_IJKL./J

one = [1; 1; 1; 0; 0; 0; 0; 0; 0];
mat1 = one*one';
matE = diag([2,2,2,1,1,1,0,0,0]);

[sigmai, cmati] = SigmaCmatNSCST3i(F,J,mateprop);
[theta1,theta2] = ThetaNS(J,mateprop);
sigmap = lam*theta1*J*one;
cpmat = lam*((theta2*J^2 + theta1*J)*mat1 - theta1*J*matE);
sigma = (sigmai + sigmap)./J;
cmat = (cmati + cpmat)./J;