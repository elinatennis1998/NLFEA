% function [tau_n1,c_n1,ee_n1,beta_n1,a_n1,s_n1,Cdev_n1] = logplast(ee_n,df,beta_n,a_n,mu,bulk,K,H,sigy)
% Tim Truster
% 05/14/12
% Kinematics for logarithmic large strain plasticity
% Based on implementation presented in deSouza plasticity book:
% Computational Methods for Plasticity - Theory and Applications

mati = mateprop(1);
matp = mateprop(3);

% Convert from vector to matrix
ee_m    = [ee_n(1) ee_n(4)/2.d0 ee_n(6)/2.d0
           ee_n(4)/2.d0 ee_n(2) ee_n(5)/2.d0
           ee_n(6)/2.d0 ee_n(5)/2.d0 ee_n(3)];

% Compute spectral decomposition
% [Vn,ee_v] = eig(ee_m);
[Vn,ee_v] = eig3(ee_m);
ee_v = diag(ee_v);
% Vn = eye(3);
% ee_v = ee_m;

% Compute Be_n = exp(2*ee_n)
be_v = diag(exp(2*diag(ee_v)));
be_n = Vn*be_v*Vn';

% Compute trial elastic strain betr_n1 = df*be_n*df'
betr_n1 = df*be_n*df';

% Compute trial logarithmic strain eetr_n1 = 1/2*ln(betr_n1)

% Compute spectral decomposition
% [Vn1,betr_v] = eig(betr_n1);
[Vn1,betr_v] = eig3(betr_n1);
betr_v = diag(betr_v);
% Vn1 = eye(3);
% betr_v = betr_n1;

% Compute Be_n = exp(2*ee_n)
eetr_v2 = diag(log(diag(betr_v)));
eetr_v = 1/2*(eetr_v2);
deetr_v = diag(1./(diag(betr_v)));
eetr_m = Vn1*eetr_v*Vn1';
eetr_n1 = [eetr_m(1,1); eetr_m(2,2); eetr_m(3,3); 2*eetr_m(1,2); 2*eetr_m(2,3); 2*eetr_m(3,1)];

if mati == 8 % Hencky model
    
PatchE = mateprop(4);
Patchv = mateprop(5);
mu = PatchE/(2*(1+Patchv));

switch matp
    
    case 1 % J2-viscoplasticity with linear isotropic-kinematic hardening

    sigy = mateprop(8);
    K = mateprop(9);
    H = mateprop(10);
    eta = mateprop(11);
    [s_n1,c_n1,ee_n1,beta_n1,a_n1] = J2RadialReturn(eetr_n1,beta_n,a_n,mu,K,H,sigy,eta,dt);

end

% ee_n1 = diag([1 1 1 2 2 2])*ee_n1;

% Compute extra 4-order tensors L and B
L = getL(Vn1,betr_v,eetr_v2,deetr_v);

vec = [1 2 3 4 5 6 4 5 6];
C_n12 = c_n1(vec,vec);
L2 = L(vec,vec);

x = betr_n1;
% B = [   x(1,1),     0,     0,     0,     0,   x(1,3),   x(1,2),     0,     0
%           0,   x(2,2),     0,   x(1,2),     0,     0,     0,   x(2,3),     0
%           0,     0,   x(3,3),     0,   x(2,3),     0,     0,     0,   x(1,3)
%       x(1,2)/2, x(1,2)/2,     0, x(1,1)/2,     0, x(2,3)/2, x(2,2)/2, x(1,3)/2,     0
%           0, x(2,3)/2, x(2,3)/2, x(1,3)/2, x(2,2)/2,     0,     0, x(3,3)/2, x(1,2)/2
%       x(1,3)/2,     0, x(1,3)/2,     0, x(1,2)/2, x(3,3)/2, x(2,3)/2,     0, x(1,1)/2
%       x(1,2)/2, x(1,2)/2,     0, x(1,1)/2,     0, x(2,3)/2, x(2,2)/2, x(1,3)/2,     0
%           0, x(2,3)/2, x(2,3)/2, x(1,3)/2, x(2,2)/2,     0,     0, x(3,3)/2, x(1,2)/2
%       x(1,3)/2,     0, x(1,3)/2,     0, x(1,2)/2, x(3,3)/2, x(2,3)/2,     0, x(1,1)/2];
B = [ 2*x(1,1),     0,     0,           x(1,2),             0,           x(1,3),           x(1,2),             0,           x(1,3)
          0, 2*x(2,2),     0,           x(1,2),           x(2,3),             0,           x(1,2),           x(2,3),             0
          0,     0, 2*x(3,3),             0,           x(2,3),           x(1,3),             0,           x(2,3),           x(1,3)
        x(1,2),   x(1,2),     0, x(1,1)/2 + x(2,2)/2,         x(1,3)/2,         x(2,3)/2, x(1,1)/2 + x(2,2)/2,         x(1,3)/2,         x(2,3)/2
          0,   x(2,3),   x(2,3),         x(1,3)/2, x(2,2)/2 + x(3,3)/2,         x(1,2)/2,         x(1,3)/2, x(2,2)/2 + x(3,3)/2,         x(1,2)/2
        x(1,3),     0,   x(1,3),         x(2,3)/2,         x(1,2)/2, x(1,1)/2 + x(3,3)/2,         x(2,3)/2,         x(1,2)/2, x(1,1)/2 + x(3,3)/2
        x(1,2),   x(1,2),     0, x(1,1)/2 + x(2,2)/2,         x(1,3)/2,         x(2,3)/2, x(1,1)/2 + x(2,2)/2,         x(1,3)/2,         x(2,3)/2
          0,   x(2,3),   x(2,3),         x(1,3)/2, x(2,2)/2 + x(3,3)/2,         x(1,2)/2,         x(1,3)/2, x(2,2)/2 + x(3,3)/2,         x(1,2)/2
        x(1,3),     0,   x(1,3),         x(2,3)/2,         x(1,2)/2, x(1,1)/2 + x(3,3)/2,         x(2,3)/2,         x(1,2)/2, x(1,1)/2 + x(3,3)/2];

c_n1 = 1/2*C_n12*L2*B;
sigma2 = s_n1;
            Smat = ...
[    sigma2(1),        0,        0,           sigma2(4)/2,                 0,           sigma2(6)/2,           sigma2(4)/2,                 0,          -sigma2(6)/2
         0,    sigma2(2),        0,           sigma2(4)/2,           sigma2(5)/2,                 0,          -sigma2(4)/2,           sigma2(5)/2,                 0
         0,        0,    sigma2(3),                 0,           sigma2(5)/2,           sigma2(6)/2,                 0,          -sigma2(5)/2,           sigma2(6)/2
   sigma2(4)/2,  sigma2(4)/2,        0, sigma2(1)/4 + sigma2(2)/4,           sigma2(6)/4,           sigma2(5)/4, sigma2(2)/4 - sigma2(1)/4,           sigma2(6)/4,          -sigma2(5)/4
         0,  sigma2(5)/2,  sigma2(5)/2,           sigma2(6)/4, sigma2(2)/4 + sigma2(3)/4,           sigma2(4)/4,          -sigma2(6)/4, sigma2(3)/4 - sigma2(2)/4,           sigma2(4)/4
   sigma2(6)/2,        0,  sigma2(6)/2,           sigma2(5)/4,           sigma2(4)/4, sigma2(1)/4 + sigma2(3)/4,           sigma2(5)/4,          -sigma2(4)/4, sigma2(1)/4 - sigma2(3)/4
   sigma2(4)/2, -sigma2(4)/2,        0, sigma2(2)/4 - sigma2(1)/4,          -sigma2(6)/4,           sigma2(5)/4, sigma2(1)/4 + sigma2(2)/4,          -sigma2(6)/4,          -sigma2(5)/4
         0,  sigma2(5)/2, -sigma2(5)/2,           sigma2(6)/4, sigma2(3)/4 - sigma2(2)/4,          -sigma2(4)/4,          -sigma2(6)/4, sigma2(2)/4 + sigma2(3)/4,          -sigma2(4)/4
  -sigma2(6)/2,        0,  sigma2(6)/2,          -sigma2(5)/4,           sigma2(4)/4, sigma2(1)/4 - sigma2(3)/4,          -sigma2(5)/4,          -sigma2(4)/4, sigma2(1)/4 + sigma2(3)/4];
c_n1 = c_n1 - 2*Smat;
c_n1 = c_n1(1:6,1:6);
s_n1 = [s_n1; zeros(3,1)];
cdev_n1 = [c_n1 zeros(6,3); zeros(3,9)];

end
