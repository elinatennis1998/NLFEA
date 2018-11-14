function [s_n1,Cdev_n1,ee_n1,beta_n1,a_n1] = J2RadialReturn(eps_n1,beta_n,a_n,mu,K,H,sigy,eta,dt)
% Tim Truster
% 04/17/2012
% Radial Return for linear isotropic/kinematic hardening with viscosity
%
% 3-D; arrangements for strains, etc is: [e11 e22 e33 e12 e23 e31]
% All tensors should be provided with all 6 components
%
% Inputs: eps_n1 = trial elastic strain at step n+1 (current value)
%         beta_n = back stress at step n
%         a_n    = equivalent plastic strain at step n
%         mu, K, H = material constants (shear, isotropic,
%                  kinematic hardening moduli)
%         sigy   = yield stress
%         eta    = viscosity
%
% Outputs: s_n1  = updated deviatoric stress at step n+1
%          Cdev_n1    = consistent algorithmic deviatoric tangent moduli
%          ee_n1   = updated elastic strain at step n+1
%          beta_n1 = updated back stress at step n+1
%          a_n1    = updated equivalent plastic strain at step n+1
%
% WARNING: cannot use both viscoplasticity and arc-length at same time

% e_n1 = zeros(6,1);
One = [1.0d0; 1.0d0; 1.0d0; 0.0d0; 0.0d0; 0.0d0];
I4 = diag([1.0d0 1.0d0 1.0d0 0.5d0 0.5d0 0.5d0]);
OneOne = (One*One');

% Constants
% R0 = 0.0d0;
% RP5 = 0.5d0;
R1 = 1.0d0;
R2 = 2.0d0;
R3 = 3.0d0;
sqr23 = sqrt(R2/R3);

% Set additional constants
Hb = K + H;
bet = K/Hb;
R2G = R2*mu;
R3G = R3*mu;

% Elastic predictors
eev = eps_n1(1) + eps_n1(2) + eps_n1(3); %volumetric strain
% press = bulk*eev;
eevd3 = eev/R3;
e_n1 = eps_n1 - eevd3*One; %trial elastic deviatoric strain
%convert engineering strain to physical component
str_n1 = R2G*I4*(e_n1); %trial elastic deviatoric stress
xtr_n1 = str_n1 - beta_n; %trial relative stress

% Check yield condition
normxtr = sqrt(   xtr_n1(1)*xtr_n1(1) +    xtr_n1(2)*xtr_n1(2) ...
             +    xtr_n1(3)*xtr_n1(3) + R2*xtr_n1(4)*xtr_n1(4) ...
             + R2*xtr_n1(5)*xtr_n1(5) + R2*xtr_n1(6)*xtr_n1(6));
ftr_n1 = normxtr - sqr23*(sigy + bet*Hb*a_n);

if ftr_n1 <= 1e-12
    
    % Elastic step
    s_n1 = str_n1;
%     sig_n1 = press*One + s_n1;
    Cdev_n1 = R2G*(I4 - R1/R3*OneOne);
%     C_n1 = bulk*OneOne + Cdev_n1;
    ee_n1 = eps_n1;
    beta_n1 = beta_n;
    a_n1 = a_n;
    
else
    
    % Plastic step
    n_n1 = xtr_n1/normxtr; %unit normal to yield surface
    dgama = ftr_n1/R2G/(eta/(R2G*dt) + (R1 + Hb/R3G)); %delta-gamma
    dH = sqr23*(R1 - bet)*Hb*dgama;
    a_n1 = a_n + sqr23*dgama;
    beta_n1 = beta_n + sqr23*dH*n_n1;
    ee_n1 = eps_n1 - dgama*diag([1 1 1 2 2 2])*n_n1;
    s_n1 = str_n1 - R2G*dgama*n_n1;
%     sig_n1 = press*One + s_n1;
    
    % Consistent algorithmic tangent
    th_n1 = 1-R2G*dgama/normxtr;
    thbar_n1 = R1/(eta/(R2G*dt) + R1+Hb/R3G) - (R1 - th_n1);
    Cdev_n1 = R2G*th_n1*(I4 - R1/R3*OneOne) - R2G*thbar_n1*(n_n1*n_n1');
%     C_n1 = bulk*OneOne + Cdev_n1;
    
end
