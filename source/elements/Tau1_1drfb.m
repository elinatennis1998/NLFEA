function [tau,intb,Mp,Kp,vol] = Tau1_1drfb(xl,D,rho,coeffm,coeffk,Nbeta,tstep)
%
% Computes tau and integral of residual free bubble for dynamics, taken
% from RFB2.m and NL_Elem5_1dRFB.m

cwave = sqrt(D/rho);
e = sqrt(Nbeta)*cwave*tstep;
h = xl(2) - xl(1);
Me = ...
[                            h/3,                            h/6,                                                h/2 - e + (2*e)/(exp(h/e) + 1) 
                             h/6,                            h/3,                                                h/2 - e + (2*e)/(exp(h/e) + 1) 
  h/2 - e + (2*e)/(exp(h/e) + 1), h/2 - e + (2*e)/(exp(h/e) + 1), (3*e + h - 3*e*exp((2*h)/e) + 4*h*exp(h/e) + h*exp((2*h)/e))/(exp(h/e) + 1)^2];
Mp = Me(3,3);
Ke = ...
[  cwave^2/h, -cwave^2/h,                                                                      0 
  -cwave^2/h,  cwave^2/h,                                                                      0 
           0,          0, (cwave^2*(e*(exp((2*h)/e) - 1) - 2*h*exp(h/e)))/(e^2*(exp(h/e) + 1)^2)];
Kp = Ke(3,3);
       
intb = 2*Me(1,3);
vol = h;

tau = (coeffm*Mp+coeffk*Kp);

tau = inv(tau);

% % hack from Harari: Trial 1
% denom = Mp+(Nbeta*tstep^2)*Kp;
% fac1 = Mp/denom;
% fac2 = (Nbeta*tstep^2)*Kp/denom;
% a = e/h;
% intb = 1 + 6*a^2*(1-cosh(1/a))/(2+cosh(1/a));
% tau = intb;
% Mp = fac1*tau;
% Kp = fac2*tau/(Nbeta*tstep^2);
% tau = coeffm*Mp + coeffk*Kp;
% tau = inv(tau);

% % hack from Harari: Trial 2
% a = e/h;
% intb_old = intb;
% intb = 1 + 6*a^2*(1-cosh(1/a))/(2+cosh(1/a));
% Mp = Mp*(intb/intb_old)^2;
% % tau = 1/intb;
% % invtau = inv(tau);
% % Kp = (invtau - Mp)/(Nbeta*tstep^2);
% Kp = (intb - Mp)/(Nbeta*tstep^2);
% tau = coeffm*Mp + coeffk*Kp;
% tau = inv(tau);

% hack from Harari: Trial 3
a = e/h;
intb = h*(1 + 6*a^2*(1-cosh(1/a))/(2+cosh(1/a)));
Kp = (intb/h)*(intb - intb^2/h)/(Nbeta*tstep^2);
Mp = (intb - (Nbeta*tstep^2)*Kp);
tau = coeffm*Mp + coeffk*Kp;
tau = inv(tau);
