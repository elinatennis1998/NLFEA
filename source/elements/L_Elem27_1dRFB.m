% Tim Truster
% 2/28/2014
% 1d reaction-diffusion
% Includes Harari's tau as an option
% Tested for linear shape functions and Harari's tau.
%
% I found that the discrete solution does oscillate around the exact
% solution when k >> 1. Adding Harari's tau damps the results closer to
% nodally exact.

% Set Material Properties

ksquar = mateprop(1);
cwave = 1;
t_on = 1;0;

switch isw %Task Switch
%%    
    case 1
        
        if ndf > 1
            
            for i = 2:ndf
                lie(i,1) = 0;
            end
            
        end
        
%%
    case 3 %Compute Stiffness and Residual
        
        k = sqrt(ksquar);
        h = xl(2) - xl(1);
        
% Me = ...
% [                                           (h*k^2)/3,                                           (h*k^2)/6,                               k^2*(h/(2*k^2) - sinh((h*k)/2)/(k^3*cosh((h*k)/2)))
%                                             (h*k^2)/6,                                           (h*k^2)/3,                               k^2*(h/(2*k^2) - sinh((h*k)/2)/(k^3*cosh((h*k)/2)))
%   k^2*(h/(2*k^2) - sinh((h*k)/2)/(k^3*cosh((h*k)/2))), k^2*(h/(2*k^2) - sinh((h*k)/2)/(k^3*cosh((h*k)/2))), (k*(h + h*exp(2*h*k)) - 3*exp(2*h*k) + 4*h*k*exp(h*k) + 3)/(k^3*(exp(h*k) + 1)^2)];
% Ke = ...
% [  1/h, -1/h,                                                         0
%   -1/h,  1/h,                                                         0
%      0,    0, -(2*h*k*exp(h*k) - exp(2*h*k) + 1)/(k^3*(exp(h*k) + 1)^2)];
Me = ...
[                                                                            (h*k^2)/3,                                                                           (h*k^2)/6,                               k^2*(h/(2*k^2) - sinh((h*k)/2)/(k^3*cosh((h*k)/2))),                        -k^2*(h/(6*k^2) + (2*(sinh((h*k)/2) - (h*k*cosh((h*k)/2))/2))/(h*k^4*sinh((h*k)/2)))
                                                                             (h*k^2)/6,                                                                           (h*k^2)/3,                               k^2*(h/(2*k^2) - sinh((h*k)/2)/(k^3*cosh((h*k)/2))),                         k^2*(h/(6*k^2) + (2*(sinh((h*k)/2) - (h*k*cosh((h*k)/2))/2))/(h*k^4*sinh((h*k)/2)))
                                   k^2*(h/(2*k^2) - sinh((h*k)/2)/(k^3*cosh((h*k)/2))),                                 k^2*(h/(2*k^2) - sinh((h*k)/2)/(k^3*cosh((h*k)/2))), (k*(h + h*exp(2*h*k)) - 3*exp(2*h*k) + 4*h*k*exp(h*k) + 3)/(k^3*(exp(h*k) + 1)^2),                                                                                                           0
  -k^2*(h/(6*k^2) + (2*(sinh((h*k)/2) - (h*k*cosh((h*k)/2))/2))/(h*k^4*sinh((h*k)/2))), k^2*(h/(6*k^2) + (2*(sinh((h*k)/2) - (h*k*cosh((h*k)/2))/2))/(h*k^4*sinh((h*k)/2))),                                                                                 0, (((h^2*cosh(h*k))/6 - (2*h^2)/3)/h - ((3*h*k*sinh(h*k))/2 - 4*cosh(h*k) + 4)/(h*k^2))/(k^2*sinh((h*k)/2)^2)];
Ke = ...
[  1/h, -1/h,                                                         0,                                                                                       0
  -1/h,  1/h,                                                         0,                                                                                       0
     0,    0, -(2*h*k*exp(h*k) - exp(2*h*k) + 1)/(k^3*(exp(h*k) + 1)^2),                                                                                       0
     0,    0,                                                         0, ((h*k^2)/2 + k*sinh((h*k)/2)*(2*sinh((h*k)/4)^2 + 1))/(k^4*sinh((h*k)/2)^2) - 4/(h*k^4)];

        if t_on
        [tau,intb,Mp,Kp,vol] = Tau1_1drfb(xl,1,ksquar,1,1,1,1);
        tau_H = intb/h;
        else
        tau_H = 0;
        end

        % Compute kinematic fields at intermediate time levels
        ulres = reshape(ul,ndf*nel,1);
        
        ElemK = (1-tau_H)*Me(1:nel,1:nel) + Ke(1:nel,1:nel);
        if nel == 4
        Mstar = ElemK(1:2,1:2) - ElemK(1:2,3:4)*inv(ElemK(3:4,3:4))*ElemK(3:4,1:2);
        end
        ElemF = - ElemK*ulres;

    case 6 %Compute Residual
        
        k = sqrt(ksquar);
        h = xl(2) - xl(1);
        
Me = ...
[                                                                            (h*k^2)/3,                                                                           (h*k^2)/6,                               k^2*(h/(2*k^2) - sinh((h*k)/2)/(k^3*cosh((h*k)/2))),                        -k^2*(h/(6*k^2) + (2*(sinh((h*k)/2) - (h*k*cosh((h*k)/2))/2))/(h*k^4*sinh((h*k)/2)))
                                                                             (h*k^2)/6,                                                                           (h*k^2)/3,                               k^2*(h/(2*k^2) - sinh((h*k)/2)/(k^3*cosh((h*k)/2))),                         k^2*(h/(6*k^2) + (2*(sinh((h*k)/2) - (h*k*cosh((h*k)/2))/2))/(h*k^4*sinh((h*k)/2)))
                                   k^2*(h/(2*k^2) - sinh((h*k)/2)/(k^3*cosh((h*k)/2))),                                 k^2*(h/(2*k^2) - sinh((h*k)/2)/(k^3*cosh((h*k)/2))), (k*(h + h*exp(2*h*k)) - 3*exp(2*h*k) + 4*h*k*exp(h*k) + 3)/(k^3*(exp(h*k) + 1)^2),                                                                                                           0
  -k^2*(h/(6*k^2) + (2*(sinh((h*k)/2) - (h*k*cosh((h*k)/2))/2))/(h*k^4*sinh((h*k)/2))), k^2*(h/(6*k^2) + (2*(sinh((h*k)/2) - (h*k*cosh((h*k)/2))/2))/(h*k^4*sinh((h*k)/2))),                                                                                 0, (((h^2*cosh(h*k))/6 - (2*h^2)/3)/h - ((3*h*k*sinh(h*k))/2 - 4*cosh(h*k) + 4)/(h*k^2))/(k^2*sinh((h*k)/2)^2)];
Ke = ...
[  1/h, -1/h,                                                         0,                                                                                       0
  -1/h,  1/h,                                                         0,                                                                                       0
     0,    0, -(2*h*k*exp(h*k) - exp(2*h*k) + 1)/(k^3*(exp(h*k) + 1)^2),                                                                                       0
     0,    0,                                                         0, ((h*k^2)/2 + k*sinh((h*k)/2)*(2*sinh((h*k)/4)^2 + 1))/(k^4*sinh((h*k)/2)^2) - 4/(h*k^4)];

        if t_on
        [tau,intb,Mp,Kp,vol] = Tau1_1drfb(xl,1,ksquar,1,1,1,1);
        tau_H = intb/h;
        else
        tau_H = 0;
        end

        % Compute kinematic fields at intermediate time levels
        ulres = reshape(ul,ndf*nel,1);
        
        ElemK = (1-tau_H)*Me(1:nel,1:nel) + Ke(1:nel,1:nel);
        ElemF = - ElemK*ulres;

    case 21 %Compute Stiffness
        
        k = sqrt(ksquar);
        h = xl(2) - xl(1);
 
Me = ...
[                                                                            (h*k^2)/3,                                                                           (h*k^2)/6,                               k^2*(h/(2*k^2) - sinh((h*k)/2)/(k^3*cosh((h*k)/2))),                        -k^2*(h/(6*k^2) + (2*(sinh((h*k)/2) - (h*k*cosh((h*k)/2))/2))/(h*k^4*sinh((h*k)/2)))
                                                                             (h*k^2)/6,                                                                           (h*k^2)/3,                               k^2*(h/(2*k^2) - sinh((h*k)/2)/(k^3*cosh((h*k)/2))),                         k^2*(h/(6*k^2) + (2*(sinh((h*k)/2) - (h*k*cosh((h*k)/2))/2))/(h*k^4*sinh((h*k)/2)))
                                   k^2*(h/(2*k^2) - sinh((h*k)/2)/(k^3*cosh((h*k)/2))),                                 k^2*(h/(2*k^2) - sinh((h*k)/2)/(k^3*cosh((h*k)/2))), (k*(h + h*exp(2*h*k)) - 3*exp(2*h*k) + 4*h*k*exp(h*k) + 3)/(k^3*(exp(h*k) + 1)^2),                                                                                                           0
  -k^2*(h/(6*k^2) + (2*(sinh((h*k)/2) - (h*k*cosh((h*k)/2))/2))/(h*k^4*sinh((h*k)/2))), k^2*(h/(6*k^2) + (2*(sinh((h*k)/2) - (h*k*cosh((h*k)/2))/2))/(h*k^4*sinh((h*k)/2))),                                                                                 0, (((h^2*cosh(h*k))/6 - (2*h^2)/3)/h - ((3*h*k*sinh(h*k))/2 - 4*cosh(h*k) + 4)/(h*k^2))/(k^2*sinh((h*k)/2)^2)];
Ke = ...
[  1/h, -1/h,                                                         0,                                                                                       0
  -1/h,  1/h,                                                         0,                                                                                       0
     0,    0, -(2*h*k*exp(h*k) - exp(2*h*k) + 1)/(k^3*(exp(h*k) + 1)^2),                                                                                       0
     0,    0,                                                         0, ((h*k^2)/2 + k*sinh((h*k)/2)*(2*sinh((h*k)/4)^2 + 1))/(k^4*sinh((h*k)/2)^2) - 4/(h*k^4)];

        if t_on
        [tau,intb,Mp,Kp,vol] = Tau1_1drfb(xl,1,ksquar,1,1,1,1);
        tau_H = intb/h;
        else
        tau_H = 0;
        end
        
        ElemK = (1-tau_H)*Me(1:nel,1:nel) + Ke(1:nel,1:nel);
        
    case 15 % Body force
        
        k = sqrt(ksquar);
        
        if iprob == 1
% %         syms x real
%         syms x x1 x2 k wave real
% %         h = xl(2) - xl(1);
%         h = x2 - x1;%xl(2) - xl(1);
%         N1 = (h/2-x)/h;
%         N2 = (h/2+x)/h;
% %         be = -(cosh(k*x) - cosh(k*h/2))/(k^2*cosh(k*h/2)); % RFB1
%         be = -(sinh(k*x) - x/(h/2)*sinh(k*h/2))/(k^2*sinh(k*h/2)); % RFB2
%         Nmat = [N1 N2 be];
% %         y =  x + (xl(2) + xl(1))/2;
%         y =  x + (x2 + x1)/2;
%         
% %         u = 32 - y^5 - (2-y)^5 + 4*sin(pi*y*wave);
%         f = 20*(y^3 + (2-y)^3) + 4*(pi*wave)^2*sin(pi*y*wave) + k^2*(32 - y^5 - (2-y)^5 + 4*sin(pi*y*wave));
% %         ElemF = double(int(Nmat'*f,x,-h/2,h/2));
%         ElemF = int(Nmat'*f,x,-h/2,h/2);

        if t_on
        [tau,intb,Mp,Kp,vol] = Tau1_1drfb(xl,1,ksquar,1,1,1,1);
        tau_H = intb/h;
        else
        tau_H = 0;
        end

        x2 = xl(2);
        x1 = xl(1);
%         ElemF = [ ((x1 - x2)*(240*x1 + 120*x2 - 80*ksquar*x1 - 40*ksquar*x2 - 60*x1*x2 + 60*ksquar*x1^2 - 24*ksquar*x1^3 + 20*ksquar*x2^2 + 5*ksquar*x1^4 - 6*ksquar*x2^3 + ksquar*x2^4 - 90*x1^2 - 30*x2^2 - 12*ksquar*x1*x2^2 - 18*ksquar*x1^2*x2 + 2*ksquar*x1*x2^3 + 4*ksquar*x1^3*x2 + 3*ksquar*x1^2*x2^2 + 40*ksquar*x1*x2 - 240))/3
%  ((x1 - x2)*(120*x1 + 240*x2 - 40*ksquar*x1 - 80*ksquar*x2 - 60*x1*x2 + 20*ksquar*x1^2 - 6*ksquar*x1^3 + 60*ksquar*x2^2 + ksquar*x1^4 - 24*ksquar*x2^3 + 5*ksquar*x2^4 - 30*x1^2 - 90*x2^2 - 18*ksquar*x1*x2^2 - 12*ksquar*x1^2*x2 + 4*ksquar*x1*x2^3 + 2*ksquar*x1^3*x2 + 3*ksquar*x1^2*x2^2 + 40*ksquar*x1*x2 - 240))/3];
        ElemF = (1-tau_H)*[ -(12*ksquar*sin(pi*wave*x1) - 12*ksquar*sin(pi*wave*x2) + 240*pi^2*wave^2*x1^2 - 240*pi^2*wave^2*x1^3 + 240*pi^2*wave^2*x2^2 + 90*pi^2*wave^2*x1^4 - 120*pi^2*wave^2*x2^3 + 30*pi^2*wave^2*x2^4 + 12*pi^2*wave^2*sin(pi*wave*x1) - 12*pi^2*wave^2*sin(pi*wave*x2) + 80*pi^2*ksquar*wave^2*x1^3 - 60*pi^2*ksquar*wave^2*x1^4 + 40*pi^2*ksquar*wave^2*x2^3 + 24*pi^2*ksquar*wave^2*x1^5 - 20*pi^2*ksquar*wave^2*x2^4 - 5*pi^2*ksquar*wave^2*x1^6 + 6*pi^2*ksquar*wave^2*x2^5 - pi^2*ksquar*wave^2*x2^6 + 360*pi^2*wave^2*x1^2*x2 - 120*pi^2*wave^2*x1^3*x2 - 12*pi^3*wave^3*x1*cos(pi*wave*x1) + 12*pi^3*wave^3*x2*cos(pi*wave*x1) - 480*pi^2*wave^2*x1*x2 - 12*pi*ksquar*wave*x1*cos(pi*wave*x1) + 12*pi*ksquar*wave*x2*cos(pi*wave*x1) - 120*pi^2*ksquar*wave^2*x1^2*x2 + 80*pi^2*ksquar*wave^2*x1^3*x2 - 30*pi^2*ksquar*wave^2*x1^4*x2 + 6*pi^2*ksquar*wave^2*x1^5*x2)/(3*pi^2*wave^2*(x1 - x2))
 -(12*ksquar*sin(pi*wave*x2) - 12*ksquar*sin(pi*wave*x1) + 240*pi^2*wave^2*x1^2 - 120*pi^2*wave^2*x1^3 + 240*pi^2*wave^2*x2^2 + 30*pi^2*wave^2*x1^4 - 240*pi^2*wave^2*x2^3 + 90*pi^2*wave^2*x2^4 - 12*pi^2*wave^2*sin(pi*wave*x1) + 12*pi^2*wave^2*sin(pi*wave*x2) + 40*pi^2*ksquar*wave^2*x1^3 - 20*pi^2*ksquar*wave^2*x1^4 + 80*pi^2*ksquar*wave^2*x2^3 + 6*pi^2*ksquar*wave^2*x1^5 - 60*pi^2*ksquar*wave^2*x2^4 - pi^2*ksquar*wave^2*x1^6 + 24*pi^2*ksquar*wave^2*x2^5 - 5*pi^2*ksquar*wave^2*x2^6 + 360*pi^2*wave^2*x1*x2^2 - 120*pi^2*wave^2*x1*x2^3 + 12*pi^3*wave^3*x1*cos(pi*wave*x2) - 12*pi^3*wave^3*x2*cos(pi*wave*x2) - 480*pi^2*wave^2*x1*x2 + 12*pi*ksquar*wave*x1*cos(pi*wave*x2) - 12*pi*ksquar*wave*x2*cos(pi*wave*x2) - 120*pi^2*ksquar*wave^2*x1*x2^2 + 80*pi^2*ksquar*wave^2*x1*x2^3 - 30*pi^2*ksquar*wave^2*x1*x2^4 + 6*pi^2*ksquar*wave^2*x1*x2^5)/(3*pi^2*wave^2*(x1 - x2))];
if nel >= 3
ElemF(3) = (80*x1^3)/3 - 40*x1^2 + 40*x2^2 - 10*x1^4 - (80*x2^3)/3 + 2*x1^5 + 10*x2^4 - 2*x2^5 + (4*sin(pi*wave*x1)*tanh((k*(x1 - x2))/2))/k + (4*sin(pi*wave*x2)*tanh((k*(x1 - x2))/2))/k + (80*x1*tanh((k*(x1 - x2))/2))/k + (80*x2*tanh((k*(x1 - x2))/2))/k + (4*cos(pi*wave*x1))/(pi*wave) - (4*cos(pi*wave*x2))/(pi*wave) - (80*x1^2*tanh((k*(x1 - x2))/2))/k + (40*x1^3*tanh((k*(x1 - x2))/2))/k - (80*x2^2*tanh((k*(x1 - x2))/2))/k - (10*x1^4*tanh((k*(x1 - x2))/2))/k + (40*x2^3*tanh((k*(x1 - x2))/2))/k - (10*x2^4*tanh((k*(x1 - x2))/2))/k;
end
if nel == 4
ElemF(4) = (8*sin(pi*wave*x1) - 8*sin(pi*wave*x2) - 4*pi*wave*x1*cos(pi*wave*x1) - 4*pi*wave*x1*cos(pi*wave*x2) + 4*pi*wave*x2*cos(pi*wave*x1) + 4*pi*wave*x2*cos(pi*wave*x2))/(pi^2*wave^2*(x1 - x2)) - (480*x1^2*sinh((k*(x1 - x2))/2) - 240*x1^3*sinh((k*(x1 - x2))/2) - 480*x2^2*sinh((k*(x1 - x2))/2) + 60*x1^4*sinh((k*(x1 - x2))/2) + 240*x2^3*sinh((k*(x1 - x2))/2) - 60*x2^4*sinh((k*(x1 - x2))/2) - 24*sin(pi*wave*x1)*sinh((k*(x1 - x2))/2) + 24*sin(pi*wave*x2)*sinh((k*(x1 - x2))/2) - 480*x1*sinh((k*(x1 - x2))/2) + 480*x2*sinh((k*(x1 - x2))/2) + 240*k*x1^2*cosh((k*(x1 - x2))/2) - 240*k*x1^3*cosh((k*(x1 - x2))/2) + 240*k*x2^2*cosh((k*(x1 - x2))/2) + 120*k*x1^4*cosh((k*(x1 - x2))/2) - 240*k*x2^3*cosh((k*(x1 - x2))/2) - 30*k*x1^5*cosh((k*(x1 - x2))/2) + 120*k*x2^4*cosh((k*(x1 - x2))/2) - 30*k*x2^5*cosh((k*(x1 - x2))/2) - 40*k^2*x1^3*sinh((k*(x1 - x2))/2) + 40*k^2*x1^4*sinh((k*(x1 - x2))/2) + 40*k^2*x2^3*sinh((k*(x1 - x2))/2) - 18*k^2*x1^5*sinh((k*(x1 - x2))/2) - 40*k^2*x2^4*sinh((k*(x1 - x2))/2) + 4*k^2*x1^6*sinh((k*(x1 - x2))/2) + 18*k^2*x2^5*sinh((k*(x1 - x2))/2) - 4*k^2*x2^6*sinh((k*(x1 - x2))/2) + 12*k*x1*sin(pi*wave*x1)*cosh((k*(x1 - x2))/2) - 12*k*x1*sin(pi*wave*x2)*cosh((k*(x1 - x2))/2) - 12*k*x2*sin(pi*wave*x1)*cosh((k*(x1 - x2))/2) + 12*k*x2*sin(pi*wave*x2)*cosh((k*(x1 - x2))/2) - 120*k^2*x1*x2^2*sinh((k*(x1 - x2))/2) + 120*k^2*x1^2*x2*sinh((k*(x1 - x2))/2) + 80*k^2*x1*x2^3*sinh((k*(x1 - x2))/2) - 80*k^2*x1^3*x2*sinh((k*(x1 - x2))/2) - 30*k^2*x1*x2^4*sinh((k*(x1 - x2))/2) + 30*k^2*x1^4*x2*sinh((k*(x1 - x2))/2) + 6*k^2*x1*x2^5*sinh((k*(x1 - x2))/2) - 6*k^2*x1^5*x2*sinh((k*(x1 - x2))/2) - 480*k*x1*x2*cosh((k*(x1 - x2))/2) + 240*k*x1*x2^2*cosh((k*(x1 - x2))/2) + 240*k*x1^2*x2*cosh((k*(x1 - x2))/2) - 120*k*x1*x2^3*cosh((k*(x1 - x2))/2) - 120*k*x1^3*x2*cosh((k*(x1 - x2))/2) + 30*k*x1*x2^4*cosh((k*(x1 - x2))/2) + 30*k*x1^4*x2*cosh((k*(x1 - x2))/2))/(3*k^2*sinh((k*(x1 - x2))/2)*(x1 - x2));
% Me = ...
% [                                                                            (h*k^2)/3,                                                                           (h*k^2)/6,                               k^2*(h/(2*k^2) - sinh((h*k)/2)/(k^3*cosh((h*k)/2))),                        -k^2*(h/(6*k^2) + (2*(sinh((h*k)/2) - (h*k*cosh((h*k)/2))/2))/(h*k^4*sinh((h*k)/2)))
%                                                                              (h*k^2)/6,                                                                           (h*k^2)/3,                               k^2*(h/(2*k^2) - sinh((h*k)/2)/(k^3*cosh((h*k)/2))),                         k^2*(h/(6*k^2) + (2*(sinh((h*k)/2) - (h*k*cosh((h*k)/2))/2))/(h*k^4*sinh((h*k)/2)))
%                                    k^2*(h/(2*k^2) - sinh((h*k)/2)/(k^3*cosh((h*k)/2))),                                 k^2*(h/(2*k^2) - sinh((h*k)/2)/(k^3*cosh((h*k)/2))), (k*(h + h*exp(2*h*k)) - 3*exp(2*h*k) + 4*h*k*exp(h*k) + 3)/(k^3*(exp(h*k) + 1)^2),                                                                                                           0
%   -k^2*(h/(6*k^2) + (2*(sinh((h*k)/2) - (h*k*cosh((h*k)/2))/2))/(h*k^4*sinh((h*k)/2))), k^2*(h/(6*k^2) + (2*(sinh((h*k)/2) - (h*k*cosh((h*k)/2))/2))/(h*k^4*sinh((h*k)/2))),                                                                                 0, (((h^2*cosh(h*k))/6 - (2*h^2)/3)/h - ((3*h*k*sinh(h*k))/2 - 4*cosh(h*k) + 4)/(h*k^2))/(k^2*sinh((h*k)/2)^2)];
% Ke = ...
% [  1/h, -1/h,                                                         0,                                                                                       0
%   -1/h,  1/h,                                                         0,                                                                                       0
%      0,    0, -(2*h*k*exp(h*k) - exp(2*h*k) + 1)/(k^3*(exp(h*k) + 1)^2),                                                                                       0
%      0,    0,                                                         0, ((h*k^2)/2 + k*sinh((h*k)/2)*(2*sinh((h*k)/4)^2 + 1))/(k^4*sinh((h*k)/2)^2) - 4/(h*k^4)];
% ElemK = (1-tau_H)*Me(1:nel,1:nel) + Ke(1:nel,1:nel);
% Fstar = ElemF(1:2) - ElemK(1:2,3:4)*inv(ElemK(3:4,3:4))*ElemF(3:4);
end


        end
        
    case 11
        
        el2el = 0;
        eprixel = 0;
        
        %Set integration number
        lint = 5;
        
        k = sqrt(ksquar);
        
        ib = 0;
        bf = 1;
        der = 1;

        % Compute kinematic fields at intermediate time levels
        ulres = reshape(ul,ndf*nel,1);
        h = xl(2) - xl(1);
        
        sw = int1d(lint);
        
        for ll = 1:lint                    

            % Evaluate 1-D basis functions at integration points
            Wgt = sw(2,ll);
            
            cint = Wgt*h/2;
            
            if iprob == 1
                x = h/2*sw(1,ll);
                y =  x + (xl(2) + xl(1))/2;
                N1 = (h/2-x)/h;
                N2 = (h/2+x)/h;
                Nmat = [N1 N2];
                B1 = -1/h;
                B2 = 1/h;
                Bmat = [B1 B2];
                if nel >= 3
                be = -(cosh(k*x) - cosh(k*h/2))/(k^2*cosh(k*h/2));
                dbe = -sinh(k*x)/(k*cosh((h*k)/2));
                Nmat = [Nmat be];
                Bmat = [Bmat dbe];
                end
                if nel == 4
                be = -(sinh(k*x) - x/(h/2)*sinh((h*k)/2))/(k^2*sinh((h*k)/2));
                dbe = -(k*cosh(k*x) - (2*sinh((h*k)/2))/h)/(k^2*sinh((h*k)/2));
                Nmat = [Nmat be];
                Bmat = [Bmat dbe];
                end
                u = Nmat*ulres;
                du = Bmat*ulres;
                ue = 32 - y^5 - (2-y)^5 + 4*sin(pi*y*wave);
                due = - 5*y^4 + 5*(2-y)^4 + 4*pi*wave*cos(pi*y*wave);
            else
                ue = 0;
                due = 0;
            end
            
            el2el   = el2el   + cint * ( (u-ue)^2 );
            eprixel = eprixel + cint * ( (du-due)^2 );

        end %je
        
        ElemE(1:2) = [el2el eprixel];
        
end