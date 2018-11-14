function [ue,duex,duey,duez] = uexact_thickcylE(x,y,z,E,v,a,b,p,alpha)
% Tim Truster
% 05/30/2014
%
% Elastic solution for thick cylinder with internal pressure, according to
% Lubliner Plasticity, page 218. Assumes origin is at (x,y) = (0,0).
% a = inner radius
% b = outer radius
% alpha = 0 for plain strain, 1-2v for closed end, -2v for open end
% p is pressure, positive in compression
% E,v = Youngs modulus, Poissons ratio

ue = zeros(4,1);
duex = ue;
duey = ue;
duez = ue;
lam = E*v/((1+v)*(1-2*v));

r = sqrt(x^2 + y^2); %radius
costheta = x/r;
sintheta = y/r;
Q = [costheta sintheta 0
    -sintheta costheta 0
    0 0 1];
ba2 = (b/a)^2;

eps_z = alpha*p/(E*(ba2 - 1));
u_r = (1+v)*p/(E*(ba2 - 1))*((1-2*v)*r+b^2/r)-v*eps_z*r;
sig_r = -p/(ba2 - 1)*(b^2/r^2 - 1);
sig_t = p/(ba2 - 1)*(b^2/r^2 + 1);
sig_z = (2*v + alpha)*p/(ba2 - 1);
eps_r = 1/E*(sig_r - v*(sig_t + sig_z));
eps_t = 1/E*(sig_t - v*(sig_r + sig_z));

epsil = [eps_r 0 0
         0 eps_t 0
         0 0 eps_z];
p = (eps_r + eps_t + eps_z)/lam;
PEX = 0;
PEY = 0;
PEZ = 0;
epsil_xy = Q'*epsil*Q;
u_xy = Q'*[u_r; 0; eps_z*z];

ue(1:3) = u_xy;
ue(4) = p;
	
duex(1) = epsil_xy(1,1);
duex(2) = epsil_xy(2,1);
duex(3) = epsil_xy(3,1);
duex(4) = PEX;

duey(1) = epsil_xy(1,2);
duey(2) = epsil_xy(2,2);
duey(3) = epsil_xy(3,2);
duey(3) = PEY;

duez(1) = epsil_xy(1,3);
duez(2) = epsil_xy(2,3);
duez(3) = epsil_xy(3,3);
duez(3) = PEZ;