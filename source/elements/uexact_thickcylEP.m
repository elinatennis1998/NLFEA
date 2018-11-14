function [ue,duex,duey,duez] = uexact_thickcylEP(x,y,z,E,v,k,a,b,c,p,alpha)
% Tim Truster
% 05/30/2014
%
% Elastic-plastic solution for thick cylinder with internal pressure, according to
% Lubliner Plasticity, page 218. Assumes origin is at (x,y) = (0,0).
% Assumes Tresca yield condition.
% a = inner radius
% b = outer radius
% alpha = 0 for plain strain, 1-2v for closed end, -2v for open end
% p is pressure, positive in compression
% E,v = Youngs modulus, Poissons ratio
% k = shear yield stress, Lubliner page 137 = SIGMAY/2

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

% Determine if yielded
pE = k*(1-1/ba2);

if pE > p % elastic everywhere

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

else % plastic somewhere
    
    if r > c % elastic
        
        bc2 = (b/c)^2;
        pb = k*(1-1/bc2 + log(c^2/a^2));
        
        % Formulas completely verified for plane strain, using either
        % displacement or stress definitions for the strains
        eps_z = alpha*pb/(E*(ba2 - 1));
        sig_r = -k*(c^2/r^2 - 1/bc2);
        sig_t = k*(c^2/r^2 + 1/bc2);
        sig_z = k*(2*v/bc2 + alpha/(ba2 - 1)*(1-1/bc2 + log(c^2/a^2)));
        eps_r = 1/E*(sig_r - v*(sig_t + sig_z));
        eps_t = 1/E*(sig_t - v*(sig_r + sig_z));
        A = 2*k*c^2*(1-v^2)/E;
        u_r = (1+v)*(1-2*v)/E*r*sig_r-v*eps_z*r + A/r;

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
        
    else % plastic
        
        bc2 = (b/c)^2;
        pb = k*(1-1/bc2 + log(c^2/a^2));
        
        eps_z = alpha*pb/(E*(ba2 - 1));
        sig_r = -k*(1-1/bc2 + log(c^2/r^2));
        sig_t = k*(1+1/bc2 - log(c^2/r^2));
        sig_z = k*(2*v*(1/bc2 - log(c^2/r^2)) + alpha/(ba2 - 1)*(1-1/bc2 + log(c^2/a^2)));
        A = 2*k*c^2*(1-v^2)/E;
        u_r = (1+v)*(1-2*v)/E*r*sig_r-v*eps_z*r + A/r;
%         eps_r = 1/E*(sig_r - v*(sig_t + sig_z)); % no longer valid
%         eps_t = 1/E*(sig_t - v*(sig_r + sig_z)); % no longer valid
        dsigr = -k*(1/(c^2/r^2)*(-2*c^2/r^3));
        eps_r = (1+v)*(1-2*v)/E*(sig_r+r*dsigr)-v*eps_z - A/r^2;
        eps_t = u_r/r;
        % confirmed 6/4/14 to give proper convergence rates

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
        
    end
    
    
end
