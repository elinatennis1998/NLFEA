function [ue,duex,duey,duez] = uexact_thicksphEP(x,y,z,E,v,sigy,a,b,c,p,alpha)
% Tim Truster
% 05/17/2015
%
% Elastic-plastic solution for thick sphere with internal pressure, according to
% Lubliner Plasticity, page 208. Assumes origin is at (x,y,z) = (0,0,0).
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

r = sqrt(x^2 + y^2 + z^2); %radius
r2 = sqrt(x^2 + y^2); %radius
% theta = atan(y/x);
% phi = atan(r2/z);
costheta = x/r2;%cos(theta);
sintheta = y/r2;%sin(theta);
cosphi = z/r;%cos(phi);
sinphi = r2/r;%sin(phi);
Q = [sinphi*costheta sinphi*sintheta cosphi
     cosphi*costheta cosphi*sintheta -sinphi
     -sintheta costheta 0];
ba3 = (b/a)^3;

% Determine if yielded
pE = 2/3*sigy*(1-1/ba3);

if pE > p % elastic everywhere

sig_r = -p/(ba3 - 1)*(b^3/r^3 - 1);
sig_t = p/(ba3 - 1)*(b^3/(2*r^3) + 1);
sig_p = sig_t;
u_r = (r/E)*((1-v)*sig_t - v*sig_r);
eps_r = 1/E*(sig_r - v*(sig_t + sig_p));
eps_t = 1/E*(sig_t - v*(sig_r + sig_p));
eps_p = 1/E*(sig_p - v*(sig_r + sig_t));

epsil = [eps_r 0 0
         0 eps_t 0
         0 0 eps_p];
pr = (eps_r + eps_t + eps_p)/lam;
PEX = 0;
PEY = 0;
PEZ = 0;
epsil_xy = Q'*epsil*Q;
u_xy = Q'*[u_r; 0; 0];

ue(1:3) = u_xy;
ue(4) = pr;
	
duex(1) = epsil_xy(1,1);
duex(2) = epsil_xy(2,1);
duex(3) = epsil_xy(3,1);
duex(4) = PEX;

duey(1) = epsil_xy(1,2);
duey(2) = epsil_xy(2,2);
duey(3) = epsil_xy(3,2);
duey(4) = PEY;

duez(1) = epsil_xy(1,3);
duez(2) = epsil_xy(2,3);
duez(3) = epsil_xy(3,3);
duez(4) = PEZ;

else % plastic somewhere
    
    if r > c % elastic
        
        bc3 = (b/c)^3;
        
        sig_r = -2/3*sigy*(c^3/r^3 - 1/bc3);
        sig_t = 2/3*sigy*(c^3/(2*r^3) + 1/bc3);
        sig_p = sig_t;
        eps_r = 1/E*(sig_r - v*(sig_t + sig_p));
        eps_t = 1/E*(sig_t - v*(sig_r + sig_p));
        eps_p = 1/E*(sig_p - v*(sig_r + sig_t));
        u_r = (r/E)*((1-v)*sig_t - v*sig_r);

        epsil = [eps_r 0 0
                 0 eps_t 0
                 0 0 eps_p];
        pr = (eps_r + eps_t + eps_p)/lam;
        PEX = 0;
        PEY = 0;
        PEZ = 0;
        epsil_xy = Q'*epsil*Q;
        u_xy = Q'*[u_r; 0; 0];

        ue(1:3) = u_xy;
        ue(4) = pr;

        duex(1) = epsil_xy(1,1);
        duex(2) = epsil_xy(2,1);
        duex(3) = epsil_xy(3,1);
        duex(4) = PEX;

        duey(1) = epsil_xy(1,2);
        duey(2) = epsil_xy(2,2);
        duey(3) = epsil_xy(3,2);
        duey(4) = PEY;

        duez(1) = epsil_xy(1,3);
        duez(2) = epsil_xy(2,3);
        duez(3) = epsil_xy(3,3);
        duez(4) = PEZ;
        
    else % plastic
        
        bc3 = (b/c)^3;
        pc = 2/3*sigy*(1-1/bc3 + log(c^3/a^3));
        
%         sig_r = -2/3*sigy*(1-1/bc3 + log(c^3/r^3));
%         sig_t = 1/3*sigy*(1+2/bc3 - 2*log(c^3/r^3));
%         sig_p = sig_t;
        u_r = sigy*r/E*((1-v)*c^3/r^3 + (1-2*v)*(2*log(r/a)-pc/sigy));
%         eps_r = 1/E*(sig_r - v*(sig_t + sig_z)); % no longer valid
%         eps_t = 1/E*(sig_t - v*(sig_r + sig_z)); % no longer valid
        eps_r = sigy/E*((1-v)*c^3/r^3 + (1-2*v)*(2*log(r/a)-pc/sigy)) ...
              + sigy*r/E*(-3*(1-v)*c^3/r^4 + (1-2*v)*(2/(r/a)*(1/a)));
        eps_t = u_r/r; % simplifies
        eps_p = u_r/r; % simplifies
        % confirmed 6/4/14 to give proper convergence rates

        epsil = [eps_r 0 0
                 0 eps_t 0
                 0 0 eps_p];
        pr = (eps_r + eps_t + eps_p)/lam;
        PEX = 0;
        PEY = 0;
        PEZ = 0;
        epsil_xy = Q'*epsil*Q;
        u_xy = Q'*[u_r; 0; 0];

        ue(1:3) = u_xy;
        ue(4) = pr;

        duex(1) = epsil_xy(1,1);
        duex(2) = epsil_xy(2,1);
        duex(3) = epsil_xy(3,1);
        duex(4) = PEX;

        duey(1) = epsil_xy(1,2);
        duey(2) = epsil_xy(2,2);
        duey(3) = epsil_xy(3,2);
        duey(4) = PEY;

        duez(1) = epsil_xy(1,3);
        duez(2) = epsil_xy(2,3);
        duez(3) = epsil_xy(3,3);
        duez(4) = PEZ;
        
    end
    
    
end
