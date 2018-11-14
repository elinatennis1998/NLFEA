% c
% c
% c     halite:
% c
% c           Actual Roters hardening function
function [props, np1, n,vec1,vec2, stress, tt, h] = mm10_h_halite(props, np1,...
    n,vec1,vec2, stress, tt)

% Load some material parameters
if props.s_type == 6
    Qbulk = props.cp_028*ones(12,1);
    c1 = props.cp_001*ones(12,1);
    c2 = props.cp_004*ones(12,1);
    c3 = props.cp_007*ones(12,1);
    c4 = props.cp_010*ones(12,1);
    c5 = props.cp_013*ones(12,1);
    c6 = props.cp_016*ones(12,1);
    c7 = props.cp_019*ones(12,1);
    c8 = props.cp_022*ones(12,1);
elseif props.s_type == 11
    Qbulk = [props.cp_028*ones(6,1);props.cp_029*ones(6,1);props.cp_030*ones(12,1)];
    c1 = [props.cp_001*ones(6,1);props.cp_002*ones(6,1);props.cp_003*ones(12,1)];
    c2 = [props.cp_004*ones(6,1);props.cp_005*ones(6,1);props.cp_006*ones(12,1)];
    c3 = [props.cp_007*ones(6,1);props.cp_008*ones(6,1);props.cp_009*ones(12,1)];
    c4 = [props.cp_010*ones(6,1);props.cp_011*ones(6,1);props.cp_012*ones(12,1)];
    c5 = [props.cp_013*ones(6,1);props.cp_014*ones(6,1);props.cp_015*ones(12,1)];
    c6 = [props.cp_016*ones(6,1);props.cp_017*ones(6,1);props.cp_018*ones(12,1)];
    c7 = [props.cp_019*ones(6,1);props.cp_020*ones(6,1);props.cp_021*ones(12,1)];
    c8 = [props.cp_022*ones(6,1);props.cp_023*ones(6,1);props.cp_024*ones(12,1)];
else
    error(' ');
end
G = props.mu_0;
b = props.burgers;
k = props.boltzman;
theta = np1.temp;
v = 0.3;%props.nu;
dt = np1.tinc;

% Compute the shear modulus using Roter's function
[G_110, G_100, G_111] = mm10_halite_Gtemp(theta);
G = [G_110*ones(6,1); G_100*ones(6,1); G_111*ones(12,1)];
% if G < 0
%     G = -G;
% else
%     G = mm10_halite_Gtemp(theta);
% end
% Load the interaction matrices for parallel and forest dislocs
[Gmat,Hmat] = mm10_halite_GH(props.s_type);

% c
%       h = zeros(1,12);

% Get dislocation density
rho = transpose(tt); % rho^a_SSD
rho_n = n.tau_tilde(1:props.num_hard); % rho^a_SSD

rs = transpose(stress*np1.ms); % tau^a
rs = sign(real(rs)).*rs;

%           [rhoF,rhoP] = mm10_rhoFP_halite(props, np1, n, tt, alpha);
rhoF = Gmat*reshape(tt(1:props.num_hard),props.num_hard,1);
rhoP = Hmat*reshape(tt(1:props.num_hard),props.num_hard,1);

tpass = c1.*G*b.*sqrt(rhoP); % (16)
ddipole = sqrt(3).*G*b/(16*pi*(1-v))./rs; % (42)
rhoM = (2*k./(c1.*c2.*c3.*G*b^3))*theta.*sqrt(rhoF.*rhoP); % (13)

slipinc = vec1(1:props.nslip);%mm10_slipinc_halite(props, np1, n, stress, tt, alpha);
%           if real(slipinc) < 0
slipinc = sign(real(slipinc)).*slipinc;
%           end
gammadot = slipinc/dt;

% Evaluate the hardening equation
%================================Ran make ddipole = 0======================
ddipole=0; %Ran made it for now
%================================Ran make ddipole = 0======================
h = rho_n + dt*(c4/b.*sqrt(rhoP).*gammadot ...
    + c6.*ddipole/b.*rhoM.*gammadot - c5.*rho.*gammadot ...
    - c7.*exp(-Qbulk/k/theta).*rs/(k*theta).*rho.^2.*gammadot.^c8); % (18)

h = transpose(h); % Don't know why he always us row vector ...
%       for alpha = 1:12
%
%           % Get dislocation density
%           rho = tt(alpha); % rho^a_SSD
%           rho_n = n.tau_tilde(alpha); % rho^a_SSD
%
%           ms = np1.ms(1:6,alpha);
%           rs = stress*ms; % tau^a
%           if real(rs) < 0
%               rs = -rs;
%           end
%
% %           [rhoF,rhoP] = mm10_rhoFP_halite(props, np1, n, tt, alpha);
%           rhoF = Gmat(alpha,:)*reshape(tt(1:12),12,1);
%           rhoP = Hmat(alpha,:)*reshape(tt(1:12),12,1);
%
%           tpass = c1*G*b*sqrt(rhoP); % (16)
%           ddipole = sqrt(3)*G*b/(16*pi*(1-v))/(rs-tpass); % (42)
%           rhoM = (2*k/(c1*c2*c3*G*b^3))*theta*sqrt(rhoF*rhoP); % (13)
%
%           slipinc = vec1(alpha);%mm10_slipinc_halite(props, np1, n, stress, tt, alpha);
%           if real(slipinc) < 0
%               slipinc = -slipinc;
%           end
%           gammadot = slipinc/dt;
%
%         % Evaluate the hardening equation
%           h(alpha) = rho_n + dt*(c4/b*sqrt(rhoP)*gammadot ...
%               + c6*ddipole/b*rhoM*gammadot - c5*rho*gammadot ...
%               - c7*exp(-Qbulk/k/theta)*rs/(k*theta)*rho^2*gammadot^c8); % (18)
%       end
end