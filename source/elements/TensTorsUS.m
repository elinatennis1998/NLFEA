function [u,du,sxy,srt] = TensTorsUS(X,Y,Z,E,v,alpha,delta)
%
% Function to calculate displacement for Tension/Torsion problem and
% stresse
% Formulas derived from Torsion_Tension.m
%
% Inputs:
% X,Y,Z = physical coordinates in reference configuration
% alpha = angle of twist
% delta = fractional elongation = stretch - 1
%
% Outputs:
% sxy = Cauchy Stress, x-y Cartesian coordinates
% srt = Cauchy stress, r-theta coordinates

% grouped constants/parameters
term1 = pi*Z*alpha*(delta + 1);
term2 = ((2*v - 1)/(4*v - (- 12*delta^2*v^2 + 8*delta^2*v - 16*delta*v^2 + 12*delta*v + 1)^(1/2) + 2*delta*v - 1));

u = [2^(1/2)*X*cos(term1)*term2^(1/2) - X - 2^(1/2)*Y*sin(term1)*term2^(1/2)
 2^(1/2)*Y*cos(term1)*term2^(1/2) - Y + 2^(1/2)*X*sin(term1)*term2^(1/2)
  Z*(delta + 1) - Z];
F =[ 2^(1/2)*cos(term1)*term2^(1/2), -2^(1/2)*sin(term1)*term2^(1/2), - 2^(1/2)*pi*Y*alpha*cos(term1)*(delta + 1)*term2^(1/2) - 2^(1/2)*pi*X*alpha*sin(term1)*(delta + 1)*term2^(1/2)
 2^(1/2)*sin(term1)*term2^(1/2),  2^(1/2)*cos(term1)*term2^(1/2),   2^(1/2)*pi*X*alpha*cos(term1)*(delta + 1)*term2^(1/2) - 2^(1/2)*pi*Y*alpha*sin(term1)*(delta + 1)*term2^(1/2)
                0,     0, delta + 1];

du = F - eye(3);

term1 = (- 12*delta^2*v^2 + 8*delta^2*v - 16*delta*v^2 + 12*delta*v + 1)^(1/2);
term2 = (4*v - term1 + 2*delta*v - 1);
term3 = 2^(1/2)*E*pi^2*alpha^2;
term4 = (2*(v + 1));
term6 = (delta + 1);
term0 = (1 - 2*v)^(1/2)*(-1/term2)^(1/2);
term8 = term0*term6;
term7 = ((2*(2*v - 1)*term6)/term2 - 1);
term9 = 2^(1/2)*E*pi*alpha;

% Formula for 1st PK stress and body force
P = ...
[0, 0, -(term9*(Y*cos(pi*Z*alpha*term6) + X*sin(pi*Z*alpha*term6))*term8)/term4
 0, 0,  (term9*(X*cos(pi*Z*alpha*term6) - Y*sin(pi*Z*alpha*term6))*term8)/term4
 [-((pi*E*Y*alpha)/(2*term6*(v + 1)) + (2*E*pi*Y*alpha*v*term7)/((v + 1)*term2))*term6, ...
   ((pi*E*X*alpha)/(2*term6*(v + 1)) + (2*E*pi*X*alpha*v*term7)/((v + 1)*term2))*term6, ...
   -term6*((E*(1/term6^2 - 1))/(2*v + 2) + (2*E*v*term7)/(term6*(v + 1)*term2))]];

if norm([X Y]) > 0
sxy = P*F';
coord_aver = [X Y Z]' + u; % mapped/deformed position on the beam of point (X,Y)
rad = sqrt((coord_aver(1))^2 + coord_aver(2)^2); %radius in deformed config, translating origin from (Rzero,0) to (0,0)
costheta = (coord_aver(1))/rad;
sintheta = coord_aver(2)/rad;
Q = [costheta sintheta 0
    -sintheta costheta 0
    0 0 1];
srt = Q*sxy*Q';
else
sxy = zeros(3,3);
srt = zeros(3,3);
end
