function [u,du,sxy,srt] = PureBendUS(X,Y,mu,psi,Ro,Ri,L,H)
%
% Function to calculate traction and body force for Pure Bending problem
% and stresses
% Formulas derived from PureBend.m
%
% Inputs:
% X,Y,Z = physical coordinates in reference configuration
% mu = material constant
% alpha = angle of arc for half of the beam
% N = outward unit normal
% Q = rotation tensor (accounts for having two materials with different
%     twists)
%
% Outputs:
% sxy = Cauchy Stress, x-y Cartesian coordinates
% srt = Cauchy stress, r-theta coordinates

% Formula for 1st PK stress and body force
phi = [...
 cos((Y*psi)/L)*(Ro^2 - (H*L)/psi + (2*L*X)/psi)^(1/2)
 sin((Y*psi)/L)*(Ro^2 - (H*L)/psi + (2*L*X)/psi)^(1/2)];
Rzero = (Ro^2- (H*L)/psi)^(1/2); % radius at mid-depth of beam, point (X,Y) = (0,0)
u = phi - [X+Rzero;Y];
F = ...
[ (L*cos((Y*psi)/L))/(psi*(Ro^2 - (H*L)/psi + (2*L*X)/psi)^(1/2)), -(psi*sin((Y*psi)/L)*(Ro^2 - (H*L)/psi + (2*L*X)/psi)^(1/2))/L
  (L*sin((Y*psi)/L))/(psi*(Ro^2 - (H*L)/psi + (2*L*X)/psi)^(1/2)),  (psi*cos((Y*psi)/L)*(Ro^2 - (H*L)/psi + (2*L*X)/psi)^(1/2))/L];du = F - eye(2);

P = ...
[ (mu*cos((Y*psi)/L)*(L^2 - Ro^2*psi^2 + H*L*psi - 2*L*X*psi))/(L*psi*((psi*Ro^2 + 2*L*X - H*L)/psi)^(1/2)),  (mu*sin((Y*psi)/L)*(L^2 - Ro^2*psi^2 + H*L*psi - 2*L*X*psi))/(L*psi*((psi*Ro^2 + 2*L*X - H*L)/psi)^(1/2))
  (mu*sin((Y*psi)/L)*(L^2 - Ro^2*psi^2 + H*L*psi - 2*L*X*psi))/(L*psi*((psi*Ro^2 + 2*L*X - H*L)/psi)^(1/2)), -(mu*cos((Y*psi)/L)*(L^2 - Ro^2*psi^2 + H*L*psi - 2*L*X*psi))/(L*psi*((psi*Ro^2 + 2*L*X - H*L)/psi)^(1/2))];

sxy = P*F';
coord_aver = [X Y]' + u; % mapped/deformed position on the beam of point (X,Y)
rad = sqrt((coord_aver(1)+Rzero)^2 + coord_aver(2)^2); %radius in deformed config, translating origin from (Rzero,0) to (0,0)
costheta = (coord_aver(1)+Rzero)/rad;
sintheta = coord_aver(2)/rad;
Q = [costheta sintheta
    -sintheta costheta];
srt = Q*sxy*Q';
    