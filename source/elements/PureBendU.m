function [u,du] = PureBendU(X,Y,psi,Ro,Ri,L,H)
%
% Function to calculate traction and body force for Pure Bending problem
% Formulas derived from PureBend.m
%
% Inputs:
% X,Y,Z = physical coordinates in reference configuration
% mu = material constant
% alpha = angle of arc for half of the beam
% N = outward unit normal
% Q = rotation tensor (accounts for having two materials with different
%     twists)

% Formula for 1st PK stress and body force
phi = [...
 cos((Y*psi)/L)*(Ro^2 - (H*L)/psi + (2*L*X)/psi)^(1/2)
 sin((Y*psi)/L)*(Ro^2 - (H*L)/psi + (2*L*X)/psi)^(1/2)];
Rzero = (Ro^2- (L*H)/(psi))^(1/2);
u = phi - [X+Rzero;Y];
F = ...
[ (L*cos((Y*psi)/L))/(psi*(Ro^2 - (H*L)/psi + (2*L*X)/psi)^(1/2)), -(psi*sin((Y*psi)/L)*(Ro^2 - (H*L)/psi + (2*L*X)/psi)^(1/2))/L
  (L*sin((Y*psi)/L))/(psi*(Ro^2 - (H*L)/psi + (2*L*X)/psi)^(1/2)),  (psi*cos((Y*psi)/L)*(Ro^2 - (H*L)/psi + (2*L*X)/psi)^(1/2))/L];
du = F - eye(2);