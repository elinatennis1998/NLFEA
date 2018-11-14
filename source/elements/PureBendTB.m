function [T,B] = PureBendTB(X,Y,mu,psi,Ro,Ri,L,H,N)
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

B = zeros(2,1);
T = B;

if psi > 0

% Formula for 1st PK stress and body force
P = ...
[ (mu*cos((Y*psi)/L)*(L^2 - Ro^2*psi^2 + H*L*psi - 2*L*X*psi))/(L*psi*((psi*Ro^2 + 2*L*X - H*L)/psi)^(1/2)),  (mu*sin((Y*psi)/L)*(L^2 - Ro^2*psi^2 + H*L*psi - 2*L*X*psi))/(L*psi*((psi*Ro^2 + 2*L*X - H*L)/psi)^(1/2))
  (mu*sin((Y*psi)/L)*(L^2 - Ro^2*psi^2 + H*L*psi - 2*L*X*psi))/(L*psi*((psi*Ro^2 + 2*L*X - H*L)/psi)^(1/2)), -(mu*cos((Y*psi)/L)*(L^2 - Ro^2*psi^2 + H*L*psi - 2*L*X*psi))/(L*psi*((psi*Ro^2 + 2*L*X - H*L)/psi)^(1/2))];
B = [...
 (mu*cos((Y*psi)/L)*(H^2*L^2*psi^2 - 4*H*L^2*X*psi^2 - 2*H*L*Ro^2*psi^3 + L^4 + 4*L^2*X^2*psi^2 + 4*L*Ro^2*X*psi^3 + Ro^4*psi^4))/(L^2*psi^2*((psi*Ro^2 + 2*L*X - H*L)/psi)^(3/2))
 (mu*sin((Y*psi)/L)*(H^2*L^2*psi^2 - 4*H*L^2*X*psi^2 - 2*H*L*Ro^2*psi^3 + L^4 + 4*L^2*X^2*psi^2 + 4*L*Ro^2*X*psi^3 + Ro^4*psi^4))/(L^2*psi^2*((psi*Ro^2 + 2*L*X - H*L)/psi)^(3/2))];
 
% % Other version to check the analytical formulas in the paper
% psi = (2*alpha);
% the= (Y*psi)/h;
% r=sqrt(2*h*X/(psi)+Ro^2-h*w/(psi));
% P2 = (mu*(h^2-(2*alpha)^2*r^2)/(2*alpha*h*r))*[cos(the) sin(the); sin(the) -cos(the)];
% k=Ro^2*(2*alpha)^4*(2*r^2-Ro^2)+h^2*(h^2+(2*alpha)^2*(w-2*X)^2);
% B2 = mu*k/(psi^2*h^2*r^3)*[cos(the); sin(the)];

T = P*N;

end