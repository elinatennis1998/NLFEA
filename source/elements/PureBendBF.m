function [T,B] = PureBendBF(X,Y,mu,lam,d,N)
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


% Formula for 1st PK stress and body force
P = ...
[ X*d*lam*(X*d + 1),                (Y*d*mu)/(X*d + 1) - X*Y*d^2*lam 
             Y*d*mu, (X*d*(lam + 2*mu + X*d*lam + X*d*mu))/(X*d + 1)];
B = [...
  - (d*mu)/(X*d + 1) - d*lam*(X*d + 1)
  0];
  
% % Other version to check the analytical formulas in the paper
% psi = (2*alpha);
% the= (Y*psi)/h;
% r=sqrt(2*h*X/(psi)+Ro^2-h*w/(psi));
% P2 = (mu*(h^2-(2*alpha)^2*r^2)/(2*alpha*h*r))*[cos(the) sin(the); sin(the) -cos(the)];
% k=Ro^2*(2*alpha)^4*(2*r^2-Ro^2)+h^2*(h^2+(2*alpha)^2*(w-2*X)^2);
% B2 = mu*k/(psi^2*h^2*r^3)*[cos(the); sin(the)];
T = P*N;
