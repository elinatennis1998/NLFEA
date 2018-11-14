function [u,du] = TensTorsU(X,Y,Z,v,alpha,delta)
%
% Function to calculate displacement for Tension/Torsion problem
% Formulas derived from Torsion_Tension.m
%
% Inputs:
% X,Y,Z = physical coordinates in reference configuration
% alpha = angle of twist
% delta = fractional elongation = stretch - 1

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