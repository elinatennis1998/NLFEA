function [T,B] = BlockJ1TB(X,Y,Z,E,v,N,Q)
%
% Function to calculate traction and body force for Tension/Torsion problem
% Formulas derived from Torsion_Tension.m
%
% Inputs:
% X,Y,Z = physical coordinates in reference configuration
% E,v = material constants
% alpha = angle of twist
% delta = fractional elongation = stretch - 1
% N = outward unit normal
% Q = rotation tensor (accounts for having two materials with different
%     twists)

mu = E/(2*(1+v));

% Formula for 1st PK stress and body force
P = ...
[ -(mu*((X - 1)^4 - 1))/(X - 1)^2,            Y*mu*(X - 1),            Z*mu*(X - 1) 
                            -Y*mu, -(X*mu*(X - 2))/(X - 1),                       0 
                            -Z*mu,                       0, -(X*mu*(X - 2))/(X - 1)];
T = P*N;
B = ...
 [2*mu*(X - 1) - (2*mu*((X - 1)^4 - 1))/(X - 1)^3
 0
 0];

% Apply rotations
T = Q*T;
B = Q*B;
