function [T,B] = TensTorsTB2(X,Y,Z,E,v,alpha,N)
%
% Function to calculate traction and body force for Torsion problem
% Formulas derived from Torsion_Tension.m
% Used to verify the general tension/torsion body force, which was found to
% be incorrect
% Inputs:
% X,Y,Z = physical coordinates in reference configuration
% E,v = material constants
% alpha = angle of twist
% delta = fractional elongation = stretch - 1
% N = outward unit normal
% Q = rotation tensor (accounts for having two materials with different
%     twists)

mu = E/(2*(1+v));
% grouped constants/parameters

% Formula for 1st PK stress and body force
P = ...
[              0,             0, -pi*alpha*mu*(Y*cos(pi*Z*alpha) + X*sin(pi*Z*alpha)) 
               0,             0,  pi*alpha*mu*(X*cos(pi*Z*alpha) - Y*sin(pi*Z*alpha)) 
  -pi*Y*alpha*mu, pi*X*alpha*mu,                                                    0];
T = P*N;
B = ...
 [ pi^2*alpha^2*mu*(X*cos(pi*Z*alpha) - Y*sin(pi*Z*alpha))
 pi^2*alpha^2*mu*(Y*cos(pi*Z*alpha) + X*sin(pi*Z*alpha))
                                                       0];
