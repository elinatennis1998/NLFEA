function [u,du] = TensTorsU2(X,Y,Z,alpha)
%
% Function to calculate displacement for Torsion problem
% Formulas derived from Torsion_Tension.m
% Used to verify the general tension/torsion body force, which was found to
% be incorrect
%
% Inputs:
% X,Y,Z = physical coordinates in reference configuration
% alpha = angle of twist
% delta = fractional elongation = stretch - 1

% grouped constants/parameters

u = [ X*cos(pi*Z*alpha) - X - Y*sin(pi*Z*alpha)
 Y*cos(pi*Z*alpha) - Y + X*sin(pi*Z*alpha)
                                         0];
F =[ cos(pi*Z*alpha), -sin(pi*Z*alpha), - pi*Y*alpha*cos(pi*Z*alpha) - pi*X*alpha*sin(pi*Z*alpha) 
  sin(pi*Z*alpha),  cos(pi*Z*alpha),   pi*X*alpha*cos(pi*Z*alpha) - pi*Y*alpha*sin(pi*Z*alpha) 
                0,                0,                                                         1];

du = F - eye(3);