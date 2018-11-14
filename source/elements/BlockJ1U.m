function [u,du] = BlockJ1U(X,Y,Z)
%
% Function to calculate displacement for Tension/Torsion problem
% Formulas derived from Torsion_Tension.m
%
% Inputs:
% X,Y,Z = physical coordinates in reference configuration
% alpha = angle of twist
% delta = fractional elongation = stretch - 1

u = [ - X - 1/(X - 1) - 1 % to translate it back to origin
 - Y - Y*(X - 1)
 - Z - Z*(X - 1)];
F =[ 1/(X - 1)^2,     0,     0 
           -Y, 1 - X,     0 
           -Z,     0, 1 - X];

du = F - eye(3);