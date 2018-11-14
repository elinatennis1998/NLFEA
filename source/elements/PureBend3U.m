function [u,du] = PureBend3U(X,Y,alpha,L)
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
-sin((X*alpha)/L)*(Y - L/alpha)
 L/alpha + cos((X*alpha)/L)*(Y - L/alpha)];
u = phi - [X;Y];
F = ...
[ -(alpha*cos((X*alpha)/L)*(Y - L/alpha))/L, -sin((X*alpha)/L)
  -(alpha*sin((X*alpha)/L)*(Y - L/alpha))/L,  cos((X*alpha)/L)];
du = F - eye(2);