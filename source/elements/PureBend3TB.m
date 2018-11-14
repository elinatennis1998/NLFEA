function [T,B] = PureBend3TB(X,Y,E,v,alpha,L,N)
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

if alpha > 0

% Formula for 1st PK stress and body force
P = ...
[ -(E*cos((X*alpha)/L)*(Y^2*alpha^2 - 2*Y^2*alpha^2*v - 2*L*Y*alpha + 2*L^2*v*log((L - Y*alpha)/L) + 4*L*Y*alpha*v))/(2*L*(L - Y*alpha)*(2*v^2 + v - 1)),  (E*v*sin((X*alpha)/L)*log((L - Y*alpha)/L))/(2*v^2 + v - 1)
  -(E*sin((X*alpha)/L)*(Y^2*alpha^2 - 2*Y^2*alpha^2*v - 2*L*Y*alpha + 2*L^2*v*log((L - Y*alpha)/L) + 4*L*Y*alpha*v))/(2*L*(L - Y*alpha)*(2*v^2 + v - 1)), -(E*v*cos((X*alpha)/L)*log((L - Y*alpha)/L))/(2*v^2 + v - 1)];

B = [...
  (E*alpha*sin((X*alpha)/L)*(2*L^2*v - Y^2*alpha^2 + 2*Y^2*alpha^2*v + 2*L*Y*alpha - 2*L^2*v*log((L - Y*alpha)/L) - 4*L*Y*alpha*v))/(2*L^2*(L - Y*alpha)*(2*v^2 + v - 1))
 -(E*alpha*cos((X*alpha)/L)*(2*L^2*v - Y^2*alpha^2 + 2*Y^2*alpha^2*v + 2*L*Y*alpha - 2*L^2*v*log((L - Y*alpha)/L) - 4*L*Y*alpha*v))/(2*L^2*(L - Y*alpha)*(2*v^2 + v - 1))];
                  
% 
% P = ...
% [ -(E*Y*alpha*cos((X*alpha)/L)*(Y*alpha - 2*L + 2*L*v))/(2*L*(L - Y*alpha)*(2*v^2 + v - 1)), -(E*Y*alpha*v*sin((X*alpha)/L)*(L - Y*alpha))/(L^2*(2*v - 1)*(v + 1))
%   -(E*Y*alpha*sin((X*alpha)/L)*(Y*alpha - 2*L + 2*L*v))/(2*L*(L - Y*alpha)*(2*v^2 + v - 1)),  (E*Y*alpha*v*cos((X*alpha)/L)*(L - Y*alpha))/(L^2*(2*v - 1)*(v + 1))];
% 
% B = [...
%   (E*alpha*sin((X*alpha)/L)*(2*L^2*v - Y^2*alpha^2 + 4*Y^2*alpha^2*v + 2*L*Y*alpha - 8*L*Y*alpha*v))/(2*L^2*(L - Y*alpha)*(2*v^2 + v - 1))
%  -(E*alpha*cos((X*alpha)/L)*(2*L^2*v - Y^2*alpha^2 + 4*Y^2*alpha^2*v + 2*L*Y*alpha - 8*L*Y*alpha*v))/(2*L^2*(L - Y*alpha)*(2*v^2 + v - 1))];
                                                          
% [  (E*Y*alpha*cos((X*alpha)/L)*(L - 2*Y*alpha*v))/(2*L^2*(2*v^2 + v - 1)), -(E*Y*alpha*sin((X*alpha)/L)*(4*L^2*v - L^2 + 2*Y^2*alpha^2*v - 4*L*Y*alpha*v))/(2*L^3*(2*v^2 + v - 1))
%    (E*Y*alpha*sin((X*alpha)/L)*(L - 2*Y*alpha*v))/(2*L^2*(2*v^2 + v - 1)),  (E*Y*alpha*cos((X*alpha)/L)*(4*L^2*v - L^2 + 2*Y^2*alpha^2*v - 4*L*Y*alpha*v))/(2*L^3*(2*v^2 + v - 1))];
% % [ -(E*cos((X*alpha)/L)*(2*L*v*log((L - Y*alpha)/L) - Y*alpha + 2*Y*alpha*v))/(2*L*(2*v^2 + v - 1)), - (E*(sin((X*alpha)/L) + (alpha*sin((X*alpha)/L)*(Y - L/alpha))/L))/(2*v + 2) - (E*alpha*v*sin((X*alpha)/L)*log((L - Y*alpha)/L)*(Y - L/alpha))/(L*(2*v - 1)*(v + 1))
% %   -(E*sin((X*alpha)/L)*(2*L*v*log((L - Y*alpha)/L) - Y*alpha + 2*Y*alpha*v))/(2*L*(2*v^2 + v - 1)),   (E*(cos((X*alpha)/L) + (alpha*cos((X*alpha)/L)*(Y - L/alpha))/L))/(2*v + 2) + (E*alpha*v*cos((X*alpha)/L)*log((L - Y*alpha)/L)*(Y - L/alpha))/(L*(2*v - 1)*(v + 1))];
% 
% B = [...
%  -(E*alpha*sin((X*alpha)/L)*(L - Y*alpha)*(L - 4*L*v + 4*Y*alpha*v))/(2*L^3*(2*v^2 + v - 1))
%   (E*alpha*cos((X*alpha)/L)*(L - Y*alpha)*(L - 4*L*v + 4*Y*alpha*v))/(2*L^3*(2*v^2 + v - 1))];
% %  -(E*alpha*sin((X*alpha)/L)*(L - Y*alpha - 4*L*v + 2*Y*alpha*v))/(2*L^2*(2*v^2 + v - 1))
% %   (E*alpha*cos((X*alpha)/L)*(L - Y*alpha - 4*L*v + 2*Y*alpha*v))/(2*L^2*(2*v^2 + v - 1))

T = P*N;

end