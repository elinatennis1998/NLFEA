function [ue,duex,duey] = uexact_BF(x,y,lam)
% 04/11/2012
% Exact solution for body force problem

ue = zeros(3,1);
duex = zeros(3,1);
duey = zeros(3,1);

ue(2) = 1.01*x*y;
ue(3) = 1.01*x*lam;
duex(2) = 1.01*y;
duex(3) = 1.01*lam;
duey(2) = 1.01*x;