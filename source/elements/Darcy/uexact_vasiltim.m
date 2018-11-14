function [ue,duex,duey] = uexact_vasiltim(x,y,mu,K,alpha,w,DS)
%
% Exact solution for Darcy-Stokes problem, Vasillev
% 05/16/2013

G = sqrt(mu*K)/alpha;
if alpha == 0
    xi = -1/2;
else
xi = (1-G)/(2*(1+G));
end
chi = (-30*xi - 17)/48;

if y < 0.5 || (nargin == 7 && DS == 0) %Darcy
    
ud = [sin(x/G+w)*exp(y/G) + w*sin(w*x)
     -cos(x/G+w)*exp(y/G)];
udx = [1/G*cos(x/G+w)*exp(y/G) + w^2*cos(w*x)
       1/G*sin(x/G+w)*exp(y/G)];
udy = 1/G*[sin(x/G+w)*exp(y/G)
      -cos(x/G+w)*exp(y/G)];
pd = mu*G/K*cos(x/G+w)*exp(y/G) + mu/K*cos(w*x);
dpd = -mu/K*[sin(x/G+w)*exp(y/G) + w*sin(w*x)
     -cos(x/G+w)*exp(y/G)];

ue = [ud; pd];
duex = [udx; dpd(1)];
duey = [udy; dpd(2)];

else %Stokes
    
us = [sin(x/G+w)*exp(y/G)
     -cos(x/G+w)*exp(y/G)];
usx = 1/G*[cos(x/G+w)*exp(y/G)
       sin(x/G+w)*exp(y/G)];
usy = 1/G*[sin(x/G+w)*exp(y/G)
      -cos(x/G+w)*exp(y/G)];
ps = y - exp(1/(2*G))*cos(w + x/G)*((2*mu)/G - (G*mu)/K) + (2*mu*y*cos(w*x))/K - 1/2;
dps =  [ (exp(1/(2*G))*sin(w + x/G)*((2*mu)/G - (G*mu)/K))/G - (2*mu*w*y*sin(w*x))/K
                                                       (2*mu*cos(w*x))/K + 1];

ue = [us; ps];
duex = [usx; dps(1)];
duey = [usy; dps(2)];

end