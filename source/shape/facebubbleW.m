function [b,db] = facebubbleW(ss)
% Tim Truster
% 10/16/2013
% face bubble function for stabilized DG 3D wedge sector

db = zeros(3,1);

r = ss(1);
s = ss(2);
t = ss(3);
twnty7 = 27.d0;
one = 1.d0;
onemt = one - t;
u = 1 - r - s;
two = 2.d0;
    
b = twnty7*r*s*u*onemt/two;
db(1) = twnty7*(s*u - r*s)*onemt/two;
db(2) = twnty7*(r*u - r*s)*onemt/two;
db(3) =-twnty7*r*s*u/two;