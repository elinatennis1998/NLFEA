function [F,J] = unl(k,rp,beta2,tn,tt,x)

len = length(x);
one = ones(len,1);
F = k^2*one - tn^2*one./(rp*one + x).^2 - beta2*tt^2*one./(rp*one + beta2*x).^2;
J = 2*tn^2*one./(rp*one + x).^3 + 2*beta2^2*tt^2*one./(rp*one + beta2*x).^3;