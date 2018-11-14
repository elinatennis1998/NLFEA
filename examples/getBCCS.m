function [NodeBC,BCLIndex] = getBCCS(n1,m1,kappa,mu)

% Function to assign BC to surface 5 of lap joint

numBC = 2*(n1+1) + 2*(m1+1);
BCLIndex = [numBC 0];
NodeBC = zeros(numBC,3);
kappamu = kappa/mu;

ind = 1;
nind = 1;
inc = 1;

for i = 1:n1+1
    x = (i-1)/n1;
    psi = -kappamu*2*pi*sin(2*pi*x);
    NodeBC(nind,:) = [ind 2 psi];
    ind = ind + inc;
    nind = nind + 1;
end

ind = 1;
inc = n1 + 1;

for i = 1:m1+1
    y = (i-1)/m1;
    psi = -kappamu*2*pi*sin(2*pi*y);
    NodeBC(nind,:) = [ind 1 psi];
    ind = ind + inc;
    nind = nind + 1;
end

ind = n1 + 1;
inc = n1 + 1;

for i = 1:m1+1
    y = (i-1)/m1;
    psi = -kappamu*2*pi*sin(2*pi*y);
    NodeBC(nind,:) = [ind 1 psi];
    ind = ind + inc;
    nind = nind + 1;
end

ind = m1*(n1+1) + 1;
inc = 1;

for i = 1:n1+1
    x = (i-1)/n1;
    psi = -kappamu*2*pi*sin(2*pi*x);
    NodeBC(nind,:) = [ind 2 psi];
    ind = ind + inc;
    nind = nind + 1;
end