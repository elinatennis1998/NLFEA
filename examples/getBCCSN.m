function [NodeBC,BCLIndex] = getBCCSN(n1,m1,n2,m2,n3,m3,n4,m4,Coordinates,kappa,mu)

% Function to assign BC to surface 5 of lap joint

numBC = (n1+1) + (m1+1) + (n2+1) + (m2+1) + (n3+1) + (m3+1) + (n4+1) + (m4+1);
BCLIndex = [numBC 0];
NodeBC = zeros(numBC,3);
kappamu = kappa/mu;
numnp1 = (n1+1)*(m1+1);
numnp2 = (n2+1)*(m2+1) + numnp1;
numnp3 = (n3+1)*(m3+1) + numnp2;

ind = 1;
nind = 1;
inc = 1;

for i = 1:n1+1
    x = Coordinates(i,1);
    psi = -kappamu*2*pi*sin(2*pi*x);
    NodeBC(nind,:) = [ind 2 psi];
    ind = ind + inc;
    nind = nind + 1;
end

ind = numnp1 + 1;
inc = 1;

for i = 1:n2+1
    x = Coordinates(numnp1+i,1);
    psi = -kappamu*2*pi*sin(2*pi*x);
    NodeBC(nind,:) = [ind 2 psi];
    ind = ind + inc;
    nind = nind + 1;
end

ind = 1;
inc = n1 + 1;

for i = 1:m1+1
    y = Coordinates((n1+1)*(i-1)+1,2);
    psi = -kappamu*2*pi*sin(2*pi*y);
    NodeBC(nind,:) = [ind 1 psi];
    ind = ind + inc;
    nind = nind + 1;
end

ind = numnp2 + 1;
inc = n3 + 1;

for i = 1:m3+1
    y = Coordinates(numnp2+(n3+1)*(i-1)+1,2);
    psi = -kappamu*2*pi*sin(2*pi*y);
    NodeBC(nind,:) = [ind 1 psi];
    ind = ind + inc;
    nind = nind + 1;
end

ind = numnp1 + n2 + 1;
inc = n2 + 1;

for i = 1:m2+1
    y = Coordinates(numnp1+(n2+1)*i,2);
    psi = -kappamu*2*pi*sin(2*pi*y);
    NodeBC(nind,:) = [ind 1 psi];
    ind = ind + inc;
    nind = nind + 1;
end

ind = numnp3 + n4 + 1;
inc = n4 + 1;

for i = 1:m4+1
    y = Coordinates(numnp3+(n4+1)*i,2);
    psi = -kappamu*2*pi*sin(2*pi*y);
    NodeBC(nind,:) = [ind 1 psi];
    ind = ind + inc;
    nind = nind + 1;
end

ind = numnp2 + m3*(n3+1) + 1;
inc = 1;

for i = 1:n3+1
    x = Coordinates(numnp2+m3*(n3+1) + i,1);
    psi = -kappamu*2*pi*sin(2*pi*x);
    NodeBC(nind,:) = [ind 2 psi];
    ind = ind + inc;
    nind = nind + 1;
end

ind = numnp3 + m4*(n4+1) + 1;
inc = 1;

for i = 1:n4+1
    x = Coordinates(numnp3+m4*(n4+1) + i,1);
    psi = -kappamu*2*pi*sin(2*pi*x);
    NodeBC(nind,:) = [ind 2 psi];
    ind = ind + inc;
    nind = nind + 1;
end