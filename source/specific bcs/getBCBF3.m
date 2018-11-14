function [NodeBC,BCLIndex] = getBCBF3(n1,m1)

% Function to assign BC to surface 5 of lap joint

numBC = 2*(n1 + 1);
BCLIndex = [numBC 0];
NodeBC = zeros(numBC,3);

ind = 1;
nind = 1;

for i = 1:n1+1
    NodeBC(nind,:) = [ind 1 0];
    nind = nind + 1;
    NodeBC(nind,:) = [ind 2 0];
    ind = ind + 1;
    nind = nind + 1;
end