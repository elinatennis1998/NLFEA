% Tim Truster
% 10/13/2013
%
% Generate mesh for Poisson RFB calculation - input file script
% Called by RFBEfuncP
% Uses triblock routine to make graded mesh toward segment edge

Coordinates = ECoordinates;
nen = 3;
nen1 = nen + 1;
iprob = 0;
shl = [0 0.74 0.26]'; % more graded
% shl = [0 0.7 0.3]'; % graded
% shl = [0 0.5 0.5]'; % uniform
x5 = Coordinates*shl;
shl = [0.74 0 0.26]'; % more graded
% shl = [0.7 0 0.3]'; % graded
% shl = [0.5 0 0.5]'; % uniform
x6 = Coordinates*shl;
xl = [1 2 3 5 6
      Coordinates [x5 x6]]';

[x,NodesOnElement,RegionOnElement,numnp,numel] = triblock('cart',minc,1,1,1,1,xl,6);

Coordinates = x';
NodesOnElement = [NodesOnElement(1:3,1:numel)'];
RegionOnElement = ones(numel,1);

NodeBC = zeros(minc*2+1,3);
nind = 1;
ind = nind;
for i = 1:minc
    NodeBC(nind,:) = [ind 1 0];
    nind = nind + 1;
    ind = ind + minc + 1 - i;
    NodeBC(nind,:) = [ind 1 0];
    nind = nind + 1;
    ind = ind + 1;
end
NodeBC(nind,:) = [ind 1 0];
numBC = (minc*2+1);
numSL = minc;

% t1 = ECoordinates(:,1)-ECoordinates(:,2);
% t1 = [t1; 0];
% [tm1, tu1] = VecNormalize(t1);
% t2 = [0; 0; 1];
% tm2 = 1;
% tu2 = t2';
% t3 = VecCrossProd(t1,t2);
% [tm3, tu3] = VecNormalize(t3);
tu3 = [1 0 0];

SurfacesL = zeros(numSL,7);
SurfacesL(1,:) = [2 1 1 1 tu3];
for i = 2:numSL;
    SurfacesL(i,:) = [1+i i 2*i-1 1 tu3];
end
         
numCL = 0;
numCn = 0;
numComp = 0;

nummat = 1;
ProbType = [numnp numel nummat 2 1 nen];
MatTypeTable = [1
                2
                0];
AlgoType = [-1; 0; 0];
OptFlag = [0 0 0 0 0 0 0]';
% MateT = [1 1 0
%          1 1 0];
MateT = [EA(1) EA(2) 1/2*EA(3)];