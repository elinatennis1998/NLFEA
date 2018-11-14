% Tim Truster
% 10/13/2013
%
% Generate mesh for Poisson RFB calculation - input file script
% Called by RFBEfuncP

Coordinates = ECoordinates';
numnp = 3;
x = Coordinates';
NodesOnElement = [1 2 3 1]';
numel = 1;
nen = 3;
nen1 = nen + 1;
iprob = 0;

% minc = 16;8;4;
celn = (2*minc+1)^2;
ndm = 2;
ncel = 2;
cel = minc^2;
epnum = [1 1 1];
epatch = [1 1 1];
gen_submesh

numel = numels;
numnp = numnps;
Coordinates = xs';
NodesOnElement = [ixs(1:3,1:numel)'];
RegionOnElement = ones(numel,1);

NodeBC = [[(1:3)'; (minc+3:3+(minc-1)+1+2*(minc-1)-1)'] ones(2*(minc-1)+3,1) zeros(2*(minc-1)+3,1)];
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
SurfacesL(1,:) = [4 1 1 1 tu3];
for i = 2:numSL-1;
    SurfacesL(i,:) = [3+i 2+i 2*i-1 1 tu3];
end
SurfacesL(numSL,:) = [2 2+minc 2*minc-1 1 tu3];
         
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