% Tim Truster
% 10/13/2013
%
% Generate mesh for Elasticity RFB calculation - input file script
% Called by RFBEfuncE

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
NodeBC = [NodeBC; [[(1:3)'; (minc+3:3+(minc-1)+1+2*(minc-1)-1)'] 2*ones(2*(minc-1)+3,1) zeros(2*(minc-1)+3,1)]];
numBC = 2*(minc*2+1);
numSL = minc;

SurfacesL = zeros(numSL,7);
SurfacesL(1,:) = [4 1 1 1 1 0 0];
for i = 2:numSL-1;
    SurfacesL(i,:) = [3+i 2+i 2*i-1 1 1 0 0];
end
SurfacesL(numSL,:) = [2 2+minc 2*minc-1 1 1 0 0];
         
numCL = 0;
numCn = 0;
numComp = 0;

nummat = 1;
ProbType = [numnp numel nummat 2 2 nen];
MatTypeTable = [1
                1
                0];
AlgoType = [-1; 0; 0];
OptFlag = [0 0 0 0 0 0 0]';

young = Emu*(2*Elam+2*Emu)/(Elam+Emu);
pois = Elam/(2*(Elam+Emu));
MateT = [young pois 1];