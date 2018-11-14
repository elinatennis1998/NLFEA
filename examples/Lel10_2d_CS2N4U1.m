% 4/29/2015
%
% Sample for: L_Elem4_2dDGb
% Notes: Analysis runs without error messages, but computed error norms
% don't match with older values; likely needs debugged or updated to new
% code standard.
%
% Tim Truster
% 09/07/2012
%
% Convergence Smooth problem, Darcy Flow
% Conforming Mesh
% Quadrilateral, Uniform, Linear, Level 1
% Use with Inter_FEA_Program

L = 1;
D = 1;
h = D/20;
p = 1;
wf = 0;
bs = 0;
n = D/h;
pu = p;
pv = p;
n1 = L/2*1*n;
m1 = D/2*1*n;
n2 = L/2*n;
m2 = D/2*n;
n3 = L/2*n;
m3 = D/2*n;
n4 = L/2*1*n;
m4 = D/2*1*n;
nel = 4;
iprob = 1;

x = zeros(0,0);
NodesOnElement = zeros(0,0);
RegionOnElement = zeros(0,0);
nen = 4;

% Block 1
xl = [1 0 0
      2 L/2 0
      4 0 D/2
      3 L/2 D/2];
type = 'cart';
rinc = n1;
sinc = m1;
node1 = 1;
elmt1 = 1;
mat = 1;
rskip = 0;
btype = 0;
[x,NodesOnElement,RegionOnElement,numnp1,numel1] = block2d(type,rinc,sinc,node1,elmt1,mat,rskip,btype,xl,nen,x,NodesOnElement,RegionOnElement);

% Block 2
xl = [1 L/2 0
      2 L 0
      4 L/2 D/2
      3 L D/2];
type = 'cart';
rinc = n2;
sinc = m2;
node1 = numnp1+1;
elmt1 = numel1+1;
mat = 2;
rskip = 0;
btype = 0;
[x,NodesOnElement,RegionOnElement,numnp2,numel2] = block2d(type,rinc,sinc,node1,elmt1,mat,rskip,btype,xl,nen,x,NodesOnElement,RegionOnElement);

% Block 3
xl = [1 0 D/2
      2 L/2 D/2
      4 0 D
      3 L/2 D];
type = 'cart';
rinc = n3;
sinc = m3;
node1 = numnp2+1;
elmt1 = numel2+1;
mat = 3;
rskip = 0;
btype = 0;
[x,NodesOnElement,RegionOnElement,numnp3,numel3] = block2d(type,rinc,sinc,node1,elmt1,mat,rskip,btype,xl,nen,x,NodesOnElement,RegionOnElement);

% Block 4
xl = [1 L/2 D/2
      2 L D/2
      4 L/2 D
      3 L D];
type = 'cart';
rinc = n4;
sinc = m4;
node1 = numnp3+1;
elmt1 = numel3+1;
mat = 4;
rskip = 0;
btype = 0;
[x,NodesOnElement,RegionOnElement,numnp,numel] = block2d(type,rinc,sinc,node1,elmt1,mat,rskip,btype,xl,nen,x,NodesOnElement,RegionOnElement);

Coordinates = x';
NodesOnElement = NodesOnElement';

numComp = 0;
numCn = 0;
numSF = 0;
numSP = 0;

% ProbType = [4 3 2];
FaceBC = 0;
numP = 1;
PatchMatTable = 1;
% kappa = 1/1;
kappa = 1;
mu = 1;
gc = 1;
rho = 1;
% MateT = [kappa mu gc rho
%          .01 mu gc rho
%          .01 mu gc rho
%          kappa mu gc rho];
MateT = [kappa mu gc rho
         kappa mu gc rho
         kappa mu gc rho
         kappa mu gc rho];
NodeLoad = 0;

nodeendR = (n1+1)*(m1+1);
nodeAR = (n1+1)*(m1)+1;
nodeAL = numnp2+1;
elemR = n1*(m1-1)+1;
elemL = numel2+1;
nodeincR = 1;
nodeincL = 1;
elemincR = 1;
elemincL = 1;
edgeR = 3;
edgeL = 1;
[SurfacesI13,numCL13] = getInterHoriz(Coordinates,nodeendR,nodeAR,nodeAL, ...
                             elemR,elemL,nodeincR,nodeincL,elemincR,elemincL,edgeR,edgeL);

nodeendR = numnp1+(n2+1)*m2+1;
nodeAR = numnp1+1;
nodeAL = (n1+1);
elemR = numel1+1;
elemL = n1;
nodeincR = n2+1;
nodeincL = n1+1;
elemincR = n2;
elemincL = n1;
edgeR = 4;
edgeL = 2;
[SurfacesI12,numCL12] = getInterVert(Coordinates,nodeendR,nodeAR,nodeAL, ...
                             elemR,elemL,nodeincR,nodeincL,elemincR,elemincL,edgeR,edgeL,numnp);

nodeendR = numnp1+(n2+1)*(m2+1);
nodeAR = numnp1+(n2+1)*(m2)+1;
nodeAL = numnp3+1;
elemR = numel1+n2*(m2-1)+1;
elemL = numel3+1;
nodeincR = 1;
nodeincL = 1;
elemincR = 1;
elemincL = 1;
edgeR = 3;
edgeL = 1;
[SurfacesI24,numCL24] = getInterHoriz(Coordinates,nodeendR,nodeAR,nodeAL, ...
                             elemR,elemL,nodeincR,nodeincL,elemincR,elemincL,edgeR,edgeL);

nodeendR = numnp3+(n4+1)*m4+1;
nodeAR = numnp3+1;
nodeAL = numnp2+(n3+1);
elemR = numel3+1;
elemL = numel2+n3;
nodeincR = n4+1;
nodeincL = n3+1;
elemincR = n4;
elemincL = n3;
edgeR = 4;
edgeL = 2;
[SurfacesI34,numCL34] = getInterVert(Coordinates,nodeendR,nodeAR,nodeAL, ...
                             elemR,elemL,nodeincR,nodeincL,elemincR,elemincL,edgeR,edgeL,numnp);
SurfacesI = [SurfacesI12
             SurfacesI13
             SurfacesI34
             SurfacesI24];
numCL = numCL13 + numCL12 + numCL34 + numCL24;
numSI = numCL;

% Boundary Conditions
[NodeBC,BCLIndex] = getBCCSN(n1,m1,n2,m2,n3,m3,n4,m4,Coordinates,kappa,mu);
NodeBC = [NodeBC; numnp1+(n2+1) 3 0];
BCLIndex(1) = BCLIndex(1) + 1;
numBC = BCLIndex(1);
MatTypeTable = [1 2 3 4
                7 7 7 7
                0 0 0 0 ];
     
% Surface Tractions
SurfacesL = 0;
numSL = 0;

[NodesOnElement,RegionOnElement,nen,numel,nummat,MatTypeTable,MateT] = ...
         FormDG(SurfacesI(:,5:8),NodesOnElement,RegionOnElement,...
         Coordinates,numCL,nen,2,numel,4,4,...
                10,0,0,MatTypeTable,MateT);

numCL = 0;
AlgoType = [-1; 1; 0];
OptFlag = [0 0 0 0 0 0 0]';
ProbType = [numnp numel nummat 2 3 nen];

expliciterr = 1;
