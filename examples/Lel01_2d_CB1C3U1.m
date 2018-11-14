% 4/29/2015
%
% Sample for: L_Elem1_2d
% Notes: shows standard L2 error convergence
%
% Tim Truster
% 12/28/2010
%
% Cantilever Beam
% Conforming Mesh
% Quadrilateral, Uniform, Linear, Level 1
% Pure-Displacement element mesh

NCR = 1;

L = 10;
D = 1;
h = D/4;
p = 1;
wf = 0;
bs = 0;
rho = 1000;
n = D/h;
pu = p;
pv = p;
nu = 2*n*L/2*pu;
nv = 2*n*pv;
n1 = 5*n;
m1 = 1*n;
iprob = 6;

x = zeros(0,0);
NodesOnElement = zeros(0,0);
RegionOnElement = zeros(0,0);
nen = 4;

% Block 1
xl = [1 0 -D
      2 L -D
      4 0 D
      3 L D];
type = 'cart';
rinc = n1;
sinc = m1;
node1 = 1;
elmt1 = 1;
mat = 1;
rskip = 0;
% btype = 5;
btype = 2;
% btype = 4;
[x,NodesOnElement,RegionOnElement,numnp,numel] = block2d(type,rinc,sinc,node1,elmt1,mat,rskip,btype,xl,nen,x,NodesOnElement,RegionOnElement);

Coordinates = x';
NodesOnElement = NodesOnElement';
nen1 = nen + 1;

% Boundary Conditions
NodeBC = [1 1 0
          (n1+1)*m1/2+1 1 0
          (n1+1)*m1/2+1 2 0
          (n1+1)*m1+1 1 0];
numBC = 4;
     
% Surface Tractions
if btype == 5 || btype == 4
    [SurfacesL,numSL] = getpressSST(n1,n1,m1,m1,0,0,2*n1*m1,0);
elseif btype == 2
    [SurfacesL,numSL] = getpressSS3(n1,n1,m1,m1,0,0,2*n1*m1,0);
end
Load = 2560;

AlgoType = [-1; 1; 0];
OptFlag = [0 1 1 0 0 0 0]';
young = 75000000;
pois = .25;%.4999;%
MateT = [young pois 1];
midloc = (n1+1)*(m1/2+1);
MatTypeTable = [1; 1; 0];

nummat = 1;
ProbType = [numnp numel nummat 2 2 nen];

expliciterr = 1;

% Errors for mesh refinement, pois = .4999;
% Standard 160  -1.6029465e+00  -2.2438190e+00  6
% Standard 640  -1.6076338e+00  -2.2486840e+00  6
% Standard 2560  -1.6279534e+00  -2.2690228e+00  6
% Clearly locks
