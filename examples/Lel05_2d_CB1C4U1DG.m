% 4/29/2015
%
% Sample for: L_Elem1_2dDG
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
btype = 0;
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
[SurfacesL,numSL] = getpressSS(n1,n1,m1,m1,0,0,n1*m1,0);
Load = 2560;

AlgoType = [-1; 1; 0];
OptFlag = [0 1 1 0 0 0 0]';
young = 75000000;
pois = .25;%.4999;%
MateT = [young pois 1];
midloc = (n1+1)*(m1/2+1);
MatTypeTable = [1; 1; 0];

nummat = 1;

% Convert to DG mesh
CGtoDGmesh
% Add interfaces between blocks; must be conforming
% top
% SurfacesI = [SurfacesI; zeros(nu1*nv1,8)];
% for elem = nu1*nv1*(nw1-1)+1:nu1*nv1*nw1
%     numCL = numCL + 1;
%     SurfacesI(numCL,:) = [0 0 0 0 elem elem+nu1*nv1 5 6];
% end
numSI = numCL;
[NodesOnElement,RegionOnElement,nen,numel,nummat,MatTypeTable,MateT] = ...
         FormDG(SurfacesI(:,5:8),NodesOnElement,RegionOnElement,...
         Coordinates,numCL,nen,2,numel,1,1,...
                5,0,0,MatTypeTable,MateT);
ProbType = [numnp numel nummat 2 2 nen];

expliciterr = 1;

% Errors for mesh refinement, pois = .4999;
% Standard 216  -1.6109594e+00  -2.2456703e+00  6
% Standard 912  -1.6429120e+00  -2.2718665e+00  6
% Standard 3744  -1.6994021e+00  -2.3281949e+00  6
% Convergence is NOT optimal for Q4 elements
