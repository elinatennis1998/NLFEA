% 4/29/2015
%
% Sample for: L_Elem4_2dVMSb
% Notes: Analysis runs without error messages, but computed error norms
% don't match with older values; likely needs debugged or updated to new
% code standard.
%
% Tim Truster
% 06/06/2011
%
% Convergence Study, Darcy Flow
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
n1 = L*n;
m1 = D*n;
nel = 4;
iprob = 1;

x = zeros(0,0);
NodesOnElement = zeros(0,0);
RegionOnElement = zeros(0,0);
nen = 4;

% Block 1
xl = [1 0 0
      2 L 0
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

numComp = 0;
numSI = 0;
numCn = 0;
numSF = 0;
numSP = 0;

% ProbType = [4 3 2];
FaceBC = 0;
numP = 1;
PatchMatTable = 1;
kappa = 1;
mu = 1;
gc = 1;
rho = 1;
MateT = [kappa mu gc rho];
NodeLoad = 0;

% Boundary Conditions
[NodeBC,BCLIndex] = getBCCS(n1,m1,kappa,mu);
NodeBC = [NodeBC; (n1+1)*n1/2+n1/2+1 3 0];
BCLIndex(1) = BCLIndex(1) + 1;
numBC = BCLIndex(1);
     
% Surface Tractions
SurfacesL = 0;
numSL = 0;

MatTypeTable = [1
                7
                0];
numCL = 0;
AlgoType = [-1; 1; 0];
OptFlag = [0 0 0 0 0 0 0]';
ProbType = [numnp numel 1 2 3 nen];
nen1 = nen + 1;

expliciterr = 1;

% NL_FEA_Program
% 100  -8.3285538e-001  -1.6085812e+000  7.5629620e-001  -1.8766982e-001  1
% 100  -Inf  -Inf 
% NL_FEA_Program
% 400  -1.4461278e+000  -2.2113668e+000  4.5250984e-001  -4.9551583e-001  1
% 400  -Inf  -Inf 
% NL_FEA_Program
% 1600  -2.0530434e+000  -2.8135822e+000  1.5051010e-001  -7.9824061e-001  1
% 1600  -Inf  -Inf 
