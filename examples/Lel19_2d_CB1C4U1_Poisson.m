% 4/29/2015
%
% Sample for: L_Elem19_2d
% Notes:
%
% Tim Truster
% 10/15/2013
%
% Cantilever Beam
% Conforming Mesh
% Quadrilateral, Uniform, Linear, Level 1
% Pure-Displacement element mesh
% 
% Used with L_Elem19_2d, which is the vectorized Poisson element. The
% entire constitutive tensor is passed to the material subroutine.

NCR = 1;

L = 10;
D = 1;
h = D/2;
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
pois = .25;
lamda = pois*young/((1+pois)*(1-2*pois));
mu = young/(2*(1+pois));
C = lamda*[1 1 0 0]'*[1 1 0 0] + mu*diag([2 2 1 0]);
A_iIjJ = CSFtoA(C,zeros(2,2),eye(2),2);
MateT = reshape(A_iIjJ,1,16);
midloc = (n1+1)*(m1/2+1);
MatTypeTable = [1; 19; 0];

nummat = 1;
numSI = 0;
ProbType = [numnp numel nummat 2 2 nen];
