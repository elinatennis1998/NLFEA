% Tim Truster
% 06/29/2015
%
% Complex 2d mesh to test CGtoDG; consists of T3 and Q4 regions

% Element refinement numbers
n1 = 8;
m1 = 8;
n2 = 8;
m2 = m1;
n3 = n1 + n2;
m3 = 4;
% Geometry
L1 = 2.0;
H1 = 2.0;
L2 = 2.0;
H2 = H1;
L3 = L1 + L2;
H3 = 1.0;

%% Generate mesh
x = zeros(0,0);
NodesOnElement = zeros(0,0);
RegionOnElement = zeros(0,0);
nen = 4;

% Block 1
xl = [1 0 0
      2 L1 0
      4 0 H1
      3 L1 H1];
type = 'cart';
rinc = n1;
sinc = m1;
node1 = 1;
elmt1 = 1;
mat = 1;
rskip = n2;
btype = 5;
[x,NodesOnElement,RegionOnElement,~,numel] = block2d(type,rinc,sinc,node1,elmt1,mat,rskip,btype,xl,nen,x,NodesOnElement,RegionOnElement);

% Block 2
xl = [1 L1 0
      2 L3 0
      4 L1 H2
      3 L3 H2];
type = 'cart';
rinc = n2;
sinc = m2;
node1 = n1+1;
elmt1 = numel+1;
mat = 1;
rskip = n1;
btype = 0;
[x,NodesOnElement,RegionOnElement,~,numel] = block2d(type,rinc,sinc,node1,elmt1,mat,rskip,btype,xl,nen,x,NodesOnElement,RegionOnElement);

% Block 3
xl = [1 0 H1
      2 L3 H2
      4 0 H1+H3
      3 L3 H2+H3];
type = 'cart';
rinc = n3;
sinc = m3;
node1 = (n3+1)*(m1+1)+1;
elmt1 = numel+1;
mat = 1;
rskip = 0;
btype = 3;
[x,NodesOnElement,RegionOnElement,numnp,numel] = block2d(type,rinc,sinc,node1,elmt1,mat,rskip,btype,xl,nen,x,NodesOnElement,RegionOnElement);

Coordinates = x';
NodesOnElement = NodesOnElement';

%% Boundary conditions
nodexm = find(abs(Coordinates(:,1)-0)<1e-9); %rollers at x=0
nodexp = find(abs(Coordinates(:,1)-L3)<1e-9); %prescribed u_x at x=2
nodeym = find(abs(Coordinates(:,2)-0)<1e-9); %rollers at y=0
NodeBC = [nodexm 1*ones(length(nodexm),1) zeros(length(nodexm),1)
          nodexp 1*ones(length(nodexp),1) .1*ones(length(nodexp),1)
          nodeym 2*ones(length(nodeym),1) zeros(length(nodeym),1)];
numBC = length(NodeBC);

%% Assign materials to particular elements; can be modified
nen1 = nen + 1;
RegionOnElement(1:2) = 3;
RegionOnElement(3:4) = 2;
RegionOnElement(7:8) = 2;
nummat = 3;
MatTypeTable = [1 2 3
                1 1 1
                0 0 0];
MateT = ones(nummat,1)*[100e3 0.25 1];
AlgoType = [0; 1; 0];
OptFlag = [0 1 1 0 0 1 1 1]';

%% Insert interface elements
% Generate CZM interface: pairs of elements and faces, duplicate the nodes,
% update the connectivities
numnpCG = numnp;
InterTypes = [0 0 0
              1 0 0
              1 1 0]; % only put CZM between the element edges between materials 1-2
DEIProgram2

% Update boundary conditions
NodeBCCG = NodeBC;
numBCCG = numBC;
[NodeBC,numBC] = UpdateNodeSet(maxel,0,RegionOnElement,ElementsOnNodeNum,...
                               ElementsOnNode,ElementsOnNodeDup,NodeBCCG,numBCCG);

% Insert DG couplers
ndm = 2;
InterDGall

%% Algorithm definitions for NL_FEA_Program
ProbType = [numnp numel nummat 2 2 nen];

stepmax = 1;
s_del_a = 1;%/21;stepmax;
mults = (s_del_a:s_del_a:s_del_a*stepmax)';
datastep = stepmax; % size of output arrays, can be larger than stepmax

itermax = 7;
Residratio = 10^-11;
reststep = 10; % dump data at steps equalting multiples of this #