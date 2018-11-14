% 05/16/2015
% Tim Truster
%
% Series of tests for DG elements after merger; shows that CGtoDG routines
% work. Q9 mesh.

clear
clc
NCR = 1; % Input file

nen = 9;
nel = 9;
% Mesh with 6x6 tiling
nu = 6;
nv = 6;
% DG insertion type
DGtype = 2;1; % 1=insert everywhere, 2=insert selectively

Coordinates = [1 0 0
             2 2 0
             4 0 1
             3 2 1];
type = 'cart';
rinc = nu;
sinc = nv;
node1 = 1;
elmt1 = 1;
mat = 1;
rskip = 0;
btype = 9;
[Coordinates,NodesOnElement,RegionOnElement,numnp,numel] = block2d(type,rinc,sinc,node1,elmt1,mat,rskip,btype,Coordinates,nen);
Coordinates = Coordinates';
NodesOnElement = NodesOnElement';

% Boundary conditions
nodexm = find(abs(Coordinates(:,1)-0)<1e-9); %rollers at x=0
nodexp = find(abs(Coordinates(:,1)-2)<1e-9); %prescribed u_x at x=2
nodeym = find(abs(Coordinates(:,2)-0)<1e-9); %rollers at y=0
NodeBC = [nodexm 1*ones(length(nodexm),1) zeros(length(nodexm),1)
          nodexp 1*ones(length(nodexp),1) .1*ones(length(nodexp),1)
          nodeym 2*ones(length(nodeym),1) zeros(length(nodeym),1)];
numBC = length(NodeBC);

nen1 = nen + 1;
RegionOnElement(1:2) = 3;
RegionOnElement(3:4) = 2;
RegionOnElement(7:8) = 2;
nummat = 3;
MatTypeTable = [1 2 3; 1 1 1; 0 0 0];
MateT = [1 1 1]'*[190e3 0.3 1];
AlgoType = [0; 1; 0];
OptFlag = [0 1 1 0 0 1 1 1]';

if DGtype == 1
    
% Convert to DG mesh
CGtoDGmesh
numSI = numCL;
[NodesOnElement,RegionOnElement,nen,numel,nummat,MatTypeTable,MateT] = ...
         FormDG(SurfacesI(:,5:8),NodesOnElement,RegionOnElement,...
         Coordinates,numCL,nen,2,numel,nummat,1,...
                5,0,0,MatTypeTable,MateT);
            
else

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

end

ProbType = [numnp numel nummat 2 2 nen];

stepmax = 1;
s_del_a = 1;%/21;stepmax;
mults = (s_del_a:s_del_a:s_del_a*stepmax)';
datastep = stepmax; % size of output arrays, can be larger than stepmax

itermax = 7;
Residratio = 10^-11;
reststep = 10; % dump data at steps equalting multiples of this #
