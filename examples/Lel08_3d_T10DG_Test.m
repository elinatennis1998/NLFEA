% 05/16/2015
% Tim Truster
% 
% Series of tests for DG elements after merger; shows that CGtoDG routines
% work. T10 mesh.

clear
clc
NCR = 1;

nen = 10;
nel = 10;
% DG insertion type
DGtype = 2;1; % 1=insert everywhere, 2=insert selectively

xl = [1 0 0 0
      2 2 0 0
      4 0 1 0
      3 2 1 0
      5 0 0 3
      6 2 0 3
      8 0 1 3
      7 2 1 3];
[Coordinates,NodesOnElement,RegionOnElement,numnp,numel] = block3d('cart',6,6,6,1,1,1,13,xl,nen);
Coordinates = Coordinates';
NodesOnElement = NodesOnElement';

% Boundary conditions
nodexm = find(abs(Coordinates(:,1)-0)<1e-9); %rollers at x=0
nodexp = find(abs(Coordinates(:,1)-2)<1e-9); %prescribed u_x at x=2
nodeym = find(abs(Coordinates(:,2)-0)<1e-9); %rollers at y=0
nodezm = find(abs(Coordinates(:,3)-0)<1e-9); %rollers at y=0
NodeBC = [nodexm 1*ones(length(nodexm),1) zeros(length(nodexm),1)
          nodexp 1*ones(length(nodexp),1) .1*ones(length(nodexp),1)
          nodeym 2*ones(length(nodeym),1) zeros(length(nodeym),1)
          nodezm 3*ones(length(nodezm),1) zeros(length(nodezm),1)];
numBC = length(NodeBC);

iprob = 0;

nummat = 3;
nen1 = nen + 1;
RegionOnElement(1:6) = 3;
RegionOnElement(13:18) = 2;
% RegionOnElement(31:36) = 2;
MatTypeTable = [1 2 3
                1 1 1
                0 0 0];
AlgoType = [0; 1; 0];
OptFlag = [0 1 1 0 0 1 1 1]';

MateT = [1 1 1]'*[190e3 0.3 1];

if DGtype == 1

% Convert to DG mesh
CGtoDGmesh3
numSI = numCL;
numelCG = numel;
[NodesOnElement,RegionOnElement,nen,numel,nummat,MatTypeTable,MateT] = ...
         FormDG(SurfacesI(:,5:8),NodesOnElement,RegionOnElement,...
         Coordinates,numCL,nen,3,numel,nummat,1,...
                5,0,0,MatTypeTable,MateT);
            
else

% Generate CZM interface: pairs of elements and faces, duplicate the nodes,
% update the connectivities
numnpCG = numnp;
InterTypes = [1 0 0
              1 1 0
              1 1 1]; % only put CZM between the element edges between materials 1-2
DEIProgram3

% Update boundary conditions
NodeBCCG = NodeBC;
numBCCG = numBC;
[NodeBC,numBC] = UpdateNodeSet(maxel,0,RegionOnElement,ElementsOnNodeNum,...
                               ElementsOnNode,ElementsOnNodeDup,NodeBCCG,numBCCG);

% Insert DG couplers
ndm = 3;
InterDGall

end

ProbType = [numnp numel nummat 3 3 nen];

stepmax = 1;
s_del_a = 1;
mults = 1;
datastep = stepmax; % size of output arrays, can be larger than stepmax

itermax = 7;
Residratio = 10^-11;
reststep = 10; % dump data at steps equalting multiples of this #
