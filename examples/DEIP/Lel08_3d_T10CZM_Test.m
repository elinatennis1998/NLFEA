% Patch test for prismatic domain of T10 elements with CZ couplers along
% interfaces between regions.
% Domain: 2x1x3 rectangle block
% Loading: Prescribed displacement of 0.1 on top edge.
%
% Last revision: 03/08/2017 TJT

clear
NCR = 1;
% clc

nen = 10;
nel = 10;

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
nodezp = find(abs(Coordinates(:,3)-3)<1e-9); %prescribed u_y at y=2
nodeym = find(abs(Coordinates(:,2)-0)<1e-9); %rollers at y=0
nodezm = find(abs(Coordinates(:,3)-0)<1e-9); %rollers at z=0
NodeBC = [nodexm 1*ones(length(nodexm),1) zeros(length(nodexm),1)
%           nodexp 1*ones(length(nodexp),1) .1*ones(length(nodexp),1)
          nodeym 2*ones(length(nodeym),1) zeros(length(nodeym),1)
          nodezp 3*ones(length(nodezp),1) .1*ones(length(nodezp),1)
          nodezm 3*ones(length(nodezm),1) zeros(length(nodezm),1)];
numBC = length(NodeBC);

nummat = 3;
RegionOnElement(1:54) = 3;
RegionOnElement(55:108) = 2;
% RegionOnElement(31:36) = 2;
MatTypeTable = [1 2 3
                1 1 1
                0 0 0];
AlgoType = [0; 1; 0];
% Output quantity flags
DHist = 1;
FHist = 1;
SHist = 1;
SEHist = 1;

MateT = [1 1 1]'*[190e3 0.3 1];

% Generate CZM interface: pairs of elements and faces, duplicate the nodes,
% update the connectivities
numnpCG = numnp;
InterTypes = [0 0 0
              1 0 0
              1 1 0]; % only put CZM between the element edges between materials 1-2
DEIProgram3

% Update boundary conditions
NodeBCCG = NodeBC;
numBCCG = numBC;
[NodeBC,numBC] = UpdateNodeSet(maxel,0,RegionOnElement,ElementsOnNodeNum,...
                               ElementsOnNode,ElementsOnNodeDup,NodeBCCG,numBCCG);

numSI = numCL;
numelCG = numel;
nen_bulk = nen;
SurfacesI = zeros(0,8);
numSI = 0;
for mat2 = 1:nummat
    for mat1 = 1:mat2
        
        matI = mat2*(mat2-1)/2 + mat1; % ID for material pair (row=mat2, col=mat1)
        if InterTypes(mat2,mat1) > 0
        numSIi = numEonF(matI);
        locF = FacetsOnInterfaceNum(matI):(FacetsOnInterfaceNum(matI+1)-1);
        facs = FacetsOnInterface(locF);
        SurfacesIi = ElementsOnFacet(facs,:);
        SurfacesI(numSI+1:numSI+numSIi,5:8) = SurfacesIi;
        numSI = numSI + numSIi;
        [NodesOnElement,RegionOnElement,nen,numel,nummat,MatTypeTable,MateT] = ...
         FormCZ(SurfacesIi,NodesOnElement,RegionOnElement,Coordinates,numSIi,nen_bulk,3,numel,nummat,6, ...
                22,0,[5000000],MatTypeTable,MateT);
        end
    end
end

ProbType = [numnp numel nummat 3 3 nen];

% plotNodeCont3(Coordinates+Node_U_V*100, Node_U_V(:,1), NodesOnElement, 2,1:54*3)
