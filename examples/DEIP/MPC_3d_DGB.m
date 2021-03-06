% Axial test of MP Periodic Boundary Conditions with CZ couplers
% Domain: 4x4x4 rectangle block
% Loading: Prescribed axial displacement.
%
% Last revision: 03/31/2017 TJT

clear
NCR = 1;
% clc

nen = 8;
nel = 8;

numc = 4;2;1;
xl = [1 0 0 0
      2 numc 0 0
      4 0 numc 0
      3 numc numc 0
      5 0 0 numc
      6 numc 0 numc
      8 0 numc numc
      7 numc numc numc];
[Coordinates,NodesOnElement,~,numnp,numel] = block3d('cart',numc,numc,numc,1,1,1,10,xl,nen);
Coordinates = Coordinates';
NodesOnElement = NodesOnElement';

% Modify regions
RegionOnElement = [1 1 2 2 1 1 2 2 3 3 4 4 3 3 4 4 ...
                   1 1 2 2 1 1 2 2 3 3 4 4 3 3 4 4 ...
                   5 5 6 6 5 5 6 6 7 7 8 8 7 7 8 8 ...
                   5 5 6 6 5 5 6 6 7 7 8 8 7 7 8 8 ]';
% RegionOnElement = [5 6 6 5 5 6 6 5 7 8 8 7 7 8 8 7 ...
%                    1 2 2 1 1 2 2 1 3 4 4 3 3 4 4 3 ...
%                    1 2 2 1 1 2 2 1 3 4 4 3 3 4 4 3 ...
%                    5 6 6 5 5 6 6 5 7 8 8 7 7 8 8 7 ...
%                    ]';
% RegionOnElement = [7 8 8 7 5 6 6 5 5 6 6 5 7 8 8 7 ...
%                    3 4 4 3 1 2 2 1 1 2 2 1 3 4 4 3 ...
%                    3 4 4 3 1 2 2 1 1 2 2 1 3 4 4 3 ...
%                    7 8 8 7 5 6 6 5 5 6 6 5 7 8 8 7 ...
%                    ]';
ix = randperm(numel);
NodesOnElement = NodesOnElement(ix,:);
RegionOnElement = RegionOnElement(ix);
nummat = 8;
MatTypeTable = [1 2 3 4 5 6 7 8
                1 1 1 1 1 1 1 1];
MateT = [100 0.25 1
         100 0.25 1
         100 0.25 1
         100 0.25 1
         100 0.25 1
         100 0.25 1
         100 0.25 1
         100 0.25 1];
AlgoType = [0; 1; 0];
     
% Boundary conditions
FIX = 1;
FRONT = (numnp-(numc+1)*(numc+1)+1); %z
RIGHT = (numc+1); %x
TOP = (numc+1)*(numc)+1; %y
TOPRIGHT = TOP + RIGHT - 1; %xy
FRONTRIGHT = FRONT + RIGHT - 1; %xz
TOPFRONT = FRONT + TOP - 1; %yz
TOPFRONTRIGHT = numnp; %yz
NodeBC = [FIX 1 0
          FIX 2 0
          FIX 3 0
%           FRONT 1 0
%           FRONT 2 0
%           FRONT 3 0.01
%           RIGHT 1 0.01
%           RIGHT 2 0
%           RIGHT 3 0
%           TOP 1 0
%           TOP 2 0.01
%           TOP 3 0
          ];
numBC = size(NodeBC,1);

% Create periodic boundary conditions
dx = 1;
dy = numc+1;
dz = dy*dy;
dxy = TOP - 1;
PBCListx = [(RIGHT:dy:TOPFRONTRIGHT)' (FIX:dy:TOPFRONT)' RIGHT*ones(dz,1) FIX*ones(dz,1)]; %yz
% remove repeats
RIGHTface = 2:dz;
PBCListx = PBCListx(RIGHTface,:);
PBCListy = zeros(dz,4); %xz
ind = 0;
TOP2 = TOP - 1;
FIX2 = 0;
for j = 1:dy
    for i = 1:dy
        TOP2 = TOP2 + 1;
        FIX2 = FIX2 + 1;
        ind = ind + 1;
        PBCListy(ind,:) = [TOP2 FIX2 TOP FIX];
    end
    TOP2 = TOP2 + TOP - 1;
    FIX2 = FIX2 + TOP - 1;
end
% remove repeats
TOPface = 1:dz;
TOPface = setdiff(TOPface,(dy:dy:dz));
TOPface = setdiff(TOPface,1);
PBCListy = PBCListy(TOPface,:);
PBCListz = [(FRONT:dx:TOPFRONTRIGHT)' (FIX:dx:TOPRIGHT)' FRONT*ones(dz,1) FIX*ones(dz,1)]; %xy
PBCListz = PBCListz(1:dz-dy,:);
% remove repeats
FRONTface = 1:dz;
FRONTface = setdiff(FRONTface,(dz-dy:dz));
FRONTface = setdiff(FRONTface,(dy:dy:dz-dy));
FRONTface = setdiff(FRONTface,1);
PBCListz = PBCListz(FRONTface,:);
PBCList = [PBCListx; PBCListy; PBCListz];
% PBCList = [PBCListx];
MPCListx = [PBCListx(:,1:2) ones(length(PBCListx),1)*[1 0 0]
            PBCListx(1,3:4) [1 0 0]];
MPCListy = [PBCListy(:,1:2) ones(length(PBCListy),1)*[0 1 0]
            PBCListy(1,3:4) [0 1 0]];
MPCListz = [PBCListz(:,1:2) ones(length(PBCListz),1)*[0 0 1]
            PBCListz(1,3:4) [0 0 1]];
MPCList = [MPCListx; MPCListy; MPCListz];


% Generate CZM interface: pairs of elements and faces, duplicate the nodes,
% update the connectivities
numnpCG = numnp;
usePBC = 2; % flag to turn on keeping PBC
InterTypes = [0 0 0 0 0 0 0 0
              1 0 0 0 0 0 0 0
              1 1 0 0 0 0 0 0
              1 1 1 0 0 0 0 0
              1 1 1 1 0 0 0 0
              1 1 1 1 1 0 0 0
              1 1 1 1 1 1 0 0
              1 1 1 1 1 1 1 0]; % only put CZM between the element edges between materials 1-2
DEIProgram3


% Insert couplers; use extended DG mesh with ghost elements so that
% couplers appear properly on the domain boundary
ndm = 3;
% CZ element stiffness
CZprop = 50000;50;1;
% MPC spring element properties
pencoeff = 1e9;
CornerXYZ = [0.000000 0.000000 0.000000
             1.000000 0.000000 0.000000
             0.000000 1.000000 0.000000];

nummatCG = nummat;
numSI = numCL;
numelCG = numel;
nen_bulk = nen;
SurfacesI = zeros(0,8);
numSI = 0;
% Arrays for new MPC links being formed
MPCListNew = zeros(0,2+ndm);
numMPCnew = 0;
CouplerNodes = zeros(0,1); % extra nodes for MPC-CZ
NodesOnLinkNew = zeros(4,numnp);
NodesOnLinknumNew = zeros(numnp,1);
% Interface region information: [materialID in MateT; mat1; mat2; regI; first coupler number; last coupler number]
RegionsOnInterface = zeros(nummat*(nummat+1)/2,6);

for mat2 = 1:nummat
    for mat1 = 1:mat2
        
        matI = GetRegionsInt(mat1,mat2);
        numSIi = numEonPBC(matI);
        if numSIi > 0 % Create PBC pairs first
            if InterTypes(mat2,mat1) > 0% then make PBC-CZ couplers
            locF = FacetsOnPBCNum(matI):(FacetsOnPBCNum(matI+1)-1);
            facs = FacetsOnPBC(locF);
            SurfacesIi = ElementsOnFacet(facs,:);
            SurfacesIi = ReverseFacets(SurfacesIi,NodesOnElement,Coordinates,numSIi,ndm);
            ElementsOnFacet(facs,:) = SurfacesIi;
            SurfacesI(numSI+1:numSI+numSIi,5:8) = SurfacesIi;
            numel_old = numel;
            [NodesOnElement,RegionOnElement,~,numnp,numel,nummat,MatTypeTable,MateT,Coordinates,...
                MPCListNew,numMPCnew,NodesOnLinkNew,NodesOnLinknumNew,CouplerNodes] = ...
             FormMPCCZ(MPCListNew,numMPCnew,NodesOnLinkNew,NodesOnLinknumNew,CouplerNodes,...
                    SurfacesIi,NodesOnElement,NodesOnElementCG,RegionOnElement,Coordinates,numSIi,nen_bulk,ndm,numnp,numel,nummat,6, ...
                    22,0,CZprop,MPCList,NodesOnLink,NodesOnLinknum,MatTypeTable,MateT);
            if numel > numel_old
            RegionsOnInterface(nummat-nummatCG,:) = [nummat mat1 mat2 matI numel_old+1 numel];
            end
            % Then form regular CZ couplers
            numSIi = numEonF(matI) - numSIi;
            locF = FacetsOnIntMinusPBCNum(matI):(FacetsOnIntMinusPBCNum(matI+1)-1);
            facs = FacetsOnIntMinusPBC(locF);
            SurfacesIi = ElementsOnFacet(facs,:);
            SurfacesI(numSI+1:numSI+numSIi,5:8) = SurfacesIi;
            numSI = numSI + numSIi;
            numel_old = numel;
            [NodesOnElement,RegionOnElement,nen,numel,nummat,MatTypeTable,MateT] = ...
             FormCZ(SurfacesIi,NodesOnElement,RegionOnElement,Coordinates,numSIi,nen_bulk,ndm,numel,nummat,6, ...
                    22,0,CZprop,MatTypeTable,MateT);
            if numel > numel_old
            RegionsOnInterface(nummat-nummatCG,:) = [nummat mat1 mat2 matI numel_old+1 numel];
            end
            else % just form PBC pairs
            locF = FacetsOnPBCNum(matI):(FacetsOnPBCNum(matI+1)-1);
            facs = FacetsOnPBC(locF);
            SurfacesIi = ElementsOnFacet(facs,:);
            SurfacesIi = ReverseFacets(SurfacesIi,NodesOnElement,Coordinates,numSIi,ndm);
            ElementsOnFacet(facs,:) = SurfacesIi;
            SurfacesI(numSI+1:numSI+numSIi,5:8) = SurfacesIi;
            numel_old = numel;
            [NodesOnElement,RegionOnElement,~,numnp,numel,nummat,MatTypeTable,MateT,Coordinates,...
                MPCListNew,numMPCnew,NodesOnLinkNew,NodesOnLinknumNew,CouplerNodes] = ...
             FormMPCCZ(MPCListNew,numMPCnew,NodesOnLinkNew,NodesOnLinknumNew,CouplerNodes,...
                    SurfacesIi,NodesOnElement,NodesOnElementCG,RegionOnElement,Coordinates,numSIi,nen_bulk,ndm,numnp,numel,nummat,6, ...
                    22,0,CZprop,MPCList,NodesOnLink,NodesOnLinknum,MatTypeTable,MateT);
            if numel > numel_old
            RegionsOnInterface(nummat-nummatCG,:) = [nummat mat1 mat2 matI numel_old+1 numel];
            end
            end
        elseif InterTypes(mat2,mat1) > 0 % just make CZ couplers
            numSIi = numEonF(matI);
            locF = FacetsOnInterfaceNum(matI):(FacetsOnInterfaceNum(matI+1)-1);
            facs = FacetsOnInterface(locF);
            SurfacesIi = ElementsOnFacet(facs,:);
            SurfacesI(numSI+1:numSI+numSIi,5:8) = SurfacesIi;
            numSI = numSI + numSIi;
            numel_old = numel;
            [NodesOnElement,RegionOnElement,nen,numel,nummat,MatTypeTable,MateT] = ...
             FormCZ(SurfacesIi,NodesOnElement,RegionOnElement,Coordinates,numSIi,nen_bulk,ndm,numel,nummat,6, ...
                    22,0,CZprop,MatTypeTable,MateT);
            if numel > numel_old
            RegionsOnInterface(nummat-nummatCG,:) = [nummat mat1 mat2 matI numel_old+1 numel];
            end
        end
    end
end
RegionsOnInterface = RegionsOnInterface(1:nummat-nummatCG,:);

[MPCList,numMPC,CouplerNodes] = CleanMPC(MPCListNew,numMPCnew,CouplerNodes,ndm);


% Node number for each group of nodes: solid, masterPBC, and Lagrange
% multipliers
NodeTypeNum = [1 numnp+1 numnp+2*ndm+1 numnp+2*ndm+numMPC+1]';
% now actually add the Lagrange multiplier nodes
[Coordinates,NodesOnElement,RegionOnElement,MatTypeTable,MateT,numnp,numel,nummat,nen] = AddMPNodes(MPCList,pencoeff,CornerXYZ,Coordinates,NodesOnElement,RegionOnElement,MatTypeTable,MateT,numnp,numel,nummat,nen);


% Applied displacements instead of forces
NodeLoad2 = [NodeTypeNum(2)+ndm 1 .01
            NodeTypeNum(2)+ndm 2 0
            NodeTypeNum(2)+ndm 3 0
            NodeTypeNum(2)+ndm+1 1 0
            NodeTypeNum(2)+ndm+1 2 .01
            NodeTypeNum(2)+ndm+1 3 0
            NodeTypeNum(2)+ndm+2 1 0
            NodeTypeNum(2)+ndm+2 2 0
            NodeTypeNum(2)+ndm+2 3 .01];
NodeBC = [NodeBC; NodeLoad2];
numBC = length(NodeBC);
% % Nodal forces allow rotation of cube
% NodeLoad = [NodeTypeNum(2)+ndm 1 1
%             NodeTypeNum(2)+ndm 2 0
%             NodeTypeNum(2)+ndm 3 0
%             NodeTypeNum(2)+ndm+1 1 0
%             NodeTypeNum(2)+ndm+1 2 1
%             NodeTypeNum(2)+ndm+1 3 0
%             NodeTypeNum(2)+ndm+2 1 0
%             NodeTypeNum(2)+ndm+2 2 0
%             NodeTypeNum(2)+ndm+2 3 1];
% numNodalF = length(NodeLoad);

ProbType = [numnp numel nummat ndm ndm nen];

% plotNodeCont3(Coordinates+Node_U_V*100, Node_U_V(:,1), NodesOnElement, 2,1:64+32*3)