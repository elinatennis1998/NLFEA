% Shear test of MP Periodic Boundary Conditions with CZ couplers
% Domain: 4x4 rectangle block
% Loading: Prescribed shear displacement.
%
% Last revision: 04/03/2017 TJT

clear
NCR = 1;
% clc

nen = 4;
nel = 4;

numc = 4;2;1;
xl = [1 0 0
      2 numc 0
      4 0 numc
      3 numc numc];
type = 'cart';
rinc = numc;
sinc = numc;
node1 = 1;
elmt1 = 1;
mat = 1;
rskip = 0;
btype = 0;
[Coordinates,NodesOnElement,~,numnp,numel] = block2d(type,rinc,sinc,node1,elmt1,mat,rskip,btype,xl,nen);
Coordinates = Coordinates';
NodesOnElement = NodesOnElement';

% Modify regions
RegionOnElement = [1 2 2 1 1 2 2 1 1 2 2 1 1 2 2 1]';
ix = [14,15,10,5,6,16,11,3,12,2,4,9,1,13,7,8];%randperm(numel);
NodesOnElement = NodesOnElement(ix,:);
RegionOnElement = RegionOnElement(ix);
nummat = 2;
MatTypeTable = [1 2
                1 1
                0 0];
MateT = [100 0.25 1
         100 0.25 1];
AlgoType = [0; 1; 0];

% Boundary conditions
FIX = 1;
RIGHT = (numc+1); %x
TOP = (numc+1)*(numc)+1; %y
TOPRIGHT = TOP + RIGHT - 1; %xy
NodeBC = [FIX 1 0
          FIX 2 0
%           RIGHT 1 0.01
%           RIGHT 2 0
%           TOP 1 0
%           TOP 2 0
          ];
numBC = size(NodeBC,1);

% Create periodic boundary conditions
dx = 1;
dy = numc+1;
dxy = TOP - 1;
PBCListx = [(RIGHT:dy:TOPRIGHT)' (FIX:dy:TOP)' RIGHT*ones(dy,1) FIX*ones(dy,1)]; %y
% remove repeats
RIGHTface = 2:dy;
PBCListx = PBCListx(RIGHTface,:);
PBCListy = [(TOP:dx:TOPRIGHT)' (FIX:dx:RIGHT)' TOP*ones(dy,1) FIX*ones(dy,1)]; %x
% remove repeats
TOPface = 1:dy;
TOPface = setdiff(TOPface,dy);
TOPface = setdiff(TOPface,1);
PBCListy = PBCListy(TOPface,:);
PBCList = [PBCListx; PBCListy];
MPCListx = [PBCListx(:,1:2) ones(size(PBCListx,1),1)*[1 0]
            PBCListx(1,3:4) [1 0]];
MPCListy = [PBCListy(:,1:2) ones(size(PBCListy,1),1)*[0 1]
            PBCListy(1,3:4) [0 1]];
MPCList = [MPCListx; MPCListy];

% Output quantity flags
DHist = 1;
FHist = 1;
SHist = 1;
SEHist = 1;


% Generate CZM interface: pairs of elements and faces, duplicate the nodes,
% update the connectivities
numnpCG = numnp;
usePBC = 2; % flag to turn on keeping PBC
InterTypes = [0 0
              1 0]; % only put CZM between the element edges between materials 1-2
DEIProgram2


% Insert couplers; use extended DG mesh with ghost elements so that
% couplers appear properly on the domain boundary
ndm = 2;
% CZ element stiffness
CZprop = 50;50000;1;
% MPC spring element properties
pencoeff = 1e9;
CornerXYZ = [0.000000 0.000000
             1.000000 0.000000];

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
            [MPCListNew,numMPCnew,NodesOnLinkNew,NodesOnLinknumNew,CouplerNodes] = ...
             FormMPC(MPCListNew,numMPCnew,NodesOnLinkNew,NodesOnLinknumNew,CouplerNodes, ...
                 SurfacesIi,NodesOnElement,NodesOnElementCG,Coordinates,numSIi,2, ...
                 MPCList,NodesOnLink,NodesOnLinknum);
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
NodeLoad2 = [NodeTypeNum(2)+ndm 1 0
            NodeTypeNum(2)+ndm 2 0.01
            NodeTypeNum(2)+ndm+1 1 0.01
            NodeTypeNum(2)+ndm+1 2 0];
NodeBC = [NodeBC; NodeLoad2];
numBC = length(NodeBC);
% % Nodal forces allow rotation of cube
% NodeLoad = [NodeTypeNum(2)+ndm 1 1
%             NodeTypeNum(2)+ndm 2 0
%             NodeTypeNum(2)+ndm+1 1 0
%             NodeTypeNum(2)+ndm+1 2 1];
% numNodalF = length(NodeLoad);

ProbType = [numnp numel nummat ndm ndm nen];

% plotMesh2(Coordinates,NodesOnElement,1,(1:size(NodesOnElement,1)),[1 1 1 1 1])
% plotNodeCont2(Coordinates+Node_U_V*100, Node_U_V(:,1), NodesOnElement, 2, 1:24, [1 1 1], 0,[3 4 6 9 0])
