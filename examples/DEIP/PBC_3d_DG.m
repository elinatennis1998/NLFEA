% Axial test of Periodic Boundary Conditions with CZ couplers, 2 grains.
% Domain: 4x4x4 rectangle block
% Loading: Prescribed shear displacement.
%
% Last revision: 03/27/2017 TJT

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
% RegionOnElement = [1 1 2 2 1 1 2 2 1 1 2 2 1 1 2 2 ...
%                    1 1 2 2 1 1 2 2 1 1 2 2 1 1 2 2 ...
%                    1 1 2 2 1 1 2 2 1 1 2 2 1 1 2 2 ...
%                    1 1 2 2 1 1 2 2 1 1 2 2 1 1 2 2 ]';
RegionOnElement = [1 2 2 1 1 2 2 1 1 2 2 1 1 2 2 1 ...
                   1 2 2 1 1 2 2 1 1 2 2 1 1 2 2 1 ...
                   1 2 2 1 1 2 2 1 1 2 2 1 1 2 2 1 ...
                   1 2 2 1 1 2 2 1 1 2 2 1 1 2 2 1 ]';
ix = randperm(numel);
% ix = [42,54,40,1,19,52,38,59,58,45,16,48,28,63,6,35,22,2,17,32,44,23,36,39,11,51,53,34,46,4,50,21,18,12,41,47,13,26,55,20,62,8,31,33,57,3,25,9,30,29,56,14,60,15,64,37,43,5,61,24,27,7,49,10];
% ix = [59,35,28,27,55,48,63,57,32,4,5,51,37,41,13,49,10,14,8,6,43,22,39,36,12,17,25,64,56,47,34,16,40,29,53,3,20,26,33,19,42,15,44,45,46,24,23,60,30,38,9,61,52,18,7,62,1,50,21,11,31,2,58,54];
% ix=[33,59,24,49,16,17,60,18,23,8,2,21,9,64,40,47,14,10,38,62,22,20,41,5,37,15,26,43,58,56,29,46,7,35,50,51,48,19,28,52,45,11,1,55,13,39,30,12,63,44,25,3,61,34,31,6,53,42,57,27,54,32,4,36];
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
          FRONT 1 0
          FRONT 2 0
          FRONT 3 0
          RIGHT 1 0.01
          RIGHT 2 0
          RIGHT 3 0
          TOP 1 0
          TOP 2 0.01
          TOP 3 0
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

% Generate CZM interface: pairs of elements and faces, duplicate the nodes,
% update the connectivities
numnpCG = numnp;
usePBC = 1; % flag to turn on keeping PBC
InterTypes = [0 0
              1 0]; % only put CZM between the element edges between materials 1-2
DEIProgram3


% Insert couplers; use extended DG mesh with ghost elements so that
% couplers appear properly on the domain boundary
numSI = numCL;
numelCG = numel;
nen_bulk = nen;
SurfacesI = zeros(0,8);
numSI = 0;
CZstiff = 50;50000;1;
PBCListNew = zeros(0,5);
numPBCnew = 0;
TieNodesNew = zeros(3,2);
NodesOnLinkNew = zeros(4,numnp);
NodesOnLinknumNew = zeros(numnp,1);
for mat2 = 1:nummat
    for mat1 = 1:mat2
        
        matI = mat2*(mat2-1)/2 + mat1; % ID for material pair (row=mat2, col=mat1)
        numSIi = numEonPBC(matI);
        if numSIi > 0 % Create PBC pairs first
            if InterTypes(mat2,mat1) > 0% then make PBC-CZ couplers
            locF = FacetsOnPBCNum(matI):(FacetsOnPBCNum(matI+1)-1);
            facs = FacetsOnPBC(locF);
            SurfacesIi = ElementsOnFacet(facs,:);
            SurfacesI(numSI+1:numSI+numSIi,5:8) = SurfacesIi;
            [NodesOnElement,RegionOnElement,~,numnp,numel,nummat,MatTypeTable,MateT,Coordinates,...
                PBCListNew,numPBCnew,NodesOnLinkNew,NodesOnLinknumNew,TieNodesNew] = ...
             FormPBCCZ(PBCListNew,numPBCnew,NodesOnLinkNew,NodesOnLinknumNew,TieNodesNew,...
                    SurfacesIi,NodesOnElement,NodesOnElementCG,RegionOnElement,Coordinates,numSIi,nen_bulk,3,numnp,numel,nummat,6, ...
                    22,0,CZstiff,PBCList,TieNodes,NodesOnLink,NodesOnLinknum,MatTypeTable,MateT);
            % Then form regular CZ couplers
            numSIi = numEonF(matI) - numSIi;
            locF = FacetsOnIntMinusPBCNum(matI):(FacetsOnIntMinusPBCNum(matI+1)-1);
            facs = FacetsOnIntMinusPBC(locF);
            SurfacesIi = ElementsOnFacet(facs,:);
            SurfacesI(numSI+1:numSI+numSIi,5:8) = SurfacesIi;
            numSI = numSI + numSIi;
            [NodesOnElement,RegionOnElement,nen,numel,nummat,MatTypeTable,MateT] = ...
             FormCZ(SurfacesIi,NodesOnElement,RegionOnElement,Coordinates,numSIi,nen_bulk,3,numel,nummat,6, ...
                    22,0,CZstiff,MatTypeTable,MateT);
            else % just form PBC pairs
            locF = FacetsOnPBCNum(matI):(FacetsOnPBCNum(matI+1)-1);
            facs = FacetsOnPBC(locF);
            SurfacesIi = ElementsOnFacet(facs,:);
            SurfacesI(numSI+1:numSI+numSIi,5:8) = SurfacesIi;
            [PBCListNew,numPBCnew,NodesOnLinkNew,NodesOnLinknumNew,TieNodesNew] = ...
            FormPBC(PBCListNew,numPBCnew,NodesOnLinkNew,NodesOnLinknumNew,TieNodesNew,SurfacesIi,...
                NodesOnElement,NodesOnElementCG,Coordinates,numSIi,3,PBCList, ...
                TieNodes,NodesOnLink,NodesOnLinknum);
            end
        elseif InterTypes(mat2,mat1) > 0 % just make CZ couplers
            numSIi = numEonF(matI);
            locF = FacetsOnInterfaceNum(matI):(FacetsOnInterfaceNum(matI+1)-1);
            facs = FacetsOnInterface(locF);
            SurfacesIi = ElementsOnFacet(facs,:);
            SurfacesI(numSI+1:numSI+numSIi,5:8) = SurfacesIi;
            numSI = numSI + numSIi;
            [NodesOnElement,RegionOnElement,nen,numel,nummat,MatTypeTable,MateT] = ...
             FormCZ(SurfacesIi,NodesOnElement,RegionOnElement,Coordinates,numSIi,nen_bulk,3,numel,nummat,6, ...
                    22,0,CZstiff,MatTypeTable,MateT);
        end
    end
end

% Update boundary conditions
NodeBCCG = NodeBC;
numBCCG = numBC;
[NodeBC,numBC] = UpdateMasterBC(TieNodes,TieNodesNew,NodeBCCG,numBCCG);

[PBCList,numPBC,TieNodes] = CleanPBC(TieNodes,3,PBCListNew,numPBCnew,TieNodesNew);

% now actually add the Lagrange multiplier nodes
[Coordinates,NodesOnElement,RegionOnElement,MatTypeTable,MateT,numnp,numel,nummat,nen] = AddLMNodes(PBCList,Coordinates,NodesOnElement,RegionOnElement,MatTypeTable,MateT,numnp,numel,nummat,nen);

ProbType = [numnp numel nummat 3 3 nen];

% plotMesh3(Coordinates,NodesOnElement,1,(1:64),(1:size(Coordinates,1))',[1 1 1 1 1])
% plotNodeCont3(Coordinates+Node_U_V*100, Node_U_V(:,1), NodesOnElement, 2,1:64+32)