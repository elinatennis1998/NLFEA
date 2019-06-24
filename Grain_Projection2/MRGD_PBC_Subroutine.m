clear
clc

% Mesh with nxn tiling
PSPS = 's'; %plane stress condition 
nen = 4;3; %number of nodes per element
nel = 4;3; %max number of nodes per element
numgh = 4;3;5;8;7;6; %number of grains along horiz. edge
numgs = 4;3;5;8;7;6; %number of grains along later. edge
numgrain = numgh*numgs; %total number of grains in RVE
bCrys = 2;3;8;1; %number of elem. along grain edge
tfact = 1;2; %tri 2, quad 1; 
numelemg = bCrys^2; %number of elements in a grain
nu = numgh*bCrys; %number of elements along x
nv = numgs*bCrys; %number of elements along y
neper_abaqus = 0;  %modification for Sundays Code
%% make the mesh of nodes and elements, centered at 0,0
Coordinates = [1 0 0
                2 4 0
                4 0 4
                3 4 4];
node1 = 1;
elmt1 = 1;
mat = 1;
rskip = 0;
btype = 0;2; %flag for meshing type (rect, triang, ...)
type = 'cart';
GrainA = ((abs(Coordinates(1,2))+abs(Coordinates(4,2)))/numgs)^2;
[Coordinates,NodesOnElement,RegionOnElement,numnp,numel] = block2d(type,nu,nv,node1,elmt1,mat,rskip,btype,Coordinates,nen);
Coordinates = Coordinates';
NodesOnElement = NodesOnElement';
NodesOnElement_PBC = NodesOnElement;

%% Map element ID onto grain ID
rsinc=nu*nv;
refine=[1 2 4 8 16 32 64]';
m=size(refine,1);
 
 for i=1:m
     if bCrys == refine(i);
        pow_r=i;
     end
 end
 startRegion=numgh^pow_r+bCrys;
 
 % stuck the material atthe bottom
 RegionOnElement(1:startRegion,1)=ones(1,startRegion);
 
 % stuck in the second material in the middle
 i=startRegion+1;
 while i < rsinc-startRegion+1
     
  RegionOnElement(i:i+bCrys*2-1,1)=2*ones(1,bCrys*2);
  RegionOnElement(i+bCrys*2:i+bCrys*4-1,1)= ones(1,bCrys*2);
    i=i+bCrys*4;   
 end
 
% the first material at the top
RegionOnElement(rsinc-startRegion+1:rsinc)=ones(1,startRegion); 
% Set up phase pattern and material properties
nummat = 2; %Number of materials for MRDG elements
nummat_PBC = nummat; %Number of materials for PBC elements
MatTypeTable = [1 2
                1 1];
MateT = [100 0.25 1
         300 0.3 1]; %Material properties

%% provide flags for desired outputs
OptFlag = [0 1 1 0 0 1 1]'; 
SEHist = 1;
%% separate the grains by duplicating nodes
% Also separeate out MRDG case
numnpCG = numnp;
InterTypes = tril(ones(nummat),-1); % only put CZM between the element edges between materials 1-2
DEIProgram2 %MRDG case 
NodesOnElement_MRDG = NodesOnElement; 
numel_MRDG = numel;
nen_MRDG = nen;
numelCG = numel;
%% Create periodic boundary conditions
FIX = 1;
RIGHT = (numgs*bCrys+1); %x
TOP = (numgh*bCrys+1)*(numgh*bCrys)+1; %y
TOPRIGHT = TOP + RIGHT - 1; %xy
NodeBC = [FIX 1 0
    FIX 2 0
    %           RIGHT 1 0.01
    %           RIGHT 2 0
    %           TOP 1 0
    %           TOP 2 0
    ];
numBC = size(NodeBC,1);

%% Create periodic boundary conditions
dx = 1;
dy = numgh*bCrys+1;
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

%% provide flags for desired outputs
OptFlag = [0 1 1 0 0 1 1]';
SEHist = 1;

%% Separate out PBC case 
numnpCG = numnp;
NodesOnElement = NodesOnElement_PBC;
numel = numelCG;
usePBC = 2; % flag to turn on keeping PBC
InterTypes = zeros(nummat,nummat); % only put CZM between the element edges between materials 1-2
DEIProgram2 %PBC case
ndm = 2;

%% Insert DG couplers along interfaces and RVE boundary, add a node for each
% grain, make some special coordinate lists
ndm = 2;
NodesOnElement = NodesOnElement_MRDG;
numel = numel_MRDG;
nen = nen_MRDG;
numnpCG = numnp;
InterTypes = tril(ones(nummat),-1); % only put CZM between the element edges between materials 1-2
ielI = 10; %MRDG elements on the interior
ielB = 9;  %PBC case on the boundary 
matepropI = [0 0 1 1]; %Meso case for interface 
matepropB = [0 1 1 0]; %Miscro case 

% InterDGallG was separated for two cases: MRDG and PBC
%% InterDGallG for MRDG
nummatCG = nummat;
numSI = numCL;
numelCG = numel;
nen_bulk = nen;
SurfacesI = zeros(0,8);
numSI = 0;
% Interface region information: [materialID in MateT; mat1; mat2; regI; first coupler number; last coupler number]
RegionsOnInterface = zeros(nummat*(nummat+1)/2,6);

% Add nodes for Grains: u, e, w
numnpMeso = (ndm+1)*nummatCG;
numnpMicro = numnp;
Coordinates = [Coordinates; zeros(numnpMeso,ndm)];
numnp = numnpMeso + numnpMicro;
% Add nodes for interfaces, holds the flux and jump
CoordinatesI = zeros((ndm+2)*nummat*(nummat+1)/2,ndm);
numnpI = 0;
% Add nodes for boundaries: u, e, w
CoordinatesB = zeros((ndm+1)*nummatCG*4,ndm);
numnpB = 0;

% Types of materials in domain
MaterTypeNum = [1 nummatCG+1 0 0]';

for mat2 = 1:nummatCG
    for mat1 = 1:mat2
        
        matI = GetRegionsInt(mat1,mat2);
        if InterTypes(mat2,mat1) > 0
        numSIi = numEonF(matI);
        locF = FacetsOnInterfaceNum(matI):(FacetsOnInterfaceNum(matI+1)-1);
        facs = FacetsOnInterface(locF);
        SurfacesIi = ElementsOnFacet(facs,:);
        SurfacesI(numSI+1:numSI+numSIi,5:8) = SurfacesIi;
        numSI = numSI + numSIi;
        numel_old = numel;
        numnpI = numnpI + (ndm+2);
        [NodesOnElement,RegionOnElement,nen,numel,nummat,MatTypeTable,MateT,CoordinatesI,numnpI] = ...
         FormDGG(SurfacesIi,NodesOnElement,RegionOnElement,Coordinates,numnpMicro,numnpMeso,numSIi,nen_bulk,ndm,numel,nummat,6, ...
                ielI,0,matepropI,MatTypeTable,MateT,CoordinatesI,numnpI,1);
        if numel > numel_old
        RegionsOnInterface(nummat-nummatCG,:) = [nummat mat1 mat2 matI numel_old+1 numel];
        end
        end
    end
end

nummat_MRDG = nummat;
MaterTypeNum(3) = nummat+1; %Additional material for interface
numel_MRDG = numel; %Number of elements for MRDG case 
%% PBC Treatment
pencoeff = 1e9;
CornerXYZ = [4.000000 0.000000
             0.000000 4.000000];
%I need to eliminate meso-scale to get PBC condition like Sunday's
%% InterDGallG for PBC with modifications
 for mat2 = 1:nummatCG
    for mat1 = mat2:mat2
        
        matI = GetRegionsInt(mat1,mat2);
        numSIi = numEonPBC(matI);
        locF = FacetsOnPBCNum(matI):(FacetsOnPBCNum(matI+1)-1);
        facs = FacetsOnPBC(locF);
        SurfacesIi = ElementsOnFacet(facs,:);
        SurfacesIi = ReverseFacets(SurfacesIi,NodesOnElement,Coordinates,numSIi,ndm);
        ElementsOnFacet(facs,:) = SurfacesIi;
        SurfacesI(numSI+1:numSI+numSIi,5:8) = SurfacesIi;
        [SurfacesIi,subFacets] = SplitFacets(SurfacesIi,NodesOnElement,Coordinates,numSIi,ndm);
        numSI = numSI + numSIi;
        for j = 1:size(subFacets,1)
            numel_old = numel;
            numSIii = subFacets(j,2) - subFacets(j,1) + 1;
            SurfacesIii = SurfacesIi(subFacets(j,1):subFacets(j,2),:);
            [NodesOnElement,RegionOnElement,nen,numel,nummat,MatTypeTable,MateT,CoordinatesB,numnpB] = ...
             FormDGG(SurfacesIii,NodesOnElement,RegionOnElement,Coordinates,numnpMicro,0,numSIii,nen_bulk,ndm,numel,nummat,6, ...
                    ielB,0,matepropB,MatTypeTable,MateT,CoordinatesB,numnpB,2);
            if numel > numel_old
            RegionsOnInterface(nummat-nummatCG,:) = [nummat mat1 mat2 matI numel_old+1 numel];
            end
        end
    end
end
RegionsOnInterface = RegionsOnInterface(1:nummat-nummatCG,:);
MaterTypeNum(4) = nummat+1;

%% Wiping out MRDG treatment on the boundary 
nen_MRDG = 8; %Nodes for Micro and solid elements 
zr = zeros(numel-numel_MRDG,nen-nen_MRDG);
NodesOnElement = [NodesOnElement(1:numel_MRDG,1:nen)
                [NodesOnElement(numel_MRDG+1:numel,1:nen_MRDG) zr]];

% Node number for each group of nodes: micro, meso, and boundary
NodeTypeNum = [1 numnpMicro+1 numnpMicro+numnpMeso+1 numnpMicro+numnpMeso+numnpB+1]';
Coordinates = [Coordinates; CoordinatesB(1:numnpB,1:ndm)];
numnp = numnp + numnpB;

%Inserting PBC nodes into the boundary elements
MPC_BCx=numnp+1;
MPC_BCy=numnp + 2;
numnp = numnp + 2;
Coordinates = [Coordinates; CornerXYZ];
% Put PBC nodes on boundary elements only
num_ins = numel-numel_MRDG; %Number of elements that need to bu sub-ed
zrt = zeros(num_ins,nen-nen_MRDG);
zrt(1:num_ins,1:2) = [MPC_BCx.*ones(num_ins,1) MPC_BCy.*ones(num_ins,1)];
NodesOnElement = [NodesOnElement(1:numel_MRDG,1:nen)
                [NodesOnElement(numel_MRDG+1:numel,1:nen_MRDG) zrt]]; 
%% RVE/macroscale BC
GrainIntegV
InterIntegV
uRVE = [0 0];
eRVE = [.02 -0.02*.25 0];
wRVE = 0;

%Need to manually insert material type for PBC condition
matPBC = zeros(nummat-nummat_MRDG,length(MateT));
matPBC(1:nummat-nummat_MRDG,1:2) = [ones(nummat-nummat_MRDG,1) ones(nummat-nummat_MRDG,1)];
matPBC(1:nummat-nummat_MRDG,3:4) = [4.*ones(nummat-nummat_MRDG,1) 4.*ones(nummat-nummat_MRDG,1)];
MateT = [MateT(1:nummat_MRDG,:)
    [matPBC]];

NodeLoad2 = [MPC_BCx 1 0
    MPC_BCx 2 0.2
    MPC_BCy 1 0
    MPC_BCy 2 0];
NodeBC = [NodeBC; NodeLoad2];
numBC = length(NodeBC);

%% Key Part: set up fluxes and jumps on the grain boundaries and RVE
TFS = 2; % FE 1; %Taylor 3; % Sachs  4; % Combo 5; % Fine Scale  
factorTS = 1;0;.5;.25;.75; % factor for Taylor vs Sachs; 1 means Taylor
factorT = factorTS;1/3;
factorS = 1-factorTS;1/3;
if TFS == 5
load('CoordinatesIT.mat');
CoordinatesI=CoordinatesIT;
end
% FormRVEDirichlet
ProbType = [numnp numel nummat ndm ndm nen];