clear
clc

% Mesh with nxn tiling
PSPS = 's'; %plane stress condition 
nen = 4;3; %number of nodes per element
nel = 4;3; %max number of nodes per element
numgh = 3;4;5;8;7;6; %number of grains along horiz. edge
numgs = 3;4;5;8;7;6; %number of grains along later. edge
numgrain = numgh*numgs; %total number of grains in RVE
bCrys = 2;1;8;3; %number of elem. along grain edge
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
%Assign individual region and material property to each grain 
%Note: alterphase is a built-in function for maerial assignment, id desired
%otherwise, must hard code
m2 = [200 0.25 1]; %material type 1
m1 = [100 0.3 1]; %material type 2
[RegionOnElement,MatTypeTable,MateT,nummat,grainG] = GeomProp(numgrain,numelemg,tfact,numgs,numgh,bCrys,m2,m1,nel,nu);
% MateT = [m1; m1; m1; m1; m2;m1; m1; m1; m1];
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
Create_PBC_2d
%% provide flags for desired outputs
OptFlag = [0 1 1 0 0 1 1]';
SEHist = 1;

%% Separate out PBC case 
numnpCG = numnp;
NodesOnElement = NodesOnElement_PBC;
numel = numelCG;
usePBC = 2; % flag to turn on keeping PBC
DEIProgram2 %PBC case

%% Insert DG couplers along interfaces and RVE boundary, add a node for each
% grain, make some special coordinate lists
ndm = 2;
NodesOnElement = NodesOnElement_MRDG;
numel = numel_MRDG;
nen = nen_MRDG;
numnpCG = numnp;
ielI = 10; %MRDG elements on the interior
ielB = 9;  %PBC case on the boundary 
matepropI = [0 0 1 1]; %Meso case for interface 
matepropB = [0 1 1 0]; %Miscro case 
pencoeff = 1e9;
CornerXYZ = [4.000000 0.000000
             0.000000 4.000000]; 
InterDGallG_PBC
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
nummatCG = nummat_PBC;
uRVE = [0 0];
eRVE = [.02 -0.02*.25 0];
wRVE = 0;
GrainIntegV
InterIntegV

%Need to manually insert material type for PBC condition
% matPBC = zeros(nummat-nummat_MRDG,length(MateT));
matPBC = zeros(nummat-nummat_MRDG,8);
matPBC(1:nummat-nummat_MRDG,1:2) = [ones(nummat-nummat_MRDG,1) ones(nummat-nummat_MRDG,1)];
matPBC(1:nummat-nummat_MRDG,3:4) = [4.*ones(nummat-nummat_MRDG,1) 4.*ones(nummat-nummat_MRDG,1)];
MateT = [MateT(1:nummat_MRDG,:)
    [matPBC]];

NodeLoad2 = [MPC_BCx 1 0
            MPC_BCx 2 -0.2
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

%% Locking Grains
locked_g = [5];
% % %     for j = 1:numgrain
% % %         locked_g(1,j) = j;
% % %     end
% %     for k = 1: numgrain 
% %         MateT(k,4) = GrainVol(1,k)/ MateT(k,3); %Grain area
% %     end 
num_locked_g = length(locked_g); 
meso_nen = 3; %Number of nodes per mose element
sub = 11; %subroutine number for a meso element
nen_bulk = 4; %Reset variables for locking case
% combo = 1; %Flag for plotting FS on CS

[NodeBC,NodesOnElement,RegionOnElement,nummat,MateT,MatTypeTable,numBC,numel] = Meso_Locking_patch(NodesOnElement,num_locked_g,NodeBC,locked_g,....
    meso_nen,grainG,nen_bulk,numelemg,GrainA,numel,numnpMicro,numnpMeso,MateT,numgrain,nummat,MatTypeTable,RegionOnElement);
ProbType = [numnp numel nummat ndm ndm nen];