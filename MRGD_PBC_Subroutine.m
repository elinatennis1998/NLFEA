clear
clc

% Mesh with nxn tiling
PSPS = 's'; %plane stress condition 
nen = 4;3; %number of nodes per element
nel = 4;3; %max number of nodes per element
numgh = 4;3;8;5;7;6; %number of grains along horiz. edge
numgs = 4;3;8;5;7;6; %number of grains along later. edge
numgrain = numgh*numgs; %total number of grains in RVE
bCrys = 2;1;8;3; %number of elem. along grain edge
tfact = 1;2; %tri 2, quad 1; 
numelemg = bCrys^2; %number of elements in a grain
nu = numgh*bCrys; %number of elements along x
nv = numgs*bCrys; %number of elements along y
%% make the mesh of nodes and elements, centered at 0,0
Coordinates = [1 0 0
                2 4 0
                4 0 4
                3 4 4];
xl = [1 0 0
    2 4 0
    4 0 4
    3 4 4];
node1 = 1;
elmt1 = 1;
mat = 1;
rskip = 0;
btype = 0;2; %flag for meshing type (rect, triang, ...)
type = 'cart';

[Coordinates,NodesOnElement,RegionOnElement,numnp,numel] = block2d(type,nu,nv,node1,elmt1,mat,rskip,btype,Coordinates,nen);
Coordinates = Coordinates';
NodesOnElement = NodesOnElement';
NodesOnElement_PBC = NodesOnElement;

%% Map element ID onto grain ID
%Assign individual region to each grain 
% for each grain, which elements belong to it
grainG = zeros(numgrain,numelemg*tfact);
grain = 0;
% for each grain, which elements belong to it
for j = 1:numgs
    for i = 1:numgh
        el = 0;
        grain = grain + 1;
        for m = 1:bCrys
            for l = 1:bCrys*tfact
                if nel == 3
                    elem = (j-1)*nu*bCrys*(j*tfact-2*(j-1))+(i-1)*bCrys*tfact; % bottom-corner element of grain
                    elem = elem + (m-1)*nu*tfact;
                    elem = elem + (l-1) + 1;
                elseif nel == 4
                    elem = (j-1)*nu*bCrys+(i-1)*bCrys;
                    elem = elem + (m-1)*nu*tfact;
                    elem = elem + (l-1) + 1;
                end
                el = el + 1;
                grainG(grain,el) = elem;
            end
        end
    end
end
% inverse map: the grain that an element belongs to
    for g = 1:numgrain
        RegionOnElement(grainG(g,:)) = g;
    end
    
%% Set up phase pattern and material properties
nummat = numgrain;
MatTypeTable = [1:numgrain; ones(1,numgrain)];
m2 = [200 0.25 1]; %material type 1
m1 = [100 0.3 1]; %material type 2
mats = [m1; m2];
alterphase = 2*rem(1:numgrain,2) - 1;
alterphase(alterphase==-1) = 2;
MateT = mats(alterphase,:);%% Output quantity flags

%% provide flags for desired outputs
OptFlag = [0 1 1 0 0 1 1]';
SEHist = 1;
%% separate the grains by duplicating nodes
numnpCG = numnp;
InterTypes = tril(ones(nummat),-1); % only put CZM between the element edges between materials 1-2
DEIProgram2
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

numnpCG = numnp;
NodesOnElement = NodesOnElement_PBC;
numel = numelCG;
usePBC = 2; % flag to turn on keeping PBC
InterTypes = zeros(16,16); % only put CZM between the element edges between materials 1-2
DEIProgram2

ndm = 2;
% MPC spring element properties
pencoeff = 1e9;
CornerXYZ = [4.000000 0.000000
             0.000000 4.000000];

nummatCG = nummat;
numSI = numCL;
numelCG = numel;
nen_bulk = nen;
SurfacesI = zeros(0,8);
numSI = 0;

%% Insert DG couplers along interfaces and RVE boundary, add a node for each
% grain, make some special coordinate lists
ndm = 2;
ielI = 10;
ielB = 10;
matepropI = [0 0 1 1];
matepropB = [0 1 1 0];
InterDGallG

%% RVE/macroscale BC
GrainIntegV
InterIntegV
uRVE = [0 0];
eRVE = [1 0 0];
wRVE = 0;
for k = 1: numgrain
    MateT(k,4) = GrainVol(1,k)/ MateT(k,3); %Grain area
end
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