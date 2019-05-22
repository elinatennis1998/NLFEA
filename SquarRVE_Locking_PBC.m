% CST solid element with 4x4 grain RVE. Coarse
% scale implemented as a material parameter

clear
clc

% Mesh with nxn tiling
PSPS = 's'; %plane stress condition 
nen = 4;3; %number of nodes per element
nel = 4;3; %max number of nodes per element
numgh = 3;5;7;6;4; %number of grains along horiz. edge
numgs = 3;5;7;6;4; %number of grains along later. edge
numgrain = numgh*numgs; %total number of grains in RVE
bCrys = 2;1; %number of elem. along grain edge
tfact = 1;2; %tri 2, quad 1; 
numelemg = bCrys^2; %number of elements in a grain
nu = numgh*bCrys; %number of elements along x
nv = numgs*bCrys; %number of elements along y

flag_PBC = 1; %Switch for periodic BC
pencoeff = 1e9;
CornerXYZ = [-4.500000 -4.500000
             -4.500000 -4.500000]; 
%% make the mesh of nodes and elements, centered at 0,0
Coordinates = [1 -4.5 -4.5
             2 4.5 -4.5
             4 -4.5 4.5
             3 4.5 4.5];
% GrainA = ((abs(Coordinates(1,2))+abs(Coordinates(4,2)))/numgs)^2;
GrainA = ((abs(Coordinates(1,2))+abs(Coordinates(4,2)))/numgs)*((abs(Coordinates(2,2))+abs(Coordinates(3,2)))/numgh);
node1 = 1;
elmt1 = 1;
mat = 1;
rskip = 0;
btype = 0;2; %flag for meshing type (rect, triang, ...)
type = 'cart';
[Coordinates,NodesOnElement,RegionOnElement,numnp,numel] = block2d(type,nu,nv,node1,elmt1,mat,rskip,btype,Coordinates,nen);
Coordinates = Coordinates';
NodesOnElement = NodesOnElement';
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
    
%     load('Locking2_EG.mat');
%     BoundXmid = Locking2_EG;
%% Set up phase pattern and material properties
nummat = numgrain;
MatTypeTable = [1:numgrain; ones(1,numgrain)];
% m1 = [100 0.25 1]; %material type 1
m2 = [100 0.3 1]; %material type 2
m1 = m2;
mats = [m1; m2];
alterphase = 2*rem(1:numgrain,2) - 1;
alterphase(alterphase==-1) = 2;
MateT = mats(alterphase,:);%% Output quantity flags
% nummat = numgrain;
% MateT = [m1; m2; m1; m2; m2; m1; m2; m1; m1; m2; m1; m2; m2; m1; m2; m1];

%% provide flags for desired outputs
OptFlag = [0 1 1 0 0 1 1]';
SEHist = 1;
%% separate the grains by duplicating nodes
numnpCG = numnp;
InterTypes = tril(ones(nummat),-1); % only put CZM between the element edges between materials 1-2
DEIProgram2

%% Insert DG couplers along interfaces and RVE boundary, add a node for each
% grain, make some special coordinate lists
ndm = 2;
ielI = 10;
ielB = 10;
matepropI = [0 0 1 1];
matepropB = [0 1 1 0];
% InterDGallG

%% RVE/macroscale BC
% GrainIntegV
% InterIntegV

%Modification of 4/16/2019 
%Periodic Boundary Conditions. Sunday's code 
if flag_PBC ==1
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
    
    % Create periodic boundary conditions
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
    
    % End of modification of 4/16/2019
else
    uRVE = [0 0];
    eRVE = [.02 -0.02*.25 0];
    wRVE = 0;
end

InterDGallG
GrainIntegV
InterIntegV

%Modification of 4/16/2019
if flag_PBC == 1
    MPC_BCx=numnp+1;
    MPC_BCy=numnp + 2;
    numnp = numnp + 2;
    Coordinates = [Coordinates; CornerXYZ];
    rem_num = (numel-numelCG);
    NodesOnElement = [NodesOnElement [zeros(numelCG,2); [MPC_BCx*ones(rem_num,1) MPC_BCy*ones(rem_num,1)]]];

end
%End of modification of 4/16/2019

%% Key Part: set up fluxes and jumps on the grain boundaries and RVE
TFS = 1; %Taylor 2; % FE 3; % Sachs 4; % Combo 5; % Fine Scale  
factorTS = 1;0;.5;.25;.75; % factor for Taylor vs Sachs; 1 means Taylor
factorT = factorTS;1/3;
factorS = 1-factorTS;1/3;
if TFS == 1
% factorTS = 1;
% TaylorSet
elseif TFS == 3
% factorTS = 0;
% SachsSet
elseif TFS == 4
% TaylorSet
% SachsSet
elseif TFS == 5
load('CoordinatesIT.mat');
CoordinatesI=CoordinatesIT;
end
if flag_PBC == 0
    FormRVEDirichlet
end
%% Zero out meso grains; for this file as it is, I don't want to compute the
% the sequence has to be retained: from low to high order of grain number

% locked_g = [9 11 13 17 19 23 25 27 31 33 37 39 41]'; %7x7 square RVE
%  locked_g = [9 10 11 12 13 17 19 23 25 27]'; %7x5 RVE
% locked_g = [1 2 3 4 5 6 10 11 15 16 20 21 22 23 24 25];
% locked_g = [5];
% num_locked_g = length(locked_g);
% meso_nen = 3;
% FSon = 0;
% sub = 11; %Sunroutine number for Meso
% %%Currently does not work for 2 or more grains in a row being locked
% [NodeBC,NodesOnElement,RegionOnElement,nummat,MateT,MatTypeTable,numBC,numel,ccell] = Meso_Locking(NodesOnElement,num_locked_g,NodeBC,locked_g,....
%     meso_nen,grainG,nen_bulk,numel,numnpMicro,numnpMeso,MateT,numgrain,nummat,MatTypeTable,RegionOnElement,sub,bCrys,tfact);
if flag_PBC == 1
        nen = nen+2; %Expanding for PBC for x and y directions
         NodeLoad = [MPC_BCx 1 0 %eRVE (1) deformtion in x
                     MPC_BCx 2 0    %eRVE (3)
                     MPC_BCy 1 2    %wRVE Shear
                     MPC_BCy 2 0];  %eRVE (2) -0.02*.25 deformation in y
        numNodalF = length(NodeLoad);
%         NodeBC = [NodeBC; NodeLoad];
        NodeBC = [NodeBC; numgh*bCrys+1 2 0];
        numBC = length(NodeBC);
end
ProbType = [numnp numel nummat 2 2 nen];

%Plotting functions:
% plotNodeCont2(Coordinates+Node_U_V, Node_U_V(:,1), NodesOnElement, 1, (1:numel));
% plotNodeCont2(Coordinates+Node_U_V, Node_U_V(:,2), NodesOnElement, 2, (1:numel));