% CST solid element with 3x3 grain RVE. Coarse
% scale implemented as a material parameter
% The grain in the middle is coarse scale only
%EG patch for new subroutine
%Created 03/2019
%Last modified 4/16/2019
clear
clc

% Mesh with 6x6 tiling
PSPS = 's'; %plane stress condition 
nen = 3;4; %number of nodes per element
nel = 3;4; %max number of nodes per element
numgh = 3;5; %number of grains along horiz. edge
numgs = 3;5; %number of grains along later. edge
numgrain = numgh*numgs; %total number of grains in RVE
bCrys = 1;2;3; %number of elem. along grain edge
tfact = 2;1; %tri 2, quad 1; 
numelemg = bCrys^2; %number of elements in a grain
nu = numgh*bCrys; %number of elements along x
nv = numgs*bCrys; %number of elements along y
FSon = 0; %Flag to utilize FS and CS simultaneously
%% make the mesh of nodes and elements, centered at 0,0
Coordinates = [1 -4.5 -4.5
             2 4.5 -4.5
             4 -4.5 4.5
             3 4.5 4.5];
node1 = 1;
elmt1 = 1;
mat = 1;
rskip = 0;
btype = 2;0; %flag for meshing type (rect, triang, ...)
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
%% Set up phase pattern and material properties
nummat = numgrain;
MatTypeTable = [1:numgrain; ones(1,numgrain)];
m1 = [100 0.25 1]; %material type 1
m2 = [200 0.3 1]; %material type 2
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
eRVE = [.02 -0.02*.25 0];
wRVE = 0;

%% Key Part: set up fluxes and jumps on the grain boundaries and RVE
TFS = 1; %Taylor 3; % Sachs 4; % Combo 5; % Fine Scale 2; % FE
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
FormRVEDirichlet

%% Zero out meso grains; for this file as it is, I don't want to compute the
% mesoscale stuff too
%     for j = 1:numgrain
%         locked_g(1,j) = j;
%     end
    locked_g = [5]';
    FSon = 0; %Flag to utilize FS and CS simultaneously
    sub = 1; %Subroutine number for Meso
%     
%     %Modification to add grain area 4/26/19
    for k = 1: numgrain 
        MateT(k,4) = GrainVol(1,k);
    end 
    num_locked_g = length(locked_g);
    meso_nen = 3;
%     
    [NodeBC,NodesOnElement,RegionOnElement,nummat,MateT,MatTypeTable,numBC,numel] = Meso_Locking_2(NodesOnElement,num_locked_g,NodeBC,locked_g,....
    meso_nen,grainG,nen_bulk,numelemg,numel,numnpMicro,numnpMeso,MateT,numgrain,nummat,MatTypeTable,RegionOnElement,FSon,bCrys,tfact);
% % 
%     [NodeBC,NodesOnElement,RegionOnElement,nummat,MateT,MatTypeTable,numBC,numel] = Meso_Locking(NodesOnElement,num_locked_g,NodeBC,locked_g,....
%     meso_nen,grainG,nen_bulk,numelemg,GrainA,numel,numnpMicro,numnpMeso,MateT,numgrain,nummat,MatTypeTable,RegionOnElement,sub,bCrys,tfact,FSon);

% MateT(34,5) = 1; %Flag = 1 to turn on CS only

ProbType = [numnp numel nummat 2 2 nen];

% micro5 = unique(reshape(NodesOnElement(grainG(5,:),1:3),2*3,1));
% l5 = length(micro5);
% NodeBC = [NodeBC
% [(numnpMicro+1:numnpMicro+4*3)'   ones(4*3,1) zeros(4*3,1)
% (numnpMicro+1:numnpMicro+4*3)' 2*ones(4*3,1) zeros(4*3,1)]
% [(numnpMicro+5*3+1:numnpMicro+numnpMeso)'   ones(4*3,1) zeros(4*3,1)
% (numnpMicro+5*3+1:numnpMicro+numnpMeso)' 2*ones(4*3,1) zeros(4*3,1)]
% [micro5 ones(l5,1) zeros(l5,1)
%  micro5 2*ones(l5,1) zeros(l5,1)]];
% numBC = length(NodeBC);
% % Correction for grain 5 of mesoscale
% MateT(13,5:8) = [1 0 0 1];
% MateT(14,5:8) = [1 0 0 1];
% MateT(16,5:8) = [0 1 1 0];
% MateT(18,5:8) = [0 1 1 0];
% MateT(5,4) = 1;
% MateT(5,5) = GrainVol(1,5);
% 
% ProbType = [numnp numel nummat 2 2 nen];
