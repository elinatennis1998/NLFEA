% CST solid element with 4x4 grain RVE. Coarse
% scale implemented as a material parameter

clear
clc

% Mesh with nxn tiling
PSPS = 's'; %plane stress condition 
nen = 3;4; %number of nodes per element
nel = 3;4; %max number of nodes per element
numgh = 3;5;7;6;4; %number of grains along horiz. edge
numgs = 3;5;7;6;4; %number of grains along later. edge
numgrain = numgh*numgs; %total number of grains in RVE
bCrys = 2;1;3; %number of elem. along grain edge
tfact = 2;1; %tri 2, quad 1; 
numelemg = bCrys^2; %number of elements in a grain
nu = numgh*bCrys; %number of elements along x
nv = numgs*bCrys; %number of elements along y
%% make the mesh of nodes and elements, centered at 0,0
Coordinates = [1 -4.5 -4.5
             2 4.5 -4.5
             4 -4.5 4.5
             3 4.5 4.5];
% GrainA = ((abs(Coordinates(1,2))+abs(Coordinates(4,2)))/numgs)^2; %For Square RVE
% GrainA = ((abs(Coordinates(1,2))+abs(Coordinates(4,2)))/numgs)*((abs(Coordinates(2,2))+abs(Coordinates(3,2)))/numgh);
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
m1 = [200 0.25 1]; %material type 1
m2 = [100 0.3 1]; %material type 2
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
InterDGallG

%% RVE/macroscale BC
GrainIntegV
InterIntegV
uRVE = [0 0];
eRVE = [.02 -0.02*.25 0];
wRVE = 0;

% load('BoundXMid_SquareL2')
% BoundXMid_SquareL2 = BoundXmid;
%% Key Part: set up fluxes and jumps on the grain boundaries and RVE
TFS = 1; %Taylor 3; % Sachs 2; % FE  4; % Combo 5; % Fine Scale  
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
% the sequence has to be retained: from low to high order of grain number

% locked_g = [9 11 13 17 19 23 25 27 31 33 37 39 41]'; %7x7 square RVE
%  locked_g = [9 10 11 12 13 17 19 23 25 27]'; %7x5 RVE %gives issues 
% locked_g = [1 2 3 4 5 6 10 11 15 16 20 21 22 23 24 25]; %gives issues
% locked_g = [1 2 3 4 5]; %gives issues 
% locked_g = [7 8 9 12 13 14 17 18 19];
locked_g = [5];
%     for j = 1:numgrain
%         locked_g(1,j) = j;
%     end
    for k = 1: numgrain 
        MateT(k,4) = GrainVol(1,k)/ MateT(k,3); %Grain area
    end 
num_locked_g = length(locked_g); 
meso_nen = 3; %Number of nodes per mose element
sub = 11; %subroutine number for a meso element
combo = 1; %Flag for plotting FS on CS

[NodeBC,NodesOnElement,RegionOnElement,nummat,MateT,MatTypeTable,numBC,numel,ccell] = Meso_Locking(NodesOnElement,num_locked_g,NodeBC,locked_g,....
    meso_nen,grainG,nen_bulk,numel,numnpMicro,numnpMeso,MateT,numgrain,nummat,MatTypeTable,RegionOnElement,sub,bCrys,tfact);
ProbType = [numnp numel nummat 2 2 nen];