% CST solid element with 4x4 grain RVE. Coarse
% scale implemented as a material parameter

clear
clc

% Mesh with nxn tiling
PSPS = 's'; %plane stress condition 
nen = 4;3; %number of nodes per element
nel = 4;3; %max number of nodes per element
numgh = 3;5;4;8;7;6; %number of grains along horiz. edge
numgs = 3;5;4;8;7;6; %number of grains along later. edge
numgrain = numgh*numgs; %total number of grains in RVE
bCrys = 1;2;8;3; %number of elem. along grain edge
tfact = 1;2; %tri 2, quad 1; 
numelemg = bCrys^2; %number of elements in a grain
nu = numgh*bCrys; %number of elements along x
nv = numgs*bCrys; %number of elements along y
%% make the mesh of nodes and elements, centered at 0,0
Coordinates = [1 0 0
             2 3 0
             4 0 3
             3 3 3];
GrainA = ((abs(Coordinates(1,2))+abs(Coordinates(4,2)))/numgs)^2; %For Square RVE
% GrainA = ((abs(Coordinates(1,2))+abs(Coordinates(4,2)))/numgs)*((abs(Coordinates(2,2))+abs(Coordinates(3,2)))/numgh);
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
%Assign individual region and material property to each grain 
%Note: alterphase is a built-in function for maerial assignment, id desired
%otherwise, must hard code
m2 = [200 0.25 1]; %material type 1
m1 = [100 0.3 1]; %material type 2
[RegionOnElement,MatTypeTable,MateT,nummat,grainG] = GeomProp(numgrain,numelemg,tfact,numgs,numgh,bCrys,m2,m1,nel,nu);
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
uRVE = [0 0];
% eRVE = [.02 -0.02*.25 0];
eRVE = [0.2 0 0];
wRVE = 0;
GrainIntegV
InterIntegV
%% Key Part: set up fluxes and jumps on the grain boundaries and RVE
TFS = 1; %Taylor 3; % Sachs 2; % FE  4; % Combo 5; % Fine Scale  
factorTS = 1;0;.5;.25;.75; % factor for Taylor vs Sachs; 1 means Taylor
factorT = factorTS;1/3;
factorS = 1-factorTS;1/3;
if TFS == 5
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
% locked_g = [5];
%     for j = 1:numgrain
%         locked_g(1,j) = j;
%     end
% %     for k = 1: numgrain 
% %         MateT(k,4) = GrainVol(1,k)/ MateT(k,3); %Grain area
% %     end 
% num_locked_g = length(locked_g); 
% meso_nen = 3; %Number of nodes per mose element
% sub = 11; %subroutine number for a meso element
% combo = 0; %Flag for plotting FS on CS
% % 
% [NodeBC,NodesOnElement,RegionOnElement,nummat,MateT,MatTypeTable,numBC,numel] = Meso_Locking_patch(NodesOnElement,num_locked_g,NodeBC,locked_g,....
%     meso_nen,grainG,nen_bulk,numelemg,GrainA,numel,numnpMicro,numnpMeso,MateT,numgrain,nummat,MatTypeTable,RegionOnElement);
ProbType = [numnp numel nummat 2 2 nen];