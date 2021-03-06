% CST solid element with 4x4 grain RVE. Coarse
% scale implemented as a material parameter

clear

% Mesh with nxn tiling
PSPS = 's'; %plane stress condition 
nen = 4;3; %number of nodes per element
nel = 4;3; %max number of nodes per element
numgh = 5;3;4;8;7;6; %number of grains along horiz. edge
numgs = 5;3;4;8;7;6; %number of grains along later. edge
numgrain = numgh*numgs; %total number of grains in RVE
bCrys = 2;1;3;8; %number of elem. along grain edge
CS_flag = 1; %Flag to collect data for post processing
if nen == 4
    tfact = 1; %quad 1
    btype = 0; %flag for meshing type (rect, triang, ...)
elseif nen == 3
    tfact = 2; %tri 2
    btype = 2;
end 
numelemg = bCrys^2; %number of elements in a grain
nu = numgh*bCrys; %number of elements along x
nv = numgs*bCrys; %number of elements along y
%% make the mesh of nodes and elements, centered at 0,0
Coordinates = [1 0 0
             2 5 0
             4 0 5
             3 5 5];
GrainA = ((abs(Coordinates(1,2))+abs(Coordinates(4,2)))/numgs)^2; %For Square RVE
% GrainA = ((abs(Coordinates(1,2))+abs(Coordinates(4,2)))/numgs)*((abs(Coordinates(2,2))+abs(Coordinates(3,2)))/numgh);
node1 = 1;
elmt1 = 1;
mat = 1;
rskip = 0;
type = 'cart';
[Coordinates,NodesOnElement,RegionOnElement,numnp,numel] = block2d(type,nu,nv,node1,elmt1,mat,rskip,btype,Coordinates,nen);
Coordinates = Coordinates';
NodesOnElement = NodesOnElement';
%% Map element ID onto grain ID
%Assign individual region and material property to each grain 
%Note: alterphase is a built-in function for maerial assignment, id desired
%otherwise, must hard code
m2 = [100 0.25 1]; %material type 1
m1 = [200 0.25 1]; %material type 2
[RegionOnElement,MatTypeTable,MateT,nummat,grainG,BoundGrains,CornerGrain] = GeomProp(numgrain,numelemg,tfact,numgs,numgh,bCrys,m2,m1,nel,nu);
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
eRVE = [0.2 0 0];
% eRVE = [1 0 0];
wRVE = 0;
GrainIntegV
InterIntegV
%% Key Part: set up fluxes and jumps on the grain boundaries and RVE
TFS = 5; % Fine Scale 2; % FE 1; %Taylor 3; % Sachs  4; % Combo 
factorTS = 1;0;.5;.25;.75; % factor for Taylor vs Sachs; 1 means Taylor
factorT = factorTS;1/3;
factorS = 1-factorTS;1/3;
if TFS == 5
    %     load('CoordinatesIT.mat');
    %     CoordinatesI=CoordinatesIT;
    load('DispFine','DispFine');
    load('GrainEps','GrainEps');
    load('GrainDisp','GrainDisp');
    load('GrainVol','GrainVol');
end
FormRVEDirichlet

%% Zero out meso grains; for this file as it is, I don't want to compute the
% the sequence has to be retained: from low to high order of grain number

locked_g = [7 8 9 12 13 14 17 18 19];
%     for j = 1:numgrain
%         locked_g(1,j) = j;
%     end
% %     for k = 1: numgrain 
%         MateT(k,4) = GrainVol(1,k)/ MateT(k,3); %Grain area
%     end 
num_locked_g = length(locked_g); 
meso_nen = 3; %Number of nodes per mose element
sub = 11; %subroutine number for a meso element
combo = 1; %Flag for plotting FS on CS
 [NodeBC,NodesOnElement,RegionOnElement,nummat,MateT,MatTypeTable,numBC,numel,ccell] = Meso_Locking_patch2(NodesOnElement,num_locked_g,NodeBC,locked_g,....
    meso_nen,grainG,nen_bulk,numelemg,GrainA,numel,numnpMicro,numnpMeso,MateT,numgrain,nummat,MatTypeTable,RegionOnElement,RegionsOnInterface);
ProbType = [numnp numel nummat 2 2 nen];