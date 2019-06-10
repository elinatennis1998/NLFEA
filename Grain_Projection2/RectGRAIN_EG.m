%MRDG Rectangular RVE pattern (2x3)
%Rectangular grains are used
%Elina Geut

clear
clc

PSPS = 's'; %plane stress condition 
nen = 4;3; %number of nodes per element
nel = 4;3; %max number of nodes per element
numgh = 2; %number of grains along horiz. edge
numgs = 3; %number of grains along later. edge
bCrys = 2;1; %number of elem. along grain edge
tfact = 1;2; %tri 2, quad 1; 
numelemg = bCrys^2; %number of elements in a grain
nu = numgh*bCrys; %number of elements along x
nv = numgs*bCrys; %number of elements along y

%% Set up inputs for Block2d (RVE meshing)
Coordinates = [1 -3 -3
               2 3 -3
               3 3 3
               4 -3 3];
rinc = nu;
sinc = nv;
node1 = 1;
elmt1 = 1;
mat = 1;
rskip = 0;
btype = 0;2; %flag for meshing type (rect, triang, ...)
type = 'cart';
[Coordinates,NodesOnElement,RegionOnElement,numnp,numel] = block2d(type,rinc,sinc,node1,elmt1,mat,rskip,btype,Coordinates,nen);
Coordinates = Coordinates';
NodesOnElement = NodesOnElement';

%% Map element ID onto grain ID
%Assign individual region to each grain 
% for each grain, which elements belong to it    
RegionOnElement = [1 1 2 2 1 1 2 2 3 4 4 5 3 4 4 5 6 6 7 7 6 6 7 7]';
% RegionOnElement = [1 2 3 4 1 2 3 4 5 6 6 7 5 6 6 7 8 8 9 9 8 8 9 9]';
numgrain = max(RegionOnElement);
%% Load material properties
nummat = numgrain;
MatTypeTable = [1:numgrain; ones(1,numgrain)];
m1 = [100 0.25 1]; %material type 1
m2 = [200 0.3 1]; %material type 2
m3 = [300 0.3 1]; %material type 2
MateT = [m1; m2; m2; m3; m1; m2; m1];
% MateT = [m1; m2; m1; m2; m2; m3; m1; m1; m2];
% MateT = [m1; m2; m2; m1; m1; m2];
%% provide flags for desired outputs
OptFlag = [0 1 1 0 0 1 1]';
SEHist = 1;

%% Duplicate nodes using DEIProgram2
numnpCG = numnp;
InterTypes = tril(ones(nummat),-1); %flag for insertion of couplers
DEIProgram2  %Converts CG mesh to DG and reassignes node numbers

%% Insert couplers using InterDGallG based on InterTypes
ndm = 2;
ielI = 10; %flag for subroutine for interface elements
ielB = 10; %flag for subroutine for boundary elements
matepropI = [0 0 1 1]; %microscale 
matepropB = [0 1 1 0]; %meso-microscale
InterDGallG %Inserts couplers as assigned by InterTypes

%% Set Boundary conditions for RVE 
GrainIntegV
InterIntegV
uRVE = [0 0];
eRVE = [.02 -0.02*.25 0];
wRVE = 0;    
%note, BCs are identical to the ones in SquareRVE2F to compare results


%% Key Part: set up fluxes and jumps on the grain boundaries and RVE
%Copied from SquareRVE2F wtitten by Dr. Trister
TFS = 1; %Taylor 2; % FE 3; % Sachs 5; % Fine Scale   4; % Combo 
factorTS = 1;0;.5;.25;.75; % factor for Taylor vs Sachs; 1 means Taylor
factorT = factorTS;1/3;
factorS = 1-factorTS;1/3;
if TFS == 1
% factorTS = 1;
TaylorSet
elseif TFS == 3
% factorTS = 0;
SachsSet
elseif TFS == 4
TaylorSet
SachsSet
elseif TFS == 5
load('CoordinatesI.mat');
end
FormRVEDirichlet


NodeBC = [NodeBC
[(numnpMicro+1:numnpMicro+numnpMeso)'   ones(numnpMeso,1) zeros(numnpMeso,1)
(numnpMicro+1:numnpMicro+numnpMeso)' 2*ones(numnpMeso,1) zeros(numnpMeso,1)]];
numBC = length(NodeBC);

ProbType = [numnp numel nummat 2 2 nen];