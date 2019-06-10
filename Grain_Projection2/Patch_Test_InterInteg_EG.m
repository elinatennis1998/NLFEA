%MRDG RVE
%Rectangular grains are used
%Patch Test for interface inegration.
%Elina Geut

clear
clc

PSPS = 's'; %plane stress condition 
nen = 3;4; %number of nodes per element
nel = 3;4; %max number of nodes per element
numgh = 2; %number of grains along horiz. edge
numgs = 3; %number of grains along later. edge
numgrain = 6;7;

bCrys = 2;1; %number of elem. along grain edge
tfact = 2;1; %tri 2, quad 1; 
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
btype = 2;0; %flag for meshing type (rect, triang, ...)
type = 'cart';
[Coordinates,NodesOnElement,RegionOnElement,numnp,numel] = block2d(type,rinc,sinc,node1,elmt1,mat,rskip,btype,Coordinates,nen);
Coordinates = Coordinates';
NodesOnElement = NodesOnElement';

%% Map element ID onto grain ID
%Assign individual region to each grain 
% for each grain, which elements belong to it   
% RegionOnElement = [1 1 2 2 1 1 2 2 3 4 4 5 3 4 4 5 6 6 7 7 6 6 7 7]';
% RegionOnElement = [1 1 2 2 1 1 2 2 3 4 4 5 3 4 4 5 6 6 5 5 6 6 5 5]';
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
    RegionOnElement(45:48,:) = 6;
    RegionOnElement(37:40,:) = 6;
    RegionOnElement(31:32,:) = 6;
    RegionOnElement(23:24,:) = 6;
% load('Patch_test_Integration.mat');
% BoundXmid = Patch_test_Integration;
%% Load material properties
nummat = numgrain;
MatTypeTable = [1:numgrain; ones(1,numgrain)];
m1 = [100 0.25 1]; %material type 1
m2 = [200 0.3 1]; %material type 2
m3 = [300 0.3 1]; %material type 2
% MateT = [m1; m2; m2; m3; m1; m2; m1]; %7
MateT = [m1; m2; m3; m1; m2; m3];       %6
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
TFS =  4; % Combo 3; % Sachs 1; %Taylor 2; % FE 5; % Fine Scale  
factorTS = .5;0;1;.25;.75; % factor for Taylor vs Sachs; 1 means Taylor
factorT = factorTS;1/3;
factorS = 1-factorTS;1/3;
if TFS == 1
% % factorTS = 1;
% TaylorSet
elseif TFS == 3
% factorTS = 0;
% SachsSet
elseif TFS == 4
% TaylorSet
% SachsSet
elseif TFS == 5
load('CoordinatesI.mat');
end
FormRVEDirichlet


NodeBC = [NodeBC
[(numnpMicro+1:numnpMicro+numnpMeso)'   ones(numnpMeso,1) zeros(numnpMeso,1)
(numnpMicro+1:numnpMicro+numnpMeso)' 2*ones(numnpMeso,1) zeros(numnpMeso,1)]];
numBC = length(NodeBC);

ProbType = [numnp numel nummat 2 2 nen];