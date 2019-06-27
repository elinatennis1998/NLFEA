%Elina Geut
%Input file for MRDG using Neper Mesh Generator 10
%2/5/2019 Created
%2/7/2019 Last Modified
%Attempt to generalize the input file 

% Plotting commands are represented below
% plotMesh2(Coordinates,NodesOnElement,1,(1:size(NodesOnElement,1)),[1 0 0 0 0])
% plotElemCont2(Coordinates,StreList(1,1:100),NodesOnElement(1:100,1:4),1,(1:size(NodesOnElement,1)),[1 0 0])
%caxis([range]); Color bar range 

clear
clc

PSPS = 's'; %plane stress condition 
meshN = 1;0; %flaf for using Neper mesh generator

%% Mesh will be generated using Neper file

if meshN == 1
    Abaqfile = 'Neper1_EG.inp';
    ndm = 2;
    numgrain = 10;
    AbaqusInputReader
    
else
    nen = 4;
    nel = 4;
    numgh = 3; %number of grains along horiz. edge
    numgs = 2; %number of grains along later. edge
    numgrain = numgh*numgs; %total number of grains in RVE
    numelemg = 2;1; %number of elements in a grain
    bCrys = 2; %number of elem. along grain edge
    nu = numgh*bCrys; %number of elements along x
    nv = numgs*bCrys; %number of elements along y
    Coordinates = [1 -3 -2
                   2 3 -2
                   3 3 2
                   4 -3 2];
    node1 = 1;
    elmt1 = 1;
    mat = 1;
    rskip = 0;
    btype = 0;2; %flag for meshing tyoe (rect, triang, ...)
    type = 'cart';
    [Coordinates,NodesOnElement,RegionOnElement,numnp,numel] = block2d(type,nu,nv,node1,elmt1,mat,rskip,btype,Coordinates,nen);
    Coordinates = Coordinates';
    NodesOnElement = NodesOnElement';

end 



%% Map element ID onto grain ID
%Assign individual region to each grain 
% inverse map: the grain that an element belongs to
if meshN == 0
    grainG = zeros(numgrain,numelemg);
    grain = 0;
    % for each grain, which elements belong to it
    for j = 1:numgs
        for i = 1:numgh
            grain = grain + 1;
            el = 0;
            for n = 1:1%bCrys, for 3d
                for m = 1:bCrys
                    for l = 1:bCrys
                        elem = (j-1)*nu*bCrys+(i-1)*bCrys; % bottom-corner element of grain
                        elem = elem + (n-1)*nu*nv;
                        elem = elem + (m-1)*nu;
                        elem = elem + (l-1) + 1;
                        el = el + 1;
                        grainG(grain,el) = elem;
                    end
                end
            end
            
        end
    end
    % inverse map: the grain that an element belongs to
    for g = 1:numgrain
        RegionOnElement(grainG(g,:)) = g;
    end
    
    %% Load material properties
    nummat = numgrain;
    MatTypeTable = [1:numgrain; ones(1,numgrain)];
    m1 = [100 0.25 1]; %material type 1
    m2 = [200 0.3 1]; %material type 2
    mats = [m1; m2];
    % automatic alternating pattern
    alterphase = 2*rem(1:numgrain,2) - 1;
    alterphase(alterphase==-1) = 2;
    MateT = mats(alterphase,:);
else
    load('RegionOnElement_N1_EG.mat');
    m1 = [100 0.25 1]; 
    m2 = [200 0.3 1];
    m3 = [300 0.25 1]; 
    m4 = [150 0.25 1];
    nummat = numgrain;
    MatTypeTable = [1:numgrain; ones(1,numgrain)];
    MateT = [m3; m2; m1; m3; m2; m1; m1; m1; m3; m4];
%     MateT = [m1; m2; m1; m2; m1; m2; m1; m2; m1; m2];
   for i = 1:numel
       numEonB(i,1) = (i);
   end 
end
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
%% Key Part: set up fluxes and jumps on the grain boundaries and RVE
%Copied from SquareRVE2F wtitten by Dr. Trister
TFS = 1; %Taylor 3; % Sachs 5; % Fine Scale  2; % FE 4; % Combo  
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
load('CoordinatesI.mat');
end
FormRVEDirichlet

ProbType = [numnp numel nummat 2 2 nen];
