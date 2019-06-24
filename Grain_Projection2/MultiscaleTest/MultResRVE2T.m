% Test problem for coarse scale solid element with 3x3 grain RVE. Coarse
% scale implemented as a material parameter
% The grain in the middle is coarse scale only

clear
clc

PSPS = 's';
nen = 4;
nel = 4;
% Mesh with 6x6 tiling
numg = 3; % number of grains along an edge
numgrain = numg*numg; % total number of grains/materials/regions
bCrys = 1;3;2;6;24;12; % number of elements along edge of a grain
numelemg = bCrys*bCrys; % number of elements per grain
nu = numg*bCrys; % number of elements along x edge
nv = numg*bCrys; % number of elements along y edge

%% make the mesh of nodes and elements, centered at 0,0
Coordinates = [1 -4.5 -4.5
             2 4.5 -4.5
             4 -4.5 4.5
             3 4.5 4.5];
type = 'cart';
rinc = nu;
sinc = nv;
node1 = 1;
elmt1 = 1;
mat = 1;
rskip = 0;
btype = 0;
[Coordinates,NodesOnElement,RegionOnElement,numnp,numel] = block2d(type,rinc,sinc,node1,elmt1,mat,rskip,btype,Coordinates,nen);
Coordinates = Coordinates';
NodesOnElement = NodesOnElement';

%% generate mapping of element ID to grain ID
grainG = zeros(numgrain,numelemg);
grain = 0;
% for each grain, which elements belong to it
for k = 1:1%numg, for 3d; this is probably from CubicGrainGen
    for j = 1:numg
        for i = 1:numg
            
            grain = grain + 1;
            el = 0;
            for n = 1:1%bCrys, for 3d
                for m = 1:bCrys
                    for l = 1:bCrys
                        elem = (k-1)*nu*nv*bCrys+(j-1)*nu*bCrys+(i-1)*bCrys; % bottom-corner element of grain
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
end
% inverse map: the grain that an element belongs to
for g = 1:numgrain
RegionOnElement(grainG(g,:)) = g;
end

%% Set up phase pattern and material properties
nummat = numgrain;
MatTypeTable = [1:numgrain; ones(1,numgrain)];
matA = [100 0.25 1];
matB = [200 0.3 1];matA;
mats = [matA; matB];
% automatic alternating pattern
alterphase = 2*rem(1:numgrain,2) - 1;
alterphase(alterphase==-1) = 2;
MateT = mats(alterphase,:);
% MateT = [matA
%          matB
%          matA
%          matB
%          matA
%          matB
%          matA
%          matB
%          matA];
% MateT = [matA
%          matA
%          matA
%          matA
%          matB
%          matA
%          matA
%          matA
%          matA];
% MateT = [matB
%          matB
%          matB
%          matB
%          matA
%          matB
%          matB
%          matB
%          matB];

%% Output quantity flags
DHist = 1;
FHist = 1;
SHist = 1;
SEHist = 1;

%% separate the grains by duplicating nodes
numnpCG = numnp;
InterTypes = tril(ones(nummat),-1); % only put CZM between the element edges between materials 1-2
DEIProgram2

%% Insert DG couplers along interfaces and RVE boundary, add a node for each
% grain, make some special coordinate lists
ndm = 2;
% InterDGarray % not sure what this is used for...
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
TFS =  1; %Taylor 3; % Sachs 4; % Combo  5; % Fine Scale 2; % FE
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
% TaylorSet
SachsSet
elseif TFS == 5
load('CoordinatesIT.mat');
CoordinatesI=CoordinatesIT;
end
FormRVEDirichlet

%% Zero out meso grains; for this file as it is, I don't want to compute the
% mesoscale stuff too
micro5 = unique(reshape(NodesOnElement(grainG(5,:),1:4),bCrys^2*4,1));
l5 = length(micro5);
NodeBC = [NodeBC
[(numnpMicro+1:numnpMicro+4*3)'   ones(4*3,1) zeros(4*3,1)
(numnpMicro+1:numnpMicro+4*3)' 2*ones(4*3,1) zeros(4*3,1)]
[(numnpMicro+5*3+1:numnpMicro+numnpMeso)'   ones(4*3,1) zeros(4*3,1)
(numnpMicro+5*3+1:numnpMicro+numnpMeso)' 2*ones(4*3,1) zeros(4*3,1)]
[micro5 ones(l5,1) zeros(l5,1)
 micro5 2*ones(l5,1) zeros(l5,1)]];
numBC = length(NodeBC);
% Correction for grain 5 of mesoscale
MateT(13,5:8) = [1 0 0 1];
MateT(14,5:8) = [1 0 0 1];
MateT(16,5:8) = [0 1 1 0];
MateT(18,5:8) = [0 1 1 0];
MateT(34,1:8) = [MateT(5,1:3) 3*3 0 0 0 0];
MatTypeTable(1:3,34) = [34 11 0]';
numel = numel + 1;
NodesOnElement(numel,1:14) = [numnpMicro+4*3+1 numnpMicro+4*3+2 numnpMicro+4*3+3 zeros(1,14-3)];
RegionOnElement(numel) = 34;
nummat = nummat + 1;

ProbType = [numnp numel nummat 2 2 nen];
