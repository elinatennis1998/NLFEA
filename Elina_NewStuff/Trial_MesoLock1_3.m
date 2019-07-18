clear

% Mesh with nxn tiling
PSPS = 's'; %plane stress condition 
nen = 4;3; %number of nodes per element
nel = 4;3; %max number of nodes per element
numgh = 3;4;5;8;7;6; %number of grains along horiz. edge
numgs = 3;4;5;8;7;6; %number of grains along later. edge
numgrain = numgh*numgs; %total number of grains in RVE
bCrys = 1;2;3;8; %number of elem. along grain edge
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
             2 3 0
             4 0 3
             3 3 3];
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
m1 = [100 0.25 1]; %material type 2
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
TFS = 1; %Taylor 2; % FE 5; % Fine Scale 3; % Sachs  4; % Combo 
factorTS = 1;0;.5;.25;.75; % factor for Taylor vs Sachs; 1 means Taylor
factorT = factorTS;1/3;
factorS = 1-factorTS;1/3;
FormRVEDirichlet

%% Zero out meso grains; for this file as it is, I don't want to compute the
% the sequence has to be retained: from low to high order of grain number
%     for j = 1:numgrain
%         locked_g(1,j) = j;
%     end
locked_g = [1 3];
num_locked_g = length(locked_g); 
meso_nen = 3; %Number of nodes per mose element
sub = 11; %subroutine number for a meso element
combo = 1; %Flag for plotting FS on CS
numelFE = numel;
for i = 1:num_locked_g
     
    micro = unique(reshape(NodesOnElement(grainG(locked_g(i),:),1:nen_bulk),numelemg*nen_bulk,1)); %Corrected
    l(i) = length(micro);
    r_g = numgrain-locked_g(i);
    NodeBC = [NodeBC
%         [(numnpMicro+1:numnpMicro+(locked_g(i)-1)*meso_nen)' ones((locked_g(i)-1)*meso_nen,1) zeros((locked_g(i)-1)*meso_nen,1)
%         (numnpMicro+1:numnpMicro+(locked_g(i)-1)*meso_nen)' 2*ones((locked_g(i)-1)*meso_nen,1) zeros((locked_g(i)-1)*meso_nen,1)]
%         [(numnpMicro+locked_g(i)*meso_nen+1:numnpMicro+numnpMeso)'   ones(r_g*meso_nen,1) zeros(r_g*meso_nen,1)
%         (numnpMicro+locked_g(i)*meso_nen+1:numnpMicro+numnpMeso)' 2*ones(r_g*meso_nen,1) zeros(r_g*meso_nen,1)]
        [micro ones(l(i),1) zeros(l(i),1)
        micro 2*ones(l(i),1) zeros(l(i),1)]];
    % Correction for grain 5 of mesoscale
    index = find(MateT(:,1) == locked_g(i));
    MateT(index(1),5:8) = [1 0 0 1];
    MateT(index(2),5:8) = [1 0 0 1];
    index2 = find(MateT(:,2) == locked_g(i));
    MateT(index2(1),5:8) = [0 1 1 0];
    MateT(index2(2),5:8) = [0 1 1 0];
    l_MT = length(MateT)+1;
    MateT(l_MT,1:8) = [MateT(locked_g(i),1:3) GrainA 0 0 0 0];
    MatTypeTable(1:3,l_MT) = [l_MT 11 0]';
    %Recised up until this point
    LockedElem = grainG(locked_g((i)),1); %Locate the corner grain in the corner element
    ConditionBC = sort(NodesOnElement(LockedElem,1:nen_bulk),2); %Find nodes for the LockedElem
    Arrange = sort(NodesOnElement(ElemInt+1:numelFE,1:nen_bulk),2);
    %                 Cond2 = sort(NodesOnElement(ElemInt+1:numel,nen_bulk+1:nen_bulk*2),2);
    Locate = find(sum(ConditionBC == Arrange,2)>=nen_bulk);
    numel = numel+1;
    NodesOnElement(numel,1:12) = [NodesOnElement(ElemInt+Locate(1),nen_bulk*2+4) NodesOnElement(ElemInt+Locate(1),nen_bulk*2+5) ...
        NodesOnElement(ElemInt+Locate(1),nen_bulk*2+6) zeros(1,12-3)];
%      numel = numel+1;
%     NodesOnElement(numel,1:12) = [NodesOnElement(ElemInt+Locate(2),nen_bulk*2+4) NodesOnElement(ElemInt+Locate(2),nen_bulk*2+5) ...
%         NodesOnElement(ElemInt+Locate(2),nen_bulk*2+6) zeros(1,12-3)];
%      numel = numel+1;
%     NodesOnElement(numel,1:12) = [numnpMicro+(locked_g(i)-1)*meso_nen+1 numnpMicro+(locked_g(i)-1)*meso_nen+2 ...
%         numnpMicro+(locked_g(i)-1)*meso_nen+3 zeros(1,12-3)];
    RegionOnElement(numel) = l_MT;
    nummat = nummat + 1;
    ccell(i) = numel;
end
     r_g = numgrain-locked_g(2); 
     NodeBC = [NodeBC
        [(numnpMicro+1:numnpMicro+(locked_g(1)-1)*meso_nen)' ones((locked_g(1)-1)*meso_nen,1) zeros((locked_g(1)-1)*meso_nen,1)
        (numnpMicro+1:numnpMicro+(locked_g(1)-1)*meso_nen)' 2*ones((locked_g(1)-1)*meso_nen,1) zeros((locked_g(1)-1)*meso_nen,1)]
        [(numnpMicro+locked_g(1)*meso_nen+1:numnpMicro+(locked_g(2)-1)*meso_nen)'   ones(1*meso_nen,1) zeros(1*meso_nen,1)
        (numnpMicro+locked_g(1)*meso_nen+1:numnpMicro+(locked_g(2)-1)*meso_nen)' 2*ones(1*meso_nen,1) zeros(1*meso_nen,1)]
        [(numnpMicro+locked_g(2)*meso_nen+1:numnpMicro+numnpMeso)'   ones(r_g*meso_nen,1) zeros(r_g*meso_nen,1)
        (numnpMicro+locked_g(2)*meso_nen+1:numnpMicro+numnpMeso)' 2*ones(r_g*meso_nen,1) zeros(r_g*meso_nen,1)]];



    numBC = length(NodeBC);
      
    ProbType = [numnp numel nummat 2 2 nen];