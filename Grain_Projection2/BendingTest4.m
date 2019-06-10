% Patch test for rectangular domain of Q4 elements with DG couplers along
% interfaces between regions. User can modify the elements belonging to
% each region.
% Domain: 2x1 rectangle
% Loading: Prescribed displacement of 0.1 on right edge.
%
% Last revision: 12/09/2015 TJT

% clear
% clc

nen = 4;
nel = 4;
% Mesh with 6x6 tiling
numg = 3;
numgrain = numg*numg;
% bCrys = 6;24;3;2;1;12;
numelemg = bCrys*bCrys;
nu = numg*bCrys;
nv = numg*bCrys;

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

grainG = zeros(numg*numg,bCrys*bCrys);
grain = 0;
for k = 1:1%numg
    for j = 1:numg
        for i = 1:numg
            
            grain = grain + 1;
            el = 0;
            for n = 1:1%bCrys
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
% grain1 = [ 1  2  3 10 11 12 19 20 21];
% grain2 = [ 4  5  6 13 14 15 22 23 24];
% grain3 = [ 7  8  9 16 17 18 25 26 27];
% grain4 = [28 29 30 37 38 39 46 47 48];
% grain5 = [31 32 33 40 41 42 49 50 51];
% grain6 = [34 35 36 43 44 45 52 53 54];
% grain7 = [55 56 57 64 65 66 73 74 75];
% grain8 = [58 59 60 67 68 69 76 77 78];
% grain9 = [61 62 63 70 71 72 79 80 81];
% grainG = [grain1
%          grain2
%          grain3
%          grain4
%          grain5
%          grain6
%          grain7
%          grain8
%          grain9];
for g = 1:numgrain
RegionOnElement(grainG(g,:)) = g;
end
% RegionOnElement(grain1) = 1;
% RegionOnElement(grain2) = 2;
% RegionOnElement(grain3) = 3;
% RegionOnElement(grain4) = 4;
% RegionOnElement(grain5) = 5;
% RegionOnElement(grain6) = 6;
% RegionOnElement(grain7) = 7;
% RegionOnElement(grain8) = 8;
% RegionOnElement(grain9) = 9;
nummat = 9;
MatTypeTable = [1 2 3 4 5 6 7 8 9; 1 1 1 1 1 1 1 1 1];
matA = [200 0.0 1];
matB = [100 0.0 1];
% MateT = [matB
%          matB
%          matB
%          matB
%          matB
%          matB
%          matB
%          matB
%          matB];
MateT = [matA
         matB
         matB
         matA
         matB
         matB
         matA
         matB
         matB];

% Boundary conditions
nodexm = find(abs(Coordinates(:,1)--4.5)<1e-9); %rollers at x=0
nodexp = find(abs(Coordinates(:,1)-4.5)<1e-9); %prescribed u_x at x=2
nodeym = find(abs(Coordinates(:,2)--4.5)<1e-9); %rollers at y=0
nodeyp = find(abs(Coordinates(:,2)-4.5)<1e-9); %rollers at y=0
nodelist = [nodexm; nodexp; nodeym; nodeyp];
numBC = length(nodelist)*2;
NodeBC = zeros(numBC,3);
bbc = 0;
for i = 1:length(nodelist)
    node = nodelist(i);
    xint = Coordinates(node,1) + 1.5; % move centerline to right
    yint = Coordinates(node,2);
    if xint < 0
    PatchE = MateT(1,1);
    Patchv = MateT(1,2);
    else
    PatchE = MateT(2,1);
    Patchv = MateT(2,2);
    end
    if abs(abs(yint)-4.5) < 1e-8 && abs(abs(xint)-0) < 1e-8
        yint
    end
    ue = uexactbb(xint,yint,PatchE,Patchv,9);
    bbc = bbc + 1;
    NodeBC(bbc,:) = [node 1 ue(1)];
    bbc = bbc + 1;
    NodeBC(bbc,:) = [node 2 ue(2)];
end
% Output quantity flags
DHist = 1;
FHist = 1;
SHist = 1;
SEHist = 1;

% Generate CZM interface: pairs of elements and faces, duplicate the nodes,
% update the connectivities
numnpCG = numnp;
InterTypes = tril(ones(nummat),0); % only put CZM between the element edges between materials 1-2
DEIProgram2

% Update boundary conditions
NodeBCCG = NodeBC;
numBCCG = numBC;
[NodeBC,numBC] = UpdateNodeSet(0,RegionOnElement,ElementsOnNodeNum,...
                               ElementsOnNode,ElementsOnNodeDup,NodeBCCG,numBCCG);
%
% Insert DG couplers
ndm = 2;
InterDGarray
CZprop = 5000000000;
% InterCZall
InterDGall

ProbType = [numnp numel nummat 2 2 nen];
