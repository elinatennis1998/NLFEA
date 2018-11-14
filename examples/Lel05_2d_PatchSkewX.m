% 4/29/2015
%
% Sample for: L_Elem1_3dDG2
% Notes:
%
% 10/15/13
% Tim Truster
% Modified version of CohShearAxialTestM3.m that utilizes L_Elem1_3DG2.m
% for the interface integration. Each square interface DG element is split
% into 2 triangles. This was verified against CohShearAxialTestM3.m for the
% same results, and thus shows that the new DG routine can handle
% interfaces in other orientations than the z-plane.

% clear
% clc
NCR = 1;

iprob = 0;
nen = 8;
nel = 8;
zL = 1;0.1;
yL = zL;0.1;
a = 0.1;0; % Note: surface must remain flat; moving one noded doesn't make it flat
b = a;0;0.1; % second node must move the same amount at node a
% The following cases give correct results with truncated sectors:
% a=b=0; yL=zL=1;
% a=b=0.1; yL=zL=1;
% a=b=0; yL=zL=0.1;
% a=0.1;b=0; yL=zL=0.1;
% a=0.1;b=0; yL=zL=1;
% a=0;b=0.1; yL=zL=0.1;

xl = [1 0 0 0
      2 1+a+b 0 0
      4 0 yL 0
      3 1+b yL 0
      5 0 0 zL
      6 1+a 0 zL
      8 0 yL zL
      7 1 yL zL];
[Coordinates,NodesOnElement,RegionOnElement,numnp,numel] = block3d('cart',2,2,1,1,1,1,10,xl,nen);
xl = [1 1+a+b 0 0
      2 2 0 0
      4 1+b yL 0
      3 2 yL 0
      5 1+a 0 zL
      6 2 0 zL
      8 1 yL zL
      7 2 yL zL];
[Coordinates,NodesOnElement,RegionOnElement,numnp,numel] = block3d('cart',2,2,1,numnp+1,numel+1,1,10,xl,nen,Coordinates,NodesOnElement,RegionOnElement);
Coordinates = Coordinates';
NodesOnElement = NodesOnElement';

NodeBC = [1 1 0; 1 2 0; 1 3 0; 4 1 0; 7 1 0; 10 1 0; 13 1 0; 16 1 0; 10 2 0];
FaceLoad = 0;
numBC = length(NodeBC);
NodeLoad = 0;

% Surface Tractions
press = 50;
shear = 5;
SurfacesL = [0 0 1 3 shear 0 0
             0 0 2 3 shear 0 0
             0 0 5 3 shear 0 0
             0 0 6 3 shear 0 0
             0 0 6 2 0 -shear 0
             0 0 8 2 0 -shear 0
             0 0 8 4 -shear 0 0
             0 0 7 4 -shear 0 0
             0 0 4 4 -shear 0 0
             0 0 3 4 -shear 0 0
             0 0 3 1 0 shear 0
             0 0 1 1 0 shear 0];
SurfacesL = [SurfacesL
             0 0 1 6 shear 0 0
             0 0 2 6 shear 0 0
             0 0 3 6 shear 0 0
             0 0 4 6 shear 0 0
             0 0 5 6 shear 0 0
             0 0 6 6 shear 0 0
             0 0 7 6 shear 0 0
             0 0 8 6 shear 0 0
             0 0 1 5 -shear 0 0
             0 0 2 5 -shear 0 0
             0 0 3 5 -shear 0 0
             0 0 4 5 -shear 0 0
             0 0 5 5 -shear 0 0
             0 0 6 5 -shear 0 0
             0 0 7 5 -shear 0 0
             0 0 8 5 -shear 0 0
             0 0 6 2 0 0 shear
             0 0 8 2 0 0 shear
             0 0 3 1 0 0 -shear
             0 0 1 1 0 0 -shear];
numSL = 32;12;

% Interface
SurfacesI = [0 0 0 0 2 5 2 1 1 2
             0 0 0 0 4 7 2 1 3 4];
% ixt = [ 3  6 12
%         6 12 15
%         6  9 15
%         9 15 18]; % Actual mesh nodes
ixt = [ 1  2 4
        2 4 5
        2  3 5
        3 5 6]; % Triangle nodes
xintt = [Coordinates(3,:)
         Coordinates(6,:)
         Coordinates(9,:)
         Coordinates(12,:)
         Coordinates(15,:)
         Coordinates(18,:)];
         
SurfacesLnp = [0 0 6 2 -press 0 0
               0 0 8 2 -press 0 0];
numCn = 0;
numComp = 0;
numCL = 2;
numSI = numCL;
numSLnp = 2;
numD = 0;

iprob = 0;

MatTypeTable = [1; 1; 0];
AlgoType = [0; 1; 0];
OptFlag = [0 0 0 0 0 0 0]';

young = 200e3;
pois = .24;
thick = 1;
MateT = [young pois thick];

[NodesOnElement,RegionOnElement,nen,numel,nummat,MatTypeTable,MateT] = ...
         FormDG(SurfacesI(:,5:8),NodesOnElement,RegionOnElement,...
         Coordinates,numCL,nen,3,numel,1,1,...
                5,0,0,MatTypeTable,MateT);
% % Calls it twice to make sure that interfaces can be appended
% [NodesOnElement,RegionOnElement,nen,numel,nummat,MatTypeTable,MateT] = ...
%          FormDG(SurfacesI(1,5:8),NodesOnElement,RegionOnElement,...
%          Coordinates,1,nen,3,numel,1,1,...
%                 9,0,0,MatTypeTable,MateT);
% [NodesOnElement,RegionOnElement,nen,numel,nummat,MatTypeTable,MateT] = ...
%          FormDG(SurfacesI(2,5:8),NodesOnElement,RegionOnElement,...
%          Coordinates,1,8,3,numel,1,1,...
%                 9,0,0,MatTypeTable,MateT);
numCL = 0;
ProbType = [numnp numel nummat 3 3 nen]; %[8 3 3 1];

% plotModelCont3(Coordinates, 0*Coordinates, ix, numel-numSI, nen, 2, 1, 1, '',(1:numnp)')