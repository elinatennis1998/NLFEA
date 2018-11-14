% 4/29/2015
%
% Sample for: NL_Elem8_3d
% Notes: shows reordering of local degrees of freedom on element, a FEAP
% feature.
%
% 4/17/13
% Tim Truster
% first sample for new program format, solves a two block "friction"
% problem where there is a linear tangential stiffness and contact in the
% normal direction. This is for verifying that multiple element types can
% be inputted into the code.

NCR = 1;

iprob = 0;
nen = 8;
nel = 8;
xl = [1 0 0 0
      2 1 0 0
      4 0 1 0
      3 1 1 0
      5 0 0 1
      6 1 0 1
      8 0 1 1
      7 1 1 1];
[Coordinates,NodesOnElement,RegionOnElement,numnp,numel] = block3d('cart',2,2,1,1,1,1,10,xl,nen);
xl = [1 1 0 0
      2 2 0 0
      4 1 1 0
      3 2 1 0
      5 1 0 1
      6 2 0 1
      8 1 1 1
      7 2 1 1];
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
SurfacesI = [0 0 0 0 2 5 2 1
             0 0 0 0 4 7 2 1];
         
SurfacesLnp = [0 0 6 2 -press 0 0
               0 0 8 2 -press 0 0];
numCn = 0;
numComp = 0;
numCL = 2;
numSI = numCL;
numSLnp = 2;
numD = 0;

[SurfacesI,ixt] = getixt(SurfacesI,numCL,NodesOnElement,8);

iprob = 0;
fmu = 0.3;
eN = 1e5;
eT = 1e8;%1e13;
% eN = 1e10;
% eT = 1e10;
sigmaI = 2.677e-3;
xin = sigmaI*sqrt(6);
eta = 2.91e-4*10^6;
R = 30.14e-3;
E = 200*10^3;
v = .24;
Estar = E/(2*(1-v^2));
Gstar = 2*(1-v)/fmu/(2-v);

MatTypeTable = [1; 1; 0; 3; 2; 1; 0];
AlgoType = [0; 1; 0];
OptFlag = [0 1 1 0 0 1 0]';

young = 200e3;
pois = .24;
thick = 1;
MateT = [young pois thick];

[NodesOnElement,RegionOnElement,nen,numel,nummat,MatTypeTable,MateT] = ...
         FormDG(SurfacesI(:,5:8),NodesOnElement,RegionOnElement,...
         Coordinates,numCL,nen,3,numel,1,1,...
                9,1,0,MatTypeTable,MateT);
% % Calls it twice to make sure that interfaces can be appended
% [NodesOnElement,RegionOnElement,nen,numel,nummat,MatTypeTable,MateT] = ...
%          FormDG(SurfacesI(1,5:8),NodesOnElement,RegionOnElement,...
%          Coordinates,1,nen,3,numel,1,1,...
%                 9,1,0,MatTypeTable,MateT);
% [NodesOnElement,RegionOnElement,nen,numel,nummat,MatTypeTable,MateT] = ...
%          FormDG(SurfacesI(2,5:8),NodesOnElement,RegionOnElement,...
%          Coordinates,1,8,3,numel,1,1,...
%                 9,1,0,MatTypeTable,MateT);
numCL = 0;
ProbType = [numnp numel nummat 3 4 nen]; %[8 3 3 1];

% plotModelCont3(Coordinates+[zeros(numnp,1) 10000*Node_U_V(:,2:3)], Node_U_V(:,3), NodesOnElement, numel, nen, 2, 1, 1, '',(1:numnp)')