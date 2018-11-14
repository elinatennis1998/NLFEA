%Input File for DG 1D
%Elina Geut
%Last modified on 11/11/2018
%Case of linear shape functions with external force applied at the most
%right node was verified to work properly. 

%The stress projection needs to be fixed

clear

%Was added 11/11/2018 now can use different order polynomials
pL = 1;                                                                    %polynomial order for Left side
pR = 1;                                                                    %polynomial order for Left side
if pL == pR                                                                %This input file is for pL=pR
    p = pR;
else
    error('use Lel57_1d2');
end
NCR = 1;                                                                   %control flag
numel = 8;                                                                 %Number of elements w/ DG
ndm = 1;                                                                   %# of dimensions
ndof = 1;
ndf = ndof;
P = 10;                                                                    %External Force applied


%Note, nodes on element is doubled in size due to the Dg elements. The
%information for Dg is taken as info for the neighbor two elements

if p == 1
    Coordinates = [0 2 4 4 7 9 9 12 15]';                                  %Nodal biunit Coordinates
    NodesOnElement = [1 2 0 0
                      2 3 0 0
                      4 5 0 0
                      5 6 0 0
                      7 8 0 0
                      8 9 0 0
                      2 3 4 5
                      5 6 7 8];
    nen = size(NodesOnElement,2);                                          %Maximum number of nodes per elem  
    NodeLoad = [9 1 P];                                                    %[node# dof mafnitude]
    lint = 2;
    
elseif p == 2
    Coordinates = [0 1 2 3 4 4 5 7 8 9 9 12 14 15 18]';
    NodesOnElement = [1 2 3 0 0 0
                      3 4 5 0 0 0
                      6 7 8 0 0 0
                      8 9 10 0 0 0
                      11 12 13 0 0 0
                      13 14 15 0 0 0
                      3 4 5 6 7 8
                      8 9 10 11 12 13];
    nen = size(NodesOnElement,2);
    NodeLoad = [15 1 P];
    lint = 3;
end
RegionOnElement = [1 1 3 3 2 2 4 5]';
numnp = size(Coordinates,1);                                               %Number of mesh nodes 
NodeBC = [1 1 0];                                                          %Assume: u(node1)=0 fixed
numBC = size(NodeBC,1);
numNodalF = 1;                                                             %Assume no exterenal forces
intfl = 1;                                                                 %Do not include BF
BodyForce = [1 5
             4 2];                                                         %Let BF act on node2, element 1
numBF = size(BodyForce,1);                                                 %Two body forces
numSL = 0;                                                                 %No surface loads
numSLnp=0;                                                                 %No non-proportional SL                                                           
numComp = 0;                                                               %Tying node command enabled                                                               
numSI = 0;                                                                 %# of surfce interfaces
numCL = numSI;
nen1 = nen + 1;

MateT = [100 1 0 0                                                         %E(1),Area(1)
         500 1 0 0                                                         %E(2),Area(2)
         200 1 0 0                                                         %E(3),Area(3)
         1 2 pL pR                                                         %DG1 element properties
         2 3 pR pL];                                                       %DG2 element properties                                                     
nummat = size(MateT,1);                                                    %Number of materials
ProbType = [numnp numel nummat ndm ndf nen];                               
AlgoType = [-1; 1; 0];                                                     %[lin.elast.sol lin el. 0]
OptFlag = [0 1 1 0 0 1 1]';

MatTypeTable = [1 3 2 4 5                                                  %Elements
                59 59 59 57 57                                             %Subroutines #59 (CG) #57 (DG)
                0 0 0 0 0];                                                %All elements are linear  
                                                                           
zeta = 0;                                                                  %initial displacement gap 
stepmax = 1;                                                               %Maximum history step

% plotMesh1(Coordinates,NodesOnElement,[1 1 1]);

% numSI = 2;
% maxmat = 2;                                                                %Max# of DG materials    
% lin = 0;                                                                   %0-> linear; ->non-linear
% iel = 57;                                                                  %insertion element subroutin case
% SurfacesI = [3 4 2 3 2 1                                                   %[left_el. right_el. cont_surf1 cont_surf2]
%              6 7 4 5 2 1];                                                 %Elements to DG for use in MATLAB                                                    

%Get the material assignments for the DG elements. Note nargin is the # of
%input arguments for FromDG
% [NodesOnElement,RegionOnElement,nen,numel,nummat,MatTypeTable,MateT] = ...
%          FormDG(SurfacesI(:,:),NodesOnElement,RegionOnElement,...
%          Coordinates,numCL,nen,ndm,numel,nummat,maxmat,...
%                 iel,lin,0,MatTypeTable,MateT);                           
            

