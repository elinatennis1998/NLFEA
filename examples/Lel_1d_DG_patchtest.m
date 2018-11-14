



clear

p = 1;                                                                     %Linear shape finctions
NCR = 1;                                                                   %control flag
nen = 4;                                                                   %Maximum number of nodes per elem
numel = 3;                                                                 %Number of elements w/ DG
ndm = 1;                                                                   %# of dimensions
ndof = 1;
ndf = ndof;
Coordinates = [0 2 2 5]';                                      %Nodal biunit Coordinates

%Note, nodes on element is doubled in size due to the Dg elements. The
%information for Dg is taken as info for the neighbor two elements

NodesOnElement = [1 2 0 0
                  3 4 0 0
                  1 2 3 4];
RegionOnElement = [1 1 2]';
numnp = size(Coordinates,1);                                                                 %Number of mesh nodes 
NodeBC = [1 1 0];                                                          %Assume: u(node1)=0 fixed
numBC = size(NodeBC,1)
numNodalF = 1;                                                             %Assume no exterenal forces
NodeLoad = [4 1 10];                                                       %[node# dof mafnitude]
intfl = 0;                                                                 %Do not include BF
BodyForce = [1 10                                                          %Let BF act on node2, element 1
             3 10];                                                        %Also BF on node3, element 
numBF = size(BodyForce,1);                                                 %Two body forces
numSL = 0;                                                                 %No surface loads
numSLnp=0;                                                                 %No non-proportional SL                                                           
numComp = 0;                                                               %Tying node command enabled                                                               
numSI = 0;                                                                 %# of surfce interfaces
numCL = numSI;
nen1 = nen + 1;

MateT = [100 1 0 0                                                         %E(1),Area(1)
         100 1 0 0                                                         %E(2),Area(2)                                                           %DG1 element properties
         2 3 2 2];                                                         %DG1 element properties                                                     
nummat = size(MateT,1);                                                    %Number of materials
ProbType = [numnp numel nummat ndm ndf nen];                               
AlgoType = [-1; 1; 0];                                                     %[lin.elast.sol lin el. 0]
OptFlag = [0 1 1 0 0 0 0]';

MatTypeTable = [1 1 2                                                      %Elements
                59 59 57                                                      %Subroutines #59 (CG) #57 (DG)
                0 0 0];                                                     %All elements are linear  
                                                                           %0-> linear; ->non-linear                                                                          
zeta = 0;                                                                  %initial displacement gap 
stepmax = 1;                                                               %Maximum history step  
