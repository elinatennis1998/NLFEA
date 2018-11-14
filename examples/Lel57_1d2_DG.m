%Elina Geut. Created 11/11/2018, last modified 11/14/2018
%The input file for 1d DG that connects two elements with different shape
%functions. The input file for same polinomials is written separatelty

clear

%Was added 11/11/2018 now can use different order polynomials
pL = 2;                                                                    
pR = 1;

if pL == pR
    error('use Lel57_1dDG input file');                          
end

NCR = 1;                                                                   %control flag
numel = 8;                                                                 %Number of elements w/ DG
ndm = 1;                                                                   %# of dimensions
ndof = 1;
ndf = ndof;
P = 10;                                                                    %External load applied 

%nel use for solids polinomial asssignments
%Note, nodes on element is doubled in size due to the Dg elements. The
%information for Dg is taken as info for the neighbor two elements
Coordinates = [0 1 2 3 4 4 7 9 9 10 12 14 15]';                            %Nodal biunit Coordinates
NodesOnElement = [1 2 3 0 0
                  3 4 5 0 0
                  6 7 0 0 0
                  7 8 0 0 0
                  9 10 11 0 0
                  11 12 13 0 0
                  3 4 5 6 7
                  7 8 9 10 11];

nen = size(NodesOnElement,2);                                              %Maximum number of nodes per elem  
NodeLoad = [13 1 P];                                                       %[node# dof mafnitude]

%Added 11/11/2018 to be able to connect different shape functions

RegionOnElement = [1 1 3 3 2 2 4 5]';
numnp = size(Coordinates,1);                                               %Number of mesh nodes 
NodeBC = [1 1 0];                                                          %Assume: u(node1)=0 fixed
numBC = size(NodeBC,1);
numNodalF = 1;                                                             %P force is applied 
intfl = 1;                                                                 %Include BF
BodyForce = [1 5                                                           %Let BF act on element 1 (q = 5)
             4 2];                                                         %Let BF act on element 4 (q = 2)
numBF = size(BodyForce,1);                                                 %Two body forces
numSL = 0;                                                                 %No surface loads
numSLnp=0;                                                                 %No non-proportional SL                                                           
numComp = 0;                                                               %Tying node command enabled                                                               
numSI = 0;                                                                 %# of surfce interfaces
numCL = numSI;
nen1 = nen + 1;

MateT = [100 3 0 0                                                         %E(1),Area(1)
         500 2 0 0                                                         %E(2),Area(2)
         200 1 0 0                                                         %E(3),Area(3)
         1 2 pL pR                                                         %DG1 element properties
         2 3 pR pL];                                                       %DG2 element properties                                                     
nummat = size(MateT,1);                                                    %Number of materials
ProbType = [numnp numel nummat ndm ndf nen];                               
AlgoType = [-1; 1; 0];                                                     %[lin.elast.sol lin el. 0]
OptFlag = [0 1 1 0 0 1 1]';                                                %Displ,Stress,Force

MatTypeTable = [1 3 2 4 5                                                  %Elements
                59 59 59 57 57                                             %Subroutines #59 (CG) #57 (DG)
                0 0 0 0 0];                                                %All elements are linear                                                                        %0-> linear; ->non-linear
                                                                           
zeta = 0;                                                                  %initial displacement gap 
stepmax = 1;                                                               %Maximum history step
