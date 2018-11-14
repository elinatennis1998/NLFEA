%Elina Geut. Last modified 11/8/2018
%Patch test for 1D subroutine. Use with L_Elem59_Elina
%Cases of p>3 have not been tested

clear
NCR = 1;
ndm = 1;
PL = 10;
p = 2;                                                                     %p=1: linear; p=2: quadratic;etc                                                                    %Length of the rod
q = 2;                                                                    %Body Force 
%This part of the code sets up the geometry/mesh for the bar at given
%conditions.

if p == 1
    Coordinates = [0 5 15 30]';                                             %Coord. -spatial coordinates of the nodes
    nen = 2;                                                               %Maximum the number of nodes per element   
    lint = 2;                                                              %# of gauss pts for integration 
    numnp = size(Coordinates,1);                                           %Total number of nodes in the mesh
    NodesOnElement = [1 2; 2 3;3 4];                                       %List of nodes connected to elements
    NodeLoad = [4 1 PL];                                                   %[node# dof mafnitude]
elseif p == 2
%     Coordinates = [0 5 10 15 20 25 30]';                                  
%     Coordinates = [0 2.5 5 10 15 22.5 30]';                              %Unequal length of elements
    Coordinates = [0 2 5 8 15 20 30]';                                   %Middle node of an elem isnt in the middle
    nen = 3;                                                              
    lint = 3;
    numnp = size(Coordinates,1); 
    NodesOnElement = [1 2 3; 3 4 5; 5 6 7];                                
    NodeLoad = [7 1 PL];                                                       
elseif p == 3
    Coordinates = [0 2 4 5 6 10 15 20 25 30]'; 
    nen = 4;
    lint = 4;
    numnp = size(Coordinates,1);
    NodesOnElement = [1 2 3 4; 4 5 6 7; 7 8 9 10];                         
    NodeLoad = [10 1 PL];                                                  
end
RegionOnElement = [1 1 1]';                                                %Region ID of an Element 
numel = 3;                                                                 %Total number of elements                 
ndf = 1;                                                                    %Number of Degrees of freedom
    
NodeBC = [1 1 0];                                                          %[node# dof magnitude]
numBC  = size(NodeBC,1);                                                   % #of BC
numNodalF = size(NodeLoad,1);                                              %How many external forces are applied
intfl = 1;                                                                 %if 1 - include BF; if 0 - dont 
BodyForce = [1 q];                                                         %Body Forces in element
numBF = size(BodyForce,1);                                                 %How many types of body forces are there
MateT = [100 1];                                                           %Material Properties
 
numSL = 0;                                                                 %Number of Surface Loads
numComp = 0;
numSI = 0;        
iprob = 0;                                                                 %Flag for case 6
nummat = 1;
nen1 = nen + 1;
ProbType = [numnp numel nummat 1 1 nen];
MatTypeTable = [1 
                59
                0];
AlgoType = [-1; 1; 0];
OptFlag = [0 1 1 0 0 1 1]'; 


