


%Quadrilateral Element 
clear
NCR = 1;
ndm = 2;
p = 1; 

nen = 4;                                                                   %Maximum the number of nodes per element 
Coordinates = [2 1 
               9 6 
               7 9
               3 7 
               10 7 
               8 10];
                                                                           %Spacial coordinates of the nodes
RegionOnElement = [1 1]';                                                  %Region ID of an Element 
NodesOnElement = [1 2 3 4
                  2 5 6 3];                                                %List of nodes connected to elements
numel = 2;                                                                 %total number of elements in the mesh
numnp = nen*numel;                                                         %total number of nodes in the mesh
ndof = 2;                 
ndf = 2;
    
NodeBC = [1 1 0
          1 2 0
          2 2 0
          5 1 0
          5 2 0];                                                          %[node# dof magnitude]

NodeLoad = [6 1 100];                                                      %[node# dof mafnitude]

numNodalF = 1;                                                             %How many external forces are applied
intfl = 0;                                                                 %if 1 - include BF; if 0 - dont 
numBF = 0;                                                                 %How many types of body forces are there
% BodyF = [1 20];  
numSL = 0;                                                                 %Number of Surface Loads
numComp = 0;
numSI = 0;        

nummat = 1;
nen1 = nen + 1;
ProbType = [numnp numel nummat 2 2 nen];
MatTypeTable = [1 
                56
                0];
AlgoType = [-1; 1; 0];
OptFlag = [0 1 1 0 0 0 0]';
%OptFlag = [0 1 1 0 0 0 1]; %will include stress
mateprop = [100 0.3 1];                                                       %E v t        
MateT = mateprop;
numBC  = length(NodeBC);    

%inputs for case 6
iprob = 0;