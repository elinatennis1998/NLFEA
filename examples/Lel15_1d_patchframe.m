% 4/29/2015
%
% Sample for: L_Elem0_1d, L_Elem7_1d
% Notes: Shows how to stack elements onto a single global material, a FEAP
% feature.
%
% Tim Truster
% 7/8/2013
% 1d bar-beam tests
% Demonstrates the combination of two elements into the same material
% Use with L_Elem0_1d.m and L_Elem7_1d.m

% nen = 2;
% Coordinates = [0 
%              5];
% NodesOnElement = [1 2 1];
% RegionOnElement = [1]';
% numnp = 2;
% numel = 1;

nen = 2;
Coordinates = [0 
             5
             9];
NodesOnElement = [1 2; 2 3];
RegionOnElement = [1 1]';
numnp = 3;
numel = 2;

NodeBC = [1 1 0; 1 2 0; 1 3 0];
FaceLoad = 0;
% BCLIndex = [length(NodeBC) 0]';
% BCLIndex = [3 2]'; %also use this with 2 elements to see rigid 2nd elem.
% NodeLoad = [2 1 7
%             2 2 4];
BCLIndex = [3 5]';
NodeLoad = [3 1 7
            3 2 4
            2 1 -2
            2 2 -15
            2 3 3];
numBC = BCLIndex(1);
numNodalF = BCLIndex(2);

intfl = 0; % flag for body force assembly
numBF = 0;
BodyForce = [1 -2];

numSL = 0;
numSLnp=0;

numCn = 0;
numComp = 0;
numCL = 0;
numSI = numCL;        

nummat = 2;
nen1 = nen + 1;
ProbType = [numnp numel nummat 1 3 nen];
MatTypeTable = [1 1
                14 15
                0 0];
AlgoType = [0; 1; 0];
OptFlag = [0 1 1 0 0 1 0]';

MateT = [100 2   %bar
         100 3]; %beam
