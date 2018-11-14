% 4/29/2015
%
% Sample for: L_Elem7_1d
% Notes:
%
% Tim Truster
% 7/8/2013
% 1d beam tests
% Verified the element stiffness matrix for these dimensions
% Use with L_Elem7_1d.m

nen = 2;
Coordinates = [0 
             5];
NodesOnElement = [1 2];
RegionOnElement = [1]';
numnp = 2;
numel = 1;

% nen = 3;
% Coordinates = [0
%              2.5
%              5];
% NodesOnElement = [1 2 3];
% RegionOnElement = [1]';
% numnp = 3;
% numel = 1;

% nen = 4;
% Coordinates = [0
%              2
%              4
%              6];
% NodesOnElement = [1 2 3 4 1];
% RegionOnElement = [1]';
% numnp = 4;
% numel = 1;

% nen = 5;
% Coordinates = [0
%              1.2
%              2.4
%              3.6
%              4.8
%              6];
% NodesOnElement = [1 2 3 4 5 1];
% numnp = 5;
% numel = 1;

NodeBC = [1 1 0; 1 2 0];
FaceLoad = 0;
% BCLIndex = [length(NodeBC) 0]';
BCLIndex = [2 1]';
numBC = BCLIndex(1);
NodeLoad = [2 1 7];
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

nummat = 1;
nen1 = nen + 1;
ProbType = [numnp numel nummat 1 2 nen];
MatTypeTable = [1
                15
                0];
AlgoType = [0; 1; 0];
OptFlag = [0 1 1 0 0 1 0]';

MateT = [100 3];
