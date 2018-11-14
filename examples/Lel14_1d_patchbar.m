% 4/29/2015
%
% Sample for: L_Elem0_1d
% Notes:
%
% Tim Truster
% 7/8/2013
% 1d rod tests
% verified p=1,2,3 for equal spacing; p>3 does not work yet
% Use with L_Elem0_1d.m

p = 2;
if p == 1
nen = 2;
Coordinates = [0
    2.5
             5];
NodesOnElement = [1 2
      2 3];
RegionOnElement = [1 1]';
numnp = 3;
numel = 2;

elseif p == 2
nen = 3;
Coordinates = [0
    1.25
             2.5
             3.75
             5];
NodesOnElement = [1 2 3
      3 4 5];
RegionOnElement = [1 1]';
numnp = 5;
numel = 2;
elseif p == 3
nen = 4;
Coordinates = [0
             5/3
             10/3
             5];
NodesOnElement = [1 2 3 4 ];
RegionOnElement = [1]';
numnp = 4;
numel = 1;
end
% nen = 5;
% Coordinates = [0
%              1.2
%              2.4
%              3.6
%              4.8
%              6];
% NodesOnElement = [1 2 3 4 5];
% RegionOnElement = [1]';
% numnp = 5;
% numel = 1;

NodeBC = [1 1 0];
FaceLoad = 0;
% BCLIndex = [length(NodeBC) 0]';
numBC = 1;
if p == 1
NodeLoad = [3 1 7];
elseif p == 2
NodeLoad = [5 1 7];
elseif p == 3
NodeLoad = [4 1 7];
end

numNodalF = 0;

intfl = 1; % flag for body force assembly
numBF = 1;
BodyForce = [1 -2];

numSL = 0;
numSLnp=0;

numCn = 0;
numComp = 0;
numCL = 0;
numSI = numCL;        

nummat = 1;
nen1 = nen + 1;
ProbType = [numnp numel nummat 1 1 nen];
MatTypeTable = [1
                14
                0];
AlgoType = [0; 1; 0];
OptFlag = [0 1 1 0 0 0 1]';

MateT = [100 3];

