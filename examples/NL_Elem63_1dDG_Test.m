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

clear
clc
NCR = 1;

nen = 4;
Coordinates = [0
             5
             5 
             10];
% Coordinates = [0
%              2.4
%              2.6
%              4.9
%              5.1 
%              7.4
%              7.6
%              10];
NodesOnElement = [1 2 0 0
      3 4 0 0
      1 2 3 4];
RegionOnElement = [1 1 2]';
numnp = 4;
numel = 3;4;
NodeBC = [1 1 0
     4 1 0.00333333333333333];
        %  4 1 0.333333333333333];
%mults= [0.1: 1:10]
numBC = 2;
% NodeLoad = [4 1 10];
% numNodalF = 1;

numSL = 0;
numSLnp=0;

numCn = 0;
numComp = 0;
numCL = 0;
numSI = numCL;        
sigmaC=6.25; %maximum stress at the interface
delC=0.1; % maximum separation in the TSL, the rato should be related to r
nummat = 2;
nen1 = nen + 1;
ProbType = [numnp numel nummat 1 1 nen];
MatTypeTable = [1 2 
                61 63
                0 1 ];
AlgoType = [0; 1; 0];
OptFlag = [0 1 1 0 0 1 1]';

MateT = [400 0.25  1   0    
         400 400 0.25  0.25 ];
stepmax=4; %maximum step
mults= [0.25:0.25:1*stepmax];
%mults=[0.1:0.1:1*stepmax]
zeta=0; % initial displacement gap
%TSL=1; % TSL= 1 for Bilinear TSL; TSL >1 for exponential functions.
    
    
    
