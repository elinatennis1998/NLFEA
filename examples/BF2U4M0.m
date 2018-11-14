% Tim Truster
% 08/01/2011
%
% Patch test, two elements for NL_Elem_2dd, axial loading, disp control

%    roller |-------.
%           |       |
%           |       |
%       pin o-------.

clear
clc
NCR = 1;

L = 1;
D = 1;
h = D;
p = 1;
wf = 0;
bs = 0;
rho = 1000;
n = 1;2;
pu = p;
pv = p;
nu = n;
nv = n;
nel = 4;
nen = 4;
iprob = 6;

Coordinates = [1 0 0
             2 L 0
             4 0 D
             3 L D];
type = 'cart';
rinc = nu;
sinc = nv;
node1 = 1;
elmt1 = 1;
mat = 1;
rskip = 0;
btype = 0;
[Coordinates,NodesOnElement,RegionOnElement,numnp,numel] = block2d(type,rinc,sinc,node1,elmt1,mat,rskip,btype,Coordinates,nen);
Coordinates = Coordinates';
NodesOnElement = NodesOnElement';

intfl = 1; % flag for body force assembly
numBF = 1;
BodyForce = [1 0 0 0];
nummat = 1;
nen1 = nen + 1;
ProbType = [numnp numel nummat 2 2 nen];
MatTypeTable = [1
                2
                1];
% ProbType = [5 3 2 1];
AlgoType = [0; 1; 0];
OptFlag = [0 1 1 0 0 1 1]';

[NodeBC,BCLIndex] = getBCBF3(nu,nv);
numBC = BCLIndex(1);
[SurfacesL,numSL] = getpressBF4(nu,nv);
PatchMatTable = [1; 1];
young = 100;
pois = .25;
% young = young*(1 + 2*pois)/(1 + pois)^2;
% pois = pois/(1 + pois);
MateT = [1 3 0 young pois];
mults = 1;
stepmax = 1;
% NodeLoad = [2,1,4.1
%             4,1,4.1];

%% flag for using script-version or function-version of element routine
use_function = 0; % script 1; % function 

% Results of study:
% Files: Pure-Displacement - BF2U4M0.m,BF2U6M0.m,BF2U9M0.m,NL_Elem2_2dM.m
%                            (Linear T3 can't reproduce the solution)
%        Basic-Mixed - BF5U4M0.m,BF5U6M0.m,BF5U9M0.m,NL_Elem5_2dM.m
%        Non-Symm-Stabilized - NL_Elem5_2dNS.m
% Material Model: S(F) = mu*(I - C^-1) + lam*(J-1)*J*C^-1;
%
% Using lint = 10 for edge integration and (100,25,100) for (Q4,T6,Q9)
% interior integration, the pure-displacement and basic-mixed elements
% reproduced the exact solution x = X; y = 1.01*X*Y + Y; z = Z to nearly
% numerical precision using 1x1 or 2x2 for Q4 and 1x1 (1x2) for Q9 (T6).
%
% With the same integration schemes (although the Q4 results were
% insensitive to the interior integration), the Non-Symm-Stabilized did not
% give the exact solution but was close to it, and the results converged as
% the mesh was refined for all three element types. (I did not evaluate the
% convergence rates).
%
% This proves that the formulation without the second-derivatives is not
% consistent for non-constant stress fields or non-zero body force.
% However, the inconsistency is of the order of the mesh resolution.
%
% Node_U_V =
%                    0                   0   0.000000000000002
%                    0                   0  40.400000000000013
%    0.000000000000000   0.000000000000000  -0.000000000000009
%   -0.000000000000000   1.010000000000000  40.400000000000006
%
% 01/08/2012 - NL_Elem5_2d4.m created from NL_Elem5_2d2.m; verified to give
% same results as NL_Elem5_2d2 for theta=log(JxX). Using BF5U4L0.m, the
% exact solution was reproduced by the symmetric stabilized mixed method to
% machine precision for 1 or 2x2 elements; thus, that method is consistent
% like the original mixed weak form. The 8-th order tensor emat was not
% properly calculated yet, as about 25 iterations are required for
% convergence, but nonetheless the residual vector (RHS) is correct and
% thus allows the proper solution to be reached. A non-symmetric element
% routine has not been created/debugged yet, but this element should also
% produce the exact solution.
% Later today, the correct epmat was used; this provided quadratic
% convergence to the exact solution.
