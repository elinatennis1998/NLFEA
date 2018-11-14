% Tim Truster
% 11/01/2011
%
% Tension test for anisotropic material, Stein CMAME190

%    roller |-------.
%           |       |
%           |       |
%       pin o-------.

clear
clc
NCR = 1;

L = 20;
D = 20;
h = 8;
p = 1;
wf = 0;
bs = 0;
rho = 1000;
n = h;
pu = p;
pv = p;
nu = n;
nv = n;
nel = 3;
nen = 3;
iprob = 1;

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
btype = 2;3;
[Coordinates,NodesOnElement,RegionOnElement,numnp,numel] = block2d(type,rinc,sinc,node1,elmt1,mat,rskip,btype,Coordinates,nen);
Coordinates = Coordinates';
NodesOnElement = NodesOnElement';

AlgoType = [0; 1; 0];
OptFlag = [0 1 1 0 0 1 0]';

[NodeBC,BCLIndex] = getBCAT3(nu,nv);
numBC = BCLIndex(1);
[SurfacesL,numSL] = getpressAT3(nu,nv,4);

%% flag for using script-version or function-version of element routine
use_function = 0; % script 1; % function 

nummat = 1;
nen1 = nen + 1;
ProbType = [numnp numel nummat 2 2 nen];
MatTypeTable = [1
                2
                1];
MateT = [6 3 0 7 0.4 0.538 -0.0685 0.0325 0.4 45];

s_del_a = 0.25;
stepmax = 2;4;
mults = (s_del_a:s_del_a:s_del_a*stepmax)';
datastep = 4; % size of output arrays, can be larger than stepmax

itermax = 7;
Residratio = 10^-11;
reststep = 10; % dump data at steps equalting multiples of this #
