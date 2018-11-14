% Example nonlinear quasi-static input file: Cook's Membrane
% Tim Truster
% 10/06/2013
%
% Part of Batch_Cooks_Membrane
% Element subroutine: NL_Elem5_2dNSCST

NCR = 1; % Input file

h = 1/2;
p = 1;
wf = 0;
bs = 0;
rho = 1000;
n = 1/h;
pu = p;
pv = p;
nu = n;
nv = n;
nel = 3;
nen = 3;
iprob = 7;

Coordinates = [1 0 0
             2 48 44
             3 48 60
             4 0 44];
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

nummat = 1;
nen1 = nen + 1;
ProbType = [numnp numel nummat 2 3 nen];
MatTypeTable = [1
                5
                1];
% ProbType = [5 3 2 1];
AlgoType = [0; 1; 0];
OptFlag = [0 1 1 0 0 1 1]';

numBC = 6;
NodeBC = [
     1     1     0
     1     2     0
     4     1     0
     4     2     0
     7     1     0
     7     2     0];
numNodalF = 0;
numSL = 2;
SurfacesL = [
     6 3 4 1 0 1.5625 0
     9 6 8 1 0 1.5625 0];
                  
PatchMatTable = [1; 1];
young = 2.405653612292506e+002;
pois = 0.499966662199283;
MateT = [2 1 0 young pois];

% Quasi-static Proportional loading
        
% Use this for ramp loading
% s_del_a = measure of increment to proportional load
% stepmax = maximum number of steps
s_del_a = 0.25;
stepmax = 2;
mults = (s_del_a:s_del_a:s_del_a*stepmax)';
datastep = 4; % size of output arrays, can be larger than stepmax

itermax = 7;
Residratio = 10^-11;
reststep = 2; % dump data at steps equaling multiples of this #


% Error estimation
expliciterr = 0;
impliciterr = 0;
implicon = 1;
subm = 1;
residflag = 0;

restartmat = 1;
printRnorm = 1;
hrstore = 0;1;
nelPn = 1; %1 for unequal-order, 0 for equal-order