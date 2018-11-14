% Tim Truster
% 03/28/2015
% Reset data for use in gensubmesh

numnp = ProbType(1);
numel = ProbType(2);
nummat = ProbType(3);
ndm = ProbType(4);
ndf = ProbType(5);
nen = ProbType(6);
residflag=1;
minc = 2; % Number of Sub-Cells in each direction
celn = (2*minc+1)^3; %(minc+1)*(minc+1)! Number of Nodes in Patch

nstar = 2; %(m+1)*(m+1)*ndf	% Number of DOFs in Patch

cel = minc^3; %(m-1)**2			% Interior nodes in Patch
ncel =(minc+1)^3;

maxel = 12;
ndfs = 3;
strwek = residflag;
strong = strwek;
x = Coordinates';
NodesOnElement = NodesOnElement';
nen1 = nen+1;
%         u = Node_U_V';
d = MateT';

gensubmesh3d

NodesOnElement = NodesOnElement';

% New mesh is:
Coordinates2 = xs';
ix2 = ixs'; % And, add the material id list to last 
mateIDs = RegionOnElement;
newmateIDs = reshape(mateIDs*ones(1,8),numel*8,1);
ix2 = [ix2 newmateIDs];
