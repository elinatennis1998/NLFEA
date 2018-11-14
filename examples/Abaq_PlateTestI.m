% Tim Truster
% 09/27/2014
%
% General fracture problem, plate under velocity. From Nyugen paper
% First test to see that pure elastic solution works, passes the patch
% test.

% clear
% clc
NCR = 1;

% Read in the .inp file
Abaqfile = 'PlateTest.inp';
AbaqusInputReader

% BCs
NodeBCholder{3}(:,3) = 1;
NodeBCholder{4}(:,3) = 1;
NodeBC = [NodeBCholder{1}; NodeBCholder{2}; NodeBCholder{3}; NodeBCholder{4}];
numBC = size(NodeBCholder{1},1) + size(NodeBCholder{2},1) + size(NodeBCholder{3},1) + size(NodeBCholder{4},1);

iprob = 0;
nen1 = nen + 1;
MateT = [1 1 1]'*[190e3 0.3 1];
MatTypeTable = [MatTypeTable; [0 0 0]];
AlgoType = [0; 1; 0];
OptFlag = [0 1 1 0 0 1 1 1]';

% Generate CZM interface: pairs of elements and faces, duplicate the nodes,
% update the connectivities
numnpCG = numnp;
InterTypes = [0 0 0
              1 0 0
              1 1 0]; % only put CZM between the element edges between materials 1-2
DEIProgram2

% Update boundary conditions
NodeBCCG = NodeBC;
numBCCG = numBC;
[NodeBC,numBC] = UpdateNodeSet(maxel,0,RegionOnElement,ElementsOnNodeNum,...
                               ElementsOnNode,ElementsOnNodeDup,NodeBCCG,numBCCG);

% Generate DG elements: create new elements in ix with format expected by
% Matlab
numSI = numCL;
numelCG = numel;
nen_bulk = nen;
for mat2 = 1:3
    for mat1 = 1:3
        
        matI = mat2*(mat2-1)/2 + mat1; % ID for material pair (row=mat2, col=mat1)
        if InterTypes(mat2,mat1) > 0
        numSIi = numEonF(matI);
        locF = FacetsOnInterfaceNum(matI):(FacetsOnInterfaceNum(matI+1)-1);
        facs = FacetsOnInterface(locF);
        SurfacesIi = ElementsOnFacet(facs,:);
        [NodesOnElement,RegionOnElement,nen,numel,nummat,MatTypeTable,MateT] = ...
         FormDG(SurfacesIi,NodesOnElement,RegionOnElement,...
         Coordinates,numCL,nen_bulk,2,numel,nummat,1,...
                5,0,0,MatTypeTable,MateT);
        end
    end
end       

ProbType = [numnp numel nummat 2 2 nen];

stepmax = 1;
s_del_a = 1;%/21;stepmax;
mults = (s_del_a:s_del_a:s_del_a*stepmax)';
datastep = stepmax; % size of output arrays, can be larger than stepmax

itermax = 7;
Residratio = 10^-11;
reststep = 10; % dump data at steps equalting multiples of this #
