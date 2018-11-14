% Test for wedge element with cohesive elements.
% Domain: two wedge elements, CZ between all faces
% Loading: Nodal displacements on top face.
%
% Last revision: 12/04/2015 TJT

clear
NCR = 1;
% clc

Coordinates = [0 0 0
             2 0 0
             0 1 0
             0 0 3
             2 0 3
             0 1 3
             2 1 0
             2 1 3];
% NodesOnElement = [1 2 3 4 5 6];
% NodesOnElement = [3 1 2 6 4 5];
NodesOnElement = [2 3 1 5 6 4
                  2 7 3 5 8 6];
RegionOnElement = [1; 1];
numel = 2;
numnp = 8;
nen = 6;
NodeBC = [1 1 0
          1 2 0
          [(1:3)'; 7] 3*ones(4,1) zeros(4,1)
          2 2 0
          4 3 0.09
          5 3 0.09
          6 3 0.09
          8 3 0.09];
numBC = length(NodeBC);

MatTypeTable = [1; 1; 0];

young = 100;
pois = .25;
thick = 1;
MateT = [young pois thick];
AlgoType = [0; 1; 0];
nummat = 1;

% Output quantity flags
DHist = 1;
FHist = 1;
SHist = 0;
SEHist = 0;

% Generate CZM interface: pairs of elements and faces, duplicate the nodes,
% update the connectivities
numnpCG = numnp;
InterTypes = [1 0
              1 1]; % only put CZM between the element edges between materials 1-2
DEIProgram3

% Update boundary conditions
NodeBCCG = NodeBC;
numBCCG = numBC;
[NodeBC,numBC] = UpdateNodeSet(maxel,0,RegionOnElement,ElementsOnNodeNum,...
                               ElementsOnNode,ElementsOnNodeDup,NodeBCCG,numBCCG);

numSI = numCL;
numelCG = numel;
nen_bulk = nen;
SurfacesI = zeros(0,8);
numSI = 0;
for mat2 = 1:nummat
    for mat1 = 1:mat2
        
        matI = mat2*(mat2-1)/2 + mat1; % ID for material pair (row=mat2, col=mat1)
        if InterTypes(mat2,mat1) > 0
        numSIi = numEonF(matI);
%         SurfacesIi = squeeze(ElementsOnFacetInt(1:numSIi,1:4,matI));
%         SurfacesI(numSI+1:numSI+numSIi,:) = squeeze(ElementsOnFacetInt(1:numSIi,1:8,matI));
        locF = FacetsOnInterfaceNum(matI):(FacetsOnInterfaceNum(matI+1)-1);
        facs = FacetsOnInterface(locF);
        SurfacesIi = ElementsOnFacet(facs,:);
        SurfacesI(numSI+1:numSI+numSIi,5:8) = SurfacesIi;
        numSI = numSI + numSIi;
        [NodesOnElement,RegionOnElement,nen,numel,nummat,MatTypeTable,MateT] = ...
         FormCZ(SurfacesIi,NodesOnElement,RegionOnElement,Coordinates,numSIi,nen_bulk,3,numel,nummat,6, ...
                22,0,[5000 3 0.2],MatTypeTable,MateT);
        end
    end
end

ProbType = [numnp numel nummat 3 3 nen]; %[8 3 3 1];
