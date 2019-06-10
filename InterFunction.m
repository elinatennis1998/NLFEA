function [NodesOnElement,RegionOnElement,Coordinates,numnp,nen,numel,nummat,MateT,MatTypeTable,RegionsOnInterface,NodeTypeNum,numMPC,MPCList...
] =InterFunction(couplertype,InterTypes,NodesOnElement,RegionOnElement,Coordinates,numnp,numel,nummat,nen,ndm,Input_data,CZprop,MateT,MatTypeTable,...
usePBC,numMPC,MPCList,...
pencoeff,CornerXYZ)

switch nargin
    case {1,2,3,4,5,6,7,8,9,10}
        error('Must supply 11 arguments')
    case 11
        iel=0;nonlin=0;mateprop=0;
        MatTypeTable = [1:nummat; ones(1,nummat); zeros(1,nummat)];
        MateT = ones(nummat,1)*[100 .25 1];
        CZprop = 0;
    case {12,13}
        error('Must supply 14 arguments')
    case 14
        usePBC = 0;
        numMPC = 0;
        MPCList = [];
        pencoeff = [];
        CornerXYZ = [];
    case {15,16,17,18}
        error('Must supply 19 arguments')
    otherwise
end

% put a select branch here based on an integer input, to change the types


numEonB = Input_data.numEonB;
numEonF = Input_data.numEonF;
ElementsOnBoundary = Input_data.ElementsOnBoundary;
numSI = Input_data.numSI;
ElementsOnFacet = Input_data.ElementsOnFacet;
ElementsOnNode = Input_data.ElementsOnNode;
ElementsOnNodeDup = Input_data.ElementsOnNodeDup;
ElementsOnNodeNum = Input_data.ElementsOnNodeNum;
numfac = Input_data.numfac;
ElementsOnNodeNum2 = Input_data.ElementsOnNodeNum2;
numinttype = Input_data.numinttype;
FacetsOnElement = Input_data.FacetsOnElement;
FacetsOnElementInt = Input_data.FacetsOnElementInt;
FacetsOnInterface = Input_data.FacetsOnInterface;
FacetsOnInterfaceNum = Input_data.FacetsOnInterfaceNum;
FacetsOnNode = Input_data.FacetsOnNode;
FacetsOnNodeCut = Input_data.FacetsOnNodeCut;
FacetsOnNodeInt = Input_data.FacetsOnNodeInt;
FacetsOnNodeNum = Input_data.FacetsOnNodeNum;
NodeCGDG = Input_data.NodeCGDG;
NodeReg = Input_data.NodeReg;
NodesOnElementCG = Input_data.NodesOnElementCG;
NodesOnElementDG = Input_data.NodesOnElementDG;
NodesOnInterface = Input_data.NodesOnInterface;
NodesOnInterfaceNum = Input_data.NodesOnInterfaceNum;
numCL = Input_data.numCL;
% arrays for multi-point constraints
if usePBC
NodesOnPBC = Input_data.NodesOnPBC;
NodesOnPBCnum = Input_data.NodesOnPBCnum;
NodesOnLink = Input_data.NodesOnLink;
NodesOnLinknum = Input_data.NodesOnLinknum;
numEonPBC = Input_data.numEonPBC;
FacetsOnPBC = Input_data.FacetsOnPBC;
FacetsOnPBCNum = Input_data.FacetsOnPBCNum;
FacetsOnIntMinusPBC = Input_data.FacetsOnIntMinusPBC;
FacetsOnIntMinusPBCNum = Input_data.FacetsOnIntMinusPBCNum;
end

switch couplertype
    case 1
        InterCZall
        NodeTypeNum = 0;
    case 2
        InterDGall
        NodeTypeNum = 0;
    case 3
        InterCZallMPC
    case 4
        InterCZallMPCsome
end