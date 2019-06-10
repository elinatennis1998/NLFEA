% Program: InterDGall
% Tim Truster
% 06/30/2017
%
% Coupler insertion, standard script
% Coupler type: DG couplers on all interfaces and intrafaces defined in
% InterTypes array
%
% Input: ndm = spatial dimension
%        remainder of arrays from DEIProgram
%
% Output: all necessary arrays for FEA_Program
%         Extra arrays:
%        -Interface region information: [materialID in MateT; mat1; mat2; regI; first coupler number; last coupler number]
%         RegionsOnInterface = zeros(nummat*(nummat+1)/2,6);

nummatCG = nummat;
numSI = numCL;
numelCG = numel;
nen_bulk = nen;
SurfacesI = zeros(0,8);
numSI = 0;
% Interface region information: [materialID in MateT; mat1; mat2; regI; first coupler number; last coupler number]
RegionsOnInterface = zeros(nummat*(nummat+1)/2,6);

% Add nodes for Grains: u, e, w
numnpMeso = (ndm+1)*nummatCG;
numnpMicro = numnp;
Coordinates = [Coordinates; zeros(numnpMeso,ndm)];
numnp = numnpMeso + numnpMicro;
% Add nodes for interfaces, holds the flux and jump
CoordinatesI = zeros((ndm+2)*nummat*(nummat+1)/2,ndm);
numnpI = 0;
% Add nodes for boundaries: u, e, w
CoordinatesB = zeros((ndm+1)*nummatCG*4,ndm);
numnpB = 0;

% Types of materials in domain
MaterTypeNum = [1 nummatCG+1 0 0]';

for mat2 = 1:nummatCG
    for mat1 = 1:mat2
        
        matI = GetRegionsInt(mat1,mat2);
        if InterTypes(mat2,mat1) > 0
        numSIi = numEonF(matI);
        locF = FacetsOnInterfaceNum(matI):(FacetsOnInterfaceNum(matI+1)-1);
        facs = FacetsOnInterface(locF);
        SurfacesIi = ElementsOnFacet(facs,:);
        SurfacesI(numSI+1:numSI+numSIi,5:8) = SurfacesIi;
        numSI = numSI + numSIi;
        numel_old = numel;
        numnpI = numnpI + (ndm+2);
        [NodesOnElement,RegionOnElement,nen,numel,nummat,MatTypeTable,MateT,CoordinatesI,numnpI] = ...
         FormDGG(SurfacesIi,NodesOnElement,RegionOnElement,Coordinates,numnpMicro,numnpMeso,numSIi,nen_bulk,ndm,numel,nummat,6, ...
                ielI,0,matepropI,MatTypeTable,MateT,CoordinatesI,numnpI,1);
        if numel > numel_old
        RegionsOnInterface(nummat-nummatCG,:) = [nummat mat1 mat2 matI numel_old+1 numel];
        end
        end
    end
end

MaterTypeNum(3) = nummat+1;

for mat2 = 1:nummatCG
    for mat1 = mat2:mat2
        
        matI = GetRegionsInt(mat1,mat2);
%         if InterTypes(mat2,mat1) > 0
        numSIi = numEonB(mat2);
        SurfacesIi = squeeze(ElementsOnBoundary(1:numSIi,1:2,mat2));
        SurfacesIi = SurfacesIi(:,[1 1 2 2]);
        SurfacesI(numSI+1:numSI+numSIi,5:8) = SurfacesIi;
        [SurfacesIi,subFacets] = SplitFacets(SurfacesIi,NodesOnElement,Coordinates,numSIi,ndm);
        numSI = numSI + numSIi;
        for j = 1:size(subFacets,1)
            numel_old = numel;
            numSIii = subFacets(j,2) - subFacets(j,1) + 1;
            SurfacesIii = SurfacesIi(subFacets(j,1):subFacets(j,2),:);
            [NodesOnElement,RegionOnElement,nen,numel,nummat,MatTypeTable,MateT,CoordinatesB,numnpB] = ...
             FormDGG(SurfacesIii,NodesOnElement,RegionOnElement,Coordinates,numnpMicro,numnpMeso,numSIii,nen_bulk,ndm,numel,nummat,6, ...
                    ielB,0,matepropB,MatTypeTable,MateT,CoordinatesB,numnpB,2);
            if numel > numel_old
            RegionsOnInterface(nummat-nummatCG,:) = [nummat mat1 mat2 matI numel_old+1 numel];
            end
        end
%         end
    end
end
RegionsOnInterface = RegionsOnInterface(1:nummat-nummatCG,:);
MaterTypeNum(4) = nummat+1;

% Node number for each group of nodes: micro, meso, and boundary
NodeTypeNum = [1 numnpMicro+1 numnpMicro+numnpMeso+1 numnpMicro+numnpMeso+numnpB+1]';
Coordinates = [Coordinates; CoordinatesB(1:numnpB,1:ndm)];
numnp = numnp + numnpB;
