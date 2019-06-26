%Elina Geut
%Created 6/26/2019
%Analog to InterDGallG but for PBC case. 
%The script inserts MRDG elements on the interior and PBC elements on the
%boundary


%% InterDGallG for MRDG
nummat_MRDG = nummat;
nummat_PBC = nummat;
numSI = numCL;
numelCG = numel;
nen_bulk = nen;
SurfacesI = zeros(0,8);
% numSI = 0;
% Interface region information: [materialID in MateT; mat1; mat2; regI; first coupler number; last coupler number]
RegionsOnInterface = zeros(nummat*(nummat+1)/2,6);

% Add nodes for Grains: u, e, w
numnpMeso = (ndm+1)*nummat_MRDG;
numnpMicro = numnp;
Coordinates = [Coordinates; zeros(numnpMeso,ndm)];
numnp = numnpMeso + numnpMicro;
% Add nodes for interfaces, holds the flux and jump
CoordinatesI = zeros((ndm+2)*nummat*(nummat+1)/2,ndm);
numnpI = 0;
% Add nodes for boundaries: u, e, w
CoordinatesB = zeros((ndm+1)*nummat_MRDG*4,ndm);
numnpB = 0;

% Types of materials in domain
MaterTypeNum = [1 nummat_MRDG+1 0 0]';

for mat2 = 1:nummat_MRDG
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
        RegionsOnInterface(nummat-nummat_MRDG,:) = [nummat mat1 mat2 matI numel_old+1 numel];
        end
        end
    end
end

nummat_MRDG = nummat;
MaterTypeNum(3) = nummat+1; %Additional material for interface
numel_MRDG = numel; %Number of elements for MRDG case 

numSI = numCL;
numelCG = numel;
nen_bulk = nen;
SurfacesI = zeros(0,8);
numSI = 0;
% Arrays for new MPC links being formed
MPCListNew = zeros(0,2+ndm);
numMPCnew = 0;
CouplerNodes = zeros(0,1); % extra nodes for MPC-CZ
NodesOnLinkNew = zeros(4,numnp);
NodesOnLinknumNew = zeros(numnp,1);
% Interface region information: [materialID in MateT; mat1; mat2; regI; first coupler number; last coupler number]
% RegionsOnInterface = zeros(nummat*(nummat+1)/2,6);

for mat2 = 1:nummat_PBC
    for mat1 = 1:mat2
        
        matI = GetRegionsInt(mat1,mat2);
        numSIi = numEonPBC(matI);
            locF = FacetsOnPBCNum(matI):(FacetsOnPBCNum(matI+1)-1);
            facs = FacetsOnPBC(locF);
            SurfacesIi = ElementsOnFacet(facs,:);
            SurfacesIi = ReverseFacets(SurfacesIi,NodesOnElement,Coordinates,numSIi,ndm);
            ElementsOnFacet(facs,:) = SurfacesIi;
            SurfacesI(numSI+1:numSI+numSIi,5:8) = SurfacesIi;
            numel_old = numel;
        [NodesOnElement,RegionOnElement,nen,numel,nummat,MatTypeTable,MateT] = ...
         FormDG(SurfacesIi,NodesOnElement,RegionOnElement,Coordinates,numSIi,nen_bulk,ndm,numel,nummat,6, ...
                9,0,[0],MatTypeTable,MateT);
        if numel > numel_old
            RegionsOnInterface(nummat-nummat_PBC,:) = [nummat mat1 mat2 matI numel_old+1 numel];
        end
    end
end

RegionsOnInterface = RegionsOnInterface(1:nummat-nummat_PBC,:);
MaterTypeNum(4) = nummat+1;
