% Program: InterCZallMPC
% Tim Truster
% 04/06/2017
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
%        -Interface region information: [materialID in MateTg; mat1; mat2; regI; first coupler number; last coupler number]
%         RegionsOnInterface = zeros(nummatg*(nummatg+1)/2,6);

NodesOnElementg=NodesOnElement;
RegionOnElementg=RegionOnElement;
neng=nen;
numelg=numel;
nummatg=nummat;
MatTypeTableg=MatTypeTable;
MateTg=MateT;

nummatgCG = nummatg;
numSI = numCL;
numelgCG = numelg;
neng_bulk = neng;
SurfacesI = zeros(0,8);
numSI = 0;
% Interface region information: [materialID in MateTg; mat1; mat2; regI; first coupler number; last coupler number]
RegionsOnInterface = zeros(nummatg*(nummatg+1)/2,6);

for mat2 = 1:nummatg
    for mat1 = 1:mat2
        
        matI = GetRegionsInt(mat1,mat2);
        if InterTypes(mat2,mat1) > 0
        numSIi = numEonF(matI);
        locF = FacetsOnInterfaceNum(matI):(FacetsOnInterfaceNum(matI+1)-1);
        facs = FacetsOnInterface(locF);
        SurfacesIi = ElementsOnFacet(facs,:);
        SurfacesI(numSI+1:numSI+numSIi,5:8) = SurfacesIi;
        numSI = numSI + numSIi;
        numelg_old = numelg;
        [NodesOnElementg,RegionOnElementg,neng,numelg,nummatg,MatTypeTableg,MateTg] = ...
         FormDG(SurfacesIi,NodesOnElementg,RegionOnElementg,Coordinates,numSIi,neng_bulk,ndm,numelg,nummatg,6, ...
                8,0,[0],MatTypeTableg,MateTg);
        if numelg > numelg_old
        RegionsOnInterface(nummatg-nummatgCG,:) = [nummatg mat1 mat2 matI numelg_old+1 numelg];
        end
        end
    end
end
RegionsOnInterface = RegionsOnInterface(1:nummatg-nummatgCG,:);
