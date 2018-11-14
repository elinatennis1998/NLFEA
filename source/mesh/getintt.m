function [NodesOnElementt,xintt,SurfacesI] = getintt(NodesOnElement,Coordinates,SurfacesI,numCL)
%
% Creates triangular segments for a conforming B8 DG mesh
NodesOnElementt = zeros(numCL*2,3);
xintt = zeros(numCL*4,3);
SurfacesI = [SurfacesI zeros(numCL,2)];

numel = size(NodesOnElement,1)-numCL;
eind = 0;
nind = 0;
for inter = 1:numCL
    nodes = NodesOnElement(numel+inter,1:4);
    eind = eind + 1;
%     NodesOnElementt(eind,:) = nodes([1 2 3]);
    NodesOnElementt(eind,:) = [1 2 3] + (inter-1)*[4 4 4];
    eind = eind + 1;
%     NodesOnElementt(eind,:) = nodes([1 3 4]);
    NodesOnElementt(eind,:) = [1 3 4] + (inter-1)*[4 4 4];
    xintt(nind+1:nind+4,:) = Coordinates(nodes,:);
    nind = nind + 4;
    SurfacesI(inter,9:10) = [eind-1 eind];
end