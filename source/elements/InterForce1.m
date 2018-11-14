% Tim Truster
% 11/26/2013
%
% Evaluate interactive force, simple model

% get values from other side
switch nreg
    case 2
        elem2 = elem - numel/2;
        ma2 = RegionOnElement(elem2);
        ElemFlag2 = NodesOnElement(elem2,1:nen);
        actnode2 = find(ElemFlag2>0);
        xl2 = zeros(ndm,nen);
        xl2(1:ndm,actnode2) = Coordinates(ElemFlag2(actnode2),1:ndm)';
[EGDOFTa2, EGDOFTi2, ELDOFTa2, ELDOFTi2] = plocal(NDOFT1, ElemFlag2, squeeze(iedof(:,:,ma2)), actnode2, nen, ndf, neq1);
        ul2 = zeros(ndf, nen);
        ul2(ELDOFTa2) = ModelDx1(EGDOFTa2)';
        ul2(ELDOFTi2) = gBC1(EGDOFTi2)';
    case 1
        elem2 = elem + numel/2;
        ma2 = RegionOnElement(elem2);
        ElemFlag2 = NodesOnElement(elem2,1:nen);
        actnode2 = find(ElemFlag2>0);
        xl2 = zeros(ndm,nen);
        xl2(1:ndm,actnode2) = Coordinates(ElemFlag2(actnode2),1:ndm)';
[EGDOFTa2, EGDOFTi2, ELDOFTa2, ELDOFTi2] = plocal(NDOFT2, ElemFlag2, squeeze(iedof(:,:,ma2)), actnode2, nen, ndf, neq2);
        ul2 = zeros(ndf, nen);
        ul2(ELDOFTa2) = ModelDx2(EGDOFTa2)';
        ul2(ELDOFTi2) = gBC2(EGDOFTi2)';
    case 0 % concurrent formation of constituent contributions
        if elem > numel_mix/2
        elem2 = elem - numel_mix/2;
        ma2 = RegionOnElement(elem2);
        ElemFlag2 = NodesOnElement(elem2,1:nen);
        actnode2 = find(ElemFlag2>0);
        xl2 = zeros(ndm,nen);
        xl2(1:ndm,actnode2) = Coordinates(ElemFlag2(actnode2),1:ndm)';
[EGDOFTa2, EGDOFTi2, ELDOFTa2, ELDOFTi2] = plocal(NDOFT, ElemFlag2, squeeze(iedof(:,:,ma2)), actnode2, nen, ndf, neq);
        ul2 = zeros(ndf, nen);
        ul2(ELDOFTa2) = ModelDx(EGDOFTa2)';
        ul2(ELDOFTi2) = gBC(EGDOFTi2)';
        else
        elem2 = elem + numel_mix/2;
        ma2 = RegionOnElement(elem2);
        ElemFlag2 = NodesOnElement(elem2,1:nen);
        actnode2 = find(ElemFlag2>0);
        xl2 = zeros(ndm,nen);
        xl2(1:ndm,actnode2) = Coordinates(ElemFlag2(actnode2),1:ndm)';
[EGDOFTa2, EGDOFTi2, ELDOFTa2, ELDOFTi2] = plocal(NDOFT, ElemFlag2, squeeze(iedof(:,:,ma2)), actnode2, nen, ndf, neq);
        ul2 = zeros(ndf, nen);
        ul2(ELDOFTa2) = ModelDx(EGDOFTa2)';
        ul2(ELDOFTi2) = gBC(EGDOFTi2)';
        end
end

% Current physical location of int pt
xint = (xl(1,1:nel)+ul(1,1:nel))*shl;
yint = (xl(2,1:nel)+ul(2,1:nel))*shl;
zint = (xl(3,1:nel)+ul(3,1:nel))*shl;

% Find current point's location in deformed config of other
% constituent
xi = POU_Coord3(xint,yint,zint,xl2(:,1:nel)+ul2(1:3,1:nel),1,nel);

% Evaluate  basis functions at integration point
shl2 = shlb(xi,nel,nel,0,0);
v1 = ul(1:3,1:nel)*shl;
v2 = ul2(1:3,1:nel)*shl2;
IntFor = JxX*kIF*(v1 - v2);