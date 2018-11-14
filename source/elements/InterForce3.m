% Tim Truster
% 02/27/2014
%
% Evaluate interactive force, simple model

IntFor = zeros(ndf,1);

% get values from other side; for 3+ constituents, will only give the
% contribution from the constituent that is numbered one below the current
% one
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
        constactiv = 1; % constituent pair is active
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
        constactiv = 1; % constituent pair is active
    case 0 % concurrent formation of constituent contributions
        if exist('numconst','var') % number of constituents is provided in the input file
            if elem <= numel_mix/numconst
            % Current element is in constituent 1
            elem2 = elem + numel_mix/numconst;
            ma2 = RegionOnElement(elem2);
            ElemFlag2 = NodesOnElement(elem2,1:nen);
            actnode2 = find(ElemFlag2>0);
            xl2 = zeros(ndm,nen);
            xl2(1:ndm,actnode2) = Coordinates(ElemFlag2(actnode2),1:ndm)';
    [EGDOFTa2, EGDOFTi2, ELDOFTa2, ELDOFTi2] = plocal(NDOFT, ElemFlag2, squeeze(iedof(:,:,ma2)), actnode2, nen, ndf, neq);
            ul2 = zeros(ndf, nen);
            ul2(ELDOFTa2) = ModelDx(EGDOFTa2)';
            ul2(ELDOFTi2) = gBC(EGDOFTi2)';
            ul2_n = zeros(ndf, nen);
            ul2_n(ELDOFTa2) = ModelDxn_1(EGDOFTa2)';
            ul2_n(ELDOFTi2) = gBC_n(EGDOFTi2)';
            constactiv = 1; % constituent pair is active
            elseif elem <= (numconst-1)*numel_mix/numconst
            % Current element is in some constituent 2 to numconst
            elem2 = elem + numel_mix/numconst;
            ma2 = RegionOnElement(elem2);
            ElemFlag2 = NodesOnElement(elem2,1:nen);
            actnode2 = find(ElemFlag2>0);
            xl2 = zeros(ndm,nen);
            xl2(1:ndm,actnode2) = Coordinates(ElemFlag2(actnode2),1:ndm)';
    [EGDOFTa2, EGDOFTi2, ELDOFTa2, ELDOFTi2] = plocal(NDOFT, ElemFlag2, squeeze(iedof(:,:,ma2)), actnode2, nen, ndf, neq);
            ul2 = zeros(ndf, nen);
            ul2(ELDOFTa2) = ModelDx(EGDOFTa2)';
            ul2(ELDOFTi2) = gBC(EGDOFTi2)';
            ul2_n = zeros(ndf, nen);
            ul2_n(ELDOFTa2) = ModelDxn_1(EGDOFTa2)';
            ul2_n(ELDOFTi2) = gBC_n(EGDOFTi2)';
            constactiv = 1; % constituent pair is active
            else
            constactiv = 0; % constituent pair is active
            end
        else % older version that only has two constituents
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
            ul2_n = zeros(ndf, nen);
            ul2_n(ELDOFTa2) = ModelDxn_1(EGDOFTa2)';
            ul2_n(ELDOFTi2) = gBC_n(EGDOFTi2)';
            constactiv = 1; % constituent pair is active
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
            ul2_n = zeros(ndf, nen);
            ul2_n(ELDOFTa2) = ModelDxn_1(EGDOFTa2)';
            ul2_n(ELDOFTi2) = gBC_n(EGDOFTi2)';
            constactiv = 1; % constituent pair is active
            end
        end
end


if constactiv % constituent pair is active
    
% Current physical location of int pt
xint = (xl(1,1:nel)+ul(1,1:nel))*shl;
yint = (xl(2,1:nel)+ul(2,1:nel))*shl;
zint = (xl(3,1:nel)+ul(3,1:nel))*shl;

% Find current point's location in deformed config of other
% constituent
xi = POU_Coord3(xint,yint,zint,xl2(:,1:nel)+ul2(1:3,1:nel),1,nel);

% Evaluate  basis functions at integration point
shl2 = shlb(xi,nel,nel,0,0);
v1 = (ul(1:3,1:nel)-ul_n(1:3,1:nel))/tstep*shl;
v2 = (ul2(1:3,1:nel)-ul2_n(1:3,1:nel))/tstep*shl2;

if exist('numconst','var') % number of constituents is provided in the input file
    constituentpair = ceil(elem/(numel_mix/numconst));
    kIF = kIFtable(constituentpair);
end

IntFor = JxX*kIF*(v1 - v2);

end


% Check for contribution from constituent numbered one greater than the
% current one.
if exist('numconst','var') % number of constituents is provided in the input file
    if elem > numel_mix/numconst
        % Current element is in some constituent 2 to numconst
        elem2 = elem - numel_mix/numconst;
        ma2 = ix(elem2,nen1);
        ElemFlag2 = ix(elem2,1:nen);
        actnode2 = find(ElemFlag2>0);
        xl2 = zeros(ndm,nen);
        xl2(1:ndm,actnode2) = Coordinates(ElemFlag2(actnode2),1:ndm)';
        [EGDOFTa2, EGDOFTi2, ELDOFTa2, ELDOFTi2] = plocal(NDOFT, ElemFlag2, squeeze(iedof(:,:,ma2)), actnode2, nen, ndf, neq);
        ul2 = zeros(ndf, nen);
        ul2(ELDOFTa2) = ModelDx(EGDOFTa2)';
        ul2(ELDOFTi2) = gBC(EGDOFTi2)';
        ul2_n = zeros(ndf, nen);
        ul2_n(ELDOFTa2) = ModelDxn_1(EGDOFTa2)';
        ul2_n(ELDOFTi2) = gBC_n(EGDOFTi2)';
        constactiv = 1; % constituent pair is active
    else
        constactiv = 0; % constituent pair is active
    end
else
    constactiv = 0; % constituent pair is active
end


if constactiv % constituent pair is active
    
% Current physical location of int pt
xint = (xl(1,1:nel)+ul(1,1:nel))*shl;
yint = (xl(2,1:nel)+ul(2,1:nel))*shl;
zint = (xl(3,1:nel)+ul(3,1:nel))*shl;

% Find current point's location in deformed config of other
% constituent
xi = POU_Coord3(xint,yint,zint,xl2(:,1:nel)+ul2(1:3,1:nel),1,nel);

% Evaluate  basis functions at integration point
shl2 = shlb(xi,nel,nel,0,0);
v1 = (ul(1:3,1:nel)-ul_n(1:3,1:nel))/tstep*shl;
v2 = (ul2(1:3,1:nel)-ul2_n(1:3,1:nel))/tstep*shl2;

constituentpair = ceil(elem/(numel_mix/numconst)) - 1;
kIF = kIFtable(constituentpair);

IntFor = IntFor + JxX*kIF*(v1 - v2);

end