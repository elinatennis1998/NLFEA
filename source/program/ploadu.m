% Tim Truster
% 12/30/2013
%
% Assemble surface loads defined by user
% 
% Flag Uiterstep is set to control:
%   Evaluate user loads at 0=step, 1=every iteration

if (Uiterstep) || (iter == 0) % only assemble if first iteration or if User desires every iteration

FcU = zeros(neq, 1);

AssemQuant = 'AssemLoadU';
isw = 7;
for FEload = 1:numUSL
    
    nodeA = SurfacesU(FEload,1);
    nodeB = SurfacesU(FEload,2);
    elem = SurfacesU(FEload,3);
    edge = SurfacesU(FEload,4);
    
    nel = nnz(NodesOnElement(elem,1:nen));
    nelP = getnelP(nel,ndm,nelP3,nelP4,nelP6,nelP9);
    nst = nel*ndf;
    
    %Extract patch nodal coordinates
    ElemFlag = NodesOnElement(elem,1:nen);
    actnode = find(ElemFlag>0);
    xl = zeros(ndm,nen);
    xl(1:ndm,actnode) = Coordinates(ElemFlag(actnode),1:ndm)';
    
    %Extract patch material properties
    ma = RegionOnElement(elem);
    mateprop = MateT(ma,:);
    iel = MatTypeTable(2,ma);
    nonlin = MatTypeTable(3,ma);

    %Compute and Assemble Patch Stiffness
%     EDOFT = LocToGlobDOF(ElemFlag, NDOFT, nel, ndf);
%     [EGDOFTa, EGDOFTi, ELDOFTa, ELDOFTi] = LocToGlobDOF2(ElemFlag, NDOFT, nel, ndf, neq);
    [EGDOFTa, EGDOFTi, ELDOFTa, ELDOFTi] = plocal(NDOFT, ElemFlag, squeeze(iedof(:,:,ma)), actnode, nen, ndf, neq);
    
    ElemRout

    %Assemble Element contribution to Model Quantity

    run(AssemQuant)
    
end

end
