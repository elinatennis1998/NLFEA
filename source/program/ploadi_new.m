% 04/30/2013

%Get Loads
Fc1_new = zeros(neq, 1);
Fc1np_new = zeros(neq, 1);
Fcount_new = 0;

AssemQuant = 'AssemLoad_new';
isw = -1;
for FEload = 1:numSL_new
    
    nodeA = SurfacesL_new(FEload,1);
    nodeB = SurfacesL_new(FEload,2);
    elem = SurfacesL_new(FEload,3);
    edge = SurfacesL_new(FEload,4);
    traction = SurfacesL_new(FEload,5:7);
    
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

% non-proportional loads
AssemQuant = 'AssemLoadnp_new';
isw = -1;
for FEload = 1:numSLnp_new
    
    nodeA = SurfacesLnp_new(FEload,1);
    nodeB = SurfacesLnp_new(FEload,2);
    elem = SurfacesLnp_new(FEload,3);
    edge = SurfacesLnp_new(FEload,4);
    traction = SurfacesLnp_new(FEload,5:7);
    
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

for i = 1:numNodalF_new
    node = NodeLoad_new(i,1);
    dir = NodeLoad_new(i,2);
    force = NodeLoad_new(i,3);
    dof = NDOFT(node,dir);
    if dof <= neq
        Fc1_new(dof) = Fc1_new(dof) + force;
        Fcount_new = Fcount_new + 1;
    end
end

for i = 1:numNodalFnp_new
    node = NodeLoadnp_new(i,1);
    dir = NodeLoadnp_new(i,2);
    force = NodeLoadnp_new(i,3);
    dof = NDOFT(node,dir);
    if dof <= neq
        Fc1np_new(dof) = Fc1np_new(dof) + force;
        Fcount_new = Fcount_new + 1;
    end
end

% Add forces to previous values
Fcount = max(Fcount,Fcount_new);
Fc1 = Fc1 + Fc1_new;
Fc1np = Fc1np + Fc1np_new;
