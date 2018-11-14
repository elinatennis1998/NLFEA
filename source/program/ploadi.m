% 04/30/2013

%Get Loads
Fc1 = zeros(neq, 1);
Fc1np = zeros(neq, 1);
Fcount = 0;

AssemQuant = 'AssemLoad';
isw = -1;
for FEload = 1:numSL
    
    nodeA = SurfacesL(FEload,1);
    nodeB = SurfacesL(FEload,2);
    elem = SurfacesL(FEload,3);
    edge = SurfacesL(FEload,4);
    traction = SurfacesL(FEload,5:7);
    
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
AssemQuant = 'AssemLoadnp';
isw = -1;
for FEload = 1:numSLnp
    
    nodeA = SurfacesLnp(FEload,1);
    nodeB = SurfacesLnp(FEload,2);
    elem = SurfacesLnp(FEload,3);
    edge = SurfacesLnp(FEload,4);
    traction = SurfacesLnp(FEload,5:7);
    
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

% Multiple tabular loads
if numSLmt > 0
AssemQuant = 'AssemLoadmt';
isw = -1;
for load_entry = 1:numSLmt
    Fc1mt{load_entry} = zeros(neq, 1);
for FEload = 1:size(SurfacesLmt{2,load_entry},1)
    
    nodeA = SurfacesLmt{2,load_entry}(FEload,1);
    nodeB = SurfacesLmt{2,load_entry}(FEload,2);
    elem = SurfacesLmt{2,load_entry}(FEload,3);
    edge = SurfacesLmt{2,load_entry}(FEload,4);
    traction = SurfacesLmt{2,load_entry}(FEload,5:7);
    
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
end

AssemQuant = 'AssemLoad';
isw = -3;
for FEload = 1:numD
    
    elem = SurfacesD(FEload,1);
    edge = SurfacesD(FEload,2);
    displac = SurfacesD(FEload,3:5);
    
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

% User surface loads
ploadu

for i = 1:numNodalF
    node = NodeLoad(i,1);
    dir = NodeLoad(i,2);
    if dir > ndf
        errmsg = ['dof ID exceeds ndf for Nodalload=' num2str(i)];
        error(errmsg)
    end
    force = NodeLoad(i,3);
    dof = NDOFT(node,dir);
    if dof <= neq
        Fc1(dof) = Fc1(dof) + force;
        Fcount = Fcount + 1;
    end
end

for i = 1:numNodalFnp
    node = NodeLoadnp(i,1);
    dir = NodeLoadnp(i,2);
    if dir > ndf
        errmsg = ['dof ID exceeds ndf for Nodalload_np=' num2str(i)];
        error(errmsg)
    end
    force = NodeLoadnp(i,3);
    dof = NDOFT(node,dir);
    if dof <= neq
        Fc1np(dof) = Fc1np(dof) + force;
        Fcount = Fcount + 1;
    end
end

ModelDx = zeros(neq,1);
ModelVx = zeros(neq,1);
ModelAx = zeros(neq,1);

if ismember(transient,transalgos)
    %Assign IC
    for i = 1:numIC
        node = NodeIC(i,1);
        dir = NodeIC(i,2);
        dis = NodeIC(i,3);
        dof = NDOFT(node,dir);
        if dof <= neq
            ModelDx(dof) = dis;
        end
    end
    %Assign Velo_IC
    for i = 1:numVIC
        node = NodeVIC(i,1);
        dir = NodeVIC(i,2);
        vel = NodeVIC(i,3);
        dof = NDOFT(node,dir);
        if dof <= neq
            ModelVx(dof) = vel;
        end
    end
end
ModelDxn_1 = ModelDx;
ModelVxn_1 = ModelVx;
ModelAxn_1 = ModelAx;

%Contact DOF

% na = 0;
% l = 1;
% CDOFT = zeros(numnp,nc);
% if twopass >= 0 %Left face multipliers
%     for conta = 1:numCn
%         nodeAL = SurfacesI(conta,3);
%         nodeBL = SurfacesI(conta,4);
%         elemL = SurfacesI(conta,5);
%         if CDOFT(nodeAL,1) == 0
%             for l = 1:nc
%                 na = na + 1;
%                 CDOFT(nodeAL,l) = na;
%             end
%         end
%         if CDOFT(nodeBL,1) == 0
%             for l = 1:nc
%                 na = na + 1;
%                 CDOFT(nodeBL,l) = na;
%             end
%         end
%     end
% end
% if abs(twopass) == 1 %right face multipliers
%     for conta = 1:numCn
%         nodeAR = SurfacesI(conta,1);
%         nodeBR = SurfacesI(conta,2);
%         elemR = SurfacesI(conta,6);
%         if CDOFT(nodeAR,1) == 0
%             for l = 1:nc
%                 na = na + 1;
%                 CDOFT(nodeAR,l) = na;
%             end
%         end
%         if CDOFT(nodeBR,1) == 0
%             for l = 1:nc
%                 na = na + 1;
%                 CDOFT(nodeBR,l) = na;
%             end
%         end
%     end
% end
% neqc = na;