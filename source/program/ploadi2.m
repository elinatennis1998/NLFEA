% 04/30/2013

%Get Loads
Fc2 = zeros(neq2, 1);
Fc2np = zeros(neq2, 1);
Fcount2 = 0;

AssemQuant = 'AssemLoad2';
isw = -1;
for load = 1:numSL2
    
    nodeA = SurfacesL2(load,1);
    nodeB = SurfacesL2(load,2);
    elem = SurfacesL2(load,3);
    edge = SurfacesL2(load,4);
    traction = SurfacesL2(load,5:7);
    
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
    [EGDOFTa, EGDOFTi, ELDOFTa, ELDOFTi] = plocal(NDOFT2, ElemFlag, squeeze(iedof(:,:,ma)), actnode, nen, ndf, neq2);
    
    ElemRout

    %Assemble Element contribution to Model Quantity

    run(AssemQuant)
    
end

% non-proportional loads
AssemQuant = 'AssemLoadnp2';
isw = -1;
for load = 1:numSLnp2
    
    nodeA = SurfacesLnp2(load,1);
    nodeB = SurfacesLnp2(load,2);
    elem = SurfacesLnp2(load,3);
    edge = SurfacesLnp2(load,4);
    traction = SurfacesLnp2(load,5:7);
    
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
    [EGDOFTa, EGDOFTi, ELDOFTa, ELDOFTi] = plocal(NDOFT2, ElemFlag, squeeze(iedof(:,:,ma)), actnode, nen, ndf, neq2);
    
    ElemRout

    %Assemble Element contribution to Model Quantity

    run(AssemQuant)
    
end

AssemQuant = 'AssemLoad2';
isw = -3;
for load = 1:numD2
    
    elem = SurfacesD2(load,1);
    edge = SurfacesD2(load,2);
    displac = SurfacesD2(load,3:5);
    
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
    [EGDOFTa, EGDOFTi, ELDOFTa, ELDOFTi] = plocal(NDOFT2, ElemFlag, squeeze(iedof(:,:,ma)), actnode, nen, ndf, neq2);
    
    ElemRout

    %Assemble Element contribution to Model Quantity

    run(AssemQuant)
    
end

for i = 1:numNodalF2
    node = NodeLoad2(i,1);
    dir = NodeLoad2(i,2);
    force = NodeLoad2(i,3);
    dof = NDOFT2(node,dir);
    if dof <= neq2
        Fc2(dof) = Fc2(dof) + force;
        Fcount2 = Fcount2 + 1;
    end
end

for i = 1:numNodalFnp2
    node = NodeLoadnp2(i,1);
    dir = NodeLoadnp2(i,2);
    force = NodeLoadnp2(i,3);
    dof = NDOFT2(node,dir);
    if dof <= neq
        Fc2np(dof) = Fc2np(dof) + force;
        Fcount2 = Fcount2 + 1;
    end
end

ModelDx2 = zeros(neq2,1);
ModelVx2 = zeros(neq2,1);
ModelAx2 = zeros(neq2,1);

if ismember(transient,transalgos)
    %Assign IC
    for i = 1:numIC2
        node = NodeIC2(i,1);
        dir = NodeIC2(i,2);
        dis = NodeIC2(i,3);
        dof = NDOFT2(node,dir);
        if dof <= neq2
            ModelDx2(dof) = dis;
        end
    end
    %Assign Velo_IC
    for i = 1:numVIC2
        node = NodeVIC2(i,1);
        dir = NodeVIC2(i,2);
        vel = NodeVIC2(i,3);
        dof = NDOFT2(node,dir);
        if dof <= neq2
            ModelVx2(dof) = vel;
        end
    end
end
ModelDx2n_1 = ModelDx2;
ModelVx2n_1 = ModelVx2;
ModelAx2n_1 = ModelAx2;

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