% 04/30/2013

%Get Loads
Fc1 = zeros(neq1, 1);
Fc1np = zeros(neq1, 1);
Fcount1 = 0;

AssemQuant = 'AssemLoad';
isw = -1;
for load = 1:numSL1
    
    nodeA = SurfacesL1(load,1);
    nodeB = SurfacesL1(load,2);
    elem = SurfacesL1(load,3);
    edge = SurfacesL1(load,4);
    traction = SurfacesL1(load,5:7);
    
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
    [EGDOFTa, EGDOFTi, ELDOFTa, ELDOFTi] = plocal(NDOFT1, ElemFlag, squeeze(iedof(:,:,ma)), actnode, nen, ndf, neq1);
    
    ElemRout

    %Assemble Element contribution to Model Quantity

    run(AssemQuant)
    
end

% non-proportional loads
AssemQuant = 'AssemLoadnp';
isw = -1;
for load = 1:numSLnp1
    
    nodeA = SurfacesLnp1(load,1);
    nodeB = SurfacesLnp1(load,2);
    elem = SurfacesLnp1(load,3);
    edge = SurfacesLnp1(load,4);
    traction = SurfacesLnp1(load,5:7);
    
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
    [EGDOFTa, EGDOFTi, ELDOFTa, ELDOFTi] = plocal(NDOFT1, ElemFlag, squeeze(iedof(:,:,ma)), actnode, nen, ndf, neq1);
    
    ElemRout

    %Assemble Element contribution to Model Quantity

    run(AssemQuant)
    
end

AssemQuant = 'AssemLoad';
isw = -3;
for load = 1:numD1
    
    elem = SurfacesD1(load,1);
    edge = SurfacesD1(load,2);
    displac = SurfacesD1(load,3:5);
    
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
    [EGDOFTa, EGDOFTi, ELDOFTa, ELDOFTi] = plocal(NDOFT1, ElemFlag, squeeze(iedof(:,:,ma)), actnode, nen, ndf, neq1);
    
    ElemRout

    %Assemble Element contribution to Model Quantity

    run(AssemQuant)
    
end

for i = 1:numNodalF1
    node = NodeLoad1(i,1);
    dir = NodeLoad1(i,2);
    force = NodeLoad1(i,3);
    dof = NDOFT1(node,dir);
    if dof <= neq1
        Fc1(dof) = Fc1(dof) + force;
        Fcount1 = Fcount1 + 1;
    end
end

for i = 1:numNodalFnp1
    node = NodeLoadnp1(i,1);
    dir = NodeLoadnp1(i,2);
    force = NodeLoadnp1(i,3);
    dof = NDOFT1(node,dir);
    if dof <= neq
        Fc1np(dof) = Fc1np(dof) + force;
        Fcount1 = Fcount1 + 1;
    end
end

ModelDx1 = zeros(neq1,1);
ModelVx1 = zeros(neq1,1);
ModelAx1 = zeros(neq1,1);

if ismember(transient,transalgos)
    %Assign IC
    for i = 1:numIC1
        node = NodeIC1(i,1);
        dir = NodeIC1(i,2);
        dis = NodeIC1(i,3);
        dof = NDOFT1(node,dir);
        if dof <= neq1
            ModelDx1(dof) = dis;
        end
    end
    %Assign Velo_IC
    for i = 1:numVIC1
        node = NodeVIC1(i,1);
        dir = NodeVIC1(i,2);
        vel = NodeVIC1(i,3);
        dof = NDOFT1(node,dir);
        if dof <= neq1
            ModelVx1(dof) = vel;
        end
    end
end
ModelDx1n_1 = ModelDx1;
ModelVx1n_1 = ModelVx1;
ModelAx1n_1 = ModelAx1;

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