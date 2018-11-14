% Assembly for one element
% 06/15/2013
    
% elem = n;
    
    %Extract patch material properties
    ma = RegionOnElement(elem);
    iel = MatTypeTable(2,ma);
    nonlin = MatTypeTable(3,ma);
    mateprop = MateT(ma,:);
    
    %Determine element size parameters
    nel = nnz(NodesOnElement(elem,1:nen));
    nelP = getnelP(nel,ndm,nelP3,nelP4,nelP6,nelP9);
    nst = nel*ndf;
    
    %Extract patch nodal coordinates
    
    ElemFlag = NodesOnElement(elem,1:nen);
    actnode = find(ElemFlag>0);
    xl = zeros(ndm,nen);
    xl(1:ndm,actnode) = Coordinates(ElemFlag(actnode),1:ndm)';

    %Compute and Assemble Patch Stiffness
%     EDOFT = LocToGlobDOF(ElemFlag, NDOFT, nel, ndf);
%     [EGDOFTa, EGDOFTi, ELDOFTa, ELDOFTi] = LocToGlobDOF2(ElemFlag, NDOFT, nel, ndf, neq);
    [EGDOFTa, EGDOFTi, ELDOFTa, ELDOFTi] = plocal(NDOFT1, ElemFlag, squeeze(iedof(:,:,ma)), actnode, nen, ndf, neq1);
    
    ElemF = zeros(nst,1);
        
    ElemRout

    %Assemble Element contribution to Model Quantity
    if(intfl1==1)
        AssemLoad
    elseif(intfl1==-1)
        AssemLoadnp
    end
    
% end
