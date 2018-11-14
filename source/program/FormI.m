% Compute interface quantities after converged solution
% isw=60 must be set outside this routine
%
% Tim Truster
% 03/25/2014

nh1 = nha;
nh2 = nhb;
nh3 = nhc;

segment = 0;

% loop over DG elements; data structure already pre-set
for elem = numel-numSI+1:numel
    
  for ma = 1:nummat
      
   if(ieFEAP(nie-2,ma) == RegionOnElement(elem))
      
    %Extract patch material properties
    iel = MatTypeTable(2,ma); %iel   = ie(nie-1,ma); same thing;
    nonlin = MatTypeTable(3,ma);
    mateprop = MateT(ma,:);
    if iscell(mateprop) % allows for crystal plasticity or other combined material listings
        mateprop = mateprop{1};
    end
    
    
%             Compute address and offset for history variables

    ht1 = 1 + ixFEAP(1,elem) + ieFEAP(nie-3,ma);
    ht2 = 1 + ixFEAP(2,elem) + ieFEAP(nie-3,ma);
    ht3 = 1 + ixFEAP(3,elem) + ieFEAP(nie-4,ma);

%             If history variables exist move into nh1,nh2

    if(ieFEAP(nie,ma) > 0) %then
        for i = 0:ieFEAP(nie,ma)-1
          hr(nha+i) = hrvec(ht1+i);
          hr(nhb+i) = hrvec(ht2+i);
        end % i
    end

%             If Element variables exist move into nh3

    if(ieFEAP(nie-5,ma) > 0) %then
        for i = 0:ieFEAP(nie-5,ma)-1
          hr(nhc+i) = hrvec(ht3+i);
        end % i
    end
    
    %Determine element size parameters
    nel = nnz(NodesOnElement(elem,1:nen));
    nelP = getnelP(nel,ndm,nelP3,nelP4,nelP6,nelP9);
    nst = nen*ndf;
    
    %Extract patch nodal coordinates
    
    ElemFlag = NodesOnElement(elem,1:nen);
    actnode = find(ElemFlag>0);
    xl = zeros(ndm,nen);
    xl(1:ndm,actnode) = Coordinates(ElemFlag(actnode),1:ndm)';

    %Compute and Assemble Patch Stiffness
%     EDOFT = LocToGlobDOF(ElemFlag, NDOFT, nel, ndf);
%     [EGDOFTa, EGDOFTi, ELDOFTa, ELDOFTi] = LocToGlobDOF2(ElemFlag, NDOFT, nel, ndf, neq);
    [EGDOFTa, EGDOFTi, ELDOFTa, ELDOFTi] = plocal(NDOFT, ElemFlag, squeeze(iedof(:,:,ma)), actnode, nen, ndf, neq);
    
    %Extract patch solution values
    ul = zeros(ndf, nen);
    ul_n = zeros(ndf, nen);
    uld = zeros(ndf, nen);
    if transient > 0
        vl = ul;
        vl_n = ul_n;
        al = ul;
        al_n = ul_n;
        ul(ELDOFTa) = ModelDx(EGDOFTa)';
        ul(ELDOFTi) = gBC(EGDOFTi)';
        ul_n(ELDOFTa) = ModelDxn_1(EGDOFTa)';
        ul_n(ELDOFTi) = gBC_n(EGDOFTi)';
        uld(ELDOFTa) = s_del_ModelDx(EGDOFTa)';
        uld(ELDOFTi) = (gBC(EGDOFTi) - gBC_n(EGDOFTi))';
        vl(ELDOFTa) = ModelVx(EGDOFTa)';
        vl_n(ELDOFTa) = ModelVxn_1(EGDOFTa)';
        al(ELDOFTa) = ModelAx(EGDOFTa)';
        al_n(ELDOFTa) = ModelAxn_1(EGDOFTa)';
    else
        ul(ELDOFTa) = ModelDx(EGDOFTa)';
        ul(ELDOFTi) = gBC(EGDOFTi)';
        ul_n(ELDOFTa) = ModelDxn_1(EGDOFTa)';
        ul_n(ELDOFTi) = gBC_n(EGDOFTi)';
        uld(ELDOFTa) = s_del_ModelDx(EGDOFTa)';
        uld(ELDOFTi) = (gBC(EGDOFTi) - gBC_n(EGDOFTi))';
    end
    
    ElemRout
    
   end %if ma
    
  end % ma
    
end % elem
