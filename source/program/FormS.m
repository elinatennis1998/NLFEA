% Assemble Quantities from Model Routine
%
% Tim Truster
% 7/2009
% UIUC

% Initialize Model-level arrays

switch isw%Task Switch
    
    case 22 %Get Mass
        Mdd11 = zeros(neqS,neqS);
        Mdf1 = zeros(neqS,nieqS);
        Mfd = zeros(nieqS,neqS);
        Mff = zeros(nieqS,nieqS);
        Fd1S = zeros(neqS,1);
        AssemQuant = 'AssemMassS';
        
end

for elem = 1:numelS
    
  for ma = 1:nummat
      
   if(ieFEAP(nie-2,ma) == RegionOnElement(elem))
      
    %Extract patch material properties
    iel = MatTypeTable(2,ma); %iel   = ie(nie-1,ma); same thing;
    nonlin = MatTypeTable(3,ma);
    mateprop = MateT(ma,:);
    if iscell(mateprop) % allows for crystal plasticity or other combined material listings
        mateprop = mateprop{1};
    end
    
    %Record time of assembly
    if Compt == 1
        tic
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
    nst = nel*ndfS;
    
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
%     if transient > 0
%         PatchVx = ul;
%         PatchAx = ul;
%         loc = 0;
%         for row = 1:nen
%             for dir = 1:ndf
%                 loc = loc + 1;
%                 grow = PDOFT(loc);
%                 if grow > neq
% %                     PatchDx(dir, row) = gBC(grow-neq);
%                     % PatchVx(dir, row) = 0;
%                     % PatchAx(dir, row) = 0;
%                 else
% %                     PatchDx(dir, row) = ModelDx(grow);
%                     PatchVx(dir, row) = ModelVx(grow);
%                     PatchAx(dir, row) = ModelAx(grow);
%                 end
%             end
%         end
%     else
%         loc = 0;
%         for row = 1:nel
%             for dir = 1:ndf
%                 loc = loc + 1;
%                 grow = EDOFT(loc);
%                 if grow > 0
%                 if grow > neq
%                     ul(dir, row) = gBC(grow-neq);
%                     ul_n(dir, row) = gBC_n(grow-neq);
%                     uld(dir, row) = gBC(grow-neq) - gBC_n(grow-neq);
%                 else
%                     ul(dir, row) = ModelDx(grow);
%                     ul_n(dir, row) = ModelDxn_1(grow);
%                     uld(dir, row) = s_del_ModelDx(grow);
%                 end
%                 end
%             end
%         end
%     end
    ul(ELDOFTa) = ModelDx(EGDOFTa)';
    ul(ELDOFTi) = gBC(EGDOFTi)';
    ul_n(ELDOFTa) = ModelDxn_1(EGDOFTa)';
    ul_n(ELDOFTi) = gBC_n(EGDOFTi)';
    uld(ELDOFTa) = s_del_ModelDx(EGDOFTa)';
    uld(ELDOFTi) = (gBC(EGDOFTi) - gBC_n(EGDOFTi))';
    
    ElemRout

    %Assemble Element contribution to Model Quantity

    run(AssemQuant)
    
    if Compt == 1
        t = toc;
        fprintf('Element %i time to assemble: %3.4f.\n',elem,t)
    end
    
   end
    
  end % ma
    
end % elem