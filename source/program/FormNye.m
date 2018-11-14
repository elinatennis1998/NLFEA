% Tim Truster
% 05/15/2014
%
% Subroutine to compute nodally-averaged fields from
% integration-point-based data by using Lie groups and algebras according
% to the Sandia Computational Mechanics paper.
        
% Initialize nodal array for storing nodal-averaged quantities
% numLie = 18; % tensors Re and Ue
LieList = zeros(numnp,numLie);
Eareas = zeros(numnp,1);

nh1 = nha;
nh2 = nhb;
nh3 = nhc;

for elem = 1:numel
    
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
    nelS = getnelP(nel,ndm,nelS3,nelS4,nelS6,nelS9);
    nst = nen*ndf;
    
    %Extract patch nodal coordinates
    
    ElemFlag = NodesOnElement(elem,1:nen);
    actnode = find(ElemFlag>0);
    xl = zeros(ndm,nen);
    xl(1:ndm,actnode) = Coordinates(ElemFlag(actnode),1:ndm)';

    %Compute and Assemble Patch Stiffness
    [EGDOFTa, EGDOFTi, ELDOFTa, ELDOFTi] = plocal(NDOFT, ElemFlag, squeeze(iedof(:,:,ma)), actnode, nen, ndf, neq);
    
    %Extract patch solution values
    ul = zeros(ndf, nen);
    ul_n = zeros(ndf, nen);
%     uld = zeros(ndf, nen);
    ul(ELDOFTa) = ModelDx(EGDOFTa)';
    ul(ELDOFTi) = gBC(EGDOFTi)';
    ul_n(ELDOFTa) = ModelDxn_1(EGDOFTa)';
    ul_n(ELDOFTi) = gBC_n(EGDOFTi)';
    uld = ul - ul_n;
    
    ElemRout

    %Assemble Element contribution to Model Quantity
    
    for k = 1:nel
        node = ElemFlag(k);
        Eareas(node) = Eareas(node) + ElemL(k,numLie+1);
    end
    
    for stres = 1:numLie
    
    for k = 1:nel
        node = ElemFlag(k);
        LieList(node, stres) = LieList(node, stres) + ElemL(k,stres)*ElemL(k,numLie+1);
    end
    
    end
    
   end
    
  end % ma
    
end % elem


% Complete post-process by dividing by areas
for stres = 1:numLie
    LieList(:,stres) = LieList(:,stres)./Eareas;
end

% % Then map back to group space
% isw = 99;
% for node = 1:numnp
%     ElemL = LieList(node,1:numLie);
%     ElemRout
%     LieList(node,1:numLie) = ElemL;
% end
% isw = 98;

% Then map back to group space
hflgu = 1;
h3flgu = 1;
isw = 99;

nh1 = nha;
nh2 = nhb;
nh3 = nhc;

for elem = 1:numel
    
  for ma = 1:nummat
      
   if(ieFEAP(nie-2,ma) == RegionOnElement(elem))
      
    %Extract patch material properties
    iel = MatTypeTable(2,ma); %iel   = ie(nie-1,ma); same thing;
    nonlin = MatTypeTable(3,ma);
    mateprop = MateT(ma,:);
    
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
    nelS = getnelP(nel,ndm,nelS3,nelS4,nelS6,nelS9);
    nst = nen*ndf;
    
    %Extract patch nodal coordinates
    
    ElemFlag = NodesOnElement(elem,1:nen);
    actnode = find(ElemFlag>0);
    xl = zeros(ndm,nen);
    xl(1:ndm,actnode) = Coordinates(ElemFlag(actnode),1:ndm)';

    %Compute and Assemble Patch Stiffness
    [EGDOFTa, EGDOFTi, ELDOFTa, ELDOFTi] = plocal(NDOFT, ElemFlag, squeeze(iedof(:,:,ma)), actnode, nen, ndf, neq);
    
    %Extract patch solution values
    ul = zeros(ndf, nen);
    ul_n = zeros(ndf, nen);
%     uld = zeros(ndf, nen);
    ul(ELDOFTa) = ModelDx(EGDOFTa)';
    ul(ELDOFTi) = gBC(EGDOFTi)';
    ul_n(ELDOFTa) = ModelDxn_1(EGDOFTa)';
    ul_n(ELDOFTi) = gBC_n(EGDOFTi)';
    uld = ul - ul_n;
    ElemL = LieList(ElemFlag,1:numLie);
    
    ElemRout

%             Position update terms 'ht1,ht2' from 'nh1,nh2' to save

    if(hflgu && ieFEAP(nie,ma) > 0) %then
      for i = 0:ieFEAP(nie,ma)-1
        temp      = hrvec(ht1+i);
        hrvec(ht1+i) = hr(nha+i);
        hr(nha+i) = temp;
        temp      = hrvec(ht2+i);
        hrvec(ht2+i) = hr(nhb+i);
        hr(nhb+i) = temp;
      end % i
    end

%             Position update terms 'ht3' from 'nh3' to save

    if(h3flgu && ieFEAP(nie-5,ma) > 0) %then
      for i = 0:ieFEAP(nie-5,ma)-1
        hrvec(ht3+i) = hr(nhc+i);
      end % i
    end
    
   end
    
  end % ma
    
end % elem
isw = 98;


LieHList(1:numLie,:,step) = LieList(:,1:numLie)';
