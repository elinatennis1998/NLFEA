% Assemble Quantities from Model Routine
%
% Tim Truster
% 7/2009
% UIUC
%
% NOTE 02/11/2014: for dynamic analyses with linear elements, the stiffness
% and mass matrix must be kept separate in order to shorten the calculation
% time. Thus, in the element routine mass is handled by isw=5 and stiffness
% by isw=3, and the% results are combined within the main program in Mdd11 
% and Kdd11.
% However, for nonlinear analyses the FEAP-standard can be used where
% K and M are combined together as the Mstar matrix. The code assumes that
% the entire static and dynamic contribution is computed by isw=3 and then
% assembled into Kdd11.
% This has been updated for transient = 1,2,5
% transient=3 still assumes separate matrices for M and K

% Initialize Model-level arrays
hflgu = 0;
h3flgu = 0;

ModelA12 = zeros(neq,1);

% if (~((isw== 3||isw==6||isw==12) && numbernonlinear==0) || (exist('initializeLinKF','var') && initializeLinKF))
% Don't assemble F or K or energy if it's a linear problem

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
    
    if (~((isw== 3||isw==6||isw==12) && nonlin==0) || (exist('initializeLinKF','var') && initializeLinKF))
    
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
%     uld = zeros(ndf, nen);
    if transient > 0
        vl = ul;
        vl_n = ul_n;
        al = ul;
        al_n = ul_n;
        ul(ELDOFTa) = ModelDx(EGDOFTa)';
        ul(ELDOFTi) = gBC(EGDOFTi)';
        ul_n(ELDOFTa) = ModelDxn_1(EGDOFTa)';
        ul_n(ELDOFTi) = gBC_n(EGDOFTi)';
        uld = ul - ul_n;
%         uld(ELDOFTa) = s_del_ModelDx(EGDOFTa)';
%         uld(ELDOFTi) = (gBC(EGDOFTi) - gBC_n(EGDOFTi))';
        vl(ELDOFTa) = ModelVx(EGDOFTa)';
        vl_n(ELDOFTa) = ModelVxn_1(EGDOFTa)';
        al(ELDOFTa) = ModelAx(EGDOFTa)';
        al_n(ELDOFTa) = ModelAxn_1(EGDOFTa)';
    else
        ul(ELDOFTa) = ModelDx(EGDOFTa)';
        ul(ELDOFTi) = gBC(EGDOFTi)';
        ul_n(ELDOFTa) = ModelDxn_1(EGDOFTa)';
        ul_n(ELDOFTi) = gBC_n(EGDOFTi)';
        uld = ul - ul_n;
%         uld(ELDOFTa) = s_del_ModelDx(EGDOFTa)';
%         uld(ELDOFTi) = (gBC(EGDOFTi) - gBC_n(EGDOFTi))';
    end
    
    ElemRout
  
    %Assemble Element contribution to Model Quantity

    AssemA12
    
    if Compt == 1
        t = toc;
        fprintf('Element %i time to assemble: %3.4f.\n',elem,t)
    end

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
    
    end % if nonlin
    
   end %if ma
    
  end % ma
    
end % elem
    
% end

hflgu = 0;
h3flgu = 0;
