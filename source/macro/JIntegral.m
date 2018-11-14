% Tim Truster
% 09/07/2014
%
% Routine to calculate J integrals for fracture mechanics
% Based off of WARP3D manual discussion, and pjint3d.f from FEAP

% Currently, 2D version only
% Method borrows from implicit error estimation macros: forms automatic
% domains similar to WARP3D by adding layers of elements around the crack
% tip node. Definition for the q-field is: 1's for nodes on interior of
% domain, 0's on the boundary nodes.

% Loop to calculate J for a series of domain sizes

if ndm == 2

nstar = JDImax; % maximum star size for J-integrals
snode = JNodes(1); % crack tip node; could be a loop over cracks in the future expansion
maxel = 8; % maximum valancy of nodes (max # of elements attached)
JstepList = zeros(ndm,nstar);

% Setup up lots of arrays
	cse = zeros(numel,1); % elements in current star
    csn = zeros(numnp,1); % nodes in current star
    c1se = zeros(numel,1); % elements in previous star
    c1sn = zeros(numnp,1); % nodes in previous star
	cseb = zeros(numel,1); % elements in boundary layer between current and previous star
    csnb = zeros(numnp,1); % nodes on boundary of current star
    c1seb = zeros(numel,1); % elements in boundary layer between previous and previous2 star
    c1snb = zeros(numnp,1); % nodes on boundary of previous star
	epatch = zeros(maxel,numnp); % elements attached to nodes in mesh
    epnum = zeros(numnp,1); % number of elements attached to nodes
	ixe = zeros(nen,1);
    cis = zeros(2,numel+1);
    cis(2,:) = nen; %REALLY: should be nel for each element.
%%
%-----------------------------------------------------
% Determine level 1 stars around coarse nodes
%-----------------------------------------------------

if nen <= 4

% Loop over elements

	for elem=1:numel

	  for k=1:nen % Loop over local Nodes

	      node = NodesOnElement(elem,k);

%   Add element to star list, increment number of elem in star
           if (node>0)

	         q = epnum(node)+1;

	         epatch(q,node) = elem;		%epatch(nel,numnp)

	         epnum(node) = q;			%epnum(numnp)

           end

	  end %k
	end %j

	elseif(nen==6)

% Loop over elements
	for elem=1:numel

	   for k=1:3 % Loop over local Nodes

	      node = NodesOnElement(elem,k);

%   Add element to star list, increment number of elem in star
           if (node>0)

	         q = epnum(node)+1;

	         epatch(q,node) = elem;		%epatch(nel,numnp)

	         epnum(node) = q;			%epnum(numnp)

           end

	  end %k
	end %j

	else

	% Loop over elements
	for elem=1:numel

	   for k=1:3 % Loop over local Nodes

	      node = NodesOnElement(elem,k);

%   Add element to star list, increment number of elem in star
          if (node>0)

	         q = epnum(node)+1;

	         epatch(q,node) = elem;		%epatch(nel,numnp)

	         epnum(node) = q;			%epnum(numnp)

          end

	   end %k

	   if(NodesOnElement(elem,9)>0)

	       node = NodesOnElement(elem,k);

%   Add element to star list, increment number of elem in star

	         q = epnum(node)+1;

	         epatch(q,node) = elem;		%epatch(nel,numnp)

	         epnum(node) = q;			%epnum(numnp)

	   end

	end %j


end
    
%%
%  -----------------------------------------------------
%   Determine level nstar stars around coarse nodes
%  -----------------------------------------------------
% Perform J integral for each one
	   
	   cstar = 0;
	   csn(1) = snode;
	   cse(1) = 0;
	   csnn = 1;
	   csen = 0;
	   csnb(1) = snode;
	   cseb(1) = 0;
	   csnbn = 1;
	   csebn = 0;

%   While nstar > cstar, add new level
       while(nstar>cstar)
           
J_int = zeros(ndm,1); % Components of J-integral in each coordinate direction; user selects crack direction on their own.

	      cstar = cstar + 1;

%     Shift elements and nodes in cstar to c_1star

	      c1snn = csnn;
	      c1sen = csen;
	      c1snbn = csnbn;
	      c1sebn = csebn;

          for i = 1:c1snn
	         c1sn(i) = csn(i);
          end
          for i = 1:c1sen
	         c1se(i) = cse(i);
          end
          for i = 1:c1snbn
	         c1snb(i) = csnb(i);
          end
          for i = 1:c1sebn
	         c1seb(i) = cseb(i);
          end

%     Loop through nodes in c1snb, shift elements into cseb

            node = c1snb(1);
	      csebn = epnum(node);
          for i = 1:csebn
	         cseb(i) = epatch(i,node);
          end

          for i = 2:c1snbn

%       Shuffle elements from nodal patch(i) into cseb
             node = c1snb(i);
             nume = epnum(node);
	         [cseb,csebn] = ShuffleIn(epatch(:,node),cseb,nume,csebn);

%     End loop through nodes in c1snb
          end

%     Combine boundary cseb into c1se and get cse

           [cse,csen,cseb,csebn] = Purge(cseb,c1se,csebn,c1sen);

%     Loop through elements in cseb, shift nodes into csnb

            elem = cseb(1);
	      csnbn = cis(2,elem);
          for i = 1:csnbn
	         csnb(i) = NodesOnElement(elem,i);
          end
	      [csnb, csnbn] = BubbleSortArray(csnb, csnbn, 1, 1);

          for i = 2:csebn

%       Shuffle nodes from element(i) into csnb
               elem = cseb(i);
               numn = cis(2,elem);
             for j = 1:numn
	            ixe(j) = NodesOnElement(elem,j);
             end
	         [ixe, numn] = BubbleSortArray(ixe, numn, 1, 1);
	         [csnb,csnbn] = ShuffleIn(ixe,csnb,numn,csnbn);

%     End loop through elements in cseb
          end

%     Combine boundary csnb into c1sn and get csn

          [csn,csnn,csnb,csnbn] = Purge(csnb,c1sn,csnbn,c1snn);
       
       
%% Set q-values
QValues = zeros(numnp,1);
DM_interior_nodes = c1sn(1:c1snn); % Nodes on interior of Domain
QValues(DM_interior_nodes) = 1;


%% Form J-integral
for elemDM = 1:csen
    
  elem = cse(elemDM);
    
  for ma = 1:nummat
      
   if(ieFEAP(nie-2,ma) == RegionOnElement(elem))
      
    %Extract patch material properties
    iel = MatTypeTable(2,ma); %iel   = ie(nie-1,ma); same thing;
    nonlin = MatTypeTable(3,ma);
    mateprop = MateT(ma,:);
    
%     if (~((isw== 3||isw==6||isw==12) && nonlin==0) || (exist('initializeLinKF','var') && initializeLinKF))
    
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

    ElemQ = QValues(ElemFlag)';
    if nel == 9 % make a smooth profile, helps the values converge faster
            ElemQ(5) = (ElemQ(1) + ElemQ(2))/2;
            ElemQ(6) = (ElemQ(2) + ElemQ(3))/2;
            ElemQ(7) = (ElemQ(3) + ElemQ(4))/2;
            ElemQ(8) = (ElemQ(1) + ElemQ(4))/2;
            ElemQ(9) = 0.50d0*(ElemQ(5) + ElemQ(6) + ElemQ(7) + ElemQ(8)) ...
                    - 0.25d0*(ElemQ(1) + ElemQ(2) + ElemQ(3) + ElemQ(4));
    end
    
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
    J_int = J_int + ElemJ;

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
    
%     end % if nonlin
    
   end %if ma
    
  end % ma
    
end % elem

    JstepList(1:ndm,cstar) = J_int;

%   End while loop on nstar 
       end

hflgu = 0;
h3flgu = 0;

end % ndm
