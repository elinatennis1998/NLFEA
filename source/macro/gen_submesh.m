%**********************************************************************
% 
%       subroutine gen_submesh(epnum,epatch,xs,ixs,x,ix,cis,ecn)
%
%...  Written by Tim Truster (Spring 2010)
%     No modifications when copied to NLFEA ver2
%...  Program to generate submesh coordinates and connectivities
%
%**********************************************************************

% 		implicit none
% 
%       integer         numnp,numel,nummat,nen,neq,ipr
%       common /cdata/  numnp,numel,nummat,nen,neq,ipr
% 
%       real*8          dm
%       integer         n,ma,mct,iel,nel
%       common /eldata/ dm,n,ma,mct,iel,nel
% 
%       integer         ior,iow
%       common /iofile/ ior,iow
% 
%       integer         nh1, nh2
%       common /hdata/  nh1, nh2
% 
%       logical debug
%       common /debugs/ debug
% 
%       integer         ndf,ndm,nen1,nst,nneq,ndl,nnlm,nadd
%       common /sdata/  ndf,ndm,nen1,nst,nneq,ndl,nnlm,nadd
% 
% 	common /counter/istep
% 	integer istep
% 			 
% 	include 'smparam.h'
% 	include 'stardata.h'
% 	include 'submeshdata.h'
% 
% %     Input Variables
%       integer ix(nen1,numel),epatch(maxel,numnp),epnum(numnp)
%       real*8 x(ndm,numnp)
% 
% %     Output Variables
%       real*8 xs(ndm,numel*celn)
% 	integer ixs(nen,numel*cel),ecn(numel*celn),cis(2,numel+1)
% 
% %     Local Variables
% 	integer i,j,k,l,node,elem,cell,mat,nele
% 	integer ixe(nen),ixc(nen+1,celn)
% 	real*8 xe(ndm,nen),xc(ndm,celn)
% 
% 	integer nloop(4,2),edge,nodeA,nodeB,elemA,elemB,numA,numB
% 	logical efound,copyn
% 	integer el2mesh(celn),nodeid,enode,nA1,nAinc,nA2,nA3
% 	integer nB1,nBinc,nB2,neleB,nume
% 	
%       logical dprt

  ixe = zeros(nen,1);
  ixc = zeros(nen+1,celn);
  xe = zeros(ndm,nen);
  xc = zeros(ndm,celn);
  nloop = zeros(4,2);
  el2mesh = zeros(celn,1);

% % Transfer coarse nodes      
% 	do k = 1,numnp
% 	   do j = 1,ndm
% 	      xs(j,k) = x(j,k)
% 	   enddo
% 	enddo
% 
% 	nodeid = numnp
% 	cell = 0
% 	cis(1,1) = 1
% 
% 	if(debug) then
% 	dprt = .true.
% 	if(dprt) then
% 	   write(24,*) 'nodeid',nodeid
% 	endif
% 	endif
% 
% % Loop over elements in domain
%       do elem = 1,numel
% 
% %   Transfer element properties to cell generator
% 	   mat = ix(nen1,elem)
% 	   if(nen.eq.9) then
% 	      if(ix(9,elem).eq.0) then
% 	         nele = 6
% 	         ncel = (2*minc+2)*(2*minc+1)/2
% 	      else
% 	         nele = 9
% 	         ncel = (2*minc+1)*(2*minc+1)
% 	      endif
% 	   elseif(nen.eq.6) then
% 	      nele = 6
% 	      ncel = (2*minc+2)*(2*minc+1)/2
% 	   elseif(nen.eq.4) then
% 	      if(ix(4,elem).eq.0) then
% 	         nele = 3
% 	         ncel = (minc+2)*(minc+1)/2 %  need to adjust this?
% 	      else
% 	         nele = 4
% 	         ncel = (minc+1)*(minc+1)
% 	      endif
% 	   else
% 	      nele = 3
% 	      ncel = (minc+2)*(minc+1)/2 %  need to adjust this?
% 	   endif
% 	   do i = 1,nele
% 	      ixe(i) = i
% 	   enddo
% 	   do i = nele+1,nen
% 	      ixe(i) = 0
% 	   enddo
% 	   do i = 1,nele
% 	      node = ix(i,elem)
% 	      do k = 1,ndm
% 	         xe(k,i) = x(k,node)
% 	      enddo
% 	   enddo
% 
% 	if(debug) then
% 	dprt = .true.
% 	if(dprt) then
% 
% 	   write(24,'(A,I4)') 'nele',nele
% 	   write(24,'(A,I4)') 'nen',nen
% 	   WRITE(24,*) 'ndm',ndm
% 	   WRITE(24,*) 'cel',cel
% 	   WRITE(24,*) 'nen1',nen1
% 	   WRITE(24,*) 'ncel',ncel
% 	   WRITE(24,*) 'celn',celn
% 	   WRITE(24,*) 'mat',mat
% 	   write(24,'(A)') ''  
% 	   write(24,'(A)') 'xe'
% 	   do i=1,nen	
% 		
% 		   WRITE(24,*) (xe(l,i),l=1,ndm)    
% 
% 	   enddo
% 	   write(24,'(A)') 'ixe'
% 	   do i=1,nen	
% 		
% 		   WRITE(24,*) ixe(i)    
% 
% 	   enddo
% 	   write(24,'(A)') ''
% 
% 	endif
% 	endif
% 
% %   Create cell nodes, node table, and ix array for element
% 	   call elemmesh2d(minc,xe,ixe,ndm,nele,nen,mat,cel,celn,ncel,
%      &	               xc,ixc)
% 
% %   Assemble cell data into submesh domain list
% 	   enode = cis(1,elem)
%          cis(1,elem+1) = ncel + enode
% 	   cis(2,elem) = nele
% 	   enode = enode - 1
% 
% 	   call pzeroi(el2mesh,celn)
% 
% 	if(debug) then
% 	dprt = .true.
% 	if(dprt) then
% 	   write(24,*) 'enode',enode
% 	endif
% 	endif

%%
% Transfer coarse nodes      
	for k = 1:numnp
       for j = 1:ndm
	      xs(j,k) = x(j,k); %#ok<AGROW>
       end
	end

	nodeid = numnp;
	cell = 0;
	cis(1,1) = 1;

% Loop over elements in domain
      for elem = 1:numel

%   Transfer element properties to cell generator
	   mat = RegionOnElement(elem);
       if(nen==9)
          if(NodesOnElement(9,elem)==0)
	         nele = 6;
	         ncel = floor((2*minc+2)*(2*minc+1)/2);
	      else
	         nele = 9;
	         ncel = floor((2*minc+1)*(2*minc+1));
          end
	   elseif(nen==6)
	      nele = 6;
	      ncel = floor((2*minc+2)*(2*minc+1)/2);
	   elseif(nen==4)
          if(NodesOnElement(4,elem)==0)
	         nele = 3;
	         ncel = floor((minc+2)*(minc+1)/2);
	      else
	         nele = 4;
	         ncel = (minc+1)*(minc+1);
          end
	   else
	      nele = 3;
	      ncel = floor((minc+2)*(minc+1)/2);
       end
       for i = 1:nele
	      ixe(i) = i;
       end
       for i = nele+1:nen
	      ixe(i) = 0;
       end
       for i = 1:nele
	      node = NodesOnElement(i,elem);
          for k = 1:ndm
	         xe(k,i) = x(k,node);
          end
       end

%   Create cell nodes, node table, and ix array for element
	   [xc,ixc] = elemmesh2d(minc,xe,ixe,ndm,nele,nen,mat,cel,celn,ncel);

%   Assemble cell data into submesh domain list
	   enode = cis(1,elem);
         cis(1,elem+1) = ncel + enode;
	   cis(2,elem) = nele;
	   enode = enode - 1;

	   el2mesh = zeros(celn,1);
    
% %     If T3 Element
% 	   if(nele.eq.3) then
% 
% 	      el2mesh(1) = NodesOnElement(1,elem)
% 	      el2mesh(minc+1) = NodesOnElement(2,elem)
% 	      el2mesh((minc+1)*minc+1) = NodesOnElement(3,elem)
% 
% 	      do i=1,3
% 	      nloop(i,1) = i
% 	      nloop(i,2) = i+1
% 	         enode = enode + 1
% 	         ecn(enode) = NodesOnElement(i,elem)
% 	      enddo
% 	      nloop(3,2) = 1
% 
% %     Loop over element edges
%             do edge = 1,3
% 
% 	         nodeA = NodesOnElement(nloop(edge,1),elem)
% 		     nodeB = NodesOnElement(nloop(edge,2),elem)
% 	         numA = epnum(nodeA)
% 	         numB = epnum(nodeB)
% 
% 	if(debug) then
% 	dprt = .false.
% 	if(dprt) then
% 	   write(24,*) 'nodeA',nodeA
% 	   write(24,*) 'nodeB',nodeB
% 	   write(24,*) 'numA',numA
% 	   write(24,*) 'numB',numB
% 	endif
% 	endif
% 
% 	         efound = .false.
% 	         copyn = .false.
% 	         k = 0
% 
% %       Search for shared edge in level 1 stars of nodes A & B
%                do 10 while((.not.efound).and.(k.lt.numA))
% 
% 	            k = k + 1
% 	            elemA = epatch(k,nodeA)
% 	if(debug) then
% 	dprt = .false.
% 	if(dprt) then
% 	   write(24,*) 'elemA',elemA
% 	endif
% 	endif
% 	            call binsearch(elemA,epatch(1,nodeB),1,numB,elemB
%      *			               ,nume,efound)
% 	 if(debug) then
% 	dprt = .false.
% 	if(dprt) then
% 	   write(24,*) 'elemB',elemB
% 	endif
% 	endif
% 	            if(efound.and.(elemB.ne.elem)) then
% 	              efound = .true.
% 	            else
% 	              efound = .false.
% 	            endif
% 
% 
% 10	         continue
% 
%       if(debug) then
% 	dprt = .false.
% 	if(dprt) then
%          write(24,*) 'efound ',efound
% 	endif
% 	endif
% 
% %       If shared, then if elem_curr > elem_old
% 	         if(efound) then
% 	            if(elem.gt.elemB) then
% 	               copyn = .true.
% 	            endif
% 	         endif
% 
%       if(debug) then
% 	dprt = .false.
% 	if(dprt) then
%          write(24,*) 'copyn ',copyn
% 	endif
% 	endif
% 
% 	         if(edge.eq.1) then
% 	             nA1 = 2
% 	             nAinc = 1
% 	             nA2 = minc
% 	         elseif(edge.eq.2) then
% 	             nA1 = 2*minc+1
% 	             nAinc = minc
% 	             nA2 = (minc-1)*(minc+1)+2
% 	         else
% 	             nA1 = (minc-1)*(minc+1)+1
% 	             nAinc = -minc-1
% 	             nA2 = minc+2
% 	         endif
% 
% 	if(debug) then
% 	dprt = .false.
% 	if(dprt) then
% 	   write(24,*) 'nA1',nA1
% 	   write(24,*) 'nAinc',nAinc
% 	   write(24,*) 'nA2',nA2
% 	   write(24,'(A)') ''
% 	endif
% 	endif
% 
% 	         if(copyn) then
% %         Copy node IDs into element table
%                   node = NodesOnElement(1,elemB)
% 	            neleB = cis(2,elemB)
% 	            k = 1
% %         Find which corner nodeA is on elemB
%                   do while((nodeA.ne.node).and.(k.lt.neleB))
% 	               k = k + 1
% 	               node = NodesOnElement(k,elemB)
% 	            enddo
% 
% 	if(debug) then
% 	dprt = .false.
% 	if(dprt) then
% 	   write(24,*) 'k',k
% 	   write(24,'(A)') ''
% 	endif
% 	endif
% 
% %         Set extraction indices for copying nodes
%                   if(neleB.eq.3) then
% 	            if(k.lt.neleB) then
% 	               if(nodeB.eq.NodesOnElement(k+1,elemB)) then
% 	                  %Forward edge
% 	                  nB1 = cis(1,elemB) + neleB + (k-1)*(minc-1)
% 	                  nBinc = 1
% 	                  nB2 = nB1 + minc - 2
% 	               else
% 	                  %Backward edge
% 	                  if(k.eq.1) then
% 	                     k = 4
% 	                  endif
% 	                  nB2 = cis(1,elemB) + neleB + (k-2)*(minc-1)
% 	                  nBinc = -1
% 	                  nB1 = nB2 + minc - 2
% 	               endif
% 	            else %k = neleB
% 	               if(nodeB.eq.NodesOnElement(1,elemB)) then
% 	                  %Forward edge
% 	                  nB1 = cis(1,elemB) + neleB + (k-1)*(minc-1)
% 	                  nBinc = 1
% 	                  nB2 = nB1 + minc - 2
% 	               else
% 	                  %Backward edge
% 	                  nB2 = cis(1,elemB) + neleB + (k-2)*(minc-1)
% 	                  nBinc = -1
% 	                  nB1 = nB2 + minc - 2
% 	               endif
% 	            endif
% 	            elseif(neleB.eq.4) then
% 	            if(k.lt.neleB) then
% 	               if(nodeB.eq.NodesOnElement(k+1,elemB)) then
% 	                  %Forward edge
% 	                  nB1 = cis(1,elemB) + neleB + (k-1)*(minc-1)
% 	                  nBinc = 1
% 	                  nB2 = nB1 + minc - 2
% 	               else
% 	                  %Backward edge
% 	                  if(k.eq.1) then
% 	                     k = 5
% 	                  endif
% 	                  nB2 = cis(1,elemB) + neleB + (k-2)*(minc-1)
% 	                  nBinc = -1
% 	                  nB1 = nB2 + minc - 2
% 	               endif
% 	            else %k = neleB
% 	               if(nodeB.eq.NodesOnElement(1,elemB)) then
% 	                  %Forward edge
% 	                  nB1 = cis(1,elemB) + neleB + (k-1)*(minc-1)
% 	                  nBinc = 1
% 	                  nB2 = nB1 + minc - 2
% 	               else
% 	                  %Backward edge
% 	                  nB2 = cis(1,elemB) + neleB + (k-2)*(minc-1)
% 	                  nBinc = -1
% 	                  nB1 = nB2 + minc - 2
% 	               endif
% 	            endif
% 	            endif
% 
% 	if(debug) then
% 	dprt = .false.
% 	if(dprt) then
% 	   write(24,*) 'nB1',nB1
% 	   write(24,*) 'nBinc',nBinc
% 	   write(24,*) 'nB2',nB2
% 	   write(24,'(A)') ''
% 	endif
% 	endif
% 
% %         Load nodes from elemB edge
%       %            if(nB1.eq.nB2) then
% 	%               node = ecn(nB1)
% 	%               el2mesh(nA1) = node
% 	%               enode = enode + 1
% 	%			   ecn(enode) = node	
% 	%            else
% 	               k = nA1
%                      do j = nB1,nB2,nBinc
% 	                  node = ecn(j)
% 	                  el2mesh(k) = node
% 	                  enode = enode + 1
% 				      ecn(enode) = node
% 	                  k = k + nAinc
% 	               enddo
% 	%            endif
% 
% %       Else, add node IDs to node table, store coordinates
%                else
% 
% 	%            if(nA1.eq.nA2) then
% 	%               nodeid = nodeid + 1
% 	%               el2mesh(nA1) = nodeid
% 	%               enode = enode + 1
% 	%			   ecn(enode) = nodeid	
% 	 %           else
%                      do j = nA1,nA2,nAinc
% 	                  nodeid = nodeid + 1
% 	                  el2mesh(j) = nodeid
% 	                  enode = enode + 1
% 				      ecn(enode) = nodeid
% 				      do i = 1,ndm
% 				         xs(i,nodeid) = xc(i,j)
% 				      enddo
% 	               enddo
% 	%            endif
% 
% 	         endif
% 
% %     End loop over edges
%             enddo
% 
% %     Store interior node IDs and coordinates
%             k = minc+1
% 	      do j = 2,minc
% 	         k = k + 1
% 	         do i = 2,minc+1-j
% 		        nodeid = nodeid + 1
% 	            k = k + 1
% 	            el2mesh(k) = nodeid
% 	            enode = enode + 1
% 				ecn(enode) = nodeid
% 				do l = 1,ndm
% 				   xs(l,nodeid) = xc(l,k)
% 				enddo
% 	         enddo
% 	         do i = minc+2-j,minc
% 	            k = k + 1
% 	         enddo
% 	         k = k + 1
% 	      enddo
% 
% 	if(debug) then
% 	dprt = .false.
% 	if(dprt) then
% 	   write(24,*) 'el2mesh'
% 	   do j = 1,celn
%          write(24,*) el2mesh(j)
% 	   enddo
% 	   write(24,'(A)') ''
% 	endif
% 	endif
% 
% %     Transfer data from element array to submesh array
%             l = 0
% 	      do k = 1,2*minc
% 	         do j = 1,2*minc-1-(k-1)*2
% 	            l = l + 1
% 	            cell = cell + 1
% 	            do i = 1,nele
% 	               node = ixc(i,l)
% 	               node = el2mesh(node)
% 	               ixs(i,cell) = node
% 	            enddo
% 	         enddo
% 	         do j = 2*minc-(k-1)*2,2*minc
% 	            l = l + 1
% 	         enddo
% 	      enddo

%%
%     If T3 Element
       if(nele==3)

	      el2mesh(1) = NodesOnElement(1,elem);
	      el2mesh(minc+1) = NodesOnElement(2,elem);
	      el2mesh((minc+1)*minc+1) = NodesOnElement(3,elem);

          for i=1:3
	      nloop(i,1) = i;
	      nloop(i,2) = i+1;
	         enode = enode + 1;
	         ecn(enode) = NodesOnElement(i,elem); %#ok<AGROW>
          end
	      nloop(3,2) = 1;

%     Loop over element edges
            for edge = 1:3

	         nodeA = NodesOnElement(nloop(edge,1),elem);
		     nodeB = NodesOnElement(nloop(edge,2),elem);
	         numA = epnum(nodeA);
	         numB = epnum(nodeB);

	         efound = 0;
	         copyn = 0;
	         k = 0;

%       Search for shared edge in level 1 stars of nodes A & B
               while((efound==0)&&(k<numA))

	            k = k + 1;
	            elemA = epatch(k,nodeA);
	            [elemB,nume,efound] = binsearch(elemA,epatch(:,nodeB),1,numB);
                if(efound&&(elemB~=elem))
	              efound = 1;
	            else
	              efound = 0;
                end


               end

%       If shared, then if elem_curr > elem_old
             if(efound)
                if(elem>elemB)
	               copyn = 1;
                end
             end

             if(edge==1)
	             nA1 = 2;
	             nAinc = 1;
	             nA2 = minc;
	         elseif(edge==2)
	             nA1 = 2*minc+1;
	             nAinc = minc;
	             nA2 = (minc-1)*(minc+1)+2;
	         else
	             nA1 = (minc-1)*(minc+1)+1;
	             nAinc = -minc-1;
	             nA2 = minc+2;
             end

             if(copyn)
%         Copy node IDs into element table
                  node = NodesOnElement(1,elemB);
	            neleB = cis(2,elemB);
	            k = 1;
%         Find which corner nodeA is on elemB
                  while((nodeA~=node)&&(k<neleB))
	               k = k + 1;
	               node = NodesOnElement(k,elemB);
                  end

%         Set extraction indices for copying nodes
               if(neleB==3)
                 if(k<neleB)
                   if(nodeB==NodesOnElement(k+1,elemB))
	                  %Forward edge
	                  nB1 = cis(1,elemB) + neleB + (k-1)*(minc-1);
	                  nBinc = 1;
	                  nB2 = nB1 + minc - 2;
	               else
	                  %Backward edge
                      if(k==1)
	                     k = 4;
                      end
	                  nB2 = cis(1,elemB) + neleB + (k-2)*(minc-1);
	                  nBinc = -1;
	                  nB1 = nB2 + minc - 2;
                   end
	            else %k = neleB
                   if(nodeB==NodesOnElement(1,elemB))
	                  %Forward edge
	                  nB1 = cis(1,elemB) + neleB + (k-1)*(minc-1);
	                  nBinc = 1;
	                  nB2 = nB1 + minc - 2;
	               else
	                  %Backward edge
	                  nB2 = cis(1,elemB) + neleB + (k-2)*(minc-1);
	                  nBinc = -1;
	                  nB1 = nB2 + minc - 2;
                   end
                end
	            elseif(neleB==4)
                  if(k<neleB)
                    if(nodeB==NodesOnElement(k+1,elemB))
	                  %Forward edge
	                  nB1 = cis(1,elemB) + neleB + (k-1)*(minc-1);
	                  nBinc = 1;
	                  nB2 = nB1 + minc - 2;
	                else
	                  %Backward edge
                      if(k==1)
	                     k = 5;
                      end
	                  nB2 = cis(1,elemB) + neleB + (k-2)*(minc-1);
	                  nBinc = -1;
	                  nB1 = nB2 + minc - 2;
                    end
	              else %k = neleB
                    if(nodeB==NodesOnElement(1,elemB))
	                  %Forward edge
	                  nB1 = cis(1,elemB) + neleB + (k-1)*(minc-1);
	                  nBinc = 1;
	                  nB2 = nB1 + minc - 2;
	                else
	                  %Backward edge
	                  nB2 = cis(1,elemB) + neleB + (k-2)*(minc-1);
	                  nBinc = -1;
	                  nB1 = nB2 + minc - 2;
                    end
                  end
               end


%         Load nodes from elemB edge
      %            if(nB1.eq.nB2) then
	%               node = ecn(nB1)
	%               el2mesh(nA1) = node
	%               enode = enode + 1
	%			   ecn(enode) = node	
	%            else
	               k = nA1;
                   for j = nB1:nBinc:nB2
	                  node = ecn(j);
	                  el2mesh(k) = node;
	                  enode = enode + 1;
				      ecn(enode) = node; %#ok<AGROW>
	                  k = k + nAinc;
                   end
	%            endif

%       Else, add node IDs to node table, store coordinates
               else

	%            if(nA1.eq.nA2) then
	%               nodeid = nodeid + 1
	%               el2mesh(nA1) = nodeid
	%               enode = enode + 1
	%			   ecn(enode) = nodeid	
	 %           else
                    for j = nA1:nAinc:nA2
	                  nodeid = nodeid + 1;
	                  el2mesh(j) = nodeid;
	                  enode = enode + 1;
				      ecn(enode) = nodeid; %#ok<AGROW>
                      for i = 1:ndm
				         xs(i,nodeid) = xc(i,j); %#ok<AGROW>
                      end
                    end
	%            endif

             end

%     End loop over edges
            end

%     Store interior node IDs and coordinates
           k = minc+1;
           for j = 2:minc
	         k = k + 1;
             for i = 2:minc+1-j
		        nodeid = nodeid + 1;
	            k = k + 1;
	            el2mesh(k) = nodeid;
	            enode = enode + 1;
				ecn(enode) = nodeid; %#ok<AGROW>
                for l = 1:ndm
				   xs(l,nodeid) = xc(l,k); %#ok<AGROW>
                end
             end
             for i = minc+2-j:minc
	            k = k + 1;
             end
	         k = k + 1;
           end

%     Transfer data from element array to submesh array
          l = 0;
          for k = 1:2*minc
             for j = 1:2*minc-1-(k-1)*2
	            l = l + 1;
	            cell = cell + 1;
                for i = 1:nele
	               node = ixc(i,l);
	               node = el2mesh(node);
	               ixs(i,cell) = node; %#ok<AGROW>
                end
             end
             for j = 2*minc-(k-1)*2:2*minc
	            l = l + 1;
             end
          end


% %     If Q4 Element
% 	   elseif(nele.eq.4) then
% 
% 	      el2mesh(1) = NodesOnElement(1,elem)
% 	      el2mesh(minc+1) = NodesOnElement(2,elem)
% 	      el2mesh((minc+1)*(minc+1)) = NodesOnElement(3,elem)
% 	      el2mesh((minc+1)*minc+1) = NodesOnElement(4,elem)
% 
% 	      do i=1,4
% 	      nloop(i,1) = i
% 	      nloop(i,2) = i+1
% 	         enode = enode + 1
% 	         ecn(enode) = NodesOnElement(i,elem)
% 	      enddo
% 	      nloop(4,2) = 1
% 
% %     Loop over element edges
%             do edge = 1,4
% 
% 	         nodeA = NodesOnElement(nloop(edge,1),elem)
% 		     nodeB = NodesOnElement(nloop(edge,2),elem)
% 	         numA = epnum(nodeA)
% 	         numB = epnum(nodeB)
% 
% 	if(debug) then
% 	dprt = .false.
% 	if(dprt) then
% 	   write(24,*) 'nodeA',nodeA
% 	   write(24,*) 'nodeB',nodeB
% 	   write(24,*) 'numA',numA
% 	   write(24,*) 'numB',numB
% 	endif
% 	endif
% 
% 	         efound = .false.
% 	         copyn = .false.
% 	         k = 0
% 
% %       Search for shared edge in level 1 stars of nodes A & B
%                do 20 while((.not.efound).and.(k.lt.numA))
% 
% 	            k = k + 1
% 	            elemA = epatch(k,nodeA)
% 	if(debug) then
% 	dprt = .false.
% 	if(dprt) then
% 	   write(24,*) 'elemA',elemA
% 	endif
% 	endif
% 	            call binsearch(elemA,epatch(1,nodeB),1,numB,elemB
%      *			               ,nume,efound)
% 	if(debug) then
% 	dprt = .false.
% 	if(dprt) then
% 	   write(24,*) 'elemB',elemB
% 	endif
% 	endif
% 	            if(efound.and.(elemB.ne.elem)) then
% 	              efound = .true.
% 	            else
% 	              efound = .false.
% 	            endif
% 
% 
% 20	         continue
% 
%       if(debug) then
% 	dprt = .false.
% 	if(dprt) then
%          write(24,*) 'efound ',efound
% 	endif
% 	endif
% 
% %       If shared, then if elem_curr > elem_old
% 	         if(efound) then
% 	            if(elem.gt.elemB) then
% 	               copyn = .true.
% 	            endif
% 	         endif
% 
%       if(debug) then
% 	dprt = .false.
% 	if(dprt) then
%          write(24,*) 'copyn ',copyn
% 	endif
% 	endif
% 
% 	         if(edge.eq.1) then
% 	             nA1 = 2
% 	             nAinc = 1
% 	             nA2 = minc
% 	         elseif(edge.eq.2) then
% 	             nA1 = 2*(minc+1)
% 	             nAinc = minc+1
% 	             nA2 = minc*(minc+1)
% 	         elseif(edge.eq.3) then
% 	             nA1 = (minc+1)*(minc+1)-1
% 	             nAinc = -1
% 	             nA2 = minc*(minc+1)+2
% 	         else
% 	             nA1 = (minc-1)*(minc+1)+1
% 	             nAinc = -minc-1
% 	             nA2 = minc+2
% 	         endif
% 
% 	if(debug) then
% 	dprt = .false.
% 	if(dprt) then
% 	   write(24,*) 'nA1',nA1
% 	   write(24,*) 'nAinc',nAinc
% 	   write(24,*) 'nA2',nA2
% 	   write(24,'(A)') ''
% 	endif
% 	endif
% 
% 	         if(copyn) then
% %         Copy node IDs into element table
%                   node = NodesOnElement(1,elemB)
% 	            neleB = cis(2,elemB)
% 	            k = 1
% %         Find which corner nodeA is on elemB
%                   do while((nodeA.ne.node).and.(k.lt.neleB))
% 	               k = k + 1
% 	               node = NodesOnElement(k,elemB)
% 	            enddo
% 
% 	if(debug) then
% 	dprt = .false.
% 	if(dprt) then
% 	   write(24,*) 'k',k
% 	   write(24,'(A)') ''
% 	endif
% 	endif
% 
% %         Set extraction indices for copying nodes
% 	            if(neleB.eq.3) then
% 	            if(k.lt.neleB) then
% 	               if(nodeB.eq.NodesOnElement(k+1,elemB)) then
% 	                  %Forward edge
% 	                  nB1 = cis(1,elemB) + neleB + (k-1)*(minc-1)
% 	                  nBinc = 1
% 	                  nB2 = nB1 + minc - 2
% 	               else
% 	                  %Backward edge
% 	                  if(k.eq.1) then
% 	                     k = 4
% 	                  endif
% 	                  nB2 = cis(1,elemB) + neleB + (k-2)*(minc-1)
% 	                  nBinc = -1
% 	                  nB1 = nB2 + minc - 2
% 	               endif
% 	            else %k = neleB
% 	               if(nodeB.eq.NodesOnElement(1,elemB)) then
% 	                  %Forward edge
% 	                  nB1 = cis(1,elemB) + neleB + (k-1)*(minc-1)
% 	                  nBinc = 1
% 	                  nB2 = nB1 + minc - 2
% 	               else
% 	                  %Backward edge
% 	                  nB2 = cis(1,elemB) + neleB + (k-2)*(minc-1)
% 	                  nBinc = -1
% 	                  nB1 = nB2 + minc - 2
% 	               endif
% 	            endif
% 	            elseif(neleB.eq.4) then
% 	            if(k.lt.neleB) then
% 	               if(nodeB.eq.NodesOnElement(k+1,elemB)) then
% 	                  %Forward edge
% 	                  nB1 = cis(1,elemB) + neleB + (k-1)*(minc-1)
% 	                  nBinc = 1
% 	                  nB2 = nB1 + minc - 2
% 	               else
% 	                  %Backward edge
% 	                  if(k.eq.1) then
% 	                     k = 5
% 	                  endif
% 	                  nB2 = cis(1,elemB) + neleB + (k-2)*(minc-1)
% 	                  nBinc = -1
% 	                  nB1 = nB2 + minc - 2
% 	               endif
% 	            else %k = neleB
% 	               if(nodeB.eq.NodesOnElement(1,elemB)) then
% 	                  %Forward edge
% 	                  nB1 = cis(1,elemB) + neleB + (k-1)*(minc-1)
% 	                  nBinc = 1
% 	                  nB2 = nB1 + minc - 2
% 	               else
% 	                  %Backward edge
% 	                  nB2 = cis(1,elemB) + neleB + (k-2)*(minc-1)
% 	                  nBinc = -1
% 	                  nB1 = nB2 + minc - 2
% 	               endif
% 	            endif
% 	            endif
% 
% 	if(debug) then
% 	dprt = .false.
% 	if(dprt) then
% 	   write(24,*) 'nB1',nB1
% 	   write(24,*) 'nBinc',nBinc
% 	   write(24,*) 'nB2',nB2
% 	   write(24,'(A)') ''
% 	endif
% 	endif
% %         Load nodes from elemB edge
%       %            if(nB1.eq.nB2) then
% 	%               node = ecn(nB1)
% 	%               el2mesh(nA1) = node
% 	%               enode = enode + 1
% 	%			   ecn(enode) = node	
% 	%            else
% 	               k = nA1
%                      do j = nB1,nB2,nBinc
% 	                  node = ecn(j)
% 	                  el2mesh(k) = node
% 	                  enode = enode + 1
% 				      ecn(enode) = node
% 	                  k = k + nAinc
% 	               enddo
% 	%            endif
% 
% %       Else, add node IDs to node table, store coordinates
%                else
% 
% 	%            if(nA1.eq.nA2) then
% 	%               nodeid = nodeid + 1
% 	%               el2mesh(nA1) = nodeid
% 	%               enode = enode + 1
% 	%			   ecn(enode) = nodeid	
% 	 %           else
%                      do j = nA1,nA2,nAinc
% 	                  nodeid = nodeid + 1
% 	                  el2mesh(j) = nodeid
% 	                  enode = enode + 1
% 				      ecn(enode) = nodeid
% 				      do i = 1,ndm
% 				         xs(i,nodeid) = xc(i,j)
% 				      enddo
% 	               enddo
% 	%            endif
% 
% 	         endif
% 
% %     End loop over edges
%             enddo
% 
% %     Store interior node IDs and coordinates
%             k = minc+1
% 	      do j = 2,minc
% 	         k = k + 1
% 	         do i = 2,minc
% 		        nodeid = nodeid + 1
% 	            k = k + 1
% 	            el2mesh(k) = nodeid
% 	            enode = enode + 1
% 				ecn(enode) = nodeid
% 				do l = 1,ndm
% 				   xs(l,nodeid) = xc(l,k)
% 				enddo
% 	         enddo
% 	         k = k + 1
% 	      enddo
% 
% 	if(debug) then
% 	dprt = .false.
% 	if(dprt) then
% 	   write(24,*) 'el2mesh'
% 	   do j = 1,ncel
%          write(24,*) el2mesh(j)
% 	   enddo
% 	   write(24,'(A)') ''
% 	endif
% 	endif
% 
% %     Transfer data from element array to submesh array
% 	      do k = 1,cel
% 	         cell = cell + 1
% 	         do j = 1,nele
% 	            node = ixc(j,k)
% 	            node = el2mesh(node)
% 	            ixs(j,cell) = node
% 	         enddo
% 	      enddo

%%
%     If Q4 Element
	   elseif(nele==4)

	      el2mesh(1) = NodesOnElement(1,elem);
	      el2mesh(minc+1) = NodesOnElement(2,elem);
	      el2mesh((minc+1)*(minc+1)) = NodesOnElement(3,elem);
	      el2mesh((minc+1)*minc+1) = NodesOnElement(4,elem);

          for i=1:4
	         nloop(i,1) = i;
	         nloop(i,2) = i+1;
	         enode = enode + 1;
	         ecn(enode) = NodesOnElement(i,elem); %#ok<AGROW>
          end
	      nloop(4,2) = 1;

%     Loop over element edges
            for edge = 1:4

	         nodeA = NodesOnElement(nloop(edge,1),elem);
		     nodeB = NodesOnElement(nloop(edge,2),elem);
	         numA = epnum(nodeA);
	         numB = epnum(nodeB);

	         efound = 0;
	         copyn = 0;
	         k = 0;

%       Search for shared edge in level 1 stars of nodes A & B
              while((efound==0)&&(k<numA))

	            k = k + 1;
	            elemA = epatch(k,nodeA);
	            [elemB,nume,efound] = binsearch(elemA,epatch(:,nodeB),1,numB);
                if(efound&&(elemB~=elem))
	              efound = 1;
	            else
	              efound = 0;
                end


              end

%       If shared, then if elem_curr > elem_old
             if(efound)
                if(elem>elemB)
	               copyn = 1;
                end
             end

             if(edge==1)
	             nA1 = 2;
	             nAinc = 1;
	             nA2 = minc;
	         elseif(edge==2)
	             nA1 = 2*(minc+1);
	             nAinc = minc+1;
	             nA2 = minc*(minc+1);
	         elseif(edge==3)
	             nA1 = (minc+1)*(minc+1)-1;
	             nAinc = -1;
	             nA2 = minc*(minc+1)+2;
	         else
	             nA1 = (minc-1)*(minc+1)+1;
	             nAinc = -minc-1;
	             nA2 = minc+2;
             end

             if(copyn)
%         Copy node IDs into element table
                  node = NodesOnElement(1,elemB);
	            neleB = cis(2,elemB);
	            k = 1;
%         Find which corner nodeA is on elemB
                  while((nodeA~=node)&&(k<neleB))
	               k = k + 1;
	               node = NodesOnElement(k,elemB);
                  end

%         Set extraction indices for copying nodes
              if(neleB==3)
                if(k<neleB)
                  if(nodeB==NodesOnElement(k+1,elemB))
	                  %Forward edge
	                  nB1 = cis(1,elemB) + neleB + (k-1)*(minc-1);
	                  nBinc = 1;
	                  nB2 = nB1 + minc - 2;
	               else
	                  %Backward edge
                      if(k==1)
	                     k = 4;
                      end
	                  nB2 = cis(1,elemB) + neleB + (k-2)*(minc-1);
	                  nBinc = -1;
	                  nB1 = nB2 + minc - 2;
                  end
	            else %k = neleB
                   if(nodeB==NodesOnElement(1,elemB))
	                  %Forward edge
	                  nB1 = cis(1,elemB) + neleB + (k-1)*(minc-1);
	                  nBinc = 1;
	                  nB2 = nB1 + minc - 2;
	               else
	                  %Backward edge
	                  nB2 = cis(1,elemB) + neleB + (k-2)*(minc-1);
	                  nBinc = -1;
	                  nB1 = nB2 + minc - 2;
                   end
                end
	          elseif(neleB==4)
                if(k<neleB)
                   if(nodeB==NodesOnElement(k+1,elemB))
	                  %Forward edge
	                  nB1 = cis(1,elemB) + neleB + (k-1)*(minc-1);
	                  nBinc = 1;
	                  nB2 = nB1 + minc - 2;
	               else
	                  %Backward edge
                      if(k==1)
	                     k = 5;
                      end
	                  nB2 = cis(1,elemB) + neleB + (k-2)*(minc-1);
	                  nBinc = -1;
	                  nB1 = nB2 + minc - 2;
                   end
	            else %k = neleB
                   if(nodeB==NodesOnElement(1,elemB))
	                  %Forward edge
	                  nB1 = cis(1,elemB) + neleB + (k-1)*(minc-1);
	                  nBinc = 1;
	                  nB2 = nB1 + minc - 2;
	               else
	                  %Backward edge
	                  nB2 = cis(1,elemB) + neleB + (k-2)*(minc-1);
	                  nBinc = -1;
	                  nB1 = nB2 + minc - 2;
                   end
                end
              end

%         Load nodes from elemB edge
      %            if(nB1.eq.nB2) then
	%               node = ecn(nB1)
	%               el2mesh(nA1) = node
	%               enode = enode + 1
	%			   ecn(enode) = node	
	%            else
	               k = nA1;
                     for j = nB1:nBinc:nB2
	                  node = ecn(j);
	                  el2mesh(k) = node;
	                  enode = enode + 1;
				      ecn(enode) = node; %#ok<AGROW>
	                  k = k + nAinc;
                     end
	%            endif

%       Else, add node IDs to node table, store coordinates
               else

	%            if(nA1.eq.nA2) then
	%               nodeid = nodeid + 1
	%               el2mesh(nA1) = nodeid
	%               enode = enode + 1
	%			   ecn(enode) = nodeid	
	 %           else
                     for j = nA1:nAinc:nA2
	                  nodeid = nodeid + 1;
	                  el2mesh(j) = nodeid;
	                  enode = enode + 1;
				      ecn(enode) = nodeid; %#ok<AGROW>
                      for i = 1:ndm
				         xs(i,nodeid) = xc(i,j); %#ok<AGROW>
                      end
                     end
	%            endif

             end

%     End loop over edges
            end

%     Store interior node IDs and coordinates
          k = minc+1;
          for j = 2:minc
	         k = k + 1;
             for i = 2:minc
		        nodeid = nodeid + 1;
	            k = k + 1;
	            el2mesh(k) = nodeid;
	            enode = enode + 1;
				ecn(enode) = nodeid; %#ok<AGROW>
                for l = 1:ndm
				   xs(l,nodeid) = xc(l,k); %#ok<AGROW>
                end
             end
	         k = k + 1;
          end

%     Transfer data from element array to submesh array
          for k = 1:cel
	         cell = cell + 1;
             for j = 1:nele
	            node = ixc(j,k);
	            node = el2mesh(node);
	            ixs(j,cell) = node; %#ok<AGROW>
             end
          end
	


% %     If T6 Element
% 	   elseif(nele.eq.6) then
% 
% 	      el2mesh(1) = NodesOnElement(1,elem)
% 	      el2mesh(2*minc+1) = NodesOnElement(2,elem)
% 	      el2mesh((2*minc+1)*2*minc+1) = NodesOnElement(3,elem)
% 	      el2mesh(minc+1) = NodesOnElement(4,elem)
% 	      el2mesh((2*minc+1)*minc+minc+1) = NodesOnElement(5,elem)
% 	      el2mesh((2*minc+1)*minc+1) = NodesOnElement(6,elem)
% 
% 	if(debug) then
% 	dprt = .true.
% 	if(dprt) then
% 	   write(24,*) 'el2mesh'
% 	   do j = 1,celn
%          write(24,*) el2mesh(j)
% 	   enddo
% 	   write(24,'(A)') ''
% 	endif
% 	endif
% 
% 	      do i=1,3
% 	      nloop(i,1) = i
% 	      nloop(i,2) = i+1
% 	         enode = enode + 1
% 	         ecn(enode) = NodesOnElement(i,elem)
% 	      enddo
% 	      do i = 4,6
% 	         enode = enode + 1
% 	         ecn(enode) = NodesOnElement(i,elem)
% 	      enddo
% 	      nloop(3,2) = 1
% 
% %     Loop over element edges
%             do edge = 1,3
% 
% 	         nodeA = NodesOnElement(nloop(edge,1),elem)
% 		     nodeB = NodesOnElement(nloop(edge,2),elem)
% 	         numA = epnum(nodeA)
% 	         numB = epnum(nodeB)
% 
% 	if(debug) then
% 	dprt = .true.
% 	if(dprt) then
% 	   write(24,*) 'nodeA',nodeA
% 	   write(24,*) 'nodeB',nodeB
% 	   write(24,*) 'numA',numA
% 	   write(24,*) 'numB',numB
% 	endif
% 	endif
% 
% 	         efound = .false.
% 	         copyn = .false.
% 	         k = 0
% 
% %       Search for shared edge in level 1 stars of nodes A & B
%                do 30 while((.not.efound).and.(k.lt.numA))
% 
% 	            k = k + 1
% 	            elemA = epatch(k,nodeA)
% 	if(debug) then
% 	dprt = .true.
% 	if(dprt) then
% 	   write(24,*) 'elemA',elemA
% 	endif
% 	endif
% 	            call binsearch(elemA,epatch(:,nodeB),1,numB,elemB
%      *			               ,nume,efound)
% 	 if(debug) then
% 	dprt = .true.
% 	if(dprt) then
% 	   write(24,*) 'elemB',elemB
% 	endif
% 	endif
% 	            if(efound.and.(elemB.ne.elem)) then
% 	              efound = .true.
% 	            else
% 	              efound = .false.
% 	            endif
% 
% 
% 30	         continue
% 
%       if(debug) then
% 	dprt = .true.
% 	if(dprt) then
%          write(24,*) 'efound ',efound
% 	endif
% 	endif
% 
% %       If shared, then if elem_curr > elem_old
% 	         if(efound) then
% 	            if(elem.gt.elemB) then
% 	               copyn = .true.
% 	            endif
% 	         endif
% 
%       if(debug) then
% 	dprt = .true.
% 	if(dprt) then
%          write(24,*) 'copyn ',copyn
% 	endif
% 	endif
% 
% 	         if(edge.eq.1) then
% 	             nA1 = 2
% 	             nAinc = 1
% 	             nA2 = 2*minc
% 	             nA3 = minc+1
% 	         elseif(edge.eq.2) then
% 	             nA1 = 4*minc+1
% 	             nAinc = 2*minc
% 	             nA2 = (2*minc-1)*(2*minc+1)+2
% 	             nA3 = (2*minc+1)*minc+minc+1
% 	         else
% 	             nA1 = (2*minc-1)*(2*minc+1)+1
% 	             nAinc = -2*minc-1
% 	             nA2 = 2*minc+2
% 	             nA3 = (2*minc+1)*minc+1
% 	         endif
% 
% 	if(debug) then
% 	dprt = .true.
% 	if(dprt) then
% 	   write(24,*) 'nA1',nA1
% 	   write(24,*) 'nAinc',nAinc
% 	   write(24,*) 'nA2',nA2
% 	   write(24,*) 'nA3',nA3
% 	   write(24,'(A)') ''
% 	endif
% 	endif
% 
% 	         if(copyn) then
% %         Copy node IDs into element table
%                   node = NodesOnElement(1,elemB)
% 	            neleB = cis(2,elemB)
% 	            k = 1
% %         Find which corner nodeA is on elemB
%                   do while((nodeA.ne.node).and.(k.lt.neleB))
% 	               k = k + 1
% 	               node = NodesOnElement(k,elemB)
% 	            enddo
% 
% 	if(debug) then
% 	dprt = .true.
% 	if(dprt) then
% 	   write(24,*) 'k',k
% 	   write(24,'(A)') ''
% 	endif
% 	endif
% 
% %         Set extraction indices for copying nodes
%                   if(neleB.eq.6) then
% 	            if(k.lt.3) then
% 	               if(nodeB.eq.NodesOnElement(k+1,elemB)) then
% 	                  %Forward edge
% 	                  nB1 = cis(1,elemB) + neleB + (k-1)*(2*minc-2)
% 	                  nBinc = 1
% 	                  nB2 = nB1 + 2*minc - 3
% 	               else
% 	                  %Backward edge
% 	                  if(k.eq.1) then
% 	                     k = 4
% 	                  endif
% 	                  nB2 = cis(1,elemB) + neleB + (k-2)*(2*minc-2)
% 	                  nBinc = -1
% 	                  nB1 = nB2 + 2*minc - 3
% 	               endif
% 	            else %k = neleB
% 	               if(nodeB.eq.NodesOnElement(1,elemB)) then
% 	                  %Forward edge
% 	                  nB1 = cis(1,elemB) + neleB + (k-1)*(2*minc-2)
% 	                  nBinc = 1
% 	                  nB2 = nB1 + 2*minc - 3
% 	               else
% 	                  %Backward edge
% 	                  nB2 = cis(1,elemB) + neleB + (k-2)*(2*minc-2)
% 	                  nBinc = -1
% 	                  nB1 = nB2 + 2*minc - 3
% 	               endif
% 	            endif
% 	            elseif(neleB.eq.9) then
% 	            if(k.lt.4) then
% 	               if(nodeB.eq.NodesOnElement(k+1,elemB)) then
% 	                  %Forward edge
% 	                  nB1 = cis(1,elemB) + neleB + (k-1)*(2*minc-2)
% 	                  nBinc = 1
% 	                  nB2 = nB1 + 2*minc - 3
% 	               else
% 	                  %Backward edge
% 	                  if(k.eq.1) then
% 	                     k = 5
% 	                  endif
% 	                  nB2 = cis(1,elemB) + neleB + (k-2)*(2*minc-2)
% 	                  nBinc = -1
% 	                  nB1 = nB2 + 2*minc - 3
% 	               endif
% 	            else %k = neleB
% 	               if(nodeB.eq.NodesOnElement(1,elemB)) then
% 	                  %Forward edge
% 	                  nB1 = cis(1,elemB) + neleB + (k-1)*(2*minc-2)
% 	                  nBinc = 1
% 	                  nB2 = nB1 + 2*minc - 3
% 	               else
% 	                  %Backward edge
% 	                  nB2 = cis(1,elemB) + neleB + (k-2)*(2*minc-2)
% 	                  nBinc = -1
% 	                  nB1 = nB2 + 2*minc - 3
% 	               endif
% 	            endif
% 	            endif
% 
% 	if(debug) then
% 	dprt = .true.
% 	if(dprt) then
% 	   write(24,*) 'nB1',nB1
% 	   write(24,*) 'nBinc',nBinc
% 	   write(24,*) 'nB2',nB2
% 	   write(24,'(A)') ''
% 	endif
% 	endif
% 
% %         Load nodes from elemB edge
%       %            if(nB1.eq.nB2) then
% 	%               node = ecn(nB1)
% 	%               el2mesh(nA1) = node
% 	%               enode = enode + 1
% 	%			   ecn(enode) = node	
% 	%            else
% 	               k = nA1
%                      do j = nB1,nB2,nBinc
% 	                  if(k.eq.nA3) then
% 	                     k = k + nAinc
% 	                  endif
% 	                  node = ecn(j)
% 	                  el2mesh(k) = node
% 	                  enode = enode + 1
% 				      ecn(enode) = node
% 	                  k = k + nAinc
% 	               enddo
% 	%            endif
% 
% %       Else, add node IDs to node table, store coordinates
%                else
% 
% 	%            if(nA1.eq.nA2) then
% 	%               nodeid = nodeid + 1
% 	%               el2mesh(nA1) = nodeid
% 	%               enode = enode + 1
% 	%			   ecn(enode) = nodeid	
% 	 %           else
%                      do j = nA1,nA2,nAinc
% 	                  if(j.ne.nA3) then
% 	                     nodeid = nodeid + 1
% 	                     el2mesh(j) = nodeid
% 	                     enode = enode + 1
% 				         ecn(enode) = nodeid
% 				         do i = 1,ndm
% 				            xs(i,nodeid) = xc(i,j)
% 				         enddo
% 	                  endif
% 	               enddo
% 	%            endif
% 
% 	         endif
% 
% 	if(debug) then
% 	dprt = .true.
% 	if(dprt) then
% 	   write(24,*) 'el2mesh'
% 	   do j = 1,celn
%          write(24,*) el2mesh(j)
% 	   enddo
% 	   write(24,'(A)') ''
% 	endif
% 	endif
% 
% %     End loop over edges
%             enddo
% 
% %     Store interior node IDs and coordinates
%             k = 2*minc+1
% 	      do j = 2,2*minc
% 	         k = k + 1
% 	         do i = 2,2*minc+1-j
% 		        nodeid = nodeid + 1
% 	            k = k + 1
% 	            el2mesh(k) = nodeid
% 	            enode = enode + 1
% 				ecn(enode) = nodeid
% 				do l = 1,ndm
% 				   xs(l,nodeid) = xc(l,k)
% 				enddo
% 	         enddo
% 	         do i = 2*minc+2-j,2*minc
% 	            k = k + 1
% 	         enddo
% 	         k = k + 1
% 	      enddo
% 
% 	if(debug) then
% 	dprt = .true.
% 	if(dprt) then
% 	   write(24,*) 'el2mesh'
% 	   do j = 1,celn
%          write(24,*) el2mesh(j)
% 	   enddo
% 	   write(24,'(A)') ''
% 	endif
% 	endif
% 
% %     Transfer data from element array to submesh array
%             l = 0
% 	      do k = 1,2*minc
% 	         do j = 1,2*minc-1-(k-1)*2
% 	            l = l + 1
% 	            cell = cell + 1
% 	            do i = 1,nele
% 	               node = ixc(i,l)
% 	               node = el2mesh(node)
% 	               ixs(i,cell) = node
% 	            enddo
% 	         enddo
% 	         do j = 2*minc-(k-1)*2,2*minc
% 	            l = l + 1
% 	         enddo
% 	      enddo

%%
%     If T6 Element
	   elseif(nele==6)

	      el2mesh(1) = NodesOnElement(1,elem);
	      el2mesh(2*minc+1) = NodesOnElement(2,elem);
	      el2mesh((2*minc+1)*2*minc+1) = NodesOnElement(3,elem);
	      el2mesh(minc+1) = NodesOnElement(4,elem);
	      el2mesh((2*minc+1)*minc+minc+1) = NodesOnElement(5,elem);
	      el2mesh((2*minc+1)*minc+1) = NodesOnElement(6,elem);

          for i=1:3
	         nloop(i,1) = i;
	         nloop(i,2) = i+1;
	         enode = enode + 1;
	         ecn(enode) = NodesOnElement(i,elem); %#ok<AGROW>
          end
          for i = 4:6
	         enode = enode + 1;
	         ecn(enode) = NodesOnElement(i,elem); %#ok<AGROW>
          end
	      nloop(3,2) = 1;

%     Loop over element edges
            for edge = 1:3

	         nodeA = NodesOnElement(nloop(edge,1),elem);
		     nodeB = NodesOnElement(nloop(edge,2),elem);
	         numA = epnum(nodeA);
	         numB = epnum(nodeB);

	         efound = 0;
	         copyn = 0;
	         k = 0;

%       Search for shared edge in level 1 stars of nodes A & B
               while((efound==0)&&(k<numA))

	            k = k + 1;
	            elemA = epatch(k,nodeA);
	            [elemB,nume,efound] = binsearch(elemA,epatch(:,nodeB),1,numB);
                if(efound&&(elemB~=elem))
	              efound = 1;
	            else
	              efound = 0;
                end


               end

%       If shared, then if elem_curr > elem_old
             if(efound)
                if(elem>elemB)
	               copyn = 1;
                end
             end

             if(edge==1)
	             nA1 = 2;
	             nAinc = 1;
	             nA2 = 2*minc;
	             nA3 = minc+1;
	         elseif(edge==2)
	             nA1 = 4*minc+1;
	             nAinc = 2*minc;
	             nA2 = (2*minc-1)*(2*minc+1)+2;
	             nA3 = (2*minc+1)*minc+minc+1;
	         else
	             nA1 = (2*minc-1)*(2*minc+1)+1;
	             nAinc = -2*minc-1;
	             nA2 = 2*minc+2;
	             nA3 = (2*minc+1)*minc+1;
             end

             if(copyn)
%         Copy node IDs into element table
                  node = NodesOnElement(1,elemB);
	            neleB = cis(2,elemB);
	            k = 1;
%         Find which corner nodeA is on elemB
                 while((nodeA~=node)&&(k<neleB))
	               k = k + 1;
	               node = NodesOnElement(k,elemB);
                 end

%         Set extraction indices for copying nodes
              if(neleB==6)
                if(k<3)
                   if(nodeB==NodesOnElement(k+1,elemB))
	                  %Forward edge
	                  nB1 = cis(1,elemB) + neleB + (k-1)*(2*minc-2);
	                  nBinc = 1;
	                  nB2 = nB1 + 2*minc - 3;
	               else
	                  %Backward edge
                      if(k==1)
	                     k = 4;
                      end
	                  nB2 = cis(1,elemB) + neleB + (k-2)*(2*minc-2);
	                  nBinc = -1;
	                  nB1 = nB2 + 2*minc - 3;
                   end
	            else %k = neleB
                   if(nodeB==NodesOnElement(1,elemB))
	                  %Forward edge
	                  nB1 = cis(1,elemB) + neleB + (k-1)*(2*minc-2);
	                  nBinc = 1;
	                  nB2 = nB1 + 2*minc - 3;
	               else
	                  %Backward edge
	                  nB2 = cis(1,elemB) + neleB + (k-2)*(2*minc-2);
	                  nBinc = -1;
	                  nB1 = nB2 + 2*minc - 3;
                   end
                end
	          elseif(neleB==9)
                if(k<4)
                   if(nodeB==NodesOnElement(k+1,elemB))
	                  %Forward edge
	                  nB1 = cis(1,elemB) + neleB + (k-1)*(2*minc-2);
	                  nBinc = 1;
	                  nB2 = nB1 + 2*minc - 3;
	               else
	                  %Backward edge
                      if(k==1)
	                     k = 5;
                      end
	                  nB2 = cis(1,elemB) + neleB + (k-2)*(2*minc-2);
	                  nBinc = -1;
	                  nB1 = nB2 + 2*minc - 3;
                   end
	            else %k = neleB
                   if(nodeB==NodesOnElement(1,elemB))
	                  %Forward edge
	                  nB1 = cis(1,elemB) + neleB + (k-1)*(2*minc-2);
	                  nBinc = 1;
	                  nB2 = nB1 + 2*minc - 3;
	               else
	                  %Backward edge
	                  nB2 = cis(1,elemB) + neleB + (k-2)*(2*minc-2);
	                  nBinc = -1;
	                  nB1 = nB2 + 2*minc - 3;
                   end
                end
              end

%         Load nodes from elemB edge
      %            if(nB1.eq.nB2) then
	%               node = ecn(nB1)
	%               el2mesh(nA1) = node
	%               enode = enode + 1
	%			   ecn(enode) = node	
	%            else
	               k = nA1;
                     for j = nB1:nBinc:nB2
                      if(k==nA3)
	                     k = k + nAinc;
                      end
	                  node = ecn(j);
	                  el2mesh(k) = node;
	                  enode = enode + 1;
				      ecn(enode) = node; %#ok<AGROW>
	                  k = k + nAinc;
                     end
	%            endif

%       Else, add node IDs to node table, store coordinates
               else

	%            if(nA1.eq.nA2) then
	%               nodeid = nodeid + 1
	%               el2mesh(nA1) = nodeid
	%               enode = enode + 1
	%			   ecn(enode) = nodeid	
	 %           else
                     for j = nA1:nAinc:nA2
                      if(j~=nA3)
	                     nodeid = nodeid + 1;
	                     el2mesh(j) = nodeid;
	                     enode = enode + 1;
				         ecn(enode) = nodeid; %#ok<AGROW>
                         for i = 1:ndm
				            xs(i,nodeid) = xc(i,j); %#ok<AGROW>
                         end
                      end
                     end
	%            endif

             end

%     End loop over edges
            end

%     Store interior node IDs and coordinates
          k = 2*minc+1;
          for j = 2:2*minc
	         k = k + 1;
             for i = 2:2*minc+1-j
		        nodeid = nodeid + 1;
	            k = k + 1;
	            el2mesh(k) = nodeid;
	            enode = enode + 1;
				ecn(enode) = nodeid; %#ok<AGROW>
                for l = 1:ndm
				   xs(l,nodeid) = xc(l,k); %#ok<AGROW>
                end
             end
             for i = 2*minc+2-j:2*minc
	            k = k + 1;
             end
	         k = k + 1;
          end

%     Transfer data from element array to submesh array
          l = 0;
          for k = 1:2*minc
             for j = 1:2*minc-1-(k-1)*2
	            l = l + 1;
	            cell = cell + 1;
                for i = 1:nele
	               node = ixc(i,l);
	               node = el2mesh(node);
	               ixs(i,cell) = node; %#ok<AGROW>
                end
             end
             for j = 2*minc-(k-1)*2:2*minc
	            l = l + 1;
             end
          end
          
          

%       %Else Q9 Element
% 	   else 
% 
% 	      el2mesh(1) = NodesOnElement(1,elem)
% 	      el2mesh(2*minc+1) = NodesOnElement(2,elem)
% 	      el2mesh((2*minc+1)*(2*minc+1)) = NodesOnElement(3,elem)
% 	      el2mesh((2*minc+1)*2*minc+1) = NodesOnElement(4,elem)
% 	      el2mesh(minc+1) = NodesOnElement(5,elem)
% 	      el2mesh((2*minc+1)*(minc+1)) = NodesOnElement(6,elem)
% 	      el2mesh((2*minc+1)*(2*minc)+minc+1) = NodesOnElement(7,elem)
% 	      el2mesh((2*minc+1)*minc+1) = NodesOnElement(8,elem)
% 	      el2mesh((2*minc+1)*minc+minc+1) = NodesOnElement(9,elem)
% 
% 	if(debug) then
% 	dprt = .true.
% 	if(dprt) then
% 	   write(24,*) 'el2mesh'
% 	   do j = 1,ncel
%          write(24,*) el2mesh(j)
% 	   enddo
% 	   write(24,'(A)') ''
% 	endif
% 	endif
% 
% 	      do i=1,4
% 	      nloop(i,1) = i
% 	      nloop(i,2) = i+1
% 	         enode = enode + 1
% 	         ecn(enode) = NodesOnElement(i,elem)
% 	      enddo
% 	      do i = 5,9
% 	         enode = enode + 1
% 	         ecn(enode) = NodesOnElement(i,elem)
% 	      enddo
% 	      nloop(4,2) = 1
% 
% %     Loop over element edges
%             do edge = 1,4
% 
% 	         nodeA = NodesOnElement(nloop(edge,1),elem)
% 		     nodeB = NodesOnElement(nloop(edge,2),elem)
% 	         numA = epnum(nodeA)
% 	         numB = epnum(nodeB)
% 
% 	if(debug) then
% 	dprt = .true.
% 	if(dprt) then
% 	   write(24,*) 'nodeA',nodeA
% 	   write(24,*) 'nodeB',nodeB
% 	   write(24,*) 'numA',numA
% 	   write(24,*) 'numB',numB
% 	endif
% 	endif
% 
% 	         efound = .false.
% 	         copyn = .false.
% 	         k = 0
% 
% %       Search for shared edge in level 1 stars of nodes A & B
%                do 40 while((.not.efound).and.(k.lt.numA))
% 
% 	            k = k + 1
% 	            elemA = epatch(k,nodeA)
% 	if(debug) then
% 	dprt = .true.
% 	if(dprt) then
% 	   write(24,*) 'elemA',elemA
% 	endif
% 	endif
% 	            call binsearch(elemA,epatch(:,nodeB),1,numB,elemB
%      *			               ,nume,efound)
% 	if(debug) then
% 	dprt = .true.
% 	if(dprt) then
% 	   write(24,*) 'elemB',elemB
% 	endif
% 	endif
% 	            if(efound.and.(elemB.ne.elem)) then
% 	              efound = .true.
% 	            else
% 	              efound = .false.
% 	            endif
% 
% 
% 40	         continue
% 
%       if(debug) then
% 	dprt = .true.
% 	if(dprt) then
%          write(24,*) 'efound ',efound
% 	endif
% 	endif
% 
% %       If shared, then if elem_curr > elem_old
% 	         if(efound) then
% 	            if(elem.gt.elemB) then
% 	               copyn = .true.
% 	            endif
% 	         endif
% 
%       if(debug) then
% 	dprt = .true.
% 	if(dprt) then
%          write(24,*) 'copyn ',copyn
% 	endif
% 	endif
% 
% 	         if(edge.eq.1) then
% 	             nA1 = 2
% 	             nAinc = 1
% 	             nA2 = 2*minc
% 	             nA3 = minc+1
% 	         elseif(edge.eq.2) then
% 	             nA1 = 2*(2*minc+1)
% 	             nAinc = 2*minc+1
% 	             nA2 = 2*minc*(2*minc+1)
% 	             nA3 = (2*minc+1)*(minc+1)
% 	         elseif(edge.eq.3) then
% 	             nA1 = (2*minc+1)*(2*minc+1)-1
% 	             nAinc = -1
% 	             nA2 = 2*minc*(2*minc+1)+2
% 	             nA3 = (2*minc+1)*(2*minc)+minc+1
% 	         else
% 	             nA1 = (2*minc-1)*(2*minc+1)+1
% 	             nAinc = -2*minc-1
% 	             nA2 = 2*minc+2
% 	             nA3 = (2*minc+1)*minc+1
% 	         endif
% 
% 	if(debug) then
% 	dprt = .true.
% 	if(dprt) then
% 	   write(24,*) 'nA1',nA1
% 	   write(24,*) 'nAinc',nAinc
% 	   write(24,*) 'nA2',nA2
% 	   write(24,*) 'nA3',nA3
% 	   write(24,'(A)') ''
% 	endif
% 	endif
% 
% 	         if(copyn) then
% %         Copy node IDs into element table
%                   node = NodesOnElement(1,elemB)
% 	            neleB = cis(2,elemB)
% 	            k = 1
% %         Find which corner nodeA is on elemB
%                   do while((nodeA.ne.node).and.(k.lt.neleB))
% 	               k = k + 1
% 	               node = NodesOnElement(k,elemB)
% 	            enddo
% 
% 	if(debug) then
% 	dprt = .false.
% 	if(dprt) then
% 	   write(24,*) 'k',k
% 	   write(24,'(A)') ''
% 	endif
% 	endif
% 
% %         Set extraction indices for copying nodes
% 	            if(neleB.eq.6) then
% 
% 	            if(k.lt.neleB) then
% 	               if(nodeB.eq.NodesOnElement(k+1,elemB)) then
% 	                  %Forward edge
% 	                  nB1 = cis(1,elemB) + neleB + (k-1)*(2*minc-2)
% 	                  nBinc = 1
% 	                  nB2 = nB1 + 2*minc - 3
% 	               else
% 	                  %Backward edge
% 	                  if(k.eq.1) then
% 	                     k = 4
% 	                  endif
% 	                  nB2 = cis(1,elemB) + neleB + (k-2)*(2*minc-2)
% 	                  nBinc = -1
% 	                  nB1 = nB2 + 2*minc - 3
% 	               endif
% 	            else %k = neleB
% 	               if(nodeB.eq.NodesOnElement(1,elemB)) then
% 	                  %Forward edge
% 	                  nB1 = cis(1,elemB) + neleB + (k-1)*(2*minc-2)
% 	                  nBinc = 1
% 	                  nB2 = nB1 + 2*minc - 3
% 	               else
% 	                  %Backward edge
% 	                  nB2 = cis(1,elemB) + neleB + (k-2)*(2*minc-2)
% 	                  nBinc = -1
% 	                  nB1 = nB2 + 2*minc - 3
% 	               endif
% 	            endif
% 	            elseif(neleB.eq.9) then
% 	            if(k.lt.4) then
% 	               if(nodeB.eq.NodesOnElement(k+1,elemB)) then
% 	                  %Forward edge
% 	                  nB1 = cis(1,elemB) + neleB + (k-1)*(2*minc-2)
% 	                  nBinc = 1
% 	                  nB2 = nB1 + 2*minc - 3
% 	               else
% 	                  %Backward edge
% 	                  if(k.eq.1) then
% 	                     k = 5
% 	                  endif
% 	                  nB2 = cis(1,elemB) + neleB + (k-2)*(2*minc-2)
% 	                  nBinc = -1
% 	                  nB1 = nB2 + 2*minc - 3
% 	               endif
% 	            else %k = neleB
% 	               if(nodeB.eq.NodesOnElement(1,elemB)) then
% 	                  %Forward edge
% 	                  nB1 = cis(1,elemB) + neleB + (k-1)*(2*minc-2)
% 	                  nBinc = 1
% 	                  nB2 = nB1 + 2*minc - 3
% 	               else
% 	                  %Backward edge
% 	                  nB2 = cis(1,elemB) + neleB + (k-2)*(2*minc-2)
% 	                  nBinc = -1
% 	                  nB1 = nB2 + 2*minc - 3
% 	               endif
% 	            endif
% 	            endif
% 
% 	if(debug) then
% 	dprt = .true.
% 	if(dprt) then
% 	   write(24,*) 'nB1',nB1
% 	   write(24,*) 'nBinc',nBinc
% 	   write(24,*) 'nB2',nB2
% 	   write(24,'(A)') ''
% 	endif
% 	endif
% %         Load nodes from elemB edge
%       %            if(nB1.eq.nB2) then
% 	%               node = ecn(nB1)
% 	%               el2mesh(nA1) = node
% 	%               enode = enode + 1
% 	%			   ecn(enode) = node	
% 	%            else
% 	               k = nA1
%                      do j = nB1,nB2,nBinc
% 	                  if(k.eq.nA3) then
% 	                     k = k + nAinc
% 	                  endif
% 	                  node = ecn(j)
% 	                  el2mesh(k) = node
% 	                  enode = enode + 1
% 				      ecn(enode) = node
% 	                  k = k + nAinc
% 	               enddo
% 	%            endif
% 
% %       Else, add node IDs to node table, store coordinates
%                else
% 
% 	%            if(nA1.eq.nA2) then
% 	%               nodeid = nodeid + 1
% 	%               el2mesh(nA1) = nodeid
% 	%               enode = enode + 1
% 	%			   ecn(enode) = nodeid	
% 	 %           else
%                      do j = nA1,nA2,nAinc
% 	                  if(j.ne.nA3) then
% 	                     nodeid = nodeid + 1
% 	                     el2mesh(j) = nodeid
% 	                     enode = enode + 1
% 				         ecn(enode) = nodeid
% 				         do i = 1,ndm
% 				            xs(i,nodeid) = xc(i,j)
% 				         enddo
% 	                  endif
% 	               enddo
% 	%            endif
% 
% 	         endif
% 
% 	if(debug) then
% 	dprt = .true.
% 	if(dprt) then
% 	   write(24,*) 'el2mesh'
% 	   do j = 1,ncel
%          write(24,*) el2mesh(j)
% 	   enddo
% 	   write(24,'(A)') ''
% 	endif
% 	endif
% 
% %     End loop over edges
%             enddo
% 
% %     Store interior node IDs and coordinates
%             k = 2*minc+1
% 	      do j = 2,2*minc
% 	         k = k + 1
% 	         do i = 2,2*minc
% 	            k = k + 1
% 	            if((j.ne.minc+1).or.(i.ne.minc+1)) then
% 		           nodeid = nodeid + 1
% 	               el2mesh(k) = nodeid
% 	               enode = enode + 1
% 				   ecn(enode) = nodeid
% 				   do l = 1,ndm
% 				      xs(l,nodeid) = xc(l,k)
% 				   enddo
% 	            endif
% 	         enddo
% 	         k = k + 1
% 	      enddo
% 
% 	if(debug) then
% 	dprt = .true.
% 	if(dprt) then
% 	   write(24,*) 'el2mesh'
% 	   do j = 1,ncel
%          write(24,*) el2mesh(j)
% 	   enddo
% 	   write(24,'(A)') ''
% 	endif
% 	endif
% 
% %     Transfer data from element array to submesh array
% 	      do k = 1,cel
% 	         cell = cell + 1
% 	         do j = 1,nele
% 	            node = ixc(j,k)
% 	            node = el2mesh(node)
% 	            ixs(j,cell) = node
% 	         enddo
% 	      enddo
% 
%          endif
% 
% 
% % End loop over elements
%       enddo
      
%%
      %Else Q9 Element
	   else 

	      el2mesh(1) = NodesOnElement(1,elem);
	      el2mesh(2*minc+1) = NodesOnElement(2,elem);
	      el2mesh((2*minc+1)*(2*minc+1)) = NodesOnElement(3,elem);
	      el2mesh((2*minc+1)*2*minc+1) = NodesOnElement(4,elem);
	      el2mesh(minc+1) = NodesOnElement(5,elem);
	      el2mesh((2*minc+1)*(minc+1)) = NodesOnElement(6,elem);
	      el2mesh((2*minc+1)*(2*minc)+minc+1) = NodesOnElement(7,elem);
	      el2mesh((2*minc+1)*minc+1) = NodesOnElement(8,elem);
	      el2mesh((2*minc+1)*minc+minc+1) = NodesOnElement(9,elem);

          for i=1:4
	         nloop(i,1) = i;
	         nloop(i,2) = i+1;
	         enode = enode + 1;
	         ecn(enode) = NodesOnElement(i,elem); %#ok<AGROW>
          end
          for i = 5:9
	         enode = enode + 1;
	         ecn(enode) = NodesOnElement(i,elem); %#ok<AGROW>
          end
	      nloop(4,2) = 1;

%     Loop over element edges
            for edge = 1:4

	         nodeA = NodesOnElement(nloop(edge,1),elem);
		     nodeB = NodesOnElement(nloop(edge,2),elem);
	         numA = epnum(nodeA);
	         numB = epnum(nodeB);

	         efound = 0;
	         copyn = 0;
	         k = 0;

%       Search for shared edge in level 1 stars of nodes A & B
               while((efound==0)&&(k<numA))

	            k = k + 1;
	            elemA = epatch(k,nodeA);
	            [elemB,nume,efound] = binsearch(elemA,epatch(:,nodeB),1,numB);
                if(efound&&(elemB~=elem))
	              efound = 1;
	            else
	              efound = 0;
                end


               end

%       If shared, then if elem_curr > elem_old
             if(efound)
                if(elem>elemB)
	               copyn = 1;
                end
             end

             if(edge==1)
	             nA1 = 2;
	             nAinc = 1;
	             nA2 = 2*minc;
	             nA3 = minc+1;
	         elseif(edge==2)
	             nA1 = 2*(2*minc+1);
	             nAinc = 2*minc+1;
	             nA2 = 2*minc*(2*minc+1);
	             nA3 = (2*minc+1)*(minc+1);
	         elseif(edge==3)
	             nA1 = (2*minc+1)*(2*minc+1)-1;
	             nAinc = -1;
	             nA2 = 2*minc*(2*minc+1)+2;
	             nA3 = (2*minc+1)*(2*minc)+minc+1;
	         else
	             nA1 = (2*minc-1)*(2*minc+1)+1;
	             nAinc = -2*minc-1;
	             nA2 = 2*minc+2;
	             nA3 = (2*minc+1)*minc+1;
             end

             if(copyn)
%         Copy node IDs into element table
                  node = NodesOnElement(1,elemB);
	            neleB = cis(2,elemB);
	            k = 1;
%         Find which corner nodeA is on elemB
                  while((nodeA~=node)&&(k<neleB))
	               k = k + 1;
	               node = NodesOnElement(k,elemB);
                  end

%         Set extraction indices for copying nodes
              if(neleB==6)

                if(k<neleB)
                   if(nodeB==NodesOnElement(k+1,elemB))
	                  %Forward edge
	                  nB1 = cis(1,elemB) + neleB + (k-1)*(2*minc-2);
	                  nBinc = 1;
	                  nB2 = nB1 + 2*minc - 3;
	               else
	                  %Backward edge
                      if(k==1)
	                     k = 4;
                      end
	                  nB2 = cis(1,elemB) + neleB + (k-2)*(2*minc-2);
	                  nBinc = -1;
	                  nB1 = nB2 + 2*minc - 3;
                   end
	            else %k = neleB
                   if(nodeB==NodesOnElement(1,elemB))
	                  %Forward edge
	                  nB1 = cis(1,elemB) + neleB + (k-1)*(2*minc-2);
	                  nBinc = 1;
	                  nB2 = nB1 + 2*minc - 3;
	               else
	                  %Backward edge
	                  nB2 = cis(1,elemB) + neleB + (k-2)*(2*minc-2);
	                  nBinc = -1;
	                  nB1 = nB2 + 2*minc - 3;
                   end
                end
	          elseif(neleB==9)
                if(k<4)
                   if(nodeB==NodesOnElement(k+1,elemB))
	                  %Forward edge
	                  nB1 = cis(1,elemB) + neleB + (k-1)*(2*minc-2);
	                  nBinc = 1;
	                  nB2 = nB1 + 2*minc - 3;
	               else
	                  %Backward edge
                      if(k==1)
	                     k = 5;
                      end
	                  nB2 = cis(1,elemB) + neleB + (k-2)*(2*minc-2);
	                  nBinc = -1;
	                  nB1 = nB2 + 2*minc - 3;
                   end
	            else %k = neleB
                   if(nodeB==NodesOnElement(1,elemB))
	                  %Forward edge
	                  nB1 = cis(1,elemB) + neleB + (k-1)*(2*minc-2);
	                  nBinc = 1;
	                  nB2 = nB1 + 2*minc - 3;
	               else
	                  %Backward edge
	                  nB2 = cis(1,elemB) + neleB + (k-2)*(2*minc-2);
	                  nBinc = -1;
	                  nB1 = nB2 + 2*minc - 3;
                   end
                end
              end

%         Load nodes from elemB edge
      %            if(nB1.eq.nB2) then
	%               node = ecn(nB1)
	%               el2mesh(nA1) = node
	%               enode = enode + 1
	%			   ecn(enode) = node	
	%            else
	               k = nA1;
                     for j = nB1:nBinc:nB2
                      if(k==nA3)
	                     k = k + nAinc;
                      end
	                  node = ecn(j);
	                  el2mesh(k) = node;
	                  enode = enode + 1;
				      ecn(enode) = node; %#ok<AGROW>
	                  k = k + nAinc;
                     end
	%            endif

%       Else, add node IDs to node table, store coordinates
               else

	%            if(nA1.eq.nA2) then
	%               nodeid = nodeid + 1
	%               el2mesh(nA1) = nodeid
	%               enode = enode + 1
	%			   ecn(enode) = nodeid	
	 %           else
                     for j = nA1:nAinc:nA2
                      if(j~=nA3)
	                     nodeid = nodeid + 1;
	                     el2mesh(j) = nodeid;
	                     enode = enode + 1;
				         ecn(enode) = nodeid; %#ok<AGROW>
                         for i = 1:ndm
				            xs(i,nodeid) = xc(i,j); %#ok<AGROW>
                         end
                      end
                     end
	%            endif

             end

%     End loop over edges
            end

%     Store interior node IDs and coordinates
          k = 2*minc+1;
          for j = 2:2*minc
	         k = k + 1;
             for i = 2:2*minc
	            k = k + 1;
                if((j~=minc+1)||(i~=minc+1))
		           nodeid = nodeid + 1;
	               el2mesh(k) = nodeid;
	               enode = enode + 1;
				   ecn(enode) = nodeid; %#ok<AGROW>
                   for l = 1:ndm
				      xs(l,nodeid) = xc(l,k); %#ok<AGROW>
                   end
                end
             end
	         k = k + 1;
          end

%     Transfer data from element array to submesh array
          for k = 1:cel
	         cell = cell + 1;
             for j = 1:nele
	            node = ixc(j,k);
	            node = el2mesh(node);
	            ixs(j,cell) = node; %#ok<AGROW>
             end
          end

       end


% End loop over elements
      end
      

	numnps = nodeid;
	numels = numel*cel;

