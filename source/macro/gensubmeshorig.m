%**********************************************************************
% 
%      subroutine sub_mesh(d,ix,x,u,iprob,strong)
%
%...  Written by Tim Truster (Fall 2009)
%     No modifications when copied to NLFEA ver2
%...  Program to compute fine scales by implicit method of solving local
%     problems on stars
%
%**********************************************************************

% 	implicit none
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

% %     Input Variables
%       integer ix(nen1,numel),iprob,strong
%       real*8 d(250,*),x(ndm,numnp),u(ndf,numnp,3)

%     Output Variables
%      real*8 
% 	include 'smparam.h'
% 	include 'stardata.h'
% 	include 'submeshdata.h'

%     Local Variables
% 	integer epatch(maxel,numnp),epnum(numnp)
% 	integer ecn(numel*celn),cis(2,numel+1)
% 
% 	integer gammab(numnp),gammah(numnp-ndm)
% 
%       integer cstar,c1sen,c1snn,c1snbn,c1sebn
% 	integer cse(numel),csn(numnp),c1se(numel),c1sn(numnp)
% 	integer cseb(numel),csnb(numnp),c1seb(numel),c1snb(numnp)
% 
% 	integer i,j,k,q,node,elem,nele
% 	integer ixe(nen)
% 
% 	logical strwek,globalm
% 
% 	real*8 xuno(ndm)
% 	integer nume,numn,locid(numel*celn)
% 	integer ind,len
% 	real    tt,etime,tary(2)
% 	
%       logical dprt,setsubmeshvar,palloc
% 
% 	include 'comblk.h'
% 	include 'pointer.h'
% 	      
%       save

	epatch = zeros(maxel,numnp);
    epnum = zeros(numnp,1);
	ecn = zeros(numel*celn,1);
    cis = zeros(2,numel+1);
	cse = zeros(numel,1);
    csn = zeros(numnp,1);
    c1se = zeros(numel,1);
    c1sn = zeros(numnp,1);
	cseb = zeros(numel,1);
    csnb = zeros(numnp,1);
    c1seb = zeros(numel,1);
    c1snb = zeros(numnp,1);
	ixe = zeros(nen,1);
	xuno = zeros(ndm,1);

% Start time
%       tt = etime(tary)

%	if(debug) then

% 	   write(24,'(A,f17.7)') 'Submesh Start = ',tt
% 	   write(24,'(A)') ''

%	endif

% 	if(strong.eq.1) then
% 	strwek = .true.
% 	else
% 	strwek = .false.
% 	endif

% 	globalm = .false. %True for global solve, false for patchwise
	globalm = 1; %True for global solve, false for patchwise
fprintf('globalm = %i, strong = %i, nstar = %i\n',globalm,strong,nstar)

% Declare submesh arrays
% 	setsubmeshvar = palloc(151,'USER1',ndm*numel*celn,2)  %xs
% 	setsubmeshvar = palloc(152,'USER2',nen*numel*cel,1)   %ixs
% 	setsubmeshvar = palloc(153,'USER3',ndfs*numel*celn,2) %us
% 	call pzeroi(ecn,numel*celn)
% 	call pzero(hr(np(151)),ndm*numel*celn)
% 	call pzeroi(mr(np(152)),nen*numel*cel)
% 	call pzero(hr(np(153)),ndfs*numel*celn)
	xs = zeros(ndm,numel*celn);  %xs
	ixs = zeros(nen,numel*cel);   %ixs
	us = zeros(ndfs,numel*celn); %us

%-----------------------------------------------------
% Determine level 1 stars around coarse nodes
%-----------------------------------------------------
% 	do j=1,numnp
% 	epnum(j) = 0
% 	do i=1,maxel
% 	epatch(i,j)=0	%Zero patch connectivity
% 	enddo
% 	enddo

% 	if(nen.le.4) then
%
% % Loop over elements
% 	do elem=1,numel
% 
% 	   do k=1,nen % Loop over local Nodes
% 
% 	      node = ix(k,elem)
% 
% %   Add element to star list, increment number of elem in star
% 	       if (node.gt.0) then
% 
% 	         q = epnum(node)+1
% 
% 	         epatch(q,node) = elem		%epatch(nel,numnp)
% 
% 	         epnum(node) = q			%epnum(numnp)
% 
% 	       endif
% 
% 	   enddo %k
%       enddo %j
% 
% 	elseif(nen.eq.6) then
% 
% % Loop over elements
% 	do elem=1,numel
% 
% 	   do k=1,3 % Loop over local Nodes
% 
% 	      node = ix(k,elem)
% 
% %   Add element to star list, increment number of elem in star
% 	       if (node.gt.0) then
% 
% 	         q = epnum(node)+1
% 
% 	         epatch(q,node) = elem		%epatch(nel,numnp)
% 
% 	         epnum(node) = q			%epnum(numnp)
% 
% 	       endif
% 
% 	   enddo %k
%       enddo %j
% 
% 	else
% 
% 	% Loop over elements
% 	do elem=1,numel
% 
% 	   do k=1,3 % Loop over local Nodes
% 
% 	      node = ix(k,elem)
% 
% %   Add element to star list, increment number of elem in star
% 	       if (node.gt.0) then
% 
% 	         q = epnum(node)+1
% 
% 	         epatch(q,node) = elem		%epatch(nel,numnp)
% 
% 	         epnum(node) = q			%epnum(numnp)
% 
% 	       endif
% 
% 	   enddo %k
% 
% 	   if(ix(9,elem).gt.0) then
% 
% 	       node = ix(4,elem)
% 
% %   Add element to star list, increment number of elem in star
% 
% 	         q = epnum(node)+1
% 
% 	         epatch(q,node) = elem		%epatch(nel,numnp)
% 
% 	         epnum(node) = q			%epnum(numnp)
% 
% 	   endif
% 
%       enddo %j
% 
% 
% 	endif
%     
% 	if(debug) then
% 	dprt = .true.
% 	if(dprt) then
% 
% 	   write(24,'(A,I4)') 'maxel',maxel
% 	   write(24,'(A,I4)') 'numnp',numnp
% 	   write(24,'(A,I4)') 'numel',numel
% 	   write(24,'(A)') ''  
% 	   WRITE(24,'(A)') 'Level 1 Stars'
% 	   write(24,'(A)') ''
% 	   write(24,'(A)') 'Node  numel   Elements'
% 	   do i=1,numnp	
% 		
% 		   WRITE(24,1001) i,epnum(i),(epatch(j,i),j=1,maxel)    
% 
% 	   enddo
% 	   write(24,'(A)') ''
%  1001 format (I4,2x,I4,2x,10(I5))
% 
% 	endif
% 	endif

%%
%-----------------------------------------------------
% Determine level 1 stars around coarse nodes
%-----------------------------------------------------

	if nen <= 4

% Loop over elements

	for elem=1:numel

	  for k=1:nen % Loop over local Nodes

	      node = NodesOnElement(k,elem);

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

	      node = NodesOnElement(k,elem);

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

	      node = NodesOnElement(k,elem);

%   Add element to star list, increment number of elem in star
          if (node>0)

	         q = epnum(node)+1;

	         epatch(q,node) = elem;		%epatch(nel,numnp)

	         epnum(node) = q;			%epnum(numnp)

          end

	   end %k

	   if(NodesOnElement(9,elem)>0)

	       node = NodesOnElement(4,elem);

%   Add element to star list, increment number of elem in star

	         q = epnum(node)+1;

	         epatch(q,node) = elem;		%epatch(nel,numnp)

	         epnum(node) = q;			%epnum(numnp)

	   end

	end %j


	end

% %-----------------------------------------------------
% % Determine nodes on domain and traction boundaries
% %-----------------------------------------------------
% 
% 	ngb = 0
% 	ngh = 0
% 
% 	do node = 1,numnp
% 	
% 	   do i = 1,ndm
% 	      xuno(i) = x(i,node)
% 	   enddo
% 
% 	   if((iprob.eq.1).or.(iprob.eq.2)) then
% 	      if((abs(xuno(1) - 0.d0).lt.10e-8).or.
%      *         (abs(xuno(1) - 10.d0).lt.10e-8)) then
% 		     ngh = ngh + 1
% 			 gammah(ngh) = node
% 			 ngb = ngb + 1
% 			 gammab(ngb) = node
% 	      elseif((abs(xuno(2) - 1.d0).lt.10e-8).or.
%      *         (abs(xuno(2) + 1.d0).lt.10e-8)) then
% 			 ngb = ngb + 1
% 			 gammab(ngb) = node
% 	      endif 
% 	   endif
% 
% 	enddo
% 
% 	if(debug) then
% 	dprt = .false.
% 	if(dprt) then
% 
% 	   write(24,'(A,I4)') 'ngb',ngb
% 	   write(24,'(A)') ''  
% 	   do i=1,ngb	
% 		
% 		   WRITE(24,'(I4)') gammab(i)    
% 
% 	   enddo
% 	   write(24,'(A)') ''
% 	   write(24,'(A,I4)') 'ngh',ngh
% 	   write(24,'(A)') ''  
% 	   do i=1,ngh	
% 		
% 		   WRITE(24,'(I4)') gammah(i)    
% 
% 	   enddo
% 	   write(24,'(A)') ''
% 
% 	endif
% 	endif
%-----------------------------------------------------
% Generate submesh across entire domain
%-----------------------------------------------------

% % Start generation time
%       tt = etime(tary)
% 
% %	if(debug) then
% 
% 	   write(24,'(A,f17.7)') 'Geom. Gen. Start = ',tt
% 	   write(24,'(A)') ''
% 
% %	endif
% 
% 	call gen_submesh(epnum,epatch,hr(np(151)),mr(np(152)),
% 	&                 x,ix,cis,ecn)
% 
% % End Generation time
%       tt = etime(tary)
% 
% %	if(debug) then
% 
% 	   write(24,'(A,f17.7)') 'Geom. Gen. End = ',tt
% 	   write(24,'(A)') ''
% 
% %	endif

% 	[xs,ixs,cis,ecn] = gen_submesh(epnum,epatch,x,ix);
	gen_submesh

