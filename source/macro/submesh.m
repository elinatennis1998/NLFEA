%**********************************************************************
% 
%      subroutine sub_mesh(d,NodesOnElement,x,u,iprob,strong)
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
%       integer NodesOnElement(nen1,numel),iprob,strong
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

sqr2 = sqrt(2);
	epatch = zeros(maxel,numnp);
    epnum = zeros(numnp,1);
	ecn = zeros(numel*celn,1);
    cis = zeros(2,numel+1);
	gammab = zeros(numnp,1);
	gammad = zeros(numnp,2);
    gammah = zeros(numnp-ndm,1);
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
	locid = zeros(numel*celn,1);

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
	globalm = 0;1; %True for global solve, false for patchwise
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

%%
%-----------------------------------------------------
% Determine nodes on domain and traction boundaries
%-----------------------------------------------------

	ngb = 0;
    ngd = zeros(2,1);
	ngh = 0;

    iel = MatTypeTable(2,1); %iel   = ie(nie-1,ma); same thing;
    nonlin = MatTypeTable(3,1);
    
    if(nonlin)
    
    if iel == 5 %NL_Elem5_2d
    
	for node = 1:numnp
	
       for i = 1:ndm
	      xuno(i) = x(i,node);
       end

	   if(iprob==1)
          if((abs(xuno(2) - 10.d0)<10e-8)&&((xuno(1))<5+10e-8))
		     ngh = ngh + 1;
			 gammah(ngh) = node;
			 ngb = ngb + 1;
			 gammab(ngb) = node;
          elseif((abs(abs(xuno(1)-5) - 5.d0)<10e-8)||(abs(abs(xuno(2)-5) - 5.d0)<10e-8))
			 ngb = ngb + 1;
			 gammab(ngb) = node;
          end
          if(abs(xuno(1) - 0.d0)<10e-8)
			 ngd(1) = ngd(1) + 1;
			 gammad(ngd(1),1) = node;
          end
          if(abs(xuno(2) - 10.d0)<10e-8)
			 ngd(1) = ngd(1) + 1;
			 gammad(ngd(1),1) = node;
          end
          if(abs(xuno(2) - 0.d0)<10e-8)
			 ngd(2) = ngd(2) + 1;
			 gammad(ngd(2),2) = node;
          end
       elseif(iprob==2)
          if((abs(xuno(1) - 1.d0)<10e-8))
		     ngh = ngh + 1;
			 gammah(ngh) = node;
			 ngb = ngb + 1;
			 gammab(ngb) = node;
          elseif(abs(abs(xuno(2)-.5) - .5)<10e-8)
			 ngb = ngb + 1;
			 gammab(ngb) = node;
          end
          if(abs(xuno(1) - 0.d0)<10e-8)
			 ngd(1) = ngd(1) + 1;
			 gammad(ngd(1),1) = node;
          end
          %%%% NEED TO ADD iprob==4
       elseif(iprob==5)
          if((abs(xuno(1) - 1.95)<10e-8))
			 ngb = ngb + 1;
			 gammab(ngb) = node;
          end
          if(abs(xuno(2) - 0.0)<10e-8)
		     ngh = ngh + 1;
			 gammah(ngh) = node;
			 ngb = ngb + 1;
			 gammab(ngb) = node;
          end
          if((abs(xuno(1) - 0.0)<10e-8))
			 ngb = ngb + 1;
			 gammab(ngb) = node;
          end
          if((abs(xuno(2) - 1.75)<10e-8)&&((xuno(1) - 1.95/3)<10e-8))
			 ngb = ngb + 1;
			 gammab(ngb) = node;
          end
          if((abs(xuno(2) - 1.75)<10e-8)&&((xuno(1) - 1.95/3)>-10e-8))
			 ngd(1) = ngd(1) + 1;
			 gammad(ngd(1),1) = node;
          end
       elseif(iprob==6)
          if(abs(xuno(1) - 0.0)<10e-8)
		     ngh = ngh + 1;
			 gammah(ngh) = node;
			 ngb = ngb + 1;
			 gammab(ngb) = node;
          end
          if((abs(xuno(2) - 1.0)<10e-8))
		     ngh = ngh + 1;
			 gammah(ngh) = node;
			 ngb = ngb + 1;
			 gammab(ngb) = node;
          end
          if(abs(xuno(1) - 1.0)<10e-8)
		     ngh = ngh + 1;
			 gammah(ngh) = node;
			 ngb = ngb + 1;
			 gammab(ngb) = node;
          end
          if((abs(xuno(2) - 0.0)<10e-8))
			 ngd(1) = ngd(1) + 1;
			 gammad(ngd(1),1) = node;
			 ngb = ngb + 1;
			 gammab(ngb) = node;
          end
       elseif(iprob==7)
           % add other boundaries too
          if(abs(xuno(1) - 48.0)<10e-8)
		     ngh = ngh + 1;
			 gammah(ngh) = node;
			 ngb = ngb + 1;
			 gammab(ngb) = node;
          end
%        elseif(iprob==7)
%           if((abs(xuno(2) + xuno(1) + sqr2)<10e-8) || ...
%                (abs(xuno(2) - xuno(1) + sqr2)<10e-8) || ...
%                (abs(xuno(2) + xuno(1) - sqr2)<10e-8) || ...
%                (abs(xuno(2) - xuno(1) - sqr2)<10e-8))
% 		     ngh = ngh + 1;
% 			 gammah(ngh) = node;
% 			 ngb = ngb + 1;
% 			 gammab(ngb) = node;
% 	      elseif((abs(xuno(2) - xuno(1))<10e-8) ...
%                 && (xuno(1)<=0.d0) || ...
%                (abs(xuno(2) + xuno(1))<10e-8) ...
%                 && (xuno(1)<=0.d0)) 
% 			 ngb = ngb + 1;
% 			 gammab(ngb) = node;
%           end 
	   end

	end
    
    end %iel == 5
    
    else %nonlin == 0
    
    if iel == 3 %L_Elem3_2d
    
	for node = 1:numnp
	
       for i = 1:ndm
	      xuno(i) = x(i,node);
       end

	   if((iprob==4)||(iprob==5))
          if((abs(xuno(1) - 0.d0)<10e-8)||(abs(xuno(1) - 10.d0)<10e-8))
		     ngh = ngh + 1;
			 gammah(ngh) = node;
			 ngb = ngb + 1;
			 gammab(ngb) = node;
	      elseif((abs(xuno(2) - 1.d0)<10e-8)||(abs(xuno(2) + 1.d0)<10e-8))
			 ngb = ngb + 1;
			 gammab(ngb) = node;
          end
       elseif(iprob==2)
          if((abs(xuno(1) - 0.d0)<10e-8)||(abs(xuno(1) - 4.d0)<10e-8))
		     ngh = ngh + 1;
			 gammah(ngh) = node;
			 ngb = ngb + 1;
			 gammab(ngb) = node;
	      elseif((abs(xuno(2) - 1.d0)<10e-8)||(abs(xuno(2) + 1.d0)<10e-8))
			 ngb = ngb + 1;
			 gammab(ngb) = node;
          end
       elseif(iprob==7)
          if((abs(xuno(2) + xuno(1) + sqr2)<10e-8) || ...
               (abs(xuno(2) - xuno(1) + sqr2)<10e-8) || ...
               (abs(xuno(2) + xuno(1) - sqr2)<10e-8) || ...
               (abs(xuno(2) - xuno(1) - sqr2)<10e-8))
		     ngh = ngh + 1;
			 gammah(ngh) = node;
			 ngb = ngb + 1;
			 gammab(ngb) = node;
	      elseif((abs(xuno(2) - xuno(1))<10e-8) ...
                && (xuno(1)<=0.d0) || ...
               (abs(xuno(2) + xuno(1))<10e-8) ...
                && (xuno(1)<=0.d0)) 
			 ngb = ngb + 1;
			 gammab(ngb) = node;
          end
       elseif(iprob==8)
          if((abs(xuno(1) - 48.0)<10e-8))
		     ngh = ngh + 1;
			 gammah(ngh) = node;
			 ngb = ngb + 1;
			 gammab(ngb) = node;
	      elseif((abs(xuno(2) - 44/48*xuno(1))<10e-8) || ...
               (abs(xuno(2) - 16/48*xuno(1) - 44)<10e-8) || ...
               (abs(xuno(1))<10e-8))
			 ngb = ngb + 1;
			 gammab(ngb) = node;
          end
	   end

	end
    
    end %iel == 3
    
    end %nonlin

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
% 	&                 x,NodesOnElement,cis,ecn)
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

% 	[xs,ixs,cis,ecn] = gen_submesh(epnum,epatch,x,NodesOnElement);
	gen_submesh

% 	if(globalm) then
    if(globalm)

% %-----------------------------------------------------
% % Solve for fine scale field in single global phase
% %-----------------------------------------------------
% 
% 	   call pzeroi(locid,numel*celn)
% 	   csnn = 0
% 	   numnpp = 0
% 	   csen = numel
% 	   do i = 1,csen
% 	      cse(i) = i
% 		  ind = cis(1,i)
% 	      len = cis(1,i+1) - ind
% 	      nele = cis(2,i)
% 	      ind = ind + nele
% %     Assign local IDs to cell nodes     	      
% 	      do j = nele+1,len
% 	         node = ecn(ind)
% 	         if(locid(node).eq.0) then
% 	            numnpp = numnpp + 1
% 	            if(node.le.numnp) then
% 	               csnn = csnn + 1
% 	            endif
% 	            locid(node) = numnpp
% 	         endif
% 	         ind = ind + 1
% 	      enddo
% 	   enddo
% 	   numelp = numels
% 
%       if(debug) then
% 	dprt = .true.
% 	if(dprt) then
% 	   write(24,*) 'numnpp',numnpp
% 	   write(24,*) 'numelp',numelp
% 	   write(24,*) 'csnn',csnn
% 	   write(24,*) 'locid'
% 	   do j = 1,numnps
%          write(24,*) locid(j)
% 	   enddo
% 	   write(24,'(A)') ''
% 	endif
% 	endif
% 
% % Start solution time
%       tt = etime(tary)
% 
% %	if(debug) then
% 
% 	   write(24,'(A,f17.7)') 'Global Solut. Start = ',tt
% 	   write(24,'(A)') ''
% 
% %	endif
% 
% 	   call starcontrol(NodesOnElement,mr(np(152)),d,hr(np(151)),
% 	&                    hr(np(153)),u,cis,ecn,cse,cseb,
%      &                    csnb,locid,gammab,gammah,strwek,
%      %                    globalm)
% 
% % End solution time
%       tt = etime(tary)
% 
% %	if(debug) then
% 
% 	   write(24,'(A,f17.7)') 'Global Solut. End = ',tt
% 	   write(24,'(A)') ''
% 
% %	endif

%%
%-----------------------------------------------------
% Solve for fine scale field in single global phase
%-----------------------------------------------------

	   locid = zeros(numel*celn,1); 

%   Assign local IDs to coarse nodes
       csnn = numnp;
       for i = 1:csnn
	      locid(i) = i;
       end
         
	   numnpp = csnn;
	   csen = numel;
       for i = 1:csen
	      cse(i) = i;
		  ind = cis(1,i);
	      len = cis(1,i+1) - ind;
	      nele = cis(2,i);
	      ind = ind + nele;
%     Assign local IDs to cell nodes
          for j = nele+1:len
	         node = ecn(ind);
             if(locid(node)==0)
	            numnpp = numnpp + 1;
                if(node<=numnp)
	               csnn = csnn + 1;
                end
	            locid(node) = numnpp;
             end
	         ind = ind + 1;
          end
       end
	   numelp = numels;

% 	   us = starcontrol(NodesOnElement,ixs,d,xs,us,u,cis,ecn,cse,cseb,csnb,locid,gammab,gammah,strong,globalm);
       starcontrol

	else

% %-----------------------------------------------------
% % Loop over coarse nodes
% % Solve for fine scale components
% %-----------------------------------------------------
% 	do 100 snode = 1,numnp
% 
% % Start patch snode time
%       tt = etime(tary)
% 
% %	if(debug) then
% 
% 	   write(24,'(A,I4,A,f17.7)') 'Start Patch ',snode,' Solut. = ',tt
% 	   write(24,'(A)') ''
% 
% %	endif
% 
% 	   if(epnum(snode).gt.0) then
% 
% %	if((snode.eq.1).or.(snode.eq.23)) then
% 
% %   Set components of partition of unity, uno
% 	   do i = 1,ndm
% 	      xuno(i) = x(i,snode)
% 	   enddo

%%
%-----------------------------------------------------
% Loop over coarse nodes
% Solve for fine scale components
%-----------------------------------------------------
    for snode = 1:numnp

       if(epnum(snode)>0)

%   Set components of partition of unity, uno
	   for i = 1:ndm
	      xuno(i) = x(i,snode);
	   end

% %  -----------------------------------------------------
% %   Determine level nstar stars around coarse nodes
% %  -----------------------------------------------------
% 	   
% 	   cstar = 0
% 	   csn(1) = snode
% 	   cse(1) = 0
% 	   csnn = 1
% 	   csen = 0
% 	   csnb(1) = snode
% 	   cseb(1) = 0
% 	   csnbn = 1
% 	   csebn = 0
% 
%       if(debug) then
% 	dprt = .false.
% 	if(dprt) then
% 	   write(24,*) 'snode',snode
% 	   write(24,*) 'csnn',csnn
% 	   write(24,*) 'csn'
% 	   do j = 1,csnn
%          write(24,*) csn(j)
% 	   enddo
% 	   write(24,*) 'csen',csen
% 	   write(24,*) 'cse'
% 	   do j = 1,csen
%          write(24,*) cse(j)
% 	   enddo
% 	   write(24,*) 'csnbn',csnbn
% 	   write(24,*) 'csnb'
% 	   do j = 1,csnbn
%          write(24,*) csnb(j)
% 	   enddo
% 	   write(24,*) 'csebn',csebn
% 	   write(24,*) 'cseb'
% 	   do j = 1,csebn
%          write(24,*) cseb(j)
% 	   enddo
% 	   write(24,'(A)') ''
% 	endif
% 	endif         
% 
% %   While nstar > cstar, add new level
% 	   do while(nstar.gt.cstar)
% 
% 	      cstar = cstar + 1
% 
% %     Shift elements and nodes in cstar to c_1star
% 
% 	      c1snn = csnn
% 	      c1sen = csen
% 	      c1snbn = csnbn
% 	      c1sebn = csebn
% 
% 	      do i = 1,c1snn
% 	         c1sn(i) = csn(i)
% 	      enddo
% 	      do i = 1,c1sen
% 	         c1se(i) = cse(i)
% 	      enddo
% 	      do i = 1,c1snbn
% 	         c1snb(i) = csnb(i)
% 	      enddo
% 	      do i = 1,c1sebn
% 	         c1seb(i) = cseb(i)
% 	      enddo
% 
%       if(debug) then
% 	dprt = .false.
% 	if(dprt) then
% 	   write(24,*) 'cstar',cstar
% 	   write(24,*) 'c1snn',c1snn
% 	   write(24,*) 'c1sn'
% 	   do j = 1,c1snn
%          write(24,*) c1sn(j)
% 	   enddo
% 	   write(24,*) 'c1sen',c1sen
% 	   write(24,*) 'c1se'
% 	   do j = 1,c1sen
%          write(24,*) c1se(j)
% 	   enddo
% 	   write(24,*) 'c1snbn',c1snbn
% 	   write(24,*) 'c1snb'
% 	   do j = 1,c1snbn
%          write(24,*) c1snb(j)
% 	   enddo
% 	   write(24,*) 'c1sebn',c1sebn
% 	   write(24,*) 'c1seb'
% 	   do j = 1,c1sebn
%          write(24,*) c1seb(j)
% 	   enddo
% 	   write(24,'(A)') ''
% 	endif
% 	endif 
% 
% %     Loop through nodes in c1snb, shift elements into cseb
% 
%             node = c1snb(1)
% 	      csebn = epnum(node)
% 	      do i = 1,csebn
% 	         cseb(i) = epatch(i,node)
% 	      enddo
% 
% 	if(debug) then
% 	dprt = .false.
% 	if(dprt) then
% 	   write(24,*) 'node',node
% 	   write(24,*) 'csebn',csebn
% 	   write(24,*) 'cseb'
% 	   do j = 1,csebn
%          write(24,*) cseb(j)
% 	   enddo
% 	   write(24,'(A)') ''
% 	endif
% 	endif
% 
% 	      do i = 2,c1snbn
% 
% %       Shuffle elements from nodal patch(i) into cseb
%                node = c1snb(i)
%                nume = epnum(node)
% 	         call shufflein(epatch(1,node),cseb,nume,csebn)
% 
% 	if(debug) then
% 	dprt = .false.
% 	if(dprt) then
% 	   write(24,*) 'node',node
% 	   write(24,*) 'csebn',csebn
% 	   write(24,*) 'cseb'
% 	   do j = 1,csebn
%          write(24,*) cseb(j)
% 	   enddo
% 	   write(24,'(A)') ''
% 	endif
% 	endif
% 
% %     End loop through nodes in c1snb
%             enddo
% 
% %     Combine boundary cseb into c1se and get cse
% 
%             call purge(cseb,c1se,cse,csebn,c1sen,csen)
% 
%       if(debug) then
% 	dprt = .false.
% 	if(dprt) then
% 	   write(24,*) 'csen',csen
% 	   write(24,*) 'cse'
% 	   do j = 1,csen
%          write(24,*) cse(j)
% 	   enddo
% 	   write(24,*) 'csebn',csebn
% 	   write(24,*) 'cseb'
% 	   do j = 1,csebn
%          write(24,*) cseb(j)
% 	   enddo
% 	   write(24,'(A)') ''
% 	endif
% 	endif
% 
% %     Loop through elements in cseb, shift nodes into csnb
% 
%             elem = cseb(1)
% 	      csnbn = cis(2,elem)
% 	      do i = 1,csnbn
% 	         csnb(i) = ix(i,elem)
% 	      enddo
% 	      call BubbleSortArray(csnb, csnbn, 1, 1)
% 
% 	if(debug) then
% 	dprt = .false.
% 	if(dprt) then
% 	   write(24,*) 'elem',elem
% 	   write(24,*) 'csnbn',csnbn
% 	   write(24,*) 'csnb'
% 	   do j = 1,csnbn
%          write(24,*) csnb(j)
% 	   enddo
% 	   write(24,'(A)') ''
% 	endif
% 	endif
% 
% 	      do i = 2,csebn
% 
% %       Shuffle nodes from element(i) into csnb
%                elem = cseb(i)
%                numn = cis(2,elem)
% 	         do j = 1,numn
% 	            ixe(j) = ix(j,elem)
% 	         enddo
% 	         call BubbleSortArray(ixe, numn, 1, 1)
% 	         call shufflein(ixe,csnb,numn,csnbn)
% 
% 	if(debug) then
% 	dprt = .false.
% 	if(dprt) then
% 	   write(24,*) 'elem',elem
% 	   write(24,*) 'csnbn',csnbn
% 	   write(24,*) 'csnb'
% 	   do j = 1,csnbn
%          write(24,*) csnb(j)
% 	   enddo
% 	   write(24,'(A)') ''
% 	endif
% 	endif
% 
% %     End loop through elements in cseb
%             enddo
% 
% %     Combine boundary csnb into c1sn and get csn
% 
%             call purge(csnb,c1sn,csn,csnbn,c1snn,csnn)
% 
%       if(debug) then
% 	dprt = .false.
% 	if(dprt) then
% 	   write(24,*) 'csnn',csnn
% 	   write(24,*) 'csn'
% 	   do j = 1,csnn
%          write(24,*) csn(j)
% 	   enddo
% 	   write(24,*) 'csnbn',csnbn
% 	   write(24,*) 'csnb'
% 	   do j = 1,csnbn
%          write(24,*) csnb(j)
% 	   enddo
% 	   write(24,'(A)') ''
% 	endif
% 	endif
% 
% %   End while loop on nstar 
% 	   enddo
% 
%       if(debug) then
% 	dprt = .true.
% 	if(dprt) then
% 	   write(24,*) 'snode',snode
% 	   write(24,*) 'csnn',csnn
% 	   write(24,*) 'csn'
% 	   do j = 1,csnn
%          write(24,*) csn(j)
% 	   enddo
% 	   write(24,*) 'csen',csen
% 	   write(24,*) 'cse'
% 	   do j = 1,csen
%          write(24,*) cse(j)
% 	   enddo
% 	   write(24,*) 'csnbn',csnbn
% 	   write(24,*) 'csnb'
% 	   do j = 1,csnbn
%          write(24,*) csnb(j)
% 	   enddo
% 	   write(24,*) 'csebn',csebn
% 	   write(24,*) 'cseb'
% 	   do j = 1,csebn
%          write(24,*) cseb(j)
% 	   enddo
% 	   write(24,'(A)') ''
% 	endif
% 	endif

%%
%  -----------------------------------------------------
%   Determine level nstar stars around coarse nodes
%  -----------------------------------------------------
	   
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
	         csnb(i) = NodesOnElement(i,elem);
          end
	      [csnb, csnbn] = BubbleSortArray(csnb, csnbn, 1, 1);

          for i = 2:csebn

%       Shuffle nodes from element(i) into csnb
               elem = cseb(i);
               numn = cis(2,elem);
             for j = 1:numn
	            ixe(j) = NodesOnElement(j,elem);
             end
	         [ixe, numn] = BubbleSortArray(ixe, numn, 1, 1);
	         [csnb,csnbn] = ShuffleIn(ixe,csnb,numn,csnbn);

%     End loop through elements in cseb
          end

%     Combine boundary csnb into c1sn and get csn

          [csn,csnn,csnb,csnbn] = Purge(csnb,c1sn,csnbn,c1snn);

%   End while loop on nstar 
       end

% %  -----------------------------------------------------
% %   Localize star problem
% %  -----------------------------------------------------
% 
%          call pzeroi(locid,numel*celn) 
% 
% %   Assign local IDs to coarse nodes
%          do i = 1,csnn
% 	      node = csn(i)
% 	      locid(node) = i
% 	   enddo
% 	   numnpp = csnn
% 
% %   Loop over elements in cstar
%          do k = 1,csen
% 	      elem = cse(k)
% 	      ind = cis(1,elem)
% 	      len = cis(1,elem+1) - ind
% 	      nele = cis(2,elem)
% 	      ind = ind + nele
% %     Assign local IDs to cell nodes     	      
% 	      do j = nele+1,len
% 	         node = ecn(ind)
% 	         if(locid(node).eq.0) then
% 	            numnpp = numnpp + 1
% 	            locid(node) = numnpp
% 	         endif
% 	         ind = ind + 1
% 	      enddo
% 
% %   End loop over elements
% 	   enddo
% 	   numelp = cel*csen
% 
%       if(debug) then
% 	dprt = .false.
% 	if(dprt) then
% 	   write(24,*) 'numnpp',numnpp
% 	   write(24,*) 'numelp',numelp
% 	   write(24,*) 'locid'
% 	   do j = 1,numnps
%          write(24,*) locid(j)
% 	   enddo
% 	   write(24,'(A)') ''
% 	endif
% 	endif

%%
%  -----------------------------------------------------
%   Localize star problem
%  -----------------------------------------------------

         locid = zeros(numel*celn,1); 

%   Assign local IDs to coarse nodes
         for i = 1:csnn
	      node = csn(i);
	      locid(node) = i;
         end
	   numnpp = csnn;

%   Loop over elements in cstar
         for k = 1:csen
	      elem = cse(k);
	      ind = cis(1,elem);
	      len = cis(1,elem+1) - ind;
	      nele = cis(2,elem);
	      ind = ind + nele;
%     Assign local IDs to cell nodes     	      
          for j = nele+1:len
	         node = ecn(ind);
             if(locid(node)==0)
	            numnpp = numnpp + 1;
	            locid(node) = numnpp;
             end
	         ind = ind + 1;
          end

%   End loop over elements
         end
	   numelp = cel*csen;


% %  -----------------------------------------------------
% %   Solve star problem
% %  -----------------------------------------------------
% 
% % Start control time
%       tt = etime(tary)
% 
% %	if(debug) then
% 
% 	   write(24,'(A,f17.7)') 'Start Starcontrol = ',tt
% 	   write(24,'(A)') ''
% 
% %	endif
% 	                 
% 	   call starcontrol(NodesOnElement,mr(np(152)),d,hr(np(151)),
% 	&                    hr(np(153)),u,cis,ecn,cse,cseb,
%      &                    csnb,locid,gammab,gammah,strwek,
%      &                    globalm)
% 
% % End control time
%       tt = etime(tary)
% 
% %	if(debug) then
% 
% 	   write(24,'(A,f17.7)') 'End Starcontrol = ',tt
% 	   write(24,'(A)') ''
% 
% %	endif

%%
%  -----------------------------------------------------
%   Solve star problem
%  -----------------------------------------------------

% 	   us = starcontrol(NodesOnElement,ixs,d,xs,us,u,cis,ecn,cse,cseb,csnb,locid,gammab,gammah,strong,globalm);

       starcontrol

% %-----------------------------------------------------
% % End loop over coarse nodes
% %-----------------------------------------------------
% 	   endif
% 
% %	endif
% 
% % End patch snode time
%       tt = etime(tary)
% 
% %	if(debug) then
% 
% 	   write(24,'(A,I4,A,f17.7)') 'End Patch ',snode,' Solut. = ',tt
% 	   write(24,'(A)') ''
% 
% %	endif
% 
% 100   continue
% 
% 	endif

%%
%-----------------------------------------------------
% End loop over coarse nodes
%-----------------------------------------------------
       end

%	endif

    end

    end

% %-----------------------------------------------------
% % Compute norms of fine scale solution
% %-----------------------------------------------------
% 
% % Start Norm time
%       tt = etime(tary)
% 
% %	if(debug) then
% 
% 	   write(24,'(A,f17.7)') 'Start Solut. Norms = ',tt
% 	   write(24,'(A)') ''
% 
% %	endif
% 
% 	call implicitnorms(NodesOnElement,mr(np(152)),hr(np(151)),hr(np(153)),cis,
% 	&                   iprob,strong)
% 
% %	endif %nen.ne.6
% 
% % End Norm time
%       tt = etime(tary)
% 
% %	if(debug) then
% 
% 	   write(24,'(A,f17.7)') 'End Solut. Norms = ',tt
% 	   write(24,'(A)') ''
% 
% %	endif
% 
%       return
% 
%       end

%%
%-----------------------------------------------------
% Compute norms of fine scale solution
%-----------------------------------------------------

% 	implicitnorms(NodesOnElement,ixs,xs,us,cis,iprob,strwek);
  implicitnorms1



