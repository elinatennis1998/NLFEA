% %**********************************************************************
% % 
%       subroutine implicitnorms(ix,ixs,xs,us,cis,iprob,strwek)
% %
% %...  Written by Tim Truster (Spring 2010)
%     No modifications when copied to NLFEA ver2
% %...  Program to compute norms of fine scale solution fields
% %
% %**********************************************************************
% 
% 
% 	implicit none
% 
%       integer         numnp,numel,nummat,nen,neq,ipr
%       common /cdata/  numnp,numel,nummat,nen,neq,ipr
% 
%       integer         ior,iow
%       common /iofile/ ior,iow
% 
%       logical debug
%       common /debugs/ debug
% 
%       integer         ndf,ndm,nen1,nst,nneq,ndl,nnlm,nadd
%       common /sdata/  ndf,ndm,nen1,nst,nneq,ndl,nnlm,nadd
% 
% 	include 'smparam.h'
% 	include 'stardata.h'
% 	include 'submeshdata.h'
% 
% %     Input Variables
%       integer ix(nen1,numel),ixs(nen,numel*cel),cis(2,numel+1)
% 	integer iprob,strwek
%       real*8  xs(ndm,numel*celn),us(ndfs,numel*celn-numnp)
% 
% %     Output Variables
% 
% %     Local Variables
% 	integer i,j,k,l,node,elem,cell,mat,nele,lint,nodeA,neleB
% 	integer ixc(nen)
% 	real*8  xc(ndm,celn),uc(ndf,nen),uf(ndfs),dux(ndfs),duy(ndfs)  
%       real*8    w,det,c1    
%       real*8    shl(3,10),shls(3,10),shg(3,10),shgs(3,10) 
% 
%       real*8    ufn,ufpnx,ufpny
%       real*8    ufl2(ndfs),ufix(ndfs),ufiy(ndfs),ufl2el(ndfs),
%      &          ufixel(ndfs),ufiyel(ndfs),ul2,uh1,dl10
% 	real*8    ufieff(numel)
%       real*8    litr,lits,celle
% 	real*8    bubble(3),sx(2,2)
% 	logical   der,bf
% 	integer   ib	
%       logical dprt 
% 
% %	include 'comblk.h'
% 	      
%       save
% 
% %	if(debug) then
% 	dprt = .true.
% 	if(dprt) then
% 	   write(24,*) 'Implicit fine scale'
% 	   write(24,*) 'node, ux, uy'
% 	   do j = 1,numnps-numnp
%          write(24,1001) j+numnp,(us(i,j),i=1,ndfs)
% 	   enddo
% 	   write(24,'(A)') ''
% 	   write(25,*) 'Implicit fine scale'
% 	   write(25,*) 'node, ux, uy'
% 	   do j = 1,numnps-numnp
%          write(25,1001) j+numnp,(us(i,j),i=1,ndfs)
% 	   enddo
% 	   write(25,'(A)') ''
%  1001 format (I7,2x,d20.14,2x,d20.14,2x,d20.14)
% 	endif
% %	endif
% 
% %.... clear the global arrays 
% 	call pzero (ufl2, ndfs)
% 	call pzero (ufix, ndfs)
% 	call pzero (ufiy, ndfs)
% 	call pzero (ufieff, numel)
% 
% %....	set shape function flags
%       ib = 0
% 	der = .false.
% 	bf = .false.
    
    ixc = zeros(nen,1);
	xc = zeros(ndm,nen);
    uc = zeros(ndf,nen);

    ufl2 = zeros(ndfs,1);
    ufix = zeros(ndfs,1);
    ufiy = zeros(ndfs,1);
	ufieff = zeros(numel,1);
    
%....	set shape function flags
    ib = 0;
	der = 0;
	bf = 0;

% %-----------------------------------------------------
% % Loop over elements in domain
% %-----------------------------------------------------
% 	do 110 elem = 1,numel
% 
% 	   mat = ix(nen1,elem)
%          nele = cis(2,elem)
%          if((nele.eq.3).or.(nele.eq.6)) then
%             neleB = 3
%          else
%             neleB = 4
%          endif
%          call IntPoint(nele,lint)
% 
%       if(debug) then
%       dprt = .true.
%       if(dprt) then
%          write(24,*) 'elem',elem
%          write(24,*) 'mat',mat
%          write(24,*) 'nele',nele
%          write(24,*) 'neleB',neleB
%          write(24,*) 'lint',lint
%          write(24,'(A)') ''
%       endif
%       endif
% %  -----------------------------------------------------
% %   Loop over cells in element
% %  -----------------------------------------------------
% 	   do 120 celle = 1,cel
% 
% %     Load cell coordinates and fine scale values
% 
%             cell = cel*(elem-1) + celle
% 	      call pzero(uc,ndfs*nen)
% 
%             do j = 1,nele
%                node = ixs(j,cell)
% 	         ixc(j) = node
% 	         nodeA = node - numnp
%                do i = 1,ndm
%                   xc(i,j) = xs(i,node)
%                enddo
% 	         if(nodeA.gt.0) then
% 	         do	i = 1,ndfs
% 	            uc(i,j) = us(i,nodeA)
% 	         enddo
% 	         else
% 	         do	i = 1,ndfs
% 	            uc(i,j) = 0.d0
% 	         enddo
% 	         endif
%             enddo
% 
%       if(debug) then
%       dprt = .true.
%       if(dprt) then
% 	   write(24,*) 'cell', cell
% 	   write(24,*) 'cell nodes'
%          do j = 1,nele
%             write(24,*) ixc(j)
%          enddo
%          write(24,*) 'cell coords'
%          do j = 1,nele
%             write(24,*) (xc(i,j),i=1,ndm)
%          enddo
%          write(24,'(A)') ''
% 	   write(24,*) 'cell values'
%          do j = 1,nele
%             write(24,*) (uc(i,j),i=1,ndfs)
%          enddo
%          write(24,'(A)') ''
%       endif
%       endif
% 
% %....	clear the element arrays
% 	      call pzero (ufl2el, ndfs)
% 	      call pzero (ufixel, ndfs)
% 	      call pzero (ufiyel, ndfs)
% 
% %    -----------------------------------------------------
% %     Loop over integration points
% %    -----------------------------------------------------
%             do 130 l=1,lint
% 
% %       Evaluate shape functions
% 
% %.... Compute Local & Global Element Shape Functions
%             if    (neleB.eq.3)then
%               call intpntt(l,lint,ib,w,litr,lits)
%               call shlt(litr,lits,nele,der,bf,shl,shls,bubble)
%               call shgt(xc,nele,shl,shls,nummat,nen,bf,der,det,shg,shgs,
%      &              bubble,sx) 
%             elseif(neleB.eq.4)then
%               call intpntq(l,lint,ib,w,litr,lits)
%               call shlq(litr,lits,nele,der,bf,shl,shls,bubble)
%               call shgq(xc,nele,shl,shls,nummat,nen,bf,der,det,shg,shgs,
%      &              bubble,sx) 
%             endif
% 	
%             c1 = det*w	
% %
% %....	clear the fine scale solutions
% %
% 	      call pzero ( uf, ndfs)
% 	      call pzero (dux, ndfs)
% 	      call pzero (duy, ndfs)
% 
% 	      do k=1,nele
% 	      do j=1,ndfs
% 
% 	      uf(j)  = uf(j)  + shg(3,k)*uc(j,k)
% 	      dux(j) = dux(j) + shg(1,k)*uc(j,k)
% 	      duy(j) = duy(j) + shg(2,k)*uc(j,k)
% 
%             enddo
%             enddo 
% 
% %	---------------------> Error Evaluation <---------------------
% 
% %....	loop over nodal vector
% 
% 	      do j=1,ndfs
% 
% 	      ufn   = c1 * ( uf(j)**2 )
% 	      ufpnx = c1 * ( dux(j)**2 )
% 	      ufpny = c1 * ( duy(j)**2 )
% 
% 	      ufl2el(j) = ufl2el(j) + ufn
% 	      ufixel(j) = ufixel(j) + ufpnx
% 	      ufiyel(j) = ufiyel(j) + ufpny
% 
% 	      enddo
% 
% %    -----------------------------------------------------
% %     End loop over integration points
% %    -----------------------------------------------------
%  130        continue
% 
% %....	add the element contribution to the global error evaluated
% 
% 	   do j = 1,ndfs
% 
% 	      ufl2(j) = ufl2(j) + ufl2el(j)
% 	      ufix(j) = ufix(j) + ufixel(j)
% 	      ufiy(j) = ufiy(j) + ufiyel(j)
% 
% 	   enddo
% 
% 	   ufieff(elem) = ufieff(elem) + ufixel(1) + ufiyel(1)
%      &                  + ufixel(2) + ufiyel(2)
% 
% %  -----------------------------------------------------
% %   End loop over cells
% %  -----------------------------------------------------
%  120     continue
% %-----------------------------------------------------
% % End loop over elements in domain
% %-----------------------------------------------------
%  110  continue
% 
% %	if(debug) then
% 	dprt = .true.
% 	if(dprt) then
% 	   write(24,*) 'Implicit fine norms^2'
% 	   write(24,*) 'elem, H1(uprime)'
% 	   do j = 1,numel
%          write(24,1002) j,ufieff(j)
% 	   enddo
% 	   write(24,'(A)') ''
% 	   write(25,*) 'Implicit fine norms^2'
% 	   write(25,*) 'elem, H1(uprime)'
% 	   do j = 1,numel
%          write(25,1002) j,ufieff(j)
% 	   enddo
% 	   write(25,'(A)') ''
%  1002 format (I7,2x,d20.14)
% 	endif
% %	endif
% 
% 	ul2 = dsqrt(ufl2(1)+ufl2(2))
% 	uh1 = dsqrt(ufix(1)+ufix(2)+ufiy(1)+ufiy(2))
% 
% %
% %....	calculate the log
% %
%       dl10  = dlog(10.d0)
%       if(ul2.gt.1.d-30)  ul2  = dlog(ul2) / dl10
%       if(uh1.gt.1.d-30)  uh1  = dlog(uh1) / dl10
% 
% 	write(11,2001)
% 	write(11,2101) numel,ul2,uh1,iprob,lint,strwek,minc,nstar
% 
% 	write(iow,2001)
% 	write(iow,2101) numel,ul2,uh1,iprob,lint,strwek,minc,nstar
%    
%  2001 format(//5x,' I M P L I C I T   F I N E   S C A L E S',// 
%      &        2x,'Elems',4x,'L2Uprime',7x,'H1Uprime' ,
%      &          4x,'Iprob',1x,'nint',1x,'strwek',1x,'minc',1x,'nstar')
%  2101	format(1x,i6,2x,e13.7,2x,e13.7,1x,i3,3x,i3,3x,i3,3x,i3,3x,i3)
% 
% 	return
% 
%       end

%%
%-----------------------------------------------------
% Loop over elements in domain
%-----------------------------------------------------
   for elem = 1:numel

	   mat = RegionOnElement(elem);
         nele = cis(2,elem);
         if((nele==3)||(nele==6))
            neleB = 3;
         else
            neleB = 4;
         end
         lint = IntPoint(nele);

%  -----------------------------------------------------
%   Loop over cells in element
%  -----------------------------------------------------
      for celle = 1:cel

%     Load cell coordinates and fine scale values

            cell = cel*(elem-1) + celle;
	        uc = zeros(ndfs,nen);

            for j = 1:nele
               node = ixs(j,cell);
	           ixc(j) = node;
% 	           nodeA = node - numnp;
               for i = 1:ndm
                  xc(i,j) = xs(i,node);
               end
%                if(nodeA>0)
                 for i = 1:ndfs
	               uc(i,j) = us(i,node);
                 end
%                else
%                  for i = 1:ndfs
% 	               uc(i,j) = 0.d0;
%                  end
%                end
            end

%....	clear the element arrays
            ufl2el = zeros(ndfs,1);
            ufixel = zeros(ndfs,1);
            ufiyel = zeros(ndfs,1);

%    -----------------------------------------------------
%     Loop over integration points
%    -----------------------------------------------------
            for l=1:lint

%       Evaluate shape functions

%.... Compute Local & Global Element Shape Functions
            if (neleB==3)
              [w,litr,lits] =  intpntt(l,lint,ib);
              [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
              [shg, shgs, det, be, sx] = shgt(xc,nel,shld,shls,nen,bf,der,be);
            else
              [w,litr,lits] =  intpntq(l,lint,ib);
              [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
              [shg, shgs, det, be, sx] = shgq(xc,nel,shld,shls,nen,bf,der,be);
            end
	
            c1 = det*w;	
%
%....	clear the fine scale solutions
%
          uf = zeros(ndfs,1);
          dux = zeros(ndfs,1);
          duy = zeros(ndfs,1);

          for k=1:nele
            for j=1:ndfs

	          uf(j)  = uf(j)  + shl(k)*uc(j,k);
	          dux(j) = dux(j) + shg(k,1)*uc(j,k);
	          duy(j) = duy(j) + shg(k,2)*uc(j,k);

            end
          end

%	---------------------> Error Evaluation <---------------------

%....	loop over nodal vector

          for j=1:ndfs

	      ufn   = c1 * ( uf(j)^2 );
	      ufpnx = c1 * ( dux(j)^2 );
	      ufpny = c1 * ( duy(j)^2 );

	      ufl2el(j) = ufl2el(j) + ufn;
	      ufixel(j) = ufixel(j) + ufpnx;
	      ufiyel(j) = ufiyel(j) + ufpny;

          end

%    -----------------------------------------------------
%     End loop over integration points
%    -----------------------------------------------------
            end

%....	add the element contribution to the global error evaluated

       for j = 1:ndfs

	      ufl2(j) = ufl2(j) + ufl2el(j);
	      ufix(j) = ufix(j) + ufixel(j);
	      ufiy(j) = ufiy(j) + ufiyel(j);

       end

	   ufieff(elem) = ufieff(elem) + ufixel(1) + ufiyel(1) + ufixel(2) + ufiyel(2);

%  -----------------------------------------------------
%   End loop over cells
%  -----------------------------------------------------
       end
%-----------------------------------------------------
% End loop over elements in domain
%-----------------------------------------------------
    end

	ul2 = sqrt(ufl2(1)+ufl2(2));
	uh1 = sqrt(ufix(1)+ufix(2)+ufiy(1)+ufiy(2));

%
%....	calculate the log
%
      dl10  = log(10.d0);
      ul2  = log(ul2) / dl10;
      uh1  = log(uh1) / dl10;
        fprintf('Implicit %i  %1.7e  %1.7e\n',numel,ul2,uh1)
        
    if ndfs == 3
        
	ul2p = sqrt(ufl2(3));
	uh1p = sqrt(ufix(3)+ufiy(3));

%
%....	calculate the log
%
      dl10  = log(10.d0);
      ul2p  = log(ul2p) / dl10;
      uh1p  = log(uh1p) / dl10;
        fprintf('Implicit %i  %1.7e  %1.7e\n',numel,ul2p,uh1p)
        
    end