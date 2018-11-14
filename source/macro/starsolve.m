%**********************************************************************
%
% 	subroutine starsolve(ix,ixs,d,xs,us,u,cis,cse,
%      &                      locid,globid,ldoft,
%      &                      gammah,strong,globalm)
%
%...  Written by Tim Truster (Spring 2010)
%     Modifications when copied to NLFEA ver2:
%        Taulist replaced by hr
%...  Program 
%
%**********************************************************************
% 
% 
%       implicit none
% 
%       integer         numnp,numel,nummat,nen,neq,ipr
%       common /cdata/  numnp,numel,nummat,nen,neq,ipr
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
%       include 'smparam.h'
%       include 'stardata.h'
%       include 'submeshdata.h'
% 	include 'starprob.h'
% 
% %     Input Variables
%       integer ix(nen1,numel),ixs(nen,numel*cel)
%       integer cis(2,numel+1),gammah(numnp-ndm),cse(numel)
%       integer locid(numel*celn),globid(numnpp),ldoft(ndfs,numnpp)
%       real*8  d(250,*),xs(ndm,numel*celn),us(ndfs,numel*celn-numnp)
%       real*8  u(ndf,numnp,3)
% 	logical strong,globalm
% 
%     Output Variables
%      real*8 
% 
%     Local Variables
%       integer i,j,k,node,elem,cell,mat,nele
%       integer ixe(nen),ixc(nen),ixl(2,nen)
%       real*8  xe(ndm,nen),xc(ndm,nen),xl(ndm,nen)
% 	real*8  ue(ndf,nen),ul(ndf,nen)
% 	real*8  temp
% 
%       integer edge,nodeA,nodeB,celli,cellj,numA,numB
%       logical efound,tractyn
% 
%       integer lint,celle,sslot,neleB,evenodd,edges(4),sslote
% 	integer edoft(ndfs*nen),numE
%       real*8  eK(nen*ndfs,nen*ndfs),eF(nen*ndfs)
%       real*8  MR,BR,MS,BS
% 	integer nA1,nAinc,nA2
% 
% 	real*8  Kdd11(neqp,neqp),Fd1(neqp)
% 	real*8  dsyswork(neqp)
% 	character ulop*1
% 	integer dsysinfo,IPIV(neqp)
% 
% 	real    tt,etime,tary(2)
% 	
%       logical dprt
% 
% 	ulop = 'U'
% 
% 	call pzero(Kdd11,neqp*neqp)
% 	call pzero(Fd1,neqp)
% 
% % Start assembly time
%       tt = etime(tary)
% 
% %	if(debug) then
% 
% 	   write(24,'(A,f17.7)') 'Start Assembly time = ',tt
% 	   write(24,'(A)') ''
% 
% %	endif

   ixe = zeros(nen,1);
   ixc = zeros(nen,1);
   ixl = zeros(2,nen);
   xe = zeros(ndm,nen);
   xc = zeros(ndm,nen);
   xl = zeros(ndm,nen);
   ue = zeros(ndf,nen);
   ul = zeros(ndf,nen);
   edges = zeros(4,1);
   edoft = zeros(ndfs*nen);
   
   eK = zeros(nen*ndfs,nen*ndfs);
   eF = zeros(nen*ndfs,1);
   KDD11 = zeros(neqp,neqp);
   FD1 = zeros(neqp,1);

% %  -----------------------------------------------------
% %   Loop over elements in star problem
% %  -----------------------------------------------------
%       do 10 k = 1,csen
% 
%          elem = cse(k)
%          mat = ix(nen1,elem)
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
% 
% %     Set flag if element is in level 1 star
%          if(globalm) then
% 	      sslot = 0
% 	      efound = .true.
% 	   else
%             j = 1
%             node = ix(j,elem)
%             if(node.eq.snode) then
%                efound = .true.
%                sslot = j
%             else
%                efound = .false.
%             endif
% 
%             do while(.not.efound.and.(j.lt.neleB))
%                j = j + 1
%                node = ix(j,elem)
%                if(node.eq.snode) then
%                   efound = .true.
%                   sslot = j
%                endif
%             enddo
% 	   endif            
% 
%       if(debug) then
%       dprt = .true.
%       if(dprt) then
%          write(24,*) 'efound in epatch: ',efound
% 	   write(24,*) 'sslot',sslot
%          write(24,'(A)') ''
%       endif
%       endif
% 
%          if(.not.efound) then	%element not in level 1 star
% %    -----------------------------------------------------
% %     Compute cell stiffness only
% %    -----------------------------------------------------

%%
%  -----------------------------------------------------
%   Loop over elements in star problem
%  -----------------------------------------------------
      for k = 1:csen

         elem = cse(k);
         mat = RegionOnElement(elem);
         iel = MatTypeTable(2,mat);
         nonlin = MatTypeTable(3,mat);
         
         nele = cis(2,elem);
         if((nele==3)||(nele==6))
            neleB = 3;
         else
            neleB = 4;
         end
         lint = IntPoint(nele);

%     Set flag if element is in level 1 star
         if(globalm)
	      sslot = 0;
	      efound = 1;
	   else
            j = 1;
            node = NodesOnElement(j,elem);
            if(node==snode)
               efound = 1;
               sslot = j;
            else
               efound = 0;
            end

            while((efound==0)&&(j<neleB))
               j = j + 1;
               node = NodesOnElement(j,elem);
               if(node==snode)
                  efound = 1;
                  sslot = j;
               end
            end
         end
%    -----------------------------------------------------
%     Compute cell stiffness only
%    -----------------------------------------------------
% 	      
% %    -----------------------------------------------------
% %     Loop over cells in element
% %    -----------------------------------------------------
%             do 20 celle = 1,cel
% 
% %       Load cell coordinates, dof IDs
% 
%                cell = cel*(elem-1) + celle
%                call pzeroi(ixc,nen)
% 
%                do j = 1,nele
%                   node = ixs(j,cell)
%                   ixc(j) = locid(node)
%                   do i = 1,ndm
%                      xc(i,j) = xs(i,node)
%                   enddo
%                enddo
% 
%       if(debug) then
%       dprt = .true.
%       if(dprt) then
% 	   write(24,*) 'cell', cell
%          write(24,*) 'cell coords'
%          do j = 1,nele
%             write(24,*) (xc(i,j),i=1,ndm)
%          enddo
%          write(24,'(A)') ''
%       endif
%       endif
% 
% %      -----------------------------------------------------
% %       Compute cell interior stiffness
% %      -----------------------------------------------------
%                call getKcell(xc,d(1,mat),ndfs,ndm,lint,nele,nen,
%      &                       nummat,eK)
% 
%       if(debug) then
%       dprt = .true.
%       if(dprt) then
%          write(24,*) 'cell stiffness'
%          do j = 1,nele*ndfs
%             write(24,1002) (eK(i,j),i=1,nele*ndfs)
%          enddo
%  1002 format (1x,d12.6,1x,d12.6,1x,d12.6,1x,d12.6,1x,d12.6,1x,d12.6,1x,
%      &        d12.6,1x,d12.6)
%          write(24,'(A)') ''
%       endif
%       endif
% 
% %       Assemble cell quantities into star quantities
% 
%                call LocToGlobDOF(ixc,ldoft,nele,nen,ndfs,edoft)
% 
%       if(debug) then
%       dprt = .true.
%       if(dprt) then
%          write(24,*) 'edoft'
%          do j = 1,nele*ndfs
%             write(24,*) edoft(j)
%          enddo
%          write(24,'(A)') ''
%       endif
%       endif
% 
% 	         call assembleK(nele,nen,ndfs,edoft,eK,Kdd11)
% 
% %    -----------------------------------------------------
% %     End loop over cells
% %    -----------------------------------------------------
%  20         continue
% 
%          else
% 
%%
%    -----------------------------------------------------
%     Loop over cells in element
%    -----------------------------------------------------
%             for celle = 1:cel
% 
% %       Load cell coordinates, dof IDs
% 
%                cell = cel*(elem-1) + celle;
%                ixc = zeros(nen,1);
% 
%                for j = 1:nele
%                   node = ixs(j,cell);
%                   ixc(j) = locid(node);
%                   for i = 1:ndm
%                      xc(i,j) = xs(i,node);
%                   end
%                end
% 
% %      -----------------------------------------------------
% %       Compute cell interior stiffness
% %      -----------------------------------------------------
%                eK = getKcell(xc,d(:,mat),xe,ue,MR,BR,MS,BS,ndfs,ndm,lint,nele,nen,nummat);
% 
% %       Assemble cell quantities into star quantities
% 
%                edoft = LocToGlobDOF(ixc, ldoft', nele, ndfs);
% 
% 	           assembleK1
% 
% %    -----------------------------------------------------
% %     End loop over cells
% %    -----------------------------------------------------
%             end
% 
%          else
% %    -----------------------------------------------------
% %     Compute cell stiffness and residual
% %    -----------------------------------------------------
%            
% %     Extract element coordinates
%             do j = 1,nele
%                node = ix(j,elem)
%                ixe(j) = node
%                do i = 1,ndm
%                   xe(i,j) = xs(i,node)
%                enddo
%                do i = 1,ndf
%                   ue(i,j) = u(i,node,1)  %%%%%????? check the indices on this
%                enddo
%             enddo
% 
%       if(debug) then
%       dprt = .true.
%       if(dprt) then
%          write(24,*) 'elem',elem
%          write(24,*) 'local nodes'
%          do j = 1,nele
%             write(24,*) ixe(j)
%          enddo
%          write(24,'(A)') ''
%          write(24,*) 'elem coords'
%          do j = 1,nele
%             write(24,*) (xe(i,j),i=1,ndm)
%          enddo
%          write(24,*) 'elem nodal values'
%          do j = 1,nele
%             write(24,*) (ue(i,j),i=1,ndf)
%          enddo
%          write(24,'(A)') ''
%       endif
%       endif               
% 	      
% %    -----------------------------------------------------
% %     Loop over cells in element
% %    -----------------------------------------------------
% %     If Triangular Element
%             if(neleB.eq.3) then
% 
%                celle = 0
% 
%                do 40 cellj = 1,minc
%                   do 30 celli = 1,2*(minc-cellj)+1
%                  
%                   celle = celle + 1
%                 
% %       Load cell coordinates, dof IDs
% 
%                   cell = cel*(elem-1) + celle
%                   call pzeroi(ixc,nen)
%    
%                   do j = 1,nele
%                      node = ixs(j,cell)
%                      ixc(j) = locid(node)
%                      do i = 1,ndm
%                         xc(i,j) = xs(i,node)
%                      enddo
%                   enddo
% 
% %       Set element-cell local coordinate map
%                   evenodd = (celli+1)/2
%                   if(celli/2.eq.evenodd) then %even celli
% 
% 	               %adjust order of cell flags for proper assembly
% 	               %due to reversal of element local coordinates
% 	               if(nele.eq.3) then
% 	                  i = ixc(1)
% 	                  ixc(1) = ixc(2)
% 	                  ixc(2) = ixc(3)
% 	                  ixc(3) = i
% 	                  do i = 1,ndm
% 	                     temp = xc(i,1)
% 	                     xc(i,1) = xc(i,2)
% 	                     xc(i,2) = xc(i,3)
% 	                     xc(i,3) = temp
% 	                  enddo
% 	               else
% 	                  i = ixc(1)
% 	                  ixc(1) = ixc(2)
% 	                  ixc(2) = ixc(3)
% 	                  ixc(3) = i
% 	                  i = ixc(4)
% 	                  ixc(4) = ixc(5)
% 	                  ixc(5) = ixc(6)
% 	                  ixc(6) = i
% 	                  do i = 1,ndm
% 	                     temp = xc(i,1)
% 	                     xc(i,1) = xc(i,2)
% 	                     xc(i,2) = xc(i,3)
% 	                     xc(i,3) = temp
% 	                     temp = xc(i,4)
% 	                     xc(i,4) = xc(i,5)
% 	                     xc(i,5) = xc(i,6)
% 	                     xc(i,6) = temp
% 	                  enddo
% 	               endif
%                      MR = -1.d0/minc
%                      BR = celli/(2.d0*minc)
%                      MS = -1.d0/minc
%                      BS = cellj/dble(minc)
%                   else
%                      MR = 1.d0/minc
%                      BR = ((celli+1)/2.d0-1.d0)/minc
%                      MS = 1.d0/minc
%                      BS = (cellj-1.d0)/minc
%                   endif
% 
%       if(debug) then
%       dprt = .true.
%       if(dprt) then
%          write(24,*) 'cell',cell
%          write(24,*) 'local nodes'
%          do j = 1,nele
%             write(24,*) ixc(j)
%          enddo
%          write(24,'(A)') ''
%          write(24,*) 'cell coords'
%          do j = 1,nele
%             write(24,*) (xc(i,j),i=1,ndm)
%          enddo
%          write(24,'(A)') ''
%          write(24,*) 'MR',MR
%          write(24,*) 'BR',BR
%          write(24,*) 'MS',MS
%          write(24,*) 'BS',BS
%          write(24,'(A)') ''
%       endif
%       endif
% 
% %      -----------------------------------------------------
% %       Compute cell interior stiffness
% %      -----------------------------------------------------
%                 call getKFcell(xc,d(1,mat),xe,ue,ndfs,ndm,lint,
%      &                      nele,nen,nummat,strong,MR,BR,MS,BS,
%      &                      sslot,eK,eF)
% 
%       if(debug) then
%       dprt = .true.
%       if(dprt) then
%          write(24,*) 'cell stiffness'
%          do j = 1,nele*ndfs
%             write(24,1002) (eK(i,j),i=1,nele*ndfs)
%          enddo
%          write(24,*) 'cell residual'
%          do j = 1,nele*ndfs
%             write(24,'(d12.6)') eF(j)
%          enddo
%          write(24,'(A)') ''
%       endif
%       endif
% 
% %       Assemble cell quantities into star quantities
% 
%                call LocToGlobDOF(ixc,ldoft,nele,nen,ndfs,edoft)
% 
%       if(debug) then
%       dprt = .true.
%       if(dprt) then
%          write(24,*) 'edoft'
%          do j = 1,nele*ndfs
%             write(24,*) edoft(j)
%          enddo
%          write(24,'(A)') ''
%       endif
%       endif
% 
% 	         call assembleKF(nele,nen,ndfs,edoft,eK,eF,Kdd11,Fd1)
% 
% %    -----------------------------------------------------
% %     End loop over cells
% %    -----------------------------------------------------
%  30               continue
%  40            continue

%%
%    -----------------------------------------------------
%     Compute cell stiffness and residual
%    -----------------------------------------------------
           
%     Extract element coordinates
            for j = 1:nele
               node = NodesOnElement(j,elem);
               ixe(j) = node;
               for i = 1:ndm
                  xe(i,j) = xs(i,node);
               end
               for i = 1:ndf
                  ue(i,j) = u(i,node);
               end
            end
            
            % Extract tau matrix for parent element
            startau
	      
%    -----------------------------------------------------
%     Loop over cells in element
%    -----------------------------------------------------
%     If Triangular Element
            if(neleB==3)

               celle = 0;

               for cellj = 1:minc
                  for celli = 1:2*(minc-cellj)+1
                 
                  celle = celle + 1;
                
%       Load cell coordinates, dof IDs

                  cell = cel*(elem-1) + celle;
                  ixc = zeros(nen,1);
   
                  for j = 1:nele
                     node = ixs(j,cell);
                     ixc(j) = locid(node);
                     for i = 1:ndm
                        xc(i,j) = xs(i,node);
                     end
                  end

%       Set element-cell local coordinate map
                  evenodd = floor((celli+1)/2);
                  if(celli/2==evenodd) %even celli

	               %adjust order of cell flags for proper assembly
	               %due to reversal of element local coordinates
                   if(nele==3)
	                  i = ixc(1);
	                  ixc(1) = ixc(2);
	                  ixc(2) = ixc(3);
	                  ixc(3) = i;
                      for i = 1:ndm
	                     temp = xc(i,1);
	                     xc(i,1) = xc(i,2);
	                     xc(i,2) = xc(i,3);
	                     xc(i,3) = temp;
                      end
	               else
	                  i = ixc(1);
	                  ixc(1) = ixc(2);
	                  ixc(2) = ixc(3);
	                  ixc(3) = i;
	                  i = ixc(4);
	                  ixc(4) = ixc(5);
	                  ixc(5) = ixc(6);
	                  ixc(6) = i;
                      for i = 1:ndm
	                     temp = xc(i,1);
	                     xc(i,1) = xc(i,2);
	                     xc(i,2) = xc(i,3);
	                     xc(i,3) = temp;
	                     temp = xc(i,4);
	                     xc(i,4) = xc(i,5);
	                     xc(i,5) = xc(i,6);
	                     xc(i,6) = temp;
                      end
                   end
                     MR = -1.d0/minc;
                     BR = celli/(2.d0*minc);
                     MS = -1.d0/minc;
                     BS = cellj/(minc);
                  else
                     MR = 1.d0/minc;
                     BR = ((celli+1)/2.d0-1.d0)/minc;
                     MS = 1.d0/minc;
                     BS = (cellj-1.d0)/minc;
                  end
                  
                  if(efound) == 0

%      -----------------------------------------------------
%       Compute cell interior stiffness
%      -----------------------------------------------------
               stargetK

%       Assemble cell quantities into star quantities

               edoft = LocToGlobDOF(ixc, ldoft', nele, ndfs);

	           assembleK1
                      
                  else

%      -----------------------------------------------------
%       Compute cell interior stiffness
%      -----------------------------------------------------
                stargetKF

%       Assemble cell quantities into star quantities

               edoft = LocToGlobDOF(ixc, ldoft', nele, ndfs);

	           assembleKF1
               
                  end

%    -----------------------------------------------------
%     End loop over cells
%    -----------------------------------------------------
                  end
               end

% %     Else Quadrilateral Element
%             else
%               
%                celle = 0
%                do 60 cellj = 1,minc
%                   do 50 celli = 1,minc
%                  
%                   celle = celle + 1
% 
% %       Load cell coordinates, dof IDs
% 
%                   cell = cel*(elem-1) + celle
%                   call pzeroi(ixc,nen)
%    
%                   do j = 1,nele
%                      node = ixs(j,cell)
%                      ixc(j) = locid(node)
%                      do i = 1,ndm
%                         xc(i,j) = xs(i,node)
%                      enddo
%                   enddo
% 
% %       Set element-cell local coordinate map
%                   %MR = 1.d0/dble(m) 
%                   MR = 1.d0/minc
%                   BR = -1.d0+(2.d0*celli-1.d0)/minc
%                   MS = 1.d0/minc
%                   BS = -1.d0+(2.d0*cellj-1.d0)/minc
% 
%       if(debug) then
%       dprt = .true.
%       if(dprt) then
%          write(24,*) 'cell',cell
%          write(24,*) 'local nodes'
%          do j = 1,nele
%             write(24,*) ixc(j)
%          enddo
%          write(24,'(A)') ''
%          write(24,*) 'cell coords'
%          do j = 1,nele
%             write(24,*) (xc(i,j),i=1,ndm)
%          enddo
%          write(24,'(A)') ''
%          write(24,*) 'MR',MR
%          write(24,*) 'BR',BR
%          write(24,*) 'MS',MS
%          write(24,*) 'BS',BS
%          write(24,'(A)') ''
%       endif
%       endif
% 
% %      -----------------------------------------------------
% %       Compute cell interior stiffness
% %      -----------------------------------------------------
%                 call getKFcell(xc,d(1,mat),xe,ue,ndfs,ndm,lint,
%      &                      nele,nen,nummat,strong,MR,BR,MS,BS,
%      &                      sslot,eK,eF)
% 
%       if(debug) then
%       dprt = .true.
%       if(dprt) then
%          write(24,*) 'cell stiffness'
%          do j = 1,nele*ndfs
%             write(24,1002) (eK(i,j),i=1,nele*ndfs)
%          enddo
%          write(24,*) 'cell residual'
%          do j = 1,nele*ndfs
%             write(24,'(d12.6)') eF(j)
%          enddo
%          write(24,'(A)') ''
%       endif
%       endif
% 
% %       Assemble cell quantities into star quantities
% 
%                call LocToGlobDOF(ixc,ldoft,nele,nen,ndfs,edoft)
% 
%       if(debug) then
%       dprt = .true.
%       if(dprt) then
%          write(24,*) 'edoft'
%          do j = 1,nele*ndfs
%             write(24,*) edoft(j)
%          enddo
%          write(24,'(A)') ''
%       endif
%       endif
% 
% 	         call assembleKF(nele,nen,ndfs,edoft,eK,eF,Kdd11,Fd1)
% 
% %    -----------------------------------------------------
% %     End loop over cells
% %    -----------------------------------------------------
%  50               continue
%  60            continue
% 
% 	      endif

%%
%     Else Quadrilateral Element
            else
              
               celle = 0;
               for cellj = 1:minc
                  for celli = 1:minc
                 
                  celle = celle + 1;

%       Load cell coordinates, dof IDs

                  cell = cel*(elem-1) + celle;
                  ixc = zeros(nen,1);
   
                  for j = 1:nele
                     node = ixs(j,cell);
                     ixc(j) = locid(node);
                     for i = 1:ndm
                        xc(i,j) = xs(i,node);
                     end
                  end

%       Set element-cell local coordinate map
                  %MR = 1.d0/dble(m) 
                  MR = 1.d0/minc;
                  BR = -1.d0+(2.d0*celli-1.d0)/minc;
                  MS = 1.d0/minc;
                  BS = -1.d0+(2.d0*cellj-1.d0)/minc;
                  
                  if(efound) == 0

%      -----------------------------------------------------
%       Compute cell interior stiffness
%      -----------------------------------------------------
               stargetK

%       Assemble cell quantities into star quantities

               edoft = LocToGlobDOF(ixc, ldoft', nele, ndfs);

	           assembleK1
                      
                  else

%      -----------------------------------------------------
%       Compute cell interior stiffness
%      -----------------------------------------------------
                stargetKF

%       Assemble cell quantities into star quantities

               edoft = LocToGlobDOF(ixc, ldoft', nele, ndfs);

	           assembleKF1
               
                  end

%    -----------------------------------------------------
%     End loop over cells
%    -----------------------------------------------------
                  end
               end

            end

         if(efound==1)
% %    -----------------------------------------------------
% %     Element is in level 1 of star: compute boun. resid.
% %    -----------------------------------------------------
% 
% %     If Triangular Element
%             if(neleB.eq.3) then
% 
% 	         if(globalm) then
% 	            numE = 3
% 	            do i = 1,numE
% 	               edges(i) = i
% 	            enddo
% 	         else
% 	            numE = 2
%                   if(sslot.eq.1) then
% 	               edges(1) = 1
% 	               edges(2) = 3
%                   elseif(sslot.eq.2) then
% 	               edges(1) = 2
% 	               edges(2) = 1
% 	            else
% 	               edges(1) = 3
% 	               edges(2) = 2
%                   endif
% 	         endif
% 
% %    -----------------------------------------------------
% %     "Loop" over element edges on interior of star
% %    -----------------------------------------------------
%                do 80 edge = 1,numE
% 
%                if(sngammah) then
% 
%                   nodeA = edges(edge)+1
% 			    if(nodeA.gt.neleB) then
% 			       nodeA = 1
% 		        endif
% 	            nodeA = ixe(nodeA)
% 			                   
% %        Check if nodeA is on the traction boundary of domain
%                   call binsearch(nodeA,gammah,1,ngh,numA
%      *			               ,numB,tractyn)
% 
% 	            if(tractyn) then
% 
% 	            nodeB = edges(edge)
% 	            nodeB = ixe(nodeB)
% 			                   
% %        Check if nodeA is on the traction boundary of domain
%                   call binsearch(nodeB,gammah,1,ngh,numA
%      *			               ,numB,tractyn)
% 
% 	            endif
% 
% 	         else
% 	            tractyn = .false.
% 	         endif
% 
%       if(debug) then
%       dprt = .true.
%       if(dprt) then
%          write(24,*) 'tractyn: ',tractyn
% 	   write(24,*) 'nodeA',nodeA
% 	     write(24,*) 'nodeB',nodeB
%          write(24,'(A)') ''
%       endif
%       endif
% 
% 	         if(strong.or.tractyn) then
% 	         if(edges(edge).eq.1) then
% 	            nA1 = -1   
% 	            nAinc = 0
% 	            nA2 = 2
% 	            ixl(1,1) = 1
% 	            ixl(1,2) = 2
% 	            ixl(1,3) = 3
% 	            ixl(2,1) = 1
% 	            ixl(2,2) = 2
% 	            ixl(2,3) = 3
% 	            if(nele.eq.6) then
% 	               ixl(1,4) = 4
% 	               ixl(1,5) = 5
% 	               ixl(1,6) = 6
% 	               ixl(2,4) = 4
% 	               ixl(2,5) = 5
% 	               ixl(2,6) = 6
% 	            endif
% 	         elseif(edges(edge).eq.2) then
% 	            nA1 = 0
% 	            nAinc = -2
% 	            nA2 = 2*minc+1
% 	            ixl(1,1) = 2
% 	            ixl(1,2) = 3
% 	            ixl(1,3) = 1
% 	            ixl(2,1) = 3
% 	            ixl(2,2) = 1
% 	            ixl(2,3) = 2
% 	            if(nele.eq.6) then
% 	               ixl(1,4) = 5
% 	               ixl(1,5) = 6
% 	               ixl(1,6) = 4
% 	               ixl(2,4) = 6
% 	               ixl(2,5) = 4
% 	               ixl(2,6) = 5
% 	            endif
% 	         else
% 	            nA1 = cel+1
% 	            nAinc = -2
% 	            nA2 = 1
% 	            ixl(1,1) = 3
% 	            ixl(1,2) = 1
% 	            ixl(1,3) = 2
% 	            ixl(2,1) = 2
% 	            ixl(2,2) = 3
% 	            ixl(2,3) = 1
% 	            if(nele.eq.6) then
% 	               ixl(1,4) = 6
% 	               ixl(1,5) = 4
% 	               ixl(1,6) = 5
% 	               ixl(2,4) = 5
% 	               ixl(2,5) = 6
% 	               ixl(2,6) = 4
% 	            endif
% 	         endif
% 	         if(globalm) then
% 	            sslote = 0
% 	         else
% 	            sslote = ixl(2,sslot)
% 	         endif
% 
% 	         do j = 1,nele
% 	            node = ixl(1,j)
%                   do i = 1,ndm
%                      xl(i,j) = xe(i,node)
%                   enddo
%                   do i = 1,ndf
%                      ul(i,j) = ue(i,node)
%                   enddo
%                enddo
% 
%       if(debug) then
%       dprt = .true.
%       if(dprt) then
%          write(24,*) 'tractyn: ',tractyn
% 	   write(24,*) 'sslote',sslote
%          write(24,*) 'local nodes'
%          do j = 1,nele
%             write(24,*) ixl(1,j)
%          enddo
%          write(24,'(A)') ''
%          write(24,*) 'elem coords'
%          do j = 1,nele
%             write(24,*) (xl(i,j),i=1,ndm)
%          enddo
%          write(24,*) 'elem nodal values'
%          do j = 1,nele
%             write(24,*) (ul(i,j),i=1,ndf)
%          enddo
%          write(24,'(A)') ''
%       endif
%       endif
% 
% %    -----------------------------------------------------
% %     Loop over cells in element
% %    -----------------------------------------------------
% 
% 	         cell = cel*(elem-1) + nA1
% 
% 	if(debug) then
% 	dprt = .true.
% 	if(dprt) then
% 	   write(24,*) 'cell',cell
% 	   write(24,*) 'nA1',nA1
% 	   write(24,*) 'nAinc',nAinc
% 	   write(24,*) 'nA2',nA2
% 	   write(24,'(A)') ''
% 	endif
% 	endif
% 
%                do 70 celli = 1,minc
% 
% %       Load cell coordinates, dof IDs
% 
%                   cell = cell + nA2 + celli*nAinc
%                   call pzeroi(ixc,nen)
%    
%                   do j = 1,nele
%                      node = ixs(ixl(1,j),cell)
%                      ixc(j) = locid(node)
%                      do i = 1,ndm
%                         xc(i,j) = xs(i,node)
%                      enddo
%                   enddo
% 
% %       Set element-cell local coordinate map
%                   MR = 1.d0/minc
%                   BR = dble(celli-1)/dble(minc)
%                   MS = 0.d0
%                   BS = 0.d0
% 
%       if(debug) then
%       dprt = .true.
%       if(dprt) then
%          write(24,*) 'cell',cell
%          write(24,*) 'local nodes'
%          do j = 1,nele
%             write(24,*) ixc(j)
%          enddo
%          write(24,'(A)') ''
%          write(24,*) 'cell coords'
%          do j = 1,nele
%             write(24,*) (xc(i,j),i=1,ndm)
%          enddo
%          write(24,'(A)') ''
%          write(24,*) 'MR',MR
%          write(24,*) 'BR',BR
%          write(24,*) 'MS',MS
%          write(24,*) 'BS',BS
%          write(24,'(A)') ''
%       endif
%       endif
% 
% %      -----------------------------------------------------
% %       Compute cell boundary residual
% %      -----------------------------------------------------
%                    call getFbcell(xc,d(1,mat),xl,ul,ndfs,ndm,
%      &                      nele,nen,nummat,strong,MR,BR,MS,BS,
%      &                      sslote,tractyn,eF)
% 
%       if(debug) then
%       dprt = .true.
%       if(dprt) then
%          write(24,*) 'cell residual'
%          do j = 1,nele*ndfs
%             write(24,'(d12.6)') eF(j)
%          enddo
%          write(24,'(A)') ''
%       endif
%       endif
% 
% %       Assemble cell quantities into star quantities
% 
%                call LocToGlobDOF(ixc,ldoft,nele,nen,ndfs,edoft)
% 
%       if(debug) then
%       dprt = .true.
%       if(dprt) then
%          write(24,*) 'edoft'
%          do j = 1,nele*ndfs
%             write(24,*) edoft(j)
%          enddo
%          write(24,'(A)') ''
%       endif
%       endif
% 
% 	           call assembleF(nele,nen,ndfs,edoft,eF,Fd1)
% 
% %    -----------------------------------------------------
% %     End loop over cells
% %    -----------------------------------------------------
%  70            continue
%                endif
% 
% %    -----------------------------------------------------
% %     End loop over element edges
% %    -----------------------------------------------------
%  80            continue
% 
%

%%
%    -----------------------------------------------------
%     Element is in level 1 of star: compute boun. resid.
%    -----------------------------------------------------

%     If Triangular Element
            if(neleB==3)

             if(globalm)
	            numE = 3;
                for i = 1:numE
	               edges(i) = i;
                end
	         else
	            numE = 2;
                  if(sslot==1)
	               edges(1) = 1;
	               edges(2) = 3;
                  elseif(sslot==2)
	               edges(1) = 2;
	               edges(2) = 1;
	            else
	               edges(1) = 3;
	               edges(2) = 2;
                  end
             end

%    -----------------------------------------------------
%     "Loop" over element edges on interior of star
%    -----------------------------------------------------
               for edge = 1:numE

               if(sngammah)

                  nodeA = edges(edge)+1;
                  if(nodeA>neleB)
			         nodeA = 1;
                  end
	              nodeA = ixe(nodeA);
			                   
%        Check if nodeA is on the traction boundary of domain
                  [numA,numB,tractyn] = binsearch(nodeA,gammah,1,ngh);

                  if(tractyn)

                    nodeB = edges(edge);
                    nodeB = ixe(nodeB);
			                   
%        Check if nodeA is on the traction boundary of domain
                    [numA,numB,tractyn] = binsearch(nodeB,gammah,1,ngh);

                  end

               else
	              tractyn = 0;
               end

             if(strong||tractyn)
             if(edges(edge)==1)
	            nA1 = -1;
	            nAinc = 0;
	            nA2 = 2;
	            ixl(1,1) = 1;
	            ixl(1,2) = 2;
	            ixl(1,3) = 3;
	            ixl(2,1) = 1;
	            ixl(2,2) = 2;
	            ixl(2,3) = 3;
                if(nele==6)
	               ixl(1,4) = 4;
	               ixl(1,5) = 5;
	               ixl(1,6) = 6;
	               ixl(2,4) = 4;
	               ixl(2,5) = 5;
	               ixl(2,6) = 6;
                end
	         elseif(edges(edge)==2)
	            nA1 = 0;
	            nAinc = -2;
	            nA2 = 2*minc+1;
	            ixl(1,1) = 2;
	            ixl(1,2) = 3;
	            ixl(1,3) = 1;
	            ixl(2,1) = 3;
	            ixl(2,2) = 1;
	            ixl(2,3) = 2;
                if(nele==6)
	               ixl(1,4) = 5;
	               ixl(1,5) = 6;
	               ixl(1,6) = 4;
	               ixl(2,4) = 6;
	               ixl(2,5) = 4;
	               ixl(2,6) = 5;
                end
	         else
	            nA1 = cel+1;
	            nAinc = -2;
	            nA2 = 1;
	            ixl(1,1) = 3;
	            ixl(1,2) = 1;
	            ixl(1,3) = 2;
	            ixl(2,1) = 2;
	            ixl(2,2) = 3;
	            ixl(2,3) = 1;
                if(nele==6)
	               ixl(1,4) = 6;
	               ixl(1,5) = 4;
	               ixl(1,6) = 5;
	               ixl(2,4) = 5;
	               ixl(2,5) = 6;
	               ixl(2,6) = 4;
                end
             end
             if(globalm)
	            sslote = 0;
	         else
	            sslote = ixl(2,sslot);
             end

             for j = 1:nele
	            node = ixl(1,j);
                  for i = 1:ndm
                     xl(i,j) = xe(i,node);
                  end
                  for i = 1:ndf
                     ul(i,j) = ue(i,node);
                  end
             end
             
%    -----------------------------------------------------
%     Loop over cells in element
%    -----------------------------------------------------

	         cell = cel*(elem-1) + nA1;

               for celli = 1:minc

%       Load cell coordinates, dof IDs

                  cell = cell + nA2 + celli*nAinc;
                  ixc = zeros(nen,1);
   
                  for j = 1:nele
                     node = ixs(ixl(1,j),cell);
                     ixc(j) = locid(node);
                     for i = 1:ndm
                        xc(i,j) = xs(i,node);
                     end
                  end

%       Set element-cell local coordinate map
                  MR = 1.d0/minc;
                  BR = (celli-1)/(minc);
                  MS = 0.d0;
                  BS = 0.d0;
                  
%      -----------------------------------------------------
%       Compute cell boundary residual
%      -----------------------------------------------------
                   stargetFb

%       Assemble cell quantities into star quantities

               edoft = LocToGlobDOF(ixc,ldoft',nele,ndfs);

	           assembleF

%    -----------------------------------------------------
%     End loop over cells
%    -----------------------------------------------------
               end
             end

%    -----------------------------------------------------
%     End loop over element edges
%    -----------------------------------------------------
               end


% %     Else Quadrilateral Element
%             else
% 
% 	         if(globalm) then
% 	            numE = 4
% 	            do i = 1,numE
% 	               edges(i) = i
% 	            enddo
% 	         else
% 	            numE = 2
%                   if(sslot.eq.1) then
% 	               edges(1) = 1
% 	               edges(2) = 4
%                   elseif(sslot.eq.2) then
% 	               edges(1) = 2
% 	               edges(2) = 1
% 	            elseif(sslot.eq.3) then
% 	               edges(1) = 3
% 	               edges(2) = 2
% 	            else
% 	               edges(1) = 4
% 	               edges(2) = 3
%                   endif
% 	         endif
% 
% %    -----------------------------------------------------
% %     "Loop" over element edges on interior of star
% %    -----------------------------------------------------
%                do 100 edge = 1,numE
% 
%                if(sngammah) then
% 
%                   nodeA = edges(edge)+1
% 			    if(nodeA.gt.neleB) then
% 			       nodeA = 1
% 		        endif
% 	            nodeA = ixe(nodeA)
% 			                   
% %        Check if nodeA is on the traction boundary of domain
%                   call binsearch(nodeA,gammah,1,ngh,numA
%      *			               ,numB,tractyn)
% 	
% 	            if(tractyn) then
% 
% 	            nodeB = edges(edge)
% 	            nodeB = ixe(nodeB)
% 			                   
% %        Check if nodeA is on the traction boundary of domain
%                   call binsearch(nodeB,gammah,1,ngh,numA
%      *			               ,numB,tractyn)
% 
% 	            endif
% 
% 	         else
% 	            tractyn = .false.
% 	         endif
% 
%       if(debug) then
%       dprt = .true.
%       if(dprt) then
%          write(24,*) 'tractyn: ',tractyn
% 	   write(24,*) 'nodeA',nodeA
% 	     write(24,*) 'nodeB',nodeB
%          write(24,'(A)') ''
%       endif
%       endif
% 
% 	         if(strong.or.tractyn) then
% 	         if(edges(edge).eq.1) then
% 	            nA1 = 0
% 	            nA2 = 1
% 	            ixl(1,1) = 1
% 	            ixl(1,2) = 2
% 	            ixl(1,3) = 3
% 	            ixl(1,4) = 4
% 	            ixl(2,1) = 1
% 	            ixl(2,2) = 2
% 	            ixl(2,3) = 3
% 	            ixl(2,4) = 4
% 	            if(nele.eq.9) then
% 	               ixl(1,5) = 5
% 	               ixl(1,6) = 6
% 	               ixl(1,7) = 7
% 	               ixl(1,8) = 8
% 	               ixl(1,9) = 9
% 	               ixl(2,5) = 5
% 	               ixl(2,6) = 6
% 	               ixl(2,7) = 7
% 	               ixl(2,8) = 8
% 	               ixl(2,9) = 9
% 	            endif
% 	         elseif(edges(edge).eq.2) then
% 	            nA1 = 0
% 	            nA2 = minc
% 	            ixl(1,1) = 2
% 	            ixl(1,2) = 3
% 	            ixl(1,3) = 4
% 	            ixl(1,4) = 1
% 	            ixl(2,1) = 4
% 	            ixl(2,2) = 1
% 	            ixl(2,3) = 2
% 	            ixl(2,4) = 3
% 	            if(nele.eq.9) then
% 	               ixl(1,5) = 6
% 	               ixl(1,6) = 7
% 	               ixl(1,7) = 8
% 	               ixl(1,8) = 5
% 	               ixl(1,9) = 9
% 	               ixl(2,5) = 8
% 	               ixl(2,6) = 5
% 	               ixl(2,7) = 6
% 	               ixl(2,8) = 7
% 	               ixl(2,9) = 9
% 	            endif
% 	         elseif(edges(edge).eq.3) then
% 	            nA1 = cel+1
% 	            nA2 = -1
% 	            ixl(1,1) = 3
% 	            ixl(1,2) = 4
% 	            ixl(1,3) = 1
% 	            ixl(1,4) = 2
% 	            ixl(2,1) = 3
% 	            ixl(2,2) = 4
% 	            ixl(2,3) = 1
% 	            ixl(2,4) = 2
% 	            if(nele.eq.9) then
% 	               ixl(1,5) = 7
% 	               ixl(1,6) = 8
% 	               ixl(1,7) = 5
% 	               ixl(1,8) = 6
% 	               ixl(1,9) = 9
% 	               ixl(2,5) = 7
% 	               ixl(2,6) = 8
% 	               ixl(2,7) = 5
% 	               ixl(2,8) = 6
% 	               ixl(2,9) = 9
% 	            endif
% 	         else
% 	            nA1 = cel+1
% 	            nA2 = -minc
% 	            ixl(1,1) = 4
% 	            ixl(1,2) = 1
% 	            ixl(1,3) = 2
% 	            ixl(1,4) = 3
% 	            ixl(2,1) = 2
% 	            ixl(2,2) = 3
% 	            ixl(2,3) = 4
% 	            ixl(2,4) = 1
% 	            if(nele.eq.9) then
% 	               ixl(1,5) = 8
% 	               ixl(1,6) = 5
% 	               ixl(1,7) = 6
% 	               ixl(1,8) = 7
% 	               ixl(1,9) = 9
% 	               ixl(2,5) = 6
% 	               ixl(2,6) = 7
% 	               ixl(2,7) = 8
% 	               ixl(2,8) = 5
% 	               ixl(2,9) = 9
% 	            endif
% 	         endif
% 	         if(globalm) then
% 	            sslote = 0
% 	         else
% 	            sslote = ixl(2,sslot)
% 	         endif
% 
% 
% 	         do j = 1,nele
% 	            node = ixl(1,j)
%                   do i = 1,ndm
%                      xl(i,j) = xe(i,node)
%                   enddo
%                   do i = 1,ndf
%                      ul(i,j) = ue(i,node)
%                   enddo
%                enddo
% 
%       if(debug) then
%       dprt = .true.
%       if(dprt) then
%          write(24,*) 'tractyn: ',tractyn
% 	   write(24,*) 'sslote',sslote
%          write(24,*) 'local nodes'
%          do j = 1,nele
%             write(24,*) ixl(1,j)
%          enddo
%          write(24,'(A)') ''
%          write(24,*) 'elem coords'
%          do j = 1,nele
%             write(24,*) (xl(i,j),i=1,ndm)
%          enddo
%          write(24,*) 'elem nodal values'
%          do j = 1,nele
%             write(24,*) (ul(i,j),i=1,ndf)
%          enddo
%          write(24,'(A)') ''
%       endif
%       endif
% 
% %    -----------------------------------------------------
% %     Loop over cells in element
% %    -----------------------------------------------------
% 
% 	         cell = cel*(elem-1) + nA1
% 
% 	if(debug) then
% 	dprt = .true.
% 	if(dprt) then
% 	   write(24,*) 'cell',cell
% 	   write(24,*) 'nA1',nA1
% 	   write(24,*) 'nA2',nA2
% 	   write(24,'(A)') ''
% 	endif
% 	endif
% 
%                do 90 celli = 1,minc
% 
% %       Load cell coordinates, dof IDs
% 
%                   cell = cell + nA2
%                   call pzeroi(ixc,nen)
%    
%                   do j = 1,nele
%                      node = ixs(ixl(1,j),cell)
%                      ixc(j) = locid(node)
%                      do i = 1,ndm
%                         xc(i,j) = xs(i,node)
%                      enddo
%                   enddo
% 
% %       Set element-cell local coordinate map
%                   MR = 1.d0/minc
%                   BR = -1.d0+(2.d0*celli-1.d0)/minc
%                   MS = 0.d0
%                   BS = -1.d0
% 
%       if(debug) then
%       dprt = .true.
%       if(dprt) then
%          write(24,*) 'cell',cell
%          write(24,*) 'local nodes'
%          do j = 1,nele
%             write(24,*) ixc(j)
%          enddo
%          write(24,'(A)') ''
%          write(24,*) 'cell coords'
%          do j = 1,nele
%             write(24,*) (xc(i,j),i=1,ndm)
%          enddo
%          write(24,'(A)') ''
%          write(24,*) 'MR',MR
%          write(24,*) 'BR',BR
%          write(24,*) 'MS',MS
%          write(24,*) 'BS',BS
%          write(24,'(A)') ''
%       endif
%       endif
% 
% %      -----------------------------------------------------
% %       Compute cell boundary residual
% %      -----------------------------------------------------
%                    call getFbcell(xc,d(1,mat),xl,ul,ndfs,ndm,
%      &                      nele,nen,nummat,strong,MR,BR,MS,BS,
%      &                      sslote,tractyn,eF)
% 
%       if(debug) then
%       dprt = .true.
%       if(dprt) then
%          write(24,*) 'cell residual'
%          do j = 1,nele*ndfs
%             write(24,'(d12.6)') eF(j)
%          enddo
%          write(24,'(A)') ''
%       endif
%       endif
% 
% %       Assemble cell quantities into star quantities
% 
%                call LocToGlobDOF(ixc,ldoft,nele,nen,ndfs,edoft)
% 
%       if(debug) then
%       dprt = .true.
%       if(dprt) then
%          write(24,*) 'edoft'
%          do j = 1,nele*ndfs
%             write(24,*) edoft(j)
%          enddo
%          write(24,'(A)') ''
%       endif
%       endif
% 
% 	           call assembleF(nele,nen,ndfs,edoft,eF,Fd1)
% 
% %    -----------------------------------------------------
% %     End loop over cells
% %    -----------------------------------------------------
%  90            continue
%                endif
% 
% %    -----------------------------------------------------
% %     End loop over element edges
% %    -----------------------------------------------------
%  100           continue
% 
% 
% 	      endif
% 
%          endif

%%
%     Else Quadrilateral Element
            else

             if(globalm)
	            numE = 4;
                for i = 1:numE
	               edges(i) = i;
                end
	         else
	            numE = 2;
                  if(sslot==1)
	               edges(1) = 1;
	               edges(2) = 4;
                  elseif(sslot==2)
	               edges(1) = 2;
	               edges(2) = 1;
	              elseif(sslot==3)
	               edges(1) = 3;
	               edges(2) = 2;
	              else
	               edges(1) = 4;
	               edges(2) = 3;
                  end
             end

%    -----------------------------------------------------
%     "Loop" over element edges on interior of star
%    -----------------------------------------------------
              for edge = 1:numE

               if(sngammah)

                  nodeA = edges(edge)+1;
                if(nodeA>neleB)
			       nodeA = 1;
                end
	            nodeA = ixe(nodeA);
			                   
%        Check if nodeA is on the traction boundary of domain
                  [numA,numB,tractyn] = binsearch(nodeA,gammah,1,ngh);
	
                if(tractyn)

	            nodeB = edges(edge);
	            nodeB = ixe(nodeB);
			                   
%        Check if nodeA is on the traction boundary of domain
                  [numA,numB,tractyn] = binsearch(nodeB,gammah,1,ngh);

                end

	           else
	            tractyn = 0;
               end

             if(strong||tractyn)
             if(edges(edge)==1)
	            nA1 = 0;
	            nA2 = 1;
	            ixl(1,1) = 1;
	            ixl(1,2) = 2;
	            ixl(1,3) = 3;
	            ixl(1,4) = 4;
	            ixl(2,1) = 1;
	            ixl(2,2) = 2;
	            ixl(2,3) = 3;
	            ixl(2,4) = 4;
                if(nele==9)
	               ixl(1,5) = 5;
	               ixl(1,6) = 6;
	               ixl(1,7) = 7;
	               ixl(1,8) = 8;
	               ixl(1,9) = 9;
	               ixl(2,5) = 5;
	               ixl(2,6) = 6;
	               ixl(2,7) = 7;
	               ixl(2,8) = 8;
	               ixl(2,9) = 9;
                end
	         elseif(edges(edge)==2)
	            nA1 = 0;
	            nA2 = minc;
	            ixl(1,1) = 2;
	            ixl(1,2) = 3;
	            ixl(1,3) = 4;
	            ixl(1,4) = 1;
	            ixl(2,1) = 4;
	            ixl(2,2) = 1;
	            ixl(2,3) = 2;
	            ixl(2,4) = 3;
                if(nele==9)
	               ixl(1,5) = 6;
	               ixl(1,6) = 7;
	               ixl(1,7) = 8;
	               ixl(1,8) = 5;
	               ixl(1,9) = 9;
	               ixl(2,5) = 8;
	               ixl(2,6) = 5;
	               ixl(2,7) = 6;
	               ixl(2,8) = 7;
	               ixl(2,9) = 9;
                end
	         elseif(edges(edge)==3)
	            nA1 = cel+1;
	            nA2 = -1;
	            ixl(1,1) = 3;
	            ixl(1,2) = 4;
	            ixl(1,3) = 1;
	            ixl(1,4) = 2;
	            ixl(2,1) = 3;
	            ixl(2,2) = 4;
	            ixl(2,3) = 1;
	            ixl(2,4) = 2;
                if(nele==9)
	               ixl(1,5) = 7;
	               ixl(1,6) = 8;
	               ixl(1,7) = 5;
	               ixl(1,8) = 6;
	               ixl(1,9) = 9;
	               ixl(2,5) = 7;
	               ixl(2,6) = 8;
	               ixl(2,7) = 5;
	               ixl(2,8) = 6;
	               ixl(2,9) = 9;
                end
	         else
	            nA1 = cel+1;
	            nA2 = -minc;
	            ixl(1,1) = 4;
	            ixl(1,2) = 1;
	            ixl(1,3) = 2;
	            ixl(1,4) = 3;
	            ixl(2,1) = 2;
	            ixl(2,2) = 3;
	            ixl(2,3) = 4;
	            ixl(2,4) = 1;
                if(nele==9)
	               ixl(1,5) = 8;
	               ixl(1,6) = 5;
	               ixl(1,7) = 6;
	               ixl(1,8) = 7;
	               ixl(1,9) = 9;
	               ixl(2,5) = 6;
	               ixl(2,6) = 7;
	               ixl(2,7) = 8;
	               ixl(2,8) = 5;
	               ixl(2,9) = 9;
                end
             end
             if(globalm)
	            sslote = 0;
	         else
	            sslote = ixl(2,sslot);
             end


            for j = 1:nele
	            node = ixl(1,j);
                  for i = 1:ndm
                     xl(i,j) = xe(i,node);
                  end
                  for i = 1:ndf
                     ul(i,j) = ue(i,node);
                  end
            end

%    -----------------------------------------------------
%     Loop over cells in element
%    -----------------------------------------------------

	         cell = cel*(elem-1) + nA1;

               for celli = 1:minc

%       Load cell coordinates, dof IDs

                  cell = cell + nA2;
                  ixc = zeros(nen,1);
   
                  for j = 1:nele
                     node = ixs(ixl(1,j),cell);
                     ixc(j) = locid(node);
                     for i = 1:ndm
                        xc(i,j) = xs(i,node);
                     end
                  end

%       Set element-cell local coordinate map
                  MR = 1.d0/minc;
                  BR = -1.d0+(2.d0*celli-1.d0)/minc;
                  MS = 0.d0;
                  BS = -1.d0;

%      -----------------------------------------------------
%       Compute cell boundary residual
%      -----------------------------------------------------
                   stargetFb

%       Assemble cell quantities into star quantities

               edoft = LocToGlobDOF(ixc,ldoft',nele,ndfs);

	           assembleF

%    -----------------------------------------------------
%     End loop over cells
%    -----------------------------------------------------
               end
             end

%    -----------------------------------------------------
%     End loop over element edges
%    -----------------------------------------------------
              end


            end

         end

% %  -----------------------------------------------------
% %   End loop over elements in star problem
% %  -----------------------------------------------------
%  10   continue
% 
% % End assembly time
%       tt = etime(tary)
% 
% %	if(debug) then
% 
% 	   write(24,'(A,f17.7)') 'End Assembly time = ',tt
% 	   write(24,'(A)') ''
% 
% %	endif
%       
% %  -----------------------------------------------------
% %   Solve star problem for fine scale components
% %  -----------------------------------------------------
% 
% 	if(debug) then
%       dprt = .true.
%       if(dprt) then
%          write(24,*) 'Fd1'
% 	   do j = 1,neqp
%             write(24,*) j,Fd1(j)
%          enddo
% 	   write(24,*) 'Kdd11'
% 	   do j = 1,neqp
%             write(24,*) (Kdd11(i,j),i=1,neqp)
%          enddo
% 	   write(24,'(A)') ''
%       endif
%       endif
% 
% % Start Factorization time
%       tt = etime(tary)
% 
% %	if(debug) then
% 
% 	   write(24,'(A,f17.7)') 'Start Factor. time = ',tt
% 	   write(24,'(A)') ''
% 
% %	endif
% 
% 	call DSYSV(ulop,neqp,1,Kdd11,neqp,IPIV,Fd1,neqp,dsyswork,neqp,
% 	&           dsysinfo)
% 
%       if(debug) then
%       dprt = .true.
%       if(dprt) then
%          write(24,*) 'info',dsysinfo
%          write(24,'(A)') ''
% 	   write(24,*) 'fine components'
% 	   do j = 1,neqp
%             write(24,*) Fd1(j)
%          enddo
% 	   write(24,'(A)') ''
%       endif
%       endif
% 
% % End Factorization time
%       tt = etime(tary)
% 
% %	if(debug) then
% 
% 	   write(24,'(A,f17.7)') 'End Factor. time = ',tt
% 	   write(24,'(A)') ''
% 
% %	endif
% 
% %  -----------------------------------------------------
% %   Globalize star problem
% %  -----------------------------------------------------
% 
% %   Assemble fine scale component into global fine scale
% 	call globalize(Fd1,globid,ldoft,us)
% 
% 	return
% 
% 	end

%%
%  -----------------------------------------------------
%   End loop over elements in star problem
%  -----------------------------------------------------
      end
      
%  -----------------------------------------------------
%   Solve star problem for fine scale components
%  -----------------------------------------------------

	usi = KDD11\FD1;

%  -----------------------------------------------------
%   Globalize star problem
%  -----------------------------------------------------

%   Assemble fine scale component into global fine scale
	globalize1 %(usi,globid,ldoft,us)
