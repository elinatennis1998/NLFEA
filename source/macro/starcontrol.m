%**********************************************************************
%
%      subroutine starcontrol(ix,ixs,d,xs,us,u,cis,ecn,cse,cseb,
%     &                       csnb,locid,gammab,gammah,strong,globalm)
%
%...  Written by Tim Truster (Fall 2009)
%     No modifications when copied to NLFEA ver2
%...  Program 
%
%**********************************************************************
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
%       integer ix(nen1,numel),ixs(nen,numel*cel),ecn(numel*celn)
%       integer cis(2,numel+1),gammab(numnp),gammah(numnp-ndm)
%       integer cse(numel),cseb(numel),csnb(numnp)
%       integer locid(numel*celn)
%       real*8  d(250,*),xs(ndm,numel*celn),us(ndfs,numel*celn-numnp)
%       real*8  u(ndf,numnp,3)
% 	logical strong,globalm

%     Output Variables
%      real*8 

%     Local Variables
%       integer i,j,k,node,elem,nele
%       integer ldoft(ndfs,numnpp),globid(numnpp)
% 
%       integer nloop(4,2),edge,nodeA,nodeB,numA,numB
%       logical efound
% 
%       integer nume,ind,len,ind1,cnode,na,ni
% 	real    tt,etime,tary(2)
% 
% %	real*8 testmA(2,2),testvB(2),work(2)
% %	character ulop*1
% %	integer info,IPIV(2)
% 	
%       logical dprt
% 
% 	call pzeroi(ldoft,ndfs*numnpp)
% 
% %   Form global ID list for visualization purposes
%       do node = 1,numnps
%          nodeA = locid(node)
%          if(nodeA.gt.0) then
%             globid(nodeA) = node
%          endif
%       enddo
% 
%       if(debug) then
%       dprt = .true.
%       if(dprt) then
%          write(24,*) 'strong ',strong
%          write(24,'(A)') ''
%       endif
%       endif
% 
%       if(globalm) then
% 	   sngammah = .true.
% 	   neiqp = csnn*ndfs
% 	else
% 
% %   Check if snode is on the traction boundary of domain
%       call binsearch(snode,gammah,1,ngh,numA
%      *			               ,numB,sngammah)

    ldoft = zeros(ndfs,numnpp);
    globid = zeros(numnpp,1);

    nloop = zeros(4,2);

%   Form global ID list for visualization purposes
      for node = 1:numnps
         nodeA = locid(node);
         if(nodeA>0)
            globid(nodeA) = node;
         end
      end

%   Constrain boundary of star
      neiqp = csnn*ndfs;
      for k = 1:csnn
         for i = 1:ndfs
            ldoft(i,k) = 1;
         end
      end

% Constrain nodes on domain Dirichlet BC
    iel = MatTypeTable(2,1); %iel   = ie(nie-1,ma); same thing;
    nonlin = MatTypeTable(3,1);
    
    if(nonlin)
    
    if iel == 5 %NL_Elem5_2d
    
      for nodeA = 1:numnpp
        node = globid(nodeA);
        xnode = xs(:,node);
        if(iprob==1)
          if(abs(xnode(1) - 0.d0)<10e-8)
             if ldoft(1,nodeA) == 0
               ldoft(1,nodeA) = 1;
               neiqp = neiqp + 1;
             end
          end
          if(abs(xnode(2) - 10.d0)<10e-8)
             if ldoft(1,nodeA) == 0
               ldoft(1,nodeA) = 1;
               neiqp = neiqp + 1;
             end
          end
          if(abs(xnode(2) - 0.d0)<10e-8)
             if ldoft(2,nodeA) == 0
               ldoft(2,nodeA) = 1;
               neiqp = neiqp + 1;
             end
          end
        elseif(iprob==3)
          if(abs(xnode(1) - 0.d0)<10e-8)
             if ldoft(1,nodeA) == 0
               ldoft(1,nodeA) = 1;
               neiqp = neiqp + 1;
             end
          end
          if(abs(xnode(1) - 2.d0)<10e-8)
             if ldoft(1,nodeA) == 0
               ldoft(1,nodeA) = 1;
               neiqp = neiqp + 1;
             end
          end
          if(abs(xnode(2) - 0.d0)<10e-8)
             if ldoft(2,nodeA) == 0
               ldoft(2,nodeA) = 1;
               neiqp = neiqp + 1;
             end
          end
          if((abs(xnode(2) - 1.d0)<10e-8)&&((xnode(1) - 1.d0)<10e-8))
             if ldoft(2,nodeA) == 0
               ldoft(2,nodeA) = 1;
               neiqp = neiqp + 1;
             end
          end
        elseif(iprob==4)
          if(abs(xnode(1) - 16.d0)<10e-8)
             if ldoft(1,nodeA) == 0
               ldoft(1,nodeA) = 1;
               neiqp = neiqp + 1;
             end
          end
          if(abs(xnode(1) + 16.d0)<10e-8)
             if ldoft(1,nodeA) == 0
               ldoft(1,nodeA) = 1;
               neiqp = neiqp + 1;
             end
          end
        elseif(iprob==5)
          if((abs(xnode(2) - 1.75)<10e-8)&&((xnode(1) - 1.95/3)>-10e-8))
             if ldoft(2,nodeA) == 0
               ldoft(2,nodeA) = 1;
               neiqp = neiqp + 1;
             end
          end
        elseif(iprob==6)
          if((abs(xnode(2) - 0.0)<10e-8))
             if ldoft(1,nodeA) == 0
               ldoft(1,nodeA) = 1;
               neiqp = neiqp + 1;
             end
             if ldoft(2,nodeA) == 0
               ldoft(2,nodeA) = 1;
               neiqp = neiqp + 1;
             end
          end
        elseif(iprob==7)
          if((abs(xnode(1) - 0.0)<10e-8))
             if ldoft(1,nodeA) == 0
               ldoft(1,nodeA) = 1;
               neiqp = neiqp + 1;
             end
             if ldoft(2,nodeA) == 0
               ldoft(2,nodeA) = 1;
               neiqp = neiqp + 1;
             end
          end
        end
      end
      
    end %iel == 5
    
    else %nonlin == 0
    
    if iel == 3 %L_Elem3_2d

      for nodeA = 1:numnpp
        node = globid(nodeA);
        xnode = xs(:,node);
        if(iprob==8)
          if((abs(xnode(1) - 0.0)<10e-8))
             if ldoft(1,nodeA) == 0
               ldoft(1,nodeA) = 1;
               neiqp = neiqp + 1;
             end
             if ldoft(2,nodeA) == 0
               ldoft(2,nodeA) = 1;
               neiqp = neiqp + 1;
             end
          end
        end
      end
        
    end %iel == 3
    
    end %nonlin
      
      if(globalm)==1
          
          sngammah = 1;
          
      else

%   Check if snode is on the traction boundary of domain
      [numA,numB,sngammah] = binsearch(snode,gammah,1,ngh);
            
%%
%  -----------------------------------------------------
%   Apply boundary conditions and assign DOF
%  -----------------------------------------------------

%   Loop over elements in cseb
      for k = 1:csebn
         elem = cseb(k);
         ind = cis(1,elem);
         len = cis(1,elem+1) - ind;
         nele = cis(2,elem);

%     If T3 Element
         if(nele==3)
            nume = 3;
            for i=1:3
               nloop(i,1) = i;
               nloop(i,2) = i+1;
            end
            nloop(3,2) = 1;

%     Loop over edges of element
            for edge = 1:nume

               nodeA = NodesOnElement(nloop(edge,1),elem);
               nodeB = NodesOnElement(nloop(edge,2),elem);

%       Search for nodes A & B in csnb
%                efound = 0;

               [numA,numB,efound] = binsearch(nodeA,csnb,1,csnbn);
               if(efound)
                  [numA,numB,efound] = binsearch(nodeB,csnb,1,csnbn);
               end

%       If both nodes A & B are in csnb, search for edge in gammab
               if(efound)

%                   efound = 0;

                  [numA,numB,efound] = binsearch(nodeA,gammab,1,ngb);
                     if(efound)
                        [numA,numB,efound] = binsearch(nodeB,gammab,1,ngb);
                     end

%       If both nodes A & B are not in gammab, constrain the edge
                  if(efound==0)
                     ind1 = ind + nele + (edge-1)*(minc-1);
                     for j = 2:minc
                        node = ecn(ind1);
                        cnode	= locid(node);
                        for i = 1:ndfs
                          if ldoft(i,cnode) == 0
                            ldoft(i,cnode) = 1;
                            neiqp = neiqp + 1;
                          end
                        end
                        ind1 = ind1 + 1;
                     end
                  end

               end

%     End loop over edges
            end


% %     If Q4 Element
%          elseif(nele.eq.4) then
%             nume = 4
%             do i=1,4
%                nloop(i,1) = i
%                nloop(i,2) = i+1
%             enddo
%             nloop(4,2) = 1
% 
% %     Loop over edges of element
%             do edge = 1,nume
% 
%                nodeA = ix(nloop(edge,1),elem)
%                nodeB = ix(nloop(edge,2),elem)
% 
%       if(debug) then
%       dprt = .false.
%       if(dprt) then
%          write(24,*) 'nodeA',nodeA
%          write(24,*) 'nodeB',nodeB
%          write(24,'(A)') ''
%       endif
%       endif
% 
% %       Search for nodes A & B in csnb
%                efound = .false.
% 
%                call binsearch(nodeA,csnb,1,csnbn,numA
%      *			               ,numB,efound)
%                if(efound) then
%                   call binsearch(nodeB,csnb,1,csnbn,numA
%      *			               ,numB,efound)
%                endif
% 
%       if(debug) then
%       dprt = .false.
%       if(dprt) then
%          write(24,*) 'efound in csnb ',efound
%          write(24,'(A)') ''
%       endif
%       endif
% 
% %       If both nodes A & B are in csnb, search for edge in gammab
%                if(efound) then
% 
%                   efound = .false.
% 
%                   call binsearch(nodeA,gammab,1,ngb,numA
%      *		   	               ,numB,efound)
%                   if(efound) then
% 	                 call binsearch(nodeB,gammab,1,ngb,numA
%      *			               ,numB,efound)
%                   endif
% 
%       if(debug) then
%       dprt = .false.
%       if(dprt) then
%          write(24,*) 'efound in gammab ',efound
%          write(24,'(A)') ''
%       endif
%       endif
% 
% %       If both nodes A & B are not in gammab, constrain the edge
%                   if(.not.efound) then
%                      ind1 = ind + nele + (edge-1)*(minc-1)
%                      do j = 2,minc
%                         node = ecn(ind1)
%                         cnode	= locid(node)
%                         do i = 1,ndfs
%                            neiqp = neiqp + 1
%                            ldoft(i,cnode) = 1
%                         enddo
%                         ind1 = ind1 + 1
%                      enddo
%                   endif
% 
%                endif
% 
% %     End loop over edges
%             enddo

%%
%     If Q4 Element
         elseif(nele==4)
            nume = 4;
            for i=1:4
               nloop(i,1) = i;
               nloop(i,2) = i+1;
            end
            nloop(4,2) = 1;

%     Loop over edges of element
            for edge = 1:nume

               nodeA = NodesOnElement(nloop(edge,1),elem);
               nodeB = NodesOnElement(nloop(edge,2),elem);

%       Search for nodes A & B in csnb
%                efound = 0;

               [numA,numB,efound] = binsearch(nodeA,csnb,1,csnbn);
               if(efound)
                  [numA,numB,efound] = binsearch(nodeB,csnb,1,csnbn);
               end

%       If both nodes A & B are in csnb, search for edge in gammab
               if(efound)

%                   efound = 0;

                  [numA,numB,efound] = binsearch(nodeA,gammab,1,ngb);
                  if(efound)
	                 [numA,numB,efound] = binsearch(nodeB,gammab,1,ngb);
                  end

%       If both nodes A & B are not in gammab, constrain the edge
                  if(efound==0)
                     ind1 = ind + nele + (edge-1)*(minc-1);
                     for j = 2:minc
                        node = ecn(ind1);
                        cnode	= locid(node);
                        for i = 1:ndfs
                          if ldoft(i,cnode) == 0
                            ldoft(i,cnode) = 1;
                            neiqp = neiqp + 1;
                          end
                        end
                        ind1 = ind1 + 1;
                     end
                  end

               end

%     End loop over edges
            end

% %     If T6 Element
%          elseif(nele.eq.6) then
%             nume = 3
%             do i=1,3
%                nloop(i,1) = i
%                nloop(i,2) = i+1
%             enddo
%             nloop(3,2) = 1
% 
% %     Loop over edges of element
%             do edge = 1,nume
% 
%                nodeA = NodesOnElement(nloop(edge,1),elem)
%                nodeB = NodesOnElement(nloop(edge,2),elem)
% 
%       if(debug) then
%       dprt = .true.
%       if(dprt) then
%          write(24,*) 'nodeA',nodeA
%          write(24,*) 'nodeB',nodeB
%          write(24,'(A)') ''
%       endif
%       endif
% 
% %       Search for nodes A & B in csnb
%                efound = .false.
% 
%                call binsearch(nodeA,csnb,1,csnbn,numA
%      *			               ,numB,efound)
%                if(efound) then
%                   call binsearch(nodeB,csnb,1,csnbn,numA
%      *			               ,numB,efound)
%                endif
% 
%       if(debug) then
%       dprt = .true.
%       if(dprt) then
%          write(24,*) 'efound in csnb ',efound
%          write(24,'(A)') ''
%       endif
%       endif
% 
% %       If both nodes A & B are in csnb, search for edge in gammab
%                if(efound) then
% 
%                   efound = .false.
% 
%                   call binsearch(nodeA,gammab,1,ngb,numA
%      *		   	               ,numB,efound)
%                   if(efound) then
% 	                 call binsearch(nodeB,gammab,1,ngb,numA
%      *			               ,numB,efound)
% 	                 if(efound) then
% 	                    call binsearch(NodesOnElement(nloop(edge,1)+4,elem),
%      *			               gammab,1,ngb,numA,numB,efound)
%                        endif
%                   endif
% 
%       if(debug) then
%       dprt = .true.
%       if(dprt) then
%          write(24,*) 'efound in gammab ',efound
%          write(24,'(A)') ''
%       endif
%       endif
% 
% %       If both nodes A & B are not in gammab, constrain the edge
%                   if(.not.efound) then
%                      ind1 = ind + nele + (edge-1)*(2*minc-2)
%                      do j = 2,2*minc
% 	                  if(j.ne.minc+1) then
%                            node = ecn(ind1)
%                            cnode	= locid(node)
%                            do i = 1,ndfs
%                               neiqp = neiqp + 1
%                               ldoft(i,cnode) = 1
%                            enddo
%                            ind1 = ind1 + 1
% 	                  endif
%                      enddo
%                   endif
% 
%                endif
% 
% %     End loop over edges
%             enddo

%%
%     If T6 Element
         elseif(nele==6)
            nume = 3;
            for i=1:3
               nloop(i,1) = i;
               nloop(i,2) = i+1;
            end
            nloop(3,2) = 1;

%     Loop over edges of element
            for edge = 1:nume

               nodeA = NodesOnElement(nloop(edge,1),elem);
               nodeB = NodesOnElement(nloop(edge,2),elem);

%       Search for nodes A & B in csnb
%                efound = 0;

               [numA,numB,efound] = binsearch(nodeA,csnb,1,csnbn);
               if(efound)
                  [numA,numB,efound] = binsearch(nodeB,csnb,1,csnbn);
               end

%       If both nodes A & B are in csnb, search for edge in gammab
               if(efound)

%                   efound = 0;

                  [numA,numB,efound] = binsearch(nodeA,gammab,1,ngb);
                  if(efound)
	                 [numA,numB,efound] = binsearch(nodeB,gammab,1,ngb);
                     if(efound)
	                    [numA,numB,efound] = binsearch(NodesOnElement(nloop(edge,1)+4,elem),gammab,1,ngb);
                     end
                  end

%       If both nodes A & B are not in gammab, constrain the edge
                  if(efound==0)
                     ind1 = ind + nele + (edge-1)*(2*minc-2);
                     for j = 2:2*minc
                       if(j~=minc+1)
                           node = ecn(ind1);
                           cnode	= locid(node);
                           for i = 1:ndfs
                          if ldoft(i,cnode) == 0
                            ldoft(i,cnode) = 1;
                            neiqp = neiqp + 1;
                          end
                           end
                           ind1 = ind1 + 1;
                       end
                     end
                  end

               end

%     End loop over edges
            end
            

% %     Else Q9 Element
%          else
%             nume = 4
%             do i=1,4
%                nloop(i,1) = i
%                nloop(i,2) = i+1
%             enddo
%             nloop(4,2) = 1
% 
% %     Loop over edges of element
%             do edge = 1,nume
% 
%                nodeA = NodesOnElement(nloop(edge,1),elem)
%                nodeB = NodesOnElement(nloop(edge,2),elem)
% 
%       if(debug) then
%       dprt = .true.
%       if(dprt) then
%          write(24,*) 'nodeA',nodeA
%          write(24,*) 'nodeB',nodeB
%          write(24,'(A)') ''
%       endif
%       endif
% 
% %       Search for nodes A & B in csnb
%                efound = .false.
% 
%                call binsearch(nodeA,csnb,1,csnbn,numA
%      *			               ,numB,efound)
%                if(efound) then
%                   call binsearch(nodeB,csnb,1,csnbn,numA
%      *			               ,numB,efound)
%                endif
% 
%       if(debug) then
%       dprt = .true.
%       if(dprt) then
%          write(24,*) 'efound in csnb ',efound
%          write(24,'(A)') ''
%       endif
%       endif
% 
% %       If both nodes A & B are in csnb, search for edge in gammab
%                if(efound) then
% 
%                   efound = .false.
% 
%                   call binsearch(nodeA,gammab,1,ngb,numA
%      *		   	               ,numB,efound)
%                   if(efound) then
% 	                 call binsearch(nodeB,gammab,1,ngb,numA
%      *			               ,numB,efound)
% 	                 if(efound) then
% 	                    call binsearch(NodesOnElement(nloop(edge,1)+4,elem),
%      *			               gammab,1,ngb,numA,numB,efound)
%                        endif
%                   endif
% 
%       if(debug) then
%       dprt = .true.
%       if(dprt) then
%          write(24,*) 'efound in gammab ',efound
%          write(24,'(A)') ''
%       endif
%       endif
% 
% %       If both nodes A & B are not in gammab, constrain the edge
%                   if(.not.efound) then
%                      ind1 = ind + nele + (edge-1)*(2*minc-2)
%                      do j = 2,2*minc
% 	                  if(j.ne.minc+1) then
%                            node = ecn(ind1)
%                            cnode	= locid(node)
%                            do i = 1,ndfs
%                               neiqp = neiqp + 1
%                               ldoft(i,cnode) = 1
%                            enddo
%                            ind1 = ind1 + 1
% 	                  endif
%                      enddo
%                   endif
% 
%                endif
% 
% %     End loop over edges
%             enddo
% 
%          endif
% 
% 
% %   End loop over elements in cseb
%       enddo
% 
% 	endif %globalm

%%
%     Else Q9 Element
         else
            nume = 4;
            for i=1:4
               nloop(i,1) = i;
               nloop(i,2) = i+1;
            end
            nloop(4,2) = 1;

%     Loop over edges of element
            for edge = 1:nume

               nodeA = NodesOnElement(nloop(edge,1),elem);
               nodeB = NodesOnElement(nloop(edge,2),elem);

%       Search for nodes A & B in csnb
%                efound = 0;

               [numA,numB,efound] = binsearch(nodeA,csnb,1,csnbn);
               if(efound)
                  [numA,numB,efound] = binsearch(nodeB,csnb,1,csnbn);
               end

%       If both nodes A & B are in csnb, search for edge in gammab
               if(efound)

%                   efound = 0;

                  [numA,numB,efound] = binsearch(nodeA,gammab,1,ngb);
                  if(efound)
	                 [numA,numB,efound] = binsearch(nodeB,gammab,1,ngb);
                     if(efound)
	                    [numA,numB,efound] = binsearch(NodesOnElement(nloop(edge,1)+4,elem),gammab,1,ngb);
                     end
                  end

%       If both nodes A & B are not in gammab, constrain the edge
                  if(efound==0)
                     ind1 = ind + nele + (edge-1)*(2*minc-2);
                     for j = 2:2*minc
                       if(j~=minc+1)
                           node = ecn(ind1);
                           cnode	= locid(node);
                           for i = 1:ndfs
                          if ldoft(i,cnode) == 0
                            ldoft(i,cnode) = 1;
                            neiqp = neiqp + 1;
                          end
                           end
                           ind1 = ind1 + 1;
                       end
                     end
                  end

               end

%     End loop over edges
            end

         end


%   End loop over elements in cseb
      end

      end %globalm

% %   Create dof table for cells
%       neqp = ndfs*numnpp - neiqp
%       na = 0
%       ni = neqp
% 
%       if(debug) then
%       dprt = .true.
%       if(dprt) then
%          write(24,*) 'neiqp',neiqp
%          write(24,*) 'neqp',neqp
% 	   write(24,*) 'sngammah ',sngammah
%          write(24,'(A)') ''
%       endif
%       endif
% 
% %   Assign boundary conditions
%       do k = 1,csnn
%          do i = 1,ndfs
%             ni = ni + 1
%             ldoft(i,k) = ni
%          enddo
%       enddo
% 
%       do k = csnn+1,numnpp
%          do i = 1,ndfs
%             if(ldoft(i,k).eq.1) then
%                ni = ni + 1
%                ldoft(i,k) = ni
%             else
%                na = na + 1
%                ldoft(i,k) = na
%             endif
%          enddo
%       enddo
% 
%       if(debug) then
%       dprt = .true.
%       if(dprt) then
%          write(24,*) 'ldoft'
%          write(24,'(A)') 'GNode  dof1  dof2  dof3'
%          do j = 1,numnpp
%             write(24,1001) globid(j),(ldoft(i,j),i=1,ndfs)
%          enddo
%          write(24,'(A)') ''
%  1001 format (I4,1x,3(I5))
%       endif
%       endif
% 
% % Start solve time
%       tt = etime(tary)
% 
% %	if(debug) then
% 
% 	   write(24,'(A,f17.7)') 'Start Solve time = ',tt
% 	   write(24,'(A)') ''
% 
% %	endif
% 
% 	call starsolve(NodesOnElement,ixs,d,xs,us,u,cis,cse,
%      &                      locid,globid,ldoft,
%      &                      gammah,strong,globalm)
% 
% % End solve time
%       tt = etime(tary)
% 
% %	if(debug) then
% 
% 	   write(24,'(A,f17.7)') 'End Solve time = ',tt
% 	   write(24,'(A)') ''
% 
% %	endif
% 
% 	return
% 
% 	end

%%
%   Create dof table for cells
      neqp = ndfs*numnpp - neiqp;
      na = 0;
      ni = neqp;

%   Assign boundary conditions
      for k = 1:csnn
         for i = 1:ndfs
            ni = ni + 1;
            ldoft(i,k) = ni;
         end
      end

      for k = csnn+1:numnpp
         for i = 1:ndfs
            if(ldoft(i,k)==1)
               ni = ni + 1;
               ldoft(i,k) = ni;
            else
               na = na + 1;
               ldoft(i,k) = na;
            end
         end
      end

% 	starsolve(NodesOnElement,ixs,d,xs,us,u,cis,cse,locid,globid,ldoft,gammah,strong,globalm)
   starsolve




