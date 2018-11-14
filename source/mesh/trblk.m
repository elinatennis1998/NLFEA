function [x,NodesOnElement,RegionOnElement] = trblk(nr,xl,ixl,x,NodesOnElement,RegionOnElement,ndm,nod1,nuel1,ma,ntyp)
% 
%       subroutine trblk(nr,xl,ixl,x,ix,ndm,nod1,nuel1,nen1,ma,ntyp,
%      &                 ctype,prt,prth)
% 
% c      * * F E A P * * A Finite Element Analysis Program
% 
% c....  Copyright (c) 1984-2008: Regents of the University of California
% c                               All rights reserved
% 
% c-----[--.----+----.----+----.-----------------------------------------]
% c     Modification log                                Date (dd/mm/year)
% c       Original version                                    01/11/2006
% c-----[--.----+----.----+----.-----------------------------------------]
% c      Purpose: Generate a triangular block of 3-node triangular elements
% 
% c      Inputs:
% c         nr        - Number elements in 1-local coordinate dir.
% c         xl(ndm,*) - Block nodal coordinate array
% c         ixl(*)    - Block nodal connection list
% c         ndm       - Spatial dimension of mesh
% c         nod1      - Initial node number for block
% c         nuel1     - Initial element number for block
% c         nen1      - Dimension of ix array
% c         ma        - Material set number for block
% c         ntyp      - Element type: 1 = 3-node; 2 = 6-node
% c         ctype     - Input coordinate types
% c         prt       - Output generated data if true
% c         prth      - Output title/header data if true
% 
% c      Outputs:
% c         x(ndm,*)  - Nodal coordinates for block
% c         ix(*)     - Element nodal connection list for block
% c-----[--.----+----.----+----.-----------------------------------------]
%       implicit  none
% 
%       include  'cdata.h'
%       include  'cdat2.h'
%       include  'iofile.h'
%       include  'trdata.h'
% 
%       logical   prt,prth, pcomp
%       character xh*6, ctype*15
%       integer   ni,nn,mct,k,i,j,i1,j1,nei,n1,n2,n3
%       integer   nr,ndm,nod1,nuel1,nen1,ma,ntyp, ixl(1),ix(nen1,1)
%       real*8    dl, rr,sn,cn, xl(3,1),x(ndm,1),xx(3),tshp(6),el(3)
% 
%       save
% 
%       data      xh/' coord'/

      ni  = nr;
      xx = zeros(3,1);
      tshp = zeros(6,1);
      el = zeros(3,1);

      for i=1:3
        j = mod(i,3) + 1;
        if((ixl(i+3)) ~= 0) %then
          for i1=1:ndm
            xl(i1,i+3) = xl(i1,i+3) - 0.5d0*(xl(i1,i) + xl(i1,j));
          end % i1
        end
        xx(i) = 0.0d0;
      end % i

%     Generate nodes

      nn = nod1;
      dl  = 1.0d0/ni;
      for i=0:ni
        el(3) = i*dl;

        for j=0:ni-i
          el(2) = j*dl;
          el(1) = 1.0d0 - el(3) - el(2);

%         Form shape functions

          for i1=1:3
            j1 = mod(i1,3) + 1;
            tshp(i1)   = el(i1);
            tshp(i1+3) = 4.0d0*el(i1)*el(j1);
          end % i1

%           if(nn.gt.numnp) then
%             write(*,*) ' Trying to generate node',nn
%             call plstop()
%           end
          for i1 = 1:ndm
            xx(i1) = 0.0d0;
            for k=1:6
              xx(i1) = xx(i1) + tshp(k)*xl(i1,k);
            end % k
          end % i1
%           if(pcomp(ctype,'pola',4) .or. pcomp(ctype,'cyli',4)) then
%             call pdegree(xx(2), sn,cn)
%             rr    = xx(1)
%             xx(1) = x0(1) + rr*cn
%             xx(2) = x0(2) + rr*sn
%           end
%           for k = 1:ndm
%             x(k,nn) = xr(k)+tr(k,1)*xx(1)+tr(k,2)*xx(2)+tr(k,3)*xx(3);
%           end % k
        for k = 1:ndm
            x(k,nn) = xx(k);
        end
          nn = nn + 1;
        end % j
      end % i

%     Generate elements

      nei = nuel1;
      n1 = nod1;
      for i=1:ntyp:ni
        for j=1:ntyp:ni-i+1
          if(ntyp == 1) %then
            n2           = n1 + ni - i + 2;
            NodesOnElement(1,nei)    = n1 + j - 1; %#ok<*AGROW>
            NodesOnElement(2,nei)    = n1 + j;
            NodesOnElement(3,nei)    = n2 + j - 1;
            RegionOnElement(nei) = ma;
            if(j < (ni-i+1)) %then
              nei = nei + 1;
              NodesOnElement(1,nei)    = n2 + j - 1;
              NodesOnElement(2,nei)    = n1 + j;
              NodesOnElement(3,nei)    = n2 + j;
              RegionOnElement(nei) = ma;
            end
          elseif(ntyp == 2) %then
            n3           = n1 + ni - i + 2;
            n2           = n3 + ni - i + 1;
            NodesOnElement(1,nei)    = n1 + j - 1;
            NodesOnElement(4,nei)    = n1 + j;
            NodesOnElement(2,nei)    = n1 + j + 1;
            NodesOnElement(6,nei)    = n3 + j - 1;
            NodesOnElement(5,nei)    = n3 + j;
            NodesOnElement(3,nei)    = n2 + j - 1;
            RegionOnElement(nei) = ma;
            if(j < (ni-i)) %then
              nei          = nei + 1;
              NodesOnElement(1,nei)    = n2 + j - 1;
              NodesOnElement(4,nei)    = n3 + j;
              NodesOnElement(2,nei)    = n1 + j + 1;
              NodesOnElement(5,nei)    = n3 + j + 1;
              NodesOnElement(6,nei)    = n2 + j;
              NodesOnElement(3,nei)    = n2 + j + 1;
              RegionOnElement(nei) = ma;
            end
          end
          nei = nei + 1;
        end % j
        n1 = n2;
      end % i

      end
