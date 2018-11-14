% c$Id:$
%       subroutine umacr6(lct,ctl)
% 
% c      * * F E A P * * A Finite Element Analysis Program
% 
% c....  Copyright (c) 1984-2014: Regents of the University of California
% c                               All rights reserved
% 
% c-----[--.----+----.----+----.-----------------------------------------]
% c     Modification log                                Date (dd/mm/year)
% c       Original version                                    01/11/2006
% c-----[--.----+----.----+----.-----------------------------------------]
% c      Purpose:  Average coordinates for mid-edge nodes of quadratic
% c                order elements.
% 
% c      Inputs:
% c         lct        - Command character parameters
% c                      'xave'rage
% c         ctl(3)     - Command numerical parameters
% 
% c      Outputs:
% c         hr(np(44)) - x(ndm,numnp): averaged coordinates
% c-----[--.----+----.----+----.-----------------------------------------]
%       implicit  none
% 
%       include  'iofile.h'
%       include  'umac1.h'
%       include  'pointer.h'
%       include  'comblk.h'
% 
%       logical   pcomp
%       character lct*15
%       real*8    ctl(3)
% 
%       save
% 
% c     Set command word
% 
%       if(pcomp(uct,'mac6',4)) then      ! Usual    form
%         uct = 'xave'                    ! Average mid-side coordinates
%       elseif(urest.eq.1) then           ! Read  restart data
% 
%       elseif(urest.eq.2) then           ! Write restart data
% 
%       else                              ! Perform user operation
%         call uxaver(mr(np(33)),hr(np(43)))
%       endif
% 
%       end
% 
function x = uxaver(ix, x,numel,ndm,nen)
% 
%       implicit   none
% 
%       include   'cdata.h'
%       include   'sdata.h'
% 
%       integer    ix(nen1,*), ne,nel, i,j, n1,n2,n3,n4,n5
%       integer    tri(3,3),qud(3,4),tet(3,6),brk(3,12),fac(5,6)
%       real*8     x(ndm,*)

tri = [1,2,4; 2,3,5; 3,1,6]';
qud = [1,2,5; 2,3,6; 3,4,7; 4,1,8]';
tet = [1,2,5; 2,3,6; 3,1,7; 1,4,8; 2,4,9; 3,4,10]';
brk = [1,2, 9; 2,3,10; 3,4,11; 4,1,12; 5,6,13; 6,7,14;
       7,8,15; 8,5,16; 1,5,17; 2,6,18; 3,7,19; 4,8,20]';
fac = [1,4,8,5,21; 2,3,7,6,22; 4,3,7,8,23;
       1,2,6,5,23; 1,2,3,4,24; 5,6,7,8,27]';

      for ne = 1:numel
        nel = 0;
        for j = 1:nen
          if(ix(j,ne) > 0) %then
            nel = j;
          end %if
        end %do ! j

        if(ndm == 2) %then
          if(nel == 6) %then  ! 6-node triangle
            for j = 1:3
              n1 = tri(1,j);
              n2 = tri(2,j);
              n3 = tri(3,j);
              for i = 1:ndm
                x(i,ix(n3,ne)) = 0.5d0*(x(i,ix(n1,ne)) + x(i,ix(n2,ne)));
              end %do ! i
            end %do ! j
          elseif(nel == 8  ||  nel == 9) %then  ! 8- or 9-node quadrilateral
            for j = 1:4
              n1 = qud(1,j);
              n2 = qud(2,j);
              n3 = qud(3,j);
              for i = 1:ndm
                x(i,ix(n3,ne)) = 0.5d0*(x(i,ix(n1,ne)) + x(i,ix(n2,ne)));
              end %do ! i
            end %do ! j
            if(nel == 9) %then
              for i = 1:ndm
                x(i,ix(9,ne)) = 0.25d0*(x(i,ix(1,ne)) + x(i,ix(2,ne)) ...
                                      + x(i,ix(3,ne)) + x(i,ix(4,ne)));
              end %do ! i
            end%if
          end%if

        elseif(ndm == 3) %then
          if(nel == 10) %then  ! 10-node tetrahedron
            for j = 1:6
              n1 = tet(1,j);
              n2 = tet(2,j);
              n3 = tet(3,j);
              for i = 1:ndm
                x(i,ix(n3,ne)) = 0.5d0*(x(i,ix(n1,ne)) + x(i,ix(n2,ne)));
              end %do ! i
            end %do ! j
          elseif(nel == 20  ||  nel ==  27) %then  ! 20- or 27-node brick
            for j = 1:12
              n1 = brk(1,j);
              n2 = brk(2,j);
              n3 = brk(3,j);
              for i = 1:ndm
                x(i,ix(n3,ne)) = 0.5d0*(x(i,ix(n1,ne)) + x(i,ix(n2,ne)));
              end %do ! i
            end %do ! j
            if(nel == 27) %then ! mid-face and central nodes
              for j = 1:4
                n1 = fac(1,j);
                n2 = fac(2,j);
                n3 = fac(3,j);
                n4 = fac(4,j);
                n5 = fac(5,j);
                for i = 1:ndm
                  x(i,ix(n5,ne)) = 0.25d0*(x(i,ix(n1,ne)) ...
                                         + x(i,ix(n2,ne)) ...
                                         + x(i,ix(n3,ne)) ...
                                         + x(i,ix(n4,ne)));
                end %do ! i
              end %do ! j
              for i = 1:ndm
                x(i,ix(27,ne)) = 0.125d0*(x(i,ix(1,ne)) + x(i,ix(2,ne)) ...
                                        + x(i,ix(3,ne)) + x(i,ix(4,ne)) ...
                                        + x(i,ix(5,ne)) + x(i,ix(6,ne)) ...
                                        + x(i,ix(7,ne)) + x(i,ix(8,ne)));
              end %do ! i
            end%if
          end%if
        end%if
      end %do ! ne

      end
