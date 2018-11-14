% c$Id:$
%       subroutine pbodyf(ix,ndf,nen1,numel,prt,prth)
% 
% c      * * F E A P * * A Finite Element Analysis Program
% 
% c....  Copyright (c) 1984-2008: Regents of the University of California
% c                               All rights reserved
% 
% c-----[--.----+----.----+----.-----------------------------------------]
% c     Modification log                                Date (dd/mm/year)
% c       Original version                                    01/11/2006
% c       1. Increase allowable dof to 30                     05/07/2007
% c-----[--.----+----.----+----.-----------------------------------------]
% c      Purpose: Compute body forces from elements and add to load
% c               vector f
% 
% c      Inputs:
% c         ix(nen1,*) - Element nodal connection list
% c         ndf        - Number dof/node
% c         nen1       - Dimension of ix array
% c         numel      - Number of elements in mesh
% c         prt        - Print results if true
% c         prth       - Print title/header information if true
% 
% c      Outputs:
% c         f(*)       - Force vector with body forces added (via pointer)
% c-----[--.----+----.----+----.-----------------------------------------]
%       implicit  none
% 
%       include  'elbody.h'
%       include  'iofile.h'
%       include  'pointer.h'
% 
%       logical   allreg
%       character text*15
%       logical   prt,prth, pcomp, tinput, errck
%       integer   i,j,n,ndf,nen1,numel
%       integer   ix(nen1,*)
%       real*8    td(31)
% 
%       save

% c     Read data from file

    for bf = 1:numBF1
        
      text = 'mate';
      td = BodyForce1(bf,:);
%       j     = min(15,ndf+1)
%       errck = tinput(text,1,td,ndf+1)
%       if(ndf.gt.14) then
%         j     = ndf - 14
%         errck = tinput(text,0,td(16),j)
%       endif

%     Set body force intensities

%       do i = 1,ndf
        bodyf = td(2:end);
%       end do ! i

%     Check for process type

      if    (strcmp(text,'mate')) 
        allreg = 0;
        j      = 0;
      elseif(strcmp(text,'regi')) 
        allreg = 0;
        j      = 1;
      elseif(strcmp(text(1:3),'all')) 
        allreg = 0;
        j      = 1;
      elseif(strcmp(text,'elem'))
        j      = 2;
      end

%     Compute and assemble body loadings

      if( j <= 1 ) 

%       Do regions or materials

        % Currently assembles only into Fc1; could use minus signs to have
        % it assemble into Fc1np;
        for elemn = 1:length(reglist1)
    
          elem = reglist1(elemn);
          if(allreg || (ixFEAP(7-j,elem) == td(1)) ) 
            isw = 15;
            FormFEnn1
          end
        end % n

%     Single element case

      elseif(j == 2) 

        elem = td(1); 
        isw = 15;
        FormFEnn1

      end

    end
