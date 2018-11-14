% c$Id:$
%       subroutine reshis(ix,nen1,numel,n1, n2)
% 
% c      * * F E A P * * A Finite Element Analysis Program
% 
% c....  Copyright (c) 1984-2008: Regents of the University of California
% c                               All rights reserved
% 
% c-----[--.----+----.----+----.-----------------------------------------]
% c     Modification log                                Date (dd/mm/year)
% c       Original version                                    01/11/2006
% c       1.                Conversion to Matlab by TJT       19/04/2013
% c       2.                Optimized, switching for -> :     30/04/2013
% c-----[--.----+----.----+----.-----------------------------------------]
% c      Purpose: Initialize t_n+1 history variables from final value
% c               of variables at t_n
% 
% c      Inputs:
% c         ix(nen1,*)  - Element connection/history pointer array
% c         nen1        - Dimension of ix array
% c         numel       - Number of elements in mesh
% c         n1          - Pointer in ix to t_n   data
% c         n2          - Pointer in ix to t_n+1 data
% 
% c      Outputs:
% c         none        - Output is retained in blank common
% c-----[--.----+----.----+----.-----------------------------------------]
%       implicit  none
% 
%       include  'hdata.h'
%       include  'pointer.h'
%       include  'comblk.h'
% 
%       integer   nen1,numel,n1,n2, n,nh,nhd, ix(nen1,numel)
% 
%       save

%     Move history variables from location 'n1' to 'n2'

      for elem = 1:numel
        nh1 = 1 + ixFEAP(hrn1,elem) - 1;
        nh2 = 1 + ixFEAP(hrn2,elem) - 1;
        if(nh2 ~= nh1) %then
          nhd = abs(nh2 - nh1);
%           for nh = 1:nhd
%             hrvec(nh+nh2) = hrvec(nh+nh1); %#ok<*SAGROW>
%           end % nh
          hrvec(nh2+1:nh2+nhd) = hrvec(nh1+1:nh1+nhd); %#ok<*SAGROW>
        end
      end % elem
      if hrstore == 1
        hrmat(:,step) = hrvec;
      end

% c     Move interface variables
% 
%       if(np(210).ne.0) then
%         call ireshis(mr(np(210)), n1,n2)
%       endif
