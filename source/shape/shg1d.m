function [shgd,shgs,xjac,bubble,shgt] = shg1d(xl,ndm,nel,shld,shls,nen,bf,der,bubble,vargin)
%
%       subroutine shp1d(s,xl,shp,ndm,nel,xjac)
% c      * * F E A P * * A Finite Element Analysis Program
% 
% c....  Copyright (c) 1984-2008: Regents of the University of California
% c                               All rights reserved
% 
% c-----[--.----+----.----+----.-----------------------------------------]
% c     Modification log                                Date (dd/mm/year)
% c       Original version                                    01/11/2006
% c-----[--.----+----.----+----.-----------------------------------------]
% c     Purpose: Compute shape functions, natural derivatives, and
% c              jacobian for 3-D beam at natural coordinate s.
% c              Linear (2 node) or quadratic (3 node) element.
% c                 o----------o                 o-----o-----o
% c                 1          2                 1     3     2
% 
% c     Inputs:
% c       s         : Natural coordinate
% c       xl(3,nel) : Nodal global coordinates
% c       ndm       : Coordinate dimension of mesh
% c       nel       : Number of nodes of element
% 
% c     Outputs:
% c       shp(2,nel): Shape functions and derivatives at s;
% c                     shp(1,1 to nel): derivatives of shape functions
% c                     shp(2,1 to nel): shape functions
% c       xjac      : Jacobian at s
% c-----[--.----+----.----+----.-----------------------------------------]

if nargin == 9
    shlt = zeros(nen,1);
else
    shlt = vargin;
end

      shgs = zeros(nen,1);
      shgt = zeros(nen,1);
% c       Convert local derivatives to global ones

%         xjac = 0.0d0
%         do j = 1,ndm
%           dxi = 0.0d0
%           do i = 1,nel
%             dxi = dxi + shp(1,i)*xl(j,i)
%           end do ! i
%           xjac = xjac + dxi*dxi
%         end do ! j
    sx = xl(:,1:nel)*shld';
    xjac = sum(sx.^2);
%         if(xjac.eq.0.d0) then
%            write(iow,3000) xjac
%            call plstop()
%         else
          xjac = sqrt(xjac);
%         endif

%         do i = 1,nel
%           shp(1,i) = shp(1,i)/xjac
%         end do ! i
    shgd = shld'/xjac;
    
        if bf == 1
%           temp = bubble(1);
%           bubble(1) = temp*xs(1,1) + bubble(2)*xs(2,1);
%           bubble(2) = temp*xs(1,2) + bubble(2)*xs(2,2);
% %           bubble(3) = bubble(3)
          bubble(2) = bubble(2)/xjac;
        end
        
      if der == 1 %%%%%%% NOT VERIFIED!!!!!!!!!!!
          
          t2 = 1/xjac^2;
          c1 = xl(:,1:nel)*shls';
          c1 = sum(c1.^2);
          t1 = -1/xjac*c1*t2;
          shgs = shld'*t1 + shls'*t2;
          
          h = xl(end)-xl(1);
          shgt = (2/h)^3*shlt;
          
      end