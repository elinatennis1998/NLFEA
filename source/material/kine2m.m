function  [fn1,detf,fi,df,fn,shpc] = kine2m(shpr,ul,uld,nel,sflag)
%       subroutine kine3m(shp,ul,f,fi,df,detf,ndf,nel,nen,lint)
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
% c      Purpose: Compute kinematic quantities for finite deformations
% 
% c      Inputs:
% c         shp(4,27,*)   - Reference configuration shape functions
% c         ul(ndf,nen,*) - Nodal displacements
% c         ndf           - Number dof/node
% c         nel           - Number nodes/element
% c         nen           - Maximum number nodes/element
% c         lint          - Number of quadrature points
% 
% c      Outputs:
% c         f(3,3,2,*)    - Deformation gradient
% c         fi(3,3,*)     - Inverse deformation gradient
% c         df(3,3,*)     - Incremental deformation gradient
% c         detf(2,*)     - Determinant of deformation gradient
% c-----[--.----+----.----+----.-----------------------------------------]
%       implicit  none
% 
%       integer   ndf,nel,nen,lint, i,j,k,l
%       real*8    detfi
%       real*8    shp(4,27,*),ul(ndf,nen,*)
%       real*8    df(3,3,*),f(3,3,2,*),fi(3,3,*),detf(2,*), cc(3)
% 
%       save

% f = eye(3,3);
% fi = zeros(3,3);
shpc = zeros(nel,2);
ad = zeros(2,2);
detf = zeros(2,1);

% c     Deformation gradient at t_n+1 : F_n+1 = I + GRAD u_n+1
% 
%       do l = 1,lint
%         do i = 1,3
%           do j = 1,3
%             f(j,i,1,l) = 0.0d0
%             df(j,i,l)  = 0.0d0
%           end do ! j
%         end do ! i
%         do i = 1,3
%           do j = 1,3
%             do k = 1,nel
%               f(i,j,1,l) = f(i,j,1,l) + ul(i,k,1)*shp(j,k,l)
%               df(i,j,l ) = df(i,j,l ) + ul(i,k,2)*shp(j,k,l)
%             end do ! k
%           end do ! j
%           f(i,i,1,l) = f(i,i,1,l) + 1.0d0
%         end do ! i
        fn1 = ul(1:2,:)*shpr+eye(2);
        df = uld(1:2,:)*shpr;

% c       Deformation gradient at t_n: F_n
% 
%         do i = 1,3
%           do j = 1,3
%             f(j,i,2,l) = f(j,i,1,l) - df(j,i,l)
%           end do ! j
%         end do ! i
        fn = fn1-df;

% c       Invert F
%     Compute adjoint to F

      ad(1,1) = fn1(2,2);
      ad(1,2) = -fn1(1,2);

      ad(2,1) = -fn1(2,1);
      ad(2,2) = fn1(1,1);

%     Compute determinant of F

      detf(1)  = fn1(1,1)*fn1(2,2)- fn1(1,2)*fn1(2,1);
      detf(2)  = fn(1,1)*fn(2,2)- fn(1,2)*fn(2,1);
%       detf = f(1,1:2)*ad(1:2,1);

      detfi = 1.d0/detf(1);
      
%     Compute F inverse

      fi(1,1) = ad(1,1)*detfi;
      fi(2,2) = ad(2,2)*detfi;
      fi(1,2) = ad(1,2)*detfi;
      fi(2,1) = ad(2,1)*detfi;

% c       Transform shape functions to current configuration
% 
%         do k = 1,nel
%           do i = 1,3
%             cc(i) = fi(1,i,l)*shp(1,k,l)
%      &            + fi(2,i,l)*shp(2,k,l)
%      &            + fi(3,i,l)*shp(3,k,l)
%           end do ! i
%           do i = 1,3
%             shp(i,k,l) = cc(i)
%           end do ! i
%         end do ! k
%       end do ! l
%         for k = 1:nel
%             for i = 1:3
%                 shpc(k,i) = fi(1,i)*shpr(k,1) + fi(2,i)*shpr(k,2) + fi(3,i)*shpr(k,3);
%             end
%         end %k
    if sflag == 1
        shpc = shpr*fi;
    end

      end