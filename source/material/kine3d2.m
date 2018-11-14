function  [f,detf,fi,shpc,shcp] = kine3d2(shpr,ul,nel,sflag,shpp,nelP)
%
%      * * F E A P * * A Finite Element Analysis Program
%
%....  Copyright (c) 1984-2002: Regents of the University of California
%                               All rights reserved
%
%....  Modified and checked by Tim Truster, UIUC,  08/20/2009
%
%-----[--.----+----.----+----.-----------------------------------------]
%      Purpose: Compute kinematic quantities for finite deformations
%
%      Inputs:
%         shpr(3,16,*)  - Reference configuration shape functions
%         xl(ndm,*)     - Nodal reference coordinates
%         ul(ndf,nen,*) - Nodal displacements
%         ndm           - Number mesh dimensions
%         ndf           - Number dof/node
%         nel           - Number nodes/element
%         nen           - Maximum number nodes/element
%         lint          - Number of quadrature points
%
%      Outputs:
%         f(3,3,2,*)    - Deformation gradient
%         fi(3,3,*)     - Inverse deformation gradient
%         df(3,3,*)     - Incremental deformation gradient
%         detf(2,*)     - Determinant of deformation gradient
%         shpc(3,16,*)  - Current configuration shape functions
%-----[--.----+----.----+----.-----------------------------------------]
%       implicit  none
% 
%       include  'pmod2d.h'
% 
%       integer   ndm,ndf,nel,nen,lint, i,j,k,l
%       real*8    detfi,temp,xx1
%       real*8    shp(3,16,*),xl(ndm,*),ul(ndf,nen,*),u0(2,*)
%       real*8    df(3,3,*),f(3,3,2,*),fi(3,3,*),detf(2,*)

f = eye(3,3);
fi = zeros(3,3);
shpc = zeros(nel,3);
shcp = zeros(nelP,3);
ad = zeros(3,3);

%     Deformation gradient at t_n+1 : F_n+1 = I + GRAD u_n+1

%       for l = 1,lint
        for i = 1:3
          for j = 1:3
%             f(i,j)  = 0.0d0
%             df(i,j,l) = 0.0d0
            for k = 1:nel
              f(i,j) = f(i,j) + ul(i,k)*shpr(k,j);
%               df(i,j,l ) = df(i,j,l ) +  ul(i,k,2)           *shp(j,k,l)
            end %k
          end %j
%           f(i,i) = f(i,i) + 1;
        end %i

%       Deformation gradient at t_n: F_n

%         f(1,1,2,l)  = f(1,1,1,l) - df(1,1,l)
%         f(2,1,2,l)  = f(2,1,1,l) - df(2,1,l)
%         f(1,2,2,l)  = f(1,2,1,l) - df(1,2,l)
%         f(2,2,2,l)  = f(2,2,1,l) - df(2,2,l)

%         f(1,3,1,l)  = 0.0d0
%         f(3,1,1,l)  = 0.0d0
% 
%         f(2,3,1,l)  = 0.0d0
%         f(3,2,1,l)  = 0.0d0
% 
%         f(1,3,2,l)  = 0.0d0
%         f(3,1,2,l)  = 0.0d0
% 
%         f(2,3,2,l)  = 0.0d0
%         f(3,2,2,l)  = 0.0d0
% 
%         df(1,3,l)   = 0.0d0
%         df(3,1,l)   = 0.0d0
% 
%         df(2,3,l)   = 0.0d0
%         df(3,2,l)   = 0.0d0

%         if(stype.eq.3) then
%           f(3,3,1,l) = 0.0d0
%           xx1      = 0.0d0
%           df(3,3,l)  = 0.0d0
%           for k = 1,nel
%             xx1      = xx1      + xl(1,k  )*shp(3,k,l)
%             f(3,3,1,l) = f(3,3,1,l) + (ul(1,k,1) - u0(1,k))*shp(3,k,l)
%             df(3,3,l)  = df(3,3,l)  +  ul(1,k,2)           *shp(3,k,l)
%           end do
%           f(3,3,1,l) = 1.d0 + f(3,3,1,l)/xx1
%           df(3,3,l)  = df(3,3,l)/xx1
%           f(3,3,2,l) = f(3,3,1,l) - df(3,3,l)
%         else
%           f(3,3,1,l) = 1.0d0
%           f(3,3,2,l) = 1.0d0
%           df(3,3,l)  = 0.0d0
%         endif

%       Invert F

%     Compute adjoint to F

      ad(1,1) = f(2,2)*f(3,3) - f(2,3)*f(3,2);
      ad(1,2) = f(3,2)*f(1,3) - f(3,3)*f(1,2);
      ad(1,3) = f(1,2)*f(2,3) - f(1,3)*f(2,2);

      ad(2,1) = f(2,3)*f(3,1) - f(2,1)*f(3,3);
      ad(2,2) = f(3,3)*f(1,1) - f(3,1)*f(1,3);
      ad(2,3) = f(1,3)*f(2,1) - f(1,1)*f(2,3);

      ad(3,1) = f(2,1)*f(3,2) - f(2,2)*f(3,1);
      ad(3,2) = f(3,1)*f(1,2) - f(3,2)*f(1,1);
      ad(3,3) = f(1,1)*f(2,2) - f(1,2)*f(2,1);

%     Compute determinant of F

      detf  = f(1,1)*ad(1,1) + f(1,2)*ad(2,1) + f(1,3)*ad(3,1);
%       detf = f(1,1:3)*ad(1:3,1);

      detfi = 1.d0/detf;
      
%     Compute F inverse

      for j = 1:3
          for i = 1:3
              fi(i,j) = ad(i,j)*detfi;
          end 
      end
%       fi = detfi*ad;

if sflag == 1

%       Transform shape functions to current configuration

        for k = 1:nel
            for i = 1:3
                shpc(k,i) = fi(1,i)*shpr(k,1) + fi(2,i)*shpr(k,2) + fi(3,i)*shpr(k,3);
            end
        end %k
        
        for k = 1:nelP
            for i = 1:3
                shcp(k,i) = fi(1,i)*shpp(k,1) + fi(2,i)*shpp(k,2) + fi(3,i)*shpp(k,3);
            end
        end %k
        
end
        
%       end do ! l