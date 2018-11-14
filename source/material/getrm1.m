function q = getrm1(r,opt)
% compute the T matrix needed for transforming rotated to unrotated
% configurations
% Q is the R in WARP3D manual
% T = [Q(1,1)^2 Q(1,2)^2 Q(1,3)^2 2*Q(1,1)*Q(1,2) 2*Q(1,2)*Q(1,3) 2*Q(1,3)*Q(1,1)
%              Q(2,1)^2 Q(2,2)^2 Q(2,3)^2 2*Q(2,1)*Q(2,2) 2*Q(2,2)*Q(2,3) 2*Q(2,3)*Q(2,1)
%              Q(3,1)^2 Q(3,2)^2 Q(3,3)^2 2*Q(3,1)*Q(3,2) 2*Q(3,2)*Q(3,3) 2*Q(3,3)*Q(3,1)
%              Q(1,1)*Q(2,1) Q(1,2)*Q(2,2) Q(1,3)*Q(2,3) Q(1,1)*Q(2,2)+Q(1,2)*Q(2,1) Q(1,2)*Q(2,3)+Q(1,3)*Q(2,2) Q(1,3)*Q(2,1)+Q(1,1)*Q(2,3)
%              Q(2,1)*Q(3,1) Q(2,2)*Q(3,2) Q(2,3)*Q(3,3) Q(2,1)*Q(3,2)+Q(2,2)*Q(3,1) Q(2,2)*Q(3,3)+Q(2,3)*Q(3,2) Q(2,3)*Q(3,1)+Q(2,1)*Q(3,3)
%              Q(3,1)*Q(1,1) Q(3,2)*Q(1,2) Q(3,3)*Q(1,3) Q(3,1)*Q(1,2)+Q(3,2)*Q(1,1) Q(3,2)*Q(1,3)+Q(3,3)*Q(1,2) Q(3,3)*Q(1,1)+Q(3,1)*Q(1,3)];
         q = zeros(6,6);
         two = 2.d0;
         
         switch opt
             case 1
% c
% c      unrotated rate of deformation vector {d} = [q] *
% c       (rotated) rate of deformation vector {D}. in tensor form:
% c
% c                 [d] = trans([R]) [D] [R]
% c
% c       both [D] and [d] are symmetric, [R] is orthogonal rotation.
% c       vector forms for {d} and {D} use engineering shear strains.
% c       vector ordering is {x,y,z,xy,yz,xz}
            q(1,1)= r(1,1)^2;
            q(1,2)= r(2,1)^2;
            q(1,3)= r(3,1)^2;
            q(1,4)= r(1,1)*r(2,1);
            q(1,5)= r(3,1)*r(2,1);
            q(1,6)= r(1,1)*r(3,1);
            q(2,1)= r(1,2)^2;
            q(2,2)= r(2,2)^2;
            q(2,3)= r(3,2)^2;
            q(2,4)= r(1,2)*r(2,2);
            q(2,5)= r(3,2)*r(2,2);
            q(2,6)= r(1,2)*r(3,2);
            q(3,1)= r(1,3)^2;
            q(3,2)= r(2,3)^2;
            q(3,3)= r(3,3)^2;
            q(3,4)= r(1,3)*r(2,3);
            q(3,5)= r(3,3)*r(2,3);
            q(3,6)= r(1,3)*r(3,3);
            q(4,1)= two*r(1,1)*r(1,2);
            q(4,2)= two*r(2,1)*r(2,2);
            q(4,3)= two*r(3,1)*r(3,2);
            q(4,4)= r(1,1)*r(2,2)+r(1,2)*r(2,1);
            q(4,5)= r(2,1)*r(3,2)+r(3,1)*r(2,2);
            q(4,6)= r(1,1)*r(3,2)+r(3,1)*r(1,2);
            q(5,1)= two*r(1,2)*r(1,3);
            q(5,2)= two*r(2,3)*r(2,2);
            q(5,3)= two*r(3,2)*r(3,3);
            q(5,4)= r(1,2)*r(2,3)+r(2,2)*r(1,3);
            q(5,5)= r(2,2)*r(3,3)+r(2,3)*r(3,2);
            q(5,6)= r(1,2)*r(3,3)+r(3,2)*r(1,3);
            q(6,1)= two*r(1,1)*r(1,3);
            q(6,2)= two*r(2,1)*r(2,3);
            q(6,3)= two*r(3,1)*r(3,3);
            q(6,4)= r(1,1)*r(2,3)+r(2,1)*r(1,3);
            q(6,5)= r(2,1)*r(3,3)+r(3,1)*r(2,3);
            q(6,6)= r(1,1)*r(3,3)+r(1,3)*r(3,1);
% c
%       if ( opt .eq. 2 ) then
% c
% c       cauchy stress {T} = [q] * (rotated) cauchy stress {t}.
% c       in tensor form:
% c
% c                 [T] = [R] [t] trans([R])
% c
% c       both [T] and [t] are symmetric, [R] is orthogonal rotation.
% c       vector ordering is {x,y,z,xy,yz,xz}. this [q] matrix
% c       is the transpose of the one above.
% c
%          do i = 1, span
             case 2
            q(1,1)= r(1,1)^2;
            q(1,2)= r(1,2)^2;
            q(1,3)= r(1,3)^2;
            q(1,4)= two*r(1,1)*r(1,2);
            q(1,5)= two*r(1,3)*r(1,2);
            q(1,6)= two*r(1,1)*r(1,3);
            q(2,1)= r(2,1)^2;
            q(2,2)= r(2,2)^2;
            q(2,3)= r(2,3)^2;
            q(2,4)= two*r(2,1)*r(2,2);
            q(2,5)= two*r(2,3)*r(2,2);
            q(2,6)= two*r(2,1)*r(2,3);
            q(3,1)= r(3,1)^2;
            q(3,2)= r(3,2)^2;
            q(3,3)= r(3,3)^2;
            q(3,4)= two*r(3,1)*r(3,2);
            q(3,5)= two*r(3,3)*r(3,2);
            q(3,6)= two*r(3,1)*r(3,3);
            q(4,1)= r(1,1)*r(2,1);
            q(4,2)= r(1,2)*r(2,2);
            q(4,3)= r(1,3)*r(2,3);
            q(4,4)= r(1,1)*r(2,2)+r(2,1)*r(1,2);
            q(4,5)= r(1,2)*r(2,3)+r(1,3)*r(2,2);
            q(4,6)= r(1,1)*r(2,3)+r(1,3)*r(2,1);
            q(5,1)= r(2,1)*r(3,1);
            q(5,2)= r(3,2)*r(2,2);
            q(5,3)= r(2,3)*r(3,3);
            q(5,4)= r(2,1)*r(3,2)+r(2,2)*r(3,1);
            q(5,5)= r(2,2)*r(3,3)+r(3,2)*r(2,3);
            q(5,6)= r(2,1)*r(3,3)+r(2,3)*r(3,1);
            q(6,1)= r(1,1)*r(3,1);
            q(6,2)= r(1,2)*r(3,2);
            q(6,3)= r(1,3)*r(3,3);
            q(6,4)= r(1,1)*r(3,2)+r(1,2)*r(3,1);
            q(6,5)= r(1,2)*r(3,3)+r(1,3)*r(3,2);
            q(6,6)= r(1,1)*r(3,3)+r(3,1)*r(1,3);
%          end do
%          return
%       end if
% c
% c
%       if ( opt .eq. 3 ) then
% c
% c       unrotated cauchy stress {t} = [q] * cauchy stress {T}.
% c       in tensor form:
% c
% c                 [t] = trans([R]) [T] [R]
% c
% c       want to use code above for opt = 2. Set rbar = trans([R])
% c       and compute [q]. We are computing the
% c
% c       both [T] and [t] are symmetric, [R] is orthogonal rotation.
% c       vector ordering is {x,y,z,xy,yz,xz}. this [q] matrix
% c       is the transpose of the one above.
% c
%          do i = 1, span
             case 3
                 rbar = zeros(3,3);
          rbar(1,1) = r(1,1);
          rbar(1,2) = r(2,1);
          rbar(1,3) = r(3,1);
          rbar(2,1) = r(1,2);
          rbar(2,2) = r(2,2);
          rbar(2,3) = r(3,2);
          rbar(3,1) = r(1,3);
          rbar(3,2) = r(2,3);
          rbar(3,3) = r(3,3);
%          end do
% c
%          do i = 1, span
           q(1,1)= rbar(1,1)^2;
           q(1,2)= rbar(1,2)^2;
           q(1,3)= rbar(1,3)^2;
           q(1,4)= two*rbar(1,1)*rbar(1,2);
           q(1,5)= two*rbar(1,3)*rbar(1,2);
           q(1,6)= two*rbar(1,1)*rbar(1,3);
           q(2,1)= rbar(2,1)^2;
           q(2,2)= rbar(2,2)^2;
           q(2,3)= rbar(2,3)^2;
           q(2,4)= two*rbar(2,1)*rbar(2,2);
           q(2,5)= two*rbar(2,3)*rbar(2,2);
           q(2,6)= two*rbar(2,1)*rbar(2,3);
           q(3,1)= rbar(3,1)^2;
           q(3,2)= rbar(3,2)^2;
           q(3,3)= rbar(3,3)^2;
           q(3,4)= two*rbar(3,1)*rbar(3,2);
           q(3,5)= two*rbar(3,3)*rbar(3,2);
           q(3,6)= two*rbar(3,1)*rbar(3,3);
           q(4,1)= rbar(1,1)*rbar(2,1);
           q(4,2)= rbar(1,2)*rbar(2,2);
           q(4,3)= rbar(1,3)*rbar(2,3);
           q(4,4)= rbar(1,1)*rbar(2,2)+rbar(2,1)*rbar(1,2);
           q(4,5)= rbar(1,2)*rbar(2,3)+rbar(1,3)*rbar(2,2);
           q(4,6)= rbar(1,1)*rbar(2,3)+rbar(1,3)*rbar(2,1);
           q(5,1)= rbar(2,1)*rbar(3,1);
           q(5,2)= rbar(3,2)*rbar(2,2);
           q(5,3)= rbar(2,3)*rbar(3,3);
           q(5,4)= rbar(2,1)*rbar(3,2)+rbar(2,2)*rbar(3,1);
           q(5,5)= rbar(2,2)*rbar(3,3)+rbar(3,2)*rbar(2,3);
           q(5,6)= rbar(2,1)*rbar(3,3)+rbar(2,3)*rbar(3,1);
           q(6,1)= rbar(1,1)*rbar(3,1);
           q(6,2)= rbar(1,2)*rbar(3,2);
           q(6,3)= rbar(1,3)*rbar(3,3);
           q(6,4)= rbar(1,1)*rbar(3,2)+rbar(1,2)*rbar(3,1);
           q(6,5)= rbar(1,2)*rbar(3,3)+rbar(1,3)*rbar(3,2);
           q(6,6)= rbar(1,1)*rbar(3,3)+rbar(3,1)*rbar(1,3);
         end