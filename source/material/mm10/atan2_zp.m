% c
% c     Helper for the above, atan2 with range 0 to 2*pi
% c
% c       from mm10_update_euler_angles.f lines 61-75
% c
function [atan2_zp] = atan2_zp(a, b)
%             implicit none
%             double precision :: atan2_zp, a, b, pi
%
pi = 4.0*atan(1.0); % 2nd definition for pi within the same file...
atan2_zp = atan2(a, b); %DATAN (fortran) = atan (matlab),
if (atan2_zp < 0.0)  % ATAN2(Y,X) (fortran) computes the arctangent of the complex number X + i Y.
    atan2_zp = atan2_zp + 2.0*pi;
end
%
%             return
end