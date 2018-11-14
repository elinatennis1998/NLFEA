function [r,ui] = rtcmp1(f)
% c     ****************************************************************
% c     *                                                              *
% c     *                      subroutine rtcmp1                       *
% c     *                                                              *
% c     *                       written by : bh                        *
% c     *                                                              *
% c     *                   last modified : 06/30/91                   *
% c     *                                                              *
% c     *     this subroutine computes the polar decompostion of the   *
% c     *     deformation gradient into the rotation tensor [R] and a  *
% c     *     deformation tensor [U] for a block of solid elements     *
% c     *                                                              *
% c     ****************************************************************
% c
% c
% c

% c
% c                       compute the inverse of the right
% c                       stretch tensor.
% c
      ui = irscp1(f);
% c
% c                       compute the rotation tensor.
% c
%       do i= 1,span
         r(1,1)= f(1,1)*ui(1)+f(1,2)*ui(2)+f(1,3)*ui(4);
         r(1,2)= f(1,1)*ui(2)+f(1,2)*ui(3)+f(1,3)*ui(5);
         r(1,3)= f(1,1)*ui(4)+f(1,2)*ui(5)+f(1,3)*ui(6);
         r(2,1)= f(2,1)*ui(1)+f(2,2)*ui(2)+f(2,3)*ui(4);
         r(2,2)= f(2,1)*ui(2)+f(2,2)*ui(3)+f(2,3)*ui(5);
         r(2,3)= f(2,1)*ui(4)+f(2,2)*ui(5)+f(2,3)*ui(6);
         r(3,1)= f(3,1)*ui(1)+f(3,2)*ui(2)+f(3,3)*ui(4);
         r(3,2)= f(3,1)*ui(2)+f(3,2)*ui(3)+f(3,3)*ui(5);
         r(3,3)= f(3,1)*ui(4)+f(3,2)*ui(5)+f(3,3)*ui(6);
%       end do
end

function ui = irscp1(f)
% c
% c     ****************************************************************
% c     *                                                              *
% c     *                      subroutine irscp1                       *
% c     *                                                              *
% c     *                       written by : bh                        *
% c     *                                                              *
% c     *                   last modified : 06/30/91                   *
% c     *                                                              *
% c     *     this subroutine computes the inverse of the right        *
% c     *     stretch tensor. the computations are for a gauss         *
% c     *     point for a block of solid elements                      *
% c     *                                                              *
% c     ****************************************************************
% c
% c
% c

one = 1.d0; two = 2.d0;
ui = zeros(6,1);
% c
% c                       ui is in symmetric upper triangular form.
% c
% c                       compute the invariants of the right
% c                       stretch tensor, the metric tensor, and
% c                       its square.
% c
      [c, cc, iu, iiu, iiiu] = ivcmp1(f);
% c
% c                       compute multipliers.
% c
%       do i = 1, span
         a2= one/(iiiu*(iu*iiu-iiiu));
         b2= iu*iiu*iiu-iiiu*(iu*iu+iiu);
         c2= -iiiu-iu*(iu*iu-two*iiu);
         d2= iu;
%       end do
% c
% c                       compute the inverse of the right
% c                       stretch tensor.
% c
%       do i = 1, span
         ui(1)= a2 * ( b2 + c2*c(1) + d2*cc(1) );
         ui(2)= a2 * (         c2*c(2) + d2*cc(2) );
         ui(3)= a2 * ( b2 + c2*c(3) + d2*cc(3) );
         ui(4)= a2 * (         c2*c(4) + d2*cc(4) );
         ui(5)= a2 * (         c2*c(5) + d2*cc(5) );
         ui(6)= a2 * ( b2 + c2*c(6) + d2*cc(6) );
%       end do
% c
%       return
end
      
function [c, cc, iu, iiu, iiiu] = ivcmp1(f)
% c     ****************************************************************
% c     *                                                              *
% c     *                      subroutine ivcmp1                       *
% c     *                                                              *
% c     *                       written by : bh                        *
% c     *                                                              *
% c     *                   last modified : 06/30/91                   *
% c     *                                 : 02/08/94                   *
% c     *                                                              *
% c     *     this subroutine computes the invariants of the right     *
% c     *     stretch tensor, the metric tensor, and its square.       *
% c     *     the computations are for a gauss point in a bblock of    *
% c     *     solid elements                                           *
% c     *                                                              *
% c     ****************************************************************
% c
% c
c = zeros(6,1);
cc = zeros(6,1);
ct = zeros(6,1);
% c
% c                       c and cc are in symmetric upper triangular
% c                       form.
% c
% c                       compute the metric tensor.
% c
%       do i = 1, span
         c(1)= f(1,1)*f(1,1)+f(2,1)*f(2,1)+f(3,1)*f(3,1);
         c(2)= f(1,1)*f(1,2)+f(2,1)*f(2,2)+f(3,1)*f(3,2);
         c(3)= f(1,2)*f(1,2)+f(2,2)*f(2,2)+f(3,2)*f(3,2);
         c(4)= f(1,1)*f(1,3)+f(2,1)*f(2,3)+f(3,1)*f(3,3);
         c(5)= f(1,2)*f(1,3)+f(2,2)*f(2,3)+f(3,2)*f(3,3);
         c(6)= f(1,3)*f(1,3)+f(2,3)*f(2,3)+f(3,3)*f(3,3);
%       end do
% c
% c                       compute the square of the metric
% c                       tensor.
% c
%       do i = 1, span
         cc(1)= c(1)*c(1)+c(2)*c(2)+c(4)*c(4);
         cc(2)= c(1)*c(2)+c(2)*c(3)+c(4)*c(5);
         cc(3)= c(2)*c(2)+c(3)*c(3)+c(5)*c(5);
         cc(4)= c(1)*c(4)+c(2)*c(5)+c(4)*c(6);
         cc(5)= c(2)*c(4)+c(3)*c(5)+c(5)*c(6);
         cc(6)= c(4)*c(4)+c(5)*c(5)+c(6)*c(6);
%       end do
% c
% c                       copy the metric tensor to stress vector
% c                       form so that principal values may be
% c                       computed.
% c
%       do i = 1, span
         ct(1)= c(1);
         ct(2)= c(3);
         ct(3)= c(6);
         ct(4)= c(2);
         ct(5)= c(5);
         ct(6)= c(4);
%       end do
% c
% c                       compute the principal values of the
% c                       metric tensor.
% c
ev = evcmp1(ct);
ct2 = [ct(1) ct(4) ct(6); ct(4) ct(2) ct(5); ct(6) ct(5) ct(3)];
ev2 = eig(ct2);
% c
% c                       set the principal values.
% c
%       do i = 1, span
         ev(1)= sqrt(ev(1));
         ev(2)= sqrt(ev(2));
         ev(3)= sqrt(ev(3));
%       end do
% c
% c
% c                       compute the invariants of the right
% c                       stretch tensor.
% c
% c
%       do i = 1, span
         iu  = ev(1)+ev(2)+ev(3);
         iiu = ev(1)*ev(2)+ev(2)*ev(3)+ev(1)*ev(3);
         iiiu= ev(1)*ev(2)*ev(3);
%       end do
% c
%       return
      end