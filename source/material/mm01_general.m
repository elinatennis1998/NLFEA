function [history1,devstr_n1] = mm01_general ...
    (history,kbar,mrts,shear_mod_n1,hprime_n1,beta,rtse,yield,lk,alpha_n)
% c     ****************************************************************
% c     *                                                              *
% c     *                 subroutine mm01_general                      *
% c     *                                                              *
% c     *                       written by : rhd                       *
% c     *                                                              *
% c     *                   last modified: 10/1/00                     *
% c     *                                                              *
% c     *      stress update procedure when one or more elements       *
% c     *      in the block undergo a temperature change over step.    *
% c     *      the update procdedure handles additional terms in       *
% c     *      isotropic and kinematic hardening due to temperature    *
% c     *      changes over the step                                   *
% c     *                                                              *
% c     ****************************************************************
% c
history1 = zeros(11,1);
devstr_n1 = zeros(6,1);
one = 1.d0; two = 2.d0; three = 3.d0; root2 = sqrt(two); twthrd = two/three;
root23 = sqrt(twthrd); four = two*two; zero = 0.d0; root2o3 = sqrt(2/9);

%       do i = 1, span
      if ( yield )
% c
% c                       compute isotropic and kinematic plastic
% c                       moduli at start and end of step.
% c
      hbari_np1 = beta * hprime_n1;
      hbark_np1 = (one - beta) * hprime_n1;
%       hbari_n   = beta * history(5);
%       hbark_n   = (one - beta) * history(5);
% c
% c                       set up terms of the quadratic equation
% c                       to solve for the plastic multiplier
% c                       lambda * deltat
% c     
% c                       a) coefficients for the 3 tensors
% c                          denoted u, v, w in writeup
% c                       b) tensor products: v dot v, v dot w,
% c                          w dot w. these are vbar, vwbar and wbar
% c                       c) t1, t2, t3 are the resulting coefficients
% c                          of terms for the quadratic equation.
% c                          t1 multiplies lambda bar ** 2,
% c                          t1 multiplies lambda bar
% c                          t3 is the constant
% c                                             
      a = ( two*shear_mod_n1 + twthrd*hbark_np1 ) / mrts;
      b = one - lk;
      d = root2o3 * hbari_np1;
% c
      vbar  = mrts * mrts;
      wbar  = history(6)^2 + history(7)^2 + history(8)^2 + ...
              two*( history(9)^2 + history(10)^2 + ...
              history(11)^2 ) ;
      vwbar = rtse(1)*history(6) + rtse(2)*history(7) ...
              + rtse(3)*history(8) + ...
              two*( rtse(4)*history(9) + rtse(5)*history(10) + ...
              rtse(6)*history(11) ) ;
% c
      t1 = a*a*vbar - two*d*d;
      t2 = -four*d*kbar - two*a*vbar - two*a*b*vwbar;
      t3 = b*b*wbar + two*b*vwbar + vbar -two*kbar*kbar;
% c
% c                       compute discriminant of the quadratic.
% c                       if it is negative, we hit a wierd case
% c                       where the material is linear elastic. issue
% c                       a warning, set platic multiplier to zero.
% c                       compute the two roots, take the smaller
% c                       root to define lambda bar = lambda * deltat
% c
% c             
      discr = t2*t2 - four*t1*t3;
      if ( discr < zero ) %then
        error('something made the material elastic,mm01_general')
%         if ( debug ) 
%      &     write(iout,9200) i, hbari_np1, hbark_np1, hbari_n, hbark_n,
%      &                   a, b, d, vbar, wbar, vwbar, t1, t2, t3,
%      &                   discr
%         lambda_deltat = zero
      else
        qroot1 = (-t2 + sqrt( discr )) / (two*t1);
        qroot2 = (-t2 - sqrt( discr )) / (two*t1);
%         lambda_deltat = min( qroot1, qroot2 );
      end %if
% c
      lambda_deltat = min( qroot1, qroot2 );
      if ( lambda_deltat < zero ) % then
        error('negative lambda_deltat')
%         write(iout,9100)
%         call die_abort
      end %if
% c     
% c                       update scalars in the history. lambda *deltat
% c                       is used by the routine (cnst1) to compute the
% c                       consistent tangent modulus. save the updated
% c                       equivalent (shear) stress to set new size
% c                       of the yield cylinder for any amount
% c                       of isotropic hardening.
% c                       updated the accumulated equivalent
% c                       plastic strain (ebarp)
% c
      history1(1) = lambda_deltat;
      history1(2) = kbar + (root2/three)*hbari_np1*lambda_deltat;
      history1(3) = history(3) + lambda_deltat * root23;
      history1(5) = hprime_n1 ;
% c
% c                 updated backstresses and compute deviators for the
% c                 updated stress state.
% c
      const1 = twthrd * hbark_np1 * lambda_deltat / mrts;
      const2 = one - lambda_deltat * ( twthrd *  hbark_np1 + ...
               two * shear_mod_n1 ) / mrts;
      b = one - lk;
% c
      history1(6:11) = lk*history(6:11)  + const1 * rtse;

% c
      devstr_n1  = history1(6:11) + b*history(6:11) + ...
                        const2 * rtse;
%       devstr_n1(i,2)  = history1(i,7) + b*history(i,7) + 
%      &                  const2 * rtse(i,2)
%       devstr_n1(i,3)  = history1(i,8) + b*history(i,8) + 
%      &                  const2 * rtse(i,3)
%       devstr_n1(i,4)  = history1(i,9) + b*history(i,9) + 
%      &                  const2 * rtse(i,4)
%       devstr_n1(i,5)  = history1(i,10) + b*history(i,10) + 
%      &                  const2 * rtse(i,5)
%       devstr_n1(i,6)  = history1(i,11) + b*history(i,11) + 
%      &                  const2 * rtse(i,6)
% c
      end
%       end do
% c
% c                       update elements that are linear elastic at this
% c                       point. note: we save the updated yield
% c                       surface size and updated back stresses at n+1
% c                       to reflect any changes due to temperature.
% c                       add the backstresses at n to rtse to define
% c                       updated deviators of the trial elastic stress
% c                       at n+1
% c
%       do i = 1, span
         if( ~yield ) 
         history1(1)    = zero;
         history1(2)    = kbar;
%          if( lnelas(i) ) history1(i,2) = zero
         history1(3)    = history(3);
         history1(5)    = hprime_n1;
         history1(6:11) = alpha_n(1:6) * lk;
         devstr_n1(1:6) = rtse(1:6) + alpha_n(1:6);
         end
%       end do