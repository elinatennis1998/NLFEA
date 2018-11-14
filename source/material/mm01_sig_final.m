function [cgn1,history1] = mm01_sig_final(cgn,history1,ym,nu,shear_mod,devstr_n1,deps,instat,eps_vol_n1)
% c     ****************************************************************
% c     *                                                              *
% c     *                    subroutine mm01_sig_final                 *
% c     *                                                              *
% c     *                       written by : rhd                       *
% c     *                                                              *
% c     *                   last modified: 10/1/00                     *
% c     *                                                              *
% c     *                                                              *
% c     ****************************************************************
% c
cgn1 = zeros(9,1);
one = 1.d0; two = 2.d0; three = 3.d0; half = one/two;
% c
% c                       compute the updated stresses from their
% c                       deviator values at state (n+1) and the
% c                       (linear elastic) mean stress contribution.
% c                       save the state variable, elastic modulus
% c                       and poisson's ratio at n+1 in history.
% c                       calculate the energy density from a 
% c                       trapezoidal numerical integration of
% c                       increments of strain and average stresses
% c
%       do i = 1, span
         sig_mean_np1 = eps_vol_n1 ...
                        *(three*ym*nu/((one+nu)* ...
                        (one-two*nu)) + two*shear_mod)/three;
         cgn1(1) = devstr_n1(1) + sig_mean_np1;
         cgn1(2) = devstr_n1(2) + sig_mean_np1;
         cgn1(3) = devstr_n1(3) + sig_mean_np1;
         cgn1(4) = devstr_n1(4);
         cgn1(5) = devstr_n1(5);
         cgn1(6) = devstr_n1(6);
         cgn1(7) = cgn(7) + half * ( ...
             deps(1) * (cgn1(1) + cgn(1)) ...
           + deps(2) * (cgn1(2) + cgn(2)) ...
           + deps(3) * (cgn1(3) + cgn(3)) ...
           + deps(4) * (cgn1(4) + cgn(4)) ...
           + deps(5) * (cgn1(5) + cgn(5)) ...
           + deps(6) * (cgn1(6) + cgn(6)) );
% c
%          iword(1)       = instat;
         history1(4)  = instat;
%       end do 