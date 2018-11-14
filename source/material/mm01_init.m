function [iostat,prior_linear,shear_mod,rtse,isothermal,lk,kbar,mrts,yf,instat,yield,eps_vol_n1,alpha_n] = ...
    mm01_init(deps,history,cgn,ym_n1,nu_n1,beta,hprime_n1,yld_n1,ym_n,nu_nn)
% c     ****************************************************************
% c     *                                                              *
% c     *                    subroutine mm01_init                      *
% c     *                                                              *
% c     *                       written by : rhd                       *
% c     *                                                              *
% c     *                   last modified: 10/1/00                     *
% c     *                                                              *
% c     *    basic set up of trial elastic state at n+1, pull terms    *
% c     *    history at n, stress deviators at n, etc. for all         *
% c     *    elements in the block                                     *
% c     *                                                              *
% c     ****************************************************************

devstr_n1_elas = zeros(6,1);

one = 1.d0; two = 2.d0; three = 3.d0; %temper_tol = 0.00001d00;
htol = 0.00001d00; root3 = sqrt(3); root2 = sqrt(2); yld_tol = 0.0000001d00;
% zero = 0.0d00;

isothermal = 1;

%       do i = 1, span

         iostat           = history(4);
         prior_linear = iostat == 3;
% c
% c                       deviatoric components of strain increment over
% c                       step.
% c 
         deps_vol    = deps(1) + deps(2) + deps(3) ;
         eps_mean    = deps_vol / three;
         de1         = deps(1) - eps_mean;
         de2         = deps(2) - eps_mean;
         de3         = deps(3) - eps_mean;
         de4         = deps(4);
         de5         = deps(5);
         de6         = deps(6);
% c
% c                      elastic components of total strain at start of
% c                      step. 
% c
         e_n         = ym_n;
         nu_n        = nu_nn;
         g_n         = e_n/two/(one+nu_n);
         een1        = (cgn(1)-nu_n*(cgn(2)+cgn(3)))/e_n;
         een2        = (cgn(2)-nu_n*(cgn(1)+cgn(3)))/e_n;
         een3        = (cgn(3)-nu_n*(cgn(1)+cgn(2)))/e_n;
         een4        = cgn(4) / g_n;
         een5        = cgn(5) / g_n;
         een6        = cgn(6) / g_n;
% c
% c                      deviatoric components of elastic strain
% c                      at start of step plus deviatoric strain
% c                      increment over the step. (volumetric term
% c                      at n+1 saved for final update operation)
% c
         eps_vol_n1 = een1 + een2 + een3 + deps_vol;
         eps_mean_n   = (een1 + een2 + een3 ) / three;
         e1           = (een1 -  eps_mean_n) + de1 ;
         e2           = (een2 -  eps_mean_n) + de2 ;
         e3           = (een3 -  eps_mean_n) + de3 ;
         e4           = een4 + de4;
         e5           = een5 + de5;
         e6           = een6 + de6;
% c
% c                       compute deviators for trial elastic
% c                       stress state. Uses deviatoric elastic strain at n
% c                       + the deviatoric strain increment over the step
% c                       and temperature dependent moduli at 
% c                       n+1. this is way to get temperature effects
% c                       on modulus and poisson's ratio properly 
% c                       included.
% c                  
         shear_mod  = ym_n1 / (two*(one+nu_n1));
% c
         devstr_n1_elas(1) = two * shear_mod * e1;
         devstr_n1_elas(2) = two * shear_mod * e2;
         devstr_n1_elas(3) = two * shear_mod * e3;
         devstr_n1_elas(4) = shear_mod * e4;
         devstr_n1_elas(5) = shear_mod * e5;
         devstr_n1_elas(6) = shear_mod * e6;
% c        
% c                       pull out backstress at start of step
% c                       at start of step) for ease of access
% c
         alpha_n =  history(6:11);   
%          alpha_n(i,2) =  history(i,7)   
%          alpha_n(i,3) =  history(i,8)   
%          alpha_n(i,4) =  history(i,9)   
%          alpha_n(i,5) =  history(i,10)   
%          alpha_n(i,6) =  history(i,11)   
% c
% c                       set isotropic and kinematic plastic hardening
% c                       moduli at start and end of step (can
% c                       be different due to temperature).
% c                       set parameter kbar - shear yield stress at
% c                       n+1 based on yield stress at n+1, current 
% c                       plastic strain and isotropic (plastic) hardening
% c                       modulus at n+1. this sets the size of yield
% c                       surface to check for yielding with trial stress
% c                       state. Also used later in stress update
% c                       procedures. 
% c  
% c                       set the lk factor for temperature dependent
% c                       kinematic hardening.
% c
         hbari_np1  = beta * hprime_n1;
         hbark_np1  = (one - beta) * hprime_n1;
         hbark_n    = (one - beta) * history(5);
         kbar    = (yld_n1 + hbari_np1*history(3))/root3;
         lk      = one;
         if ( abs( hbark_n ) > htol ) 
             lk = hbark_np1 / hbark_n;
         end
% c          
% c                       compute deviators of relative trial elastic
% c                       stresses at n+1.
% c                       note use of the backstress at n updated
% c                       to the temperature at n+1. compute norm
% c                       of relative trial stress and evaluate yield
% c                       criterion.
% c
         rtse = devstr_n1_elas - alpha_n*lk;
%          rtse(i,2) = devstr_n1_elas(2) - alpha_n(i,2)*lk(i)
%          rtse(i,3) = devstr_n1_elas(3) - alpha_n(i,3)*lk(i)
%          rtse(i,4) = devstr_n1_elas(4) - alpha_n(i,4)*lk(i)
%          rtse(i,5) = devstr_n1_elas(5) - alpha_n(i,5)*lk(i)
%          rtse(i,6) = devstr_n1_elas(6) - alpha_n(i,6)*lk(i)
% c
         mrts = sqrt( rtse(1)^2+rtse(2)^2+rtse(3)^2 + ...
                    two*(rtse(4)^2+rtse(5)^2+rtse(6)^2) );
         yf = mrts - root2 * kbar;
% c
% c                      the relative trial (deviator) stress used
% c                      in the actual update procedures does not
% c                      reflect the back stress at n updated for
% c                      the temperature at n+1.
% c
         rtse = devstr_n1_elas - alpha_n;
%          rtse(i,2) = devstr_n1_elas(2) - alpha_n(i,2)
%          rtse(i,3) = devstr_n1_elas(3) - alpha_n(i,3)
%          rtse(i,4) = devstr_n1_elas(4) - alpha_n(i,4)
%          rtse(i,5) = devstr_n1_elas(5) - alpha_n(i,5)
%          rtse(i,6) = devstr_n1_elas(6) - alpha_n(i,6)
         mrts = sqrt( rtse(1)^2+rtse(2)^2+rtse(3)^2 + ...
                    two*(rtse(4)^2+rtse(5)^2+rtse(6)^2) );
% c
% c                      update the isothermal flag to reflect the
% c                      status of this element
% c
%          if ( abs( dtemps(i) ) .gt. temper_tol )  isothermal = .false.
% c
% c                      set various logical flags based on element
% c                      status determined by the trial elastic state.
% c
% c                      if the lnelas flag is true, the element response
% c                      must always be linear elastic no matter what.
% c
% c                      if this gauss point for element is yielding,
% c                      set flags.
% c
% c                      state variable:
% c                         = 1,  point is actively yielding
% c                         = 3,  point is not actively yielding
% c
% c
         instat = 3;
         yield  = 0;
%          if( lnelas(i) ) cycle
% c
         if( yf >= yld_tol*root2*kbar ) % then
            yield = 1;
            instat = 1;
         end %if
% c
%       end do