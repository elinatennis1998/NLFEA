function [cgn1,history1,rtse,yield] = mm01(ym_n1,nu_n1,beta,hprime_n1,yld_n1,cgn,deps,history,ym_n,nu_n)
% c     ****************************************************************
% c     *                                                              *
% c     *                      subroutine mm01                         *
% c     *                                                              *
% c     *                       written by : rhd                       *
% c     *                                                              *
% c     *                   last modified : 6/21/12                    *
% c     *                                                              *
% c     *     this subroutine performs the recovery of the             *
% c     *     unrotated cauchy stress for the rate-independent         *
% c     *     mises material model for a block of similar, non-        *
% c     *     conflicting elements. mixed isotropic-kinematic          *
% c     *     hardening is supported. the material uniaxial stress-    *
% c     *     strain curve is bilinear. material response can be       *
% c     *     temperature dependent (e, nu, yield, hardening modulus)  *
% c     *     temperature dependent thermal expansion coefficients do  *
% c     *     not come into play here for stress updates - thermal     *
% c     *     strain increments are computed and subtracted out        *
% c     *     before we get here                                       *
% c     *                                                              *
% c     ****************************************************************
% c
% c    
% c
% c          update the stresses and internal state variables at gauss
% c          point "gpn" for "span" elements in the current block.
% c          values marked with (*) are updated by this routine
% c
% c          span     -- number of elements in this block
% c          felem    -- actual element number for first element in blk
% c          gpn      -- gauss point number being processed on
% c                      this call for elements in the block
% c          step     -- current load step number
% c          iter     -- current (global) netwon iteration number (>1)
% c          ym_n1    -- young's modulus for elements
% c                      in the block at temperature for step n+1
% c          nu_n1    -- poisson's ratio for elements
% c                      in the block at temperature for step n+1
% c          beta     -- isotropic/kinematic fractional factor for elements
% c                      in the block
% c         hprime_n1 -- plastic hardening modulus for elements
% c                      in the block at temperature for step n+1
% c          lnelas   -- logical flag indicating element response
% c                      is always linear elastic and temperature
% c                      independent
% c          yld_n1   -- uniaxial yield stress for elements
% c                      in the block at temperature for step n+1
% c          cgn      -- cartesian stresses at step n for elements
% c                      in the block
% c    (*)   cgn1     -- cartesian stresses at end of this iteration for
% c                      load step n+1 (output) for elements in the block
% c          deps     -- cartesian increments of "mechanical" strain
% c                      between n and n+1 for elements in the block
% c          history  -- history data at step n for elements in the block
% c    (*)   history1 -- history data at step n+1 for elements in the block
% c    (*)   rtse     -- "relative" trial elstic stress state for elements
% c                      in the block (output). computed and used here and
% c                      later by consistent tangent generator (routine cnst1).
% c          dtemps   -- temperature change over step at this gp for
% c                      each element in block
% c          ym_n     -- young's modulus for elements
% c                      in the block at temperature for step n
% c          nu_n     -- poisson's ratio for elements
% c                      in the block at temperature for step n
% c
% c      Notes:
% c      ------
% c
% c        o   the variable mxvl is defined as a parameter in the
% c            include file param_def. it sets the maximum possible
% c            number of elements in a block allowed in analyses.
% c            most arrays are sized based on this variable
% c
% c        o   most parameters are vectors of length mxvl
% c
% c        o   the stresses for elements in the block (cgn, cgn1)
% c            are arrays sized mxvl by *. The ordering of terms is
% c            xx, yy, zz, xy, yz, xz, energy density, plastic
% c            work density, acummulated (incremental) plastic
% c            strain
% c
% c        o   the strain increment (deps) for elements in the block
% c            is an array sized mxvl by *. The ordering of terms is
% c            xx, yy, zz, gam-xy, gam-yz, gam-xz. The thermal
% c            strain contribution has been subtracted before this
% c            routine is called.
% c
% c        o   the rtse array (mxvl by 6) stores the trial elastic stresses
% c            at n+1. They are computed here and returned for use
% c            later by the consistent tangent routine.
% c
% c        o   the ordering of terms in the history vectors is
% c            described in mm01_set_history

% c
% c                       do the basic setup of the trial elastic stress
% c                       state, deviators at n, pull backstress at n from
% c                       history, etc. to compute elastic trial
% c                       state we use e, nu at n+1. evaluate the
% c                       material state as elastic or currently plastic.
% c         
            [iostat,prior_linear,shear_mod_n1,rtse,isothermal,lk,kbar,mrts,yf,instat,yield,eps_vol_n1,alpha_n] = ...
    mm01_init(deps,history,cgn,ym_n1,nu_n1,beta,hprime_n1,yld_n1,ym_n,nu_n);
% c
% c                       compute updated yield surface size, updated
% c                       backstresses and save in history. compute
% c                       deviators of updated stress stateat n+1. if all
% c                       elements at this gauss point have no temperature
% c                       change over the step, use the faster
% c                       isothermal update procedure. these update
% c                       routines use their own loops over span and skip
% c                       linear elastic elements 
% c         
            [history1,devstr_n1] = mm01_general ...
    (history,kbar,mrts,shear_mod_n1,hprime_n1,beta,rtse,yield,lk,alpha_n);

% c
% c                       compute the total updated stresses from their
% c                       deviator values at state (n+1) and the
% c                       mean stress (linear elastic) contribution.
% c                       save the state variable, the elastic
% c                       modulus and poisson's ratio at n+1 in
% c                       updated history vector.
% c                       calculate the energy density from a 
% c                       trapezoidal numerical integration of
% c                       increments of strain and average stresses
% c
         [cgn1,history1] = mm01_sig_final(cgn,history1,ym_n1,nu_n1,shear_mod_n1,devstr_n1,deps,instat,eps_vol_n1);
% c
% c                       update the plastic
% c                       work density for elements in the block. 
% c
         [cgn1] = mm01_plastic_work(cgn,cgn1,yield,deps,nu_n1,ym_n1,shear_mod_n1);