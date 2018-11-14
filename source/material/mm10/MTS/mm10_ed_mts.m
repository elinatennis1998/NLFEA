% c
% c
% c     MTS:
% c
% c           Derivative of hardening fn wrt strain
function [props, np1, n, stress, tt, ed] = mm10_ed_mts(props, np1,...
    n, stress, tt, ed)
%       use mm10_defs
%       implicit none
% c
%       type(crystal_props) :: props
%       type(crystal_state) :: np1, n
%       double precision, dimension(6) :: stress
%       double precision :: tt
%       double precision, dimension(6) :: ed
% c
%       double precision :: mm10_slipinc
% c
%       double precision :: dslip, lnv, lny, dgc, ty, tv, mnp0, m0np, sc
%       double precision, dimension(6) :: d_mod, dydd, dvdd
%       integer :: s
% c
% c     Form a bunch of simple constitutive things
% c
% c
      d_mod = np1.D;
      d_mod(4:6) = 0.5 * d_mod(4:6);
% c
      dgc = np1.dg / np1.tinc;
% c
% c     Form d_tauy/d_deltad by steps
% c
      lny = log(props.eps_dot_0_y/dgc);
      ty = props.boltzman*np1.temp...
          /(np1.mu_harden*(props.burgers^3)*props.G_0_y)*lny;
      dydd = 2.0*props.tau_hat_y/(3.0*np1.dg^2.0*props.q_y*props.p_y*...
          lny)*(1.0-ty^(1.0/props.q_y))^(1.0/props.p_y-1.0)*...
          ty^(1.0/props.q_y)*d_mod;
% 
% c
% c     Form d_tauv/d_deltad by steps
% c
      lnv = log(props.eps_dot_0_v/dgc);
      tv = props.boltzman*np1.temp...
          /(np1.mu_harden*(props.burgers^3)*props.G_0_v)*lnv;
      dvdd = 2.0*props.tau_hat_v/(3.0*np1.dg^2.0*props.q_v*props.p_v*...
          lnv)*(1.0-tv^(1.0/props.q_v))^(1.0/props.p_v-1.0)*...
          tv^(1.0/props.q_v)*d_mod;
% 
% c
% c     Form a couple more common components
% c
      mnp0 = np1.mu_harden / props.mu_0;
      sc = tt/mnp0 - props.tau_a/mnp0 - np1.tau_y;
% c
% c     Glue everything together
% c
      edmat = (dydd'*(props.voche_m*(1/np1.tau_v+np1.tau_l'/sc^2.0).*...
            (1.0-sc/np1.tau_v+np1.tau_l'/sc).^(props.voche_m-1.0))...
            + dvdd'*props.voche_m/(np1.tau_v)^2.0*sc*...
            (1.0-sc/np1.tau_v+np1.tau_l'/sc).^(props.voche_m-1.0)...
            + d_mod'*2.0/(3.0*np1.dg^2.0)*(1.0-sc/np1.tau_v+np1.tau_l'/sc)...
            .^(props.voche_m)).*...
            (ones(6,1)*abs(np1.dg * ((stress*np1.ms/tt).^(props.rate_n).*sign(stress*np1.ms))));
      ed = sum(edmat(:,1:props.nslip),2);
      
      ed = props.theta_0 * mnp0 * ed + mnp0*dydd';
% 
%       return
end