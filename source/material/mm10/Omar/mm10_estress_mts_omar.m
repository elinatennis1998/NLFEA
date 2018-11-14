% c
% c
% c     MTS:
% c
% c           Derivative of hardening fn wrt stress
function [props, np1, n, stress, tt, et] = mm10_estress_mts_omar(props,...
    np1, n, stress, tt)
%       use mm10_defs
%       implicit none
% c
%       type(crystal_props) :: props
%       type(crystal_state) :: np1, n
%       double precision, dimension(6) :: stress
%       double precision :: tt
%       double precision, dimension(6) :: et
% c
%       double precision :: mm10_rs
% c
%       double precision :: rs, cta, ct
%       integer :: i
      cta = (props.mu_0/...
          np1.mu_harden)*tt - (props.mu_0/np1.mu_harden)*props.tau_a -...
          np1.tau_y;
      ct = 1.0 - cta/np1.tau_v;
      
      teff_omar = stress*np1.ms-np1.backstress_omar;
      etmat = (ones(6,1)*((ct+np1.tau_l'/cta).^(props.voche_m).*abs(teff_omar).^(props.rate_n-2.0).*(teff_omar))).*np1.ms;
      et = sum(etmat(:,1:props.nslip),2);
      
      et =  props.theta_0 * ...
          (np1.mu_harden/props.mu_0)*et*props.rate_n*...
          np1.dg/tt^props.rate_n;
      et = et';
% 
% c
%       return
end