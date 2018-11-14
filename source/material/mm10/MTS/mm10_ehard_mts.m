% c
% c
% c     MTS:
% c
% c           Derivative of hardening fn wrt hardening
function [props, np1, n, stress, tt, etau] = mm10_ehard_mts(props,...
    np1, n, stress, tt)
%       use mm10_defs
%       implicit none
% c
%       type(crystal_props) :: props
%       type(crystal_state) :: np1, n
%       double precision, dimension(6) :: stress
%       double precision :: tt, etau
% c
%       double precision :: mm10_slipinc
% c
%       double precision :: ur, s, ct, cta
%       integer :: i
% c
      cta = (props.mu_0/...
          np1.mu_harden)*tt - (props.mu_0/np1.mu_harden)*props.tau_a -...
          np1.tau_y;
      ct = 1.0 - cta/np1.tau_v;
% c
      ur = np1.mu_harden / props.mu_0;
% c
      rss = stress*np1.ms;
      rss2 = rss.*sign(real(rss));
       emat = (props.voche_m*(1/np1.tau_v+np1.tau_l'/cta^2.0).*...
            (ct+np1.tau_l'/cta).^(-1.0) + ur*props.rate_n/tt).*...
            (ct+np1.tau_l'/cta).^(props.voche_m).*(np1.dg * ((rss2/tt).^(props.rate_n).*sign(rss)));
      emat = emat.*sign(real(emat));
      etau = sum(emat(:,1:props.nslip),2);
% c
      etau = -props.theta_0*etau;
      if imag(etau) ~= 0 % Catch any negative values of (ct + np1.tau_l(i)/cta)
          etau = NaN; % Mark's code sets them to NaN, which causes an increment to the line search in the outer function
      end
      
      etau = 1.0 - etau;
%       return
end