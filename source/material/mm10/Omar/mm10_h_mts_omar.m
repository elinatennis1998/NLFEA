% c
% c
% c     MTS:
% c
% c           Actual MTS hardening function
function [props, np1, n, stress, tt, h] = mm10_h_mts(props, np1,...
    n, stress, tt)
%       use mm10_defs
%       implicit none
% c
%       type(crystal_props) :: props
%       type(crystal_state) :: np1, n
%       double precision, dimension(6) :: stress
%       double precision :: tt
%       double precision :: h
%       integer :: i
% c
%       double precision :: mm10_slipinc
%       double precision :: ct, cta
% c
      cta = (props.mu_0/...
          np1.mu_harden)*tt - (props.mu_0/np1.mu_harden)*props.tau_a -...
          np1.tau_y; % (3.12.7) solved for tau_bar := cta
      ct = 1.0 - cta/np1.tau_v;
      
      teff_omar = stress*np1.ms-np1.backstress_omar;
      rss = teff_omar;
      rss2 = rss.*sign(real(rss));
       hmat = ((ct + np1.tau_l'/cta).^(props.voche_m).*(np1.dg * (rss2/tt).^(props.rate_n).*sign(real(rss))));
      hmat = hmat.*sign(real(hmat));
       h = sum(hmat(:,1:props.nslip),2); % without the theta_0 factor, this is (3.12.13) when voche_m=1.0
%       if imag(h) ~= 0 % Catch any negative values of (ct + np1.tau_l(i)/cta)
%           h = NaN; % Mark's code sets them to NaN, which causes an increment to the line search in the outer function
%       end
% c
      h = props.tau_a*(1.0 - np1.mu_harden/n.mu_harden) + ...
          (np1.mu_harden/props.mu_0)*(np1.tau_y - n.tau_y) + ...
          (np1.mu_harden/n.mu_harden)*n.tau_tilde(1) + props.theta_0 * ...
          (np1.mu_harden/props.mu_0)*h; % seems like a convoluted solution of (3.12.7) in rate form, where tau_tilde is initialized in mm10_setup_mts
%       return
end