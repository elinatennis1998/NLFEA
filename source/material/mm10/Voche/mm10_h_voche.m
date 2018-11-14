% c
% c *****************************************************************************
% c *                                                                           *
% c *         Built in hardening routines                                       *
% c *                                                                           *
% c *****************************************************************************
% c
% c
% c           Actual voche law hardening function
function [props, np1, n, stress, tt, h] = ...
    mm10_h_voche(props, np1, n, stress, tt)
%       use mm10_defs
% 
%       implicit none
% c
%       type(crystal_props) :: props
%       type(crystal_state) :: np1, n
%       double precision, dimension(6) :: stress
%       double precision :: tt
%       double precision :: h, dg
%       integer :: i
% c
%       double precision :: mm10_slipinc
% c
       h = 0.0;
      for i = 1:props.nslip
        h = h + (1.0-(tt-props.tau_y)/props.tau_v+np1.tau_l(i)...
            /(tt-props.tau_y))^(props.voche_m)*...
            abs(mm10_slipinc(props.rate_n,np1.dg, np1.ms(:,i), stress, tt));
      end
      h = n.tau_tilde(1) + props.theta_0*h;
% c
%       return
end