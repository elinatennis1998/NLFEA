% c
% c *****************************************************************************
% c *                                                                           *
% c *         Built in hardening routines                                       *
% c *                                                                           *
% c *****************************************************************************
% c
% c
% c           Derivative of hardening fn wrt hardening
function [props, np1, n, stress, tt, etau] = ...
    mm10_ehard_voche(props, np1, n, stress, tt)
%       use mm10_defs
%       implicit none
% c
%       type(crystal_props) :: props
%       type(crystal_state) :: np1, n
%       double precision, dimension(6) :: stress
%       double precision :: tt, etau
% c
%       double precision :: mm10_slipinc
%       integer :: i
% c
       etau = 0.0;
      for i = 1:props.nslip
        etau = etau + (props.voche_m*(1.0/props.tau_v+np1.tau_l(i)...
            /(tt-props.tau_y)^2.0)*(1.0-(tt-props.tau_y)...
            /props.tau_v+np1.tau_l(i)/(tt-props.tau_y))^(-1.0)...
            +props.rate_n/tt)*(1.0-(tt-props.tau_y)/props.tau_v+...
            np1.tau_l(i)/(tt-props.tau_y))^(props.voche_m)*...
            abs(mm10_slipinc(props.rate_n,np1.dg, np1.ms(:,i), stress, tt));
      end
% 
      etau = -props.theta_0*etau;
      etau = 1.0 - etau;
% c
%       return
end