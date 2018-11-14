% c
% c *****************************************************************************
% c *                                                                           *
% c *         Built in hardening routines                                       *
% c *                                                                           *
% c *****************************************************************************
% c
% c
% c           Derivative of hardening fn wrt stress
function [props, np1, n, stress, tt, et] = ...
    mm10_estress_voche(props, np1, n, stress, tt)
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
%       double precision :: rs
%       integer :: i
% c
       et = zeros(6,1);
% c
      for i = 1:props.nslip
        rs = mm10_rs(np1, stress, i);
        et = et + (1.0-(tt-props.tau_y)/props.tau_v+np1.tau_l(i)...
            /(tt-props.tau_y))^(props.voche_m)*abs(rs)^...
            (props.rate_n-2)*rs*np1.ms(1:6,i);
      end
% c
      et = props.theta_0*np1.dg*props.rate_n/tt^(props.rate_n)*et;
% c
%       return
end