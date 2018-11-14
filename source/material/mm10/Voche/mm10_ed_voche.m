% c
% c *****************************************************************************
% c *                                                                           *
% c *         Built in hardening routines                                       *
% c *                                                                           *
% c *****************************************************************************
% c
% c
% c           Derivative of hardening fn wrt strain
function [props, np1, n, stress, tt, ed] = ...
    mm10_ed_voche(props, np1, n, stress, tt)
%       use mm10_defs
%       implicit none
% c
%       type(crystal_props) :: props
%       type(crystal_state) :: np1, n
%       double precision, dimension(6) :: stress
%       double precision :: tt
%       double precision, dimension(6) :: ed
% c
%       double precision :: h
%       double precision, dimension(6) :: d_mod
% c
%       ed = 0.0;
% c
      [props, np1, n, stress, tt, h] = ...
          mm10_h_voche(props, np1, n, stress, tt);
      d_mod = np1.D;
      d_mod(4:6) = 0.5 * d_mod(4:6);
% c
      ed = 2.0*(h - n.tau_tilde)/(3.0*np1.dg^2.0)*d_mod;
% c
%       return
end