% c
% c *****************************************************************************
% c *                                                                           *
% c *         Built in hardening routines                                       *
% c *                                                                           *
% c *****************************************************************************
% c
% c     Simple voche:
% c
% c           Initialize history
function [props, tau_tilde, uhist] = mm10_init_voche(props, tau_tilde, uhist)
%       use mm10_defs
% init_hard = ?;
% tau_y = 0;
%       implicit integer(a-z)
% c
%       type(crystal_props) :: props
%       double precision :: tau_tilde
%       double precision, dimension(max_uhard) :: uhist
%       double precision :: init_hard
init_hard = 0.1;
% c
      tau_tilde = props.tau_y+init_hard;
      % Nothing with the user history
%       return
end