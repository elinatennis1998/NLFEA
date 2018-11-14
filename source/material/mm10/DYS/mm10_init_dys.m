% c
% c
% c     DYS:
% c
% c           Initialize history
function [props, tau_tilde, uhist] = mm10_init_dys(props, tau_tilde, uhist)
%       use mm10_defs
%       implicit integer(a-z)
% c
%       type(crystal_props) :: props
%       double precision :: tau_tilde
%       double precision, dimension(max_uhard) :: uhist
% c
      tau_tilde(1) = props.c1; % Initial densities, page 5 left column
% c
%       return
end