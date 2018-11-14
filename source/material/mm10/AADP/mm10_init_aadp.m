% c
% c
% c     ORNL:
% c
% c           Initialize history
function [props, tau_tilde, uhist] = mm10_init_aadp(props, tau_tilde, uhist)
%       use mm10_defs
%       implicit integer(a-z)
% c
%       type(crystal_props) :: props
%       double precision :: tau_tilde
%       double precision, dimension(max_uhard) :: uhist
% c
      tau_tilde(1:12) = 4.165e10;4.165e12; % Initial densities for edges
      tau_tilde(13:24) = 4.165e10;4.165e12; % Initial densities for screws
% c
%       return
end