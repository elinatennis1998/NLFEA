% c
% c
% c     MTS:
% c
% c           Initialize history
function [props, tau_tilde, uhist] = mm10_init_mts_omar(props, tau_tilde, uhist)
%       use mm10_defs
%       implicit integer(a-z)
% c
%       type(crystal_props) :: props
%       double precision :: tau_tilde
%       double precision, dimension(max_uhard) :: uhist
% c
      tau_tilde = -1.0; % This only works because these are actually flags
      uhist(1) = -1.0;
      uhist(2) = -1.0;
% c
%       return
end