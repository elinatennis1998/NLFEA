% c
% c
% c     ORNL:
% c
% c           Initialize history
function [props, tau_tilde, uhist] = mm10_init_ornl(props, tau_tilde, uhist)
%       use mm10_defs
%       implicit integer(a-z)
% c
%       type(crystal_props) :: props
%       double precision :: tau_tilde
%       double precision, dimension(max_uhard) :: uhist
% c
      if props.theta_0 > 100
          tau_tilde(1:props.num_hard) = props.theta_0; % user initialized
      else
      tau_tilde(1:props.num_hard) = 1e8; % Initial densities for edges
      end
% c
%       return
end