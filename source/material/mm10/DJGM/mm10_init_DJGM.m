% c
% c
% c     MRR:
% c
% c           Initialize history
function [props, tau_tilde, uhist] = mm10_init_DJGM(props, tau_tilde, uhist)
%       use mm10_defs
%       implicit integer(a-z)
% c
%       type(crystal_props) :: props
%       double precision :: tau_tilde
%       double precision, dimension(max_uhard) :: uhist
% c
      if props.k_0 == -100
          if props.s_type == 9
              tau_tilde(1:props.num_hard) = [ones(3,1)*props.G_0_y;ones(3,1)*props.eps_dot_0_y]; % Initial g_0 (MPa)
          elseif props.s_type == 10
              tau_tilde(1:props.num_hard) = [ones(3,1)*props.G_0_y;ones(3,1)*props.eps_dot_0_y;ones(12,1)*props.G_0_v]; % Initial g_0 (MPa)
          else
              error('values for tau_tilde initial not specified for this slip system')
          end
      else
          if props.s_type == 9
              if props.theta_0 > 100*1e6
                  tau_tilde(1:props.num_hard) = props.theta_0; % user initialized
              elseif props.theta_0 == -1
                tau_tilde(1:props.num_hard) = 1e6*[3;3;3;2.4;2.4;2.4]*100; % Initial g_0 (MPa)

              else
              tau_tilde(1:props.num_hard) = 1e6*[3;3;3;2.4;2.4;2.4]*100; % Initial g_0 (MPa)
              end
          elseif props.s_type == 10
              if props.theta_0 > 100*1e6
                  tau_tilde(1:props.num_hard) = props.theta_0; % user initialized
              else
              tau_tilde(1:props.num_hard) = 1e6*[3;3;3;2.4;2.4;2.4;3*3*ones(12,1)]*100; % Initial g_0 (MPa)
              end
          else
              error('values for tau_tilde initial not specified for this slip system')
          end
      end
% c
%       return
end