% c
% c
% c     AADP:
% c
% c           Actual AADP hardening function
function [props, np1, n, stress, tt, h] = mm10_h_aadp(props, np1,...
    n, stress, tt)
%       use mm10_defs
%       implicit none
% c
%       type(crystal_props) :: props
%       type(crystal_state) :: np1, n
%       double precision, dimension(6) :: stress
%       double precision :: tt
%       double precision :: h
%       integer :: i
% c
%       double precision :: mm10_slipinc
%       double precision :: ct, cta
% c
      % Set some properties
      dt = np1.tinc;
      R_e = props.c2;
      R_s = props.c4;
      
      h = zeros(1,24);
      for i = 1:24
          if i <= 12

              % Get dislocation density
              rho_e = tt(i); % rho^a_edge
              rho_n = n.tau_tilde(i); % rho^a_edge_n-1
              rho_s = tt(i+12); % rho^a_edge
              
              alpha = i;
              ms = np1.ms(1:6,alpha);
              rs = stress*ms; % tau^a
              
              % Compute two dependencies
              vbar_e = mm10_vbare_aadp(props, np1, n, rs, tt, alpha,dt); % (A.11)
%               lbar_e = mm10_lbare_aadp(props, np1, n, rs, tt, i); % (A.9)
              vbar_s = mm10_vbars_aadp(props, np1, n, rs, tt, alpha,dt); % (A.12)
              lbar_s = mm10_lbars_aadp(props, np1, n, rs, tt, alpha); % (A.10)
              
              % Evaluate the equation
              h(i) = rho_n + dt*(2.d0*rho_s*vbar_s/lbar_s - rho_e^2.d0*R_e*vbar_e); % (A.7)
              
          else

              % Get dislocation density
              rho_s = tt(i); % rho^a_screw
              rho_n = n.tau_tilde(i); % rho^a_screw_n-1
              rho_e = tt(i-12); % rho^a_edge
              
              alpha = i - 12;
              ms = np1.ms(1:6,alpha);
              rs = stress*ms; % tau^a
              
              % Compute two dependencies
              vbar_e = mm10_vbare_aadp(props, np1, n, rs, tt, alpha,dt); % (A.11)
              lbar_e = mm10_lbare_aadp(props, np1, n, rs, tt, alpha); % (A.9)
              vbar_s = mm10_vbars_aadp(props, np1, n, rs, tt, alpha,dt); % (A.12)
%               lbar_s = mm10_lbars_aadp(props, np1, n, rs, tt, i); % (A.10)
              
              % Evaluate the equation
              h(i) = rho_n + dt*(2.d0*rho_e*vbar_e/lbar_e - rho_s^2.d0*R_s*vbar_s); % (A.7)
              
          end
      end
end