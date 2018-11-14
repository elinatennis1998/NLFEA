% c
% c
% c     AADP:
% c
% c           Derivative of hardening fn wrt hardening
function [props, np1, n, stress, tt, etau] = mm10_ehard_aadp(props,...
    np1, n, stress, tt)
%       use mm10_defs
%       implicit none
% c
%       type(crystal_props) :: props
%       type(crystal_state) :: np1, n
%       double precision, dimension(6) :: stress
%       double precision :: tt, etau
% c
%       double precision :: mm10_slipinc
% c
%       double precision :: ur, s, ct, cta
%       integer :: i
% c
      % Set some properties
      dt = np1.tinc;
      R_e = props.c2;
      R_s = props.c4;
      
      etau = zeros(24,24);
      
      for j = 1:24
          if j <= 12
              
              beta = j;
              es = 1; % edge dislocation
              
      for i = 1:24
          if i <= 12

              % Get dislocation density
              rho_e = tt(i); % rho^a_edge
              rho_s = tt(i+12); % rho^a_edge
              
              alpha = i;
              ms = np1.ms(1:6,alpha);
              rs = stress*ms; % tau^a
              
              % Compute two dependencies
              vbar_e = mm10_vbare_aadp(props, np1, n, rs, tt, alpha,dt); % (A.11)
%               lbar_e = mm10_lbare_aadp(props, np1, n, rs, tt, i); % (A.9)
              vbar_s = mm10_vbars_aadp(props, np1, n, rs, tt, alpha,dt); % (A.12)
              lbar_s = mm10_lbars_aadp(props, np1, n, rs, tt, alpha); % (A.10)
              
              dvedr = mm10_dvedr_aadp(props,np1,n,rs,tt,vbar_e,alpha,beta,es,dt);
              dvsdr = mm10_dvsdr_aadp(props,np1,n,rs,tt,vbar_s,alpha,beta,es,dt);
%               dledr = mm10_dledr_aadp(props,np1,rs,tt,alpha,beta,es);
              dlsdr = mm10_dlsdr_aadp(props,np1,n,rs,tt,alpha,beta,es);
              
              if i == j
                  delta_ij = 1.d0;
              else
                  delta_ij = 0.d0;
              end
              
              % Evaluate the equation
              etau(i,j) = delta_ij - ...
                  dt*(2.d0*rho_s*dvsdr/lbar_s -2.d0*rho_s*vbar_s/(lbar_s^2.d0)*dlsdr ...
                      - 2.d0*rho_e*delta_ij*R_e*vbar_e - rho_e^2.d0*R_e*dvedr); % (A.7)
              
          else

              % Get dislocation density
              rho_s = tt(i); % rho^a_screw
              rho_e = tt(i-12); % rho^a_edge
              
              alpha = i - 12;
              ms = np1.ms(1:6,alpha);
              rs = stress*ms; % tau^a
              
              % Compute two dependencies
              vbar_e = mm10_vbare_aadp(props, np1, n, rs, tt, alpha,dt); % (A.11)
              lbar_e = mm10_lbare_aadp(props, np1, n, rs, tt, alpha); % (A.9)
              vbar_s = mm10_vbars_aadp(props, np1, n, rs, tt, alpha,dt); % (A.12)
%               lbar_s = mm10_lbars_aadp(props, np1, n, rs, tt, i); % (A.10)
              
              dvedr = mm10_dvedr_aadp(props,np1,n,rs,tt,vbar_e,alpha,beta,es,dt);
              dvsdr = mm10_dvsdr_aadp(props,np1,n,rs,tt,vbar_s,alpha,beta,es,dt);
              dledr = mm10_dledr_aadp(props,np1,n,rs,tt,alpha,beta,es);
%               dlsdr = mm10_dlsdr_aadp(props,np1,n,rs,tt,alpha,beta,es);
              
              if i == j
                  delta_ij = 1.d0;
              else
                  delta_ij = 0.d0;
              end
              
              % Evaluate the equation
              etau(i,j) = delta_ij - ...
                  dt*(2.d0*rho_e*dvedr/lbar_e -2.d0*rho_e*vbar_e/(lbar_e^2.d0)*dledr ...
                      - 2.d0*rho_s*delta_ij*R_s*vbar_s - rho_s^2.d0*R_s*dvsdr); % (A.7)
              
          end
      end
      
          else % screw dislocation
              
              beta = j - 12;
              es = 2; % screw dislocation
              
      for i = 1:24
          if i <= 12

              % Get dislocation density
              rho_e = tt(i); % rho^a_edge
              rho_s = tt(i+12); % rho^a_edge
              
              alpha = i;
              ms = np1.ms(1:6,alpha);
              rs = stress*ms; % tau^a
              
              % Compute two dependencies
              vbar_e = mm10_vbare_aadp(props, np1, n, rs, tt, alpha,dt); % (A.11)
%               lbar_e = mm10_lbare_aadp(props, np1, n, rs, tt, i); % (A.9)
              vbar_s = mm10_vbars_aadp(props, np1, n, rs, tt, alpha,dt); % (A.12)
              lbar_s = mm10_lbars_aadp(props, np1, n, rs, tt, alpha); % (A.10)
              
              dvedr = mm10_dvedr_aadp(props,np1,n,rs,tt,vbar_e,alpha,beta,es,dt);
              dvsdr = mm10_dvsdr_aadp(props,np1,n,rs,tt,vbar_s,alpha,beta,es,dt);
%               dledr = mm10_dledr_aadp(props,np1,rs,tt,alpha,beta,es);
              dlsdr = mm10_dlsdr_aadp(props,np1,n,rs,tt,alpha,beta,es);
              
              if i == j
                  delta_ij = 1.d0;
              else
                  delta_ij = 0.d0;
              end
              
              % Evaluate the equation
              etau(i,j) = delta_ij - ...
                  dt*(2.d0*rho_s*dvsdr/lbar_s -2.d0*rho_s*vbar_s/(lbar_s^2.d0)*dlsdr ...
                      - 2.d0*rho_e*delta_ij*R_e*vbar_e - rho_e^2.d0*R_e*dvedr); % (A.7)
              
          else

              % Get dislocation density
              rho_s = tt(i); % rho^a_screw
              rho_e = tt(i-12); % rho^a_edge
              
              alpha = i - 12;
              ms = np1.ms(1:6,alpha);
              rs = stress*ms; % tau^a
              
              % Compute two dependencies
              vbar_e = mm10_vbare_aadp(props, np1, n, rs, tt, alpha,dt); % (A.11)
              lbar_e = mm10_lbare_aadp(props, np1, n, rs, tt, alpha); % (A.9)
              vbar_s = mm10_vbars_aadp(props, np1, n, rs, tt, alpha,dt); % (A.12)
%               lbar_s = mm10_lbars_aadp(props, np1, n, rs, tt, i); % (A.10)
              
              dvedr = mm10_dvedr_aadp(props,np1,n,rs,tt,vbar_e,alpha,beta,es,dt);
              dvsdr = mm10_dvsdr_aadp(props,np1,n,rs,tt,vbar_s,alpha,beta,es,dt);
              dledr = mm10_dledr_aadp(props,np1,n,rs,tt,alpha,beta,es);
%               dlsdr = mm10_dlsdr_aadp(props,np1,n,rs,tt,alpha,beta,es);
              
              if i == j
                  delta_ij = 1.d0;
              else
                  delta_ij = 0.d0;
              end
              
              % Evaluate the equation
              etau(i,j) = delta_ij - ...
                  dt*(2.d0*rho_e*dvedr/lbar_e -2.d0*rho_e*vbar_e/(lbar_e^2.d0)*dledr ...
                      - 2.d0*rho_s*delta_ij*R_s*vbar_s - rho_s^2.d0*R_s*dvsdr); % (A.7)
              
          end
      end
              
          end
          
      end % j
      
%       etau = 1.0 - etau;
%       return
end