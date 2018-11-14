% c
% c
% c     AADP:
% c
% c           Derivative of hardening fn wrt stress
function [props, np1, n, stress, tt, et] = mm10_estress_aadp(props,...
    np1, n, stress, tt)
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
% c
%       double precision :: rs, cta, ct
%       integer :: i
      
      dt = np1.tinc;
      R_e = props.c2;
      R_s = props.c4;
      
      et = zeros(24,6);
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
              
              dvedt = mm10_dvedt_aadp(props,np1,n,rs,tt,vbar_e,alpha,dt);
              dvsdt = mm10_dvsdt_aadp(props,np1,n,rs,tt,vbar_s,alpha,dt);  
              
              dtdstress = ms;
              
              % Evaluate the equation
              et(i,1:6) = dt*(2.d0*rho_s*dvsdt/lbar_s - rho_e^2.d0*R_e*dvedt)*dtdstress; % (A.7)
              
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
              
              dvedt = mm10_dvedt_aadp(props,np1,n,rs,tt,vbar_e,alpha,dt);
              dvsdt = mm10_dvsdt_aadp(props,np1,n,rs,tt,vbar_s,alpha,dt);  
              
              dtdstress = ms;
              
              % Evaluate the equation
              et(i,1:6) = dt*(2.d0*rho_e*dvedt/lbar_e - rho_s^2.d0*R_s*dvsdt)*dtdstress; % (A.7)
              
          end
      end
% 
% c
%       return
end