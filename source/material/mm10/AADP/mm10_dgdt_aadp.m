% c
% c     ****************************************************************
% c     *                                                              *
% c     *                 subroutine mm10_dgdt_mts                     *
% c     *                                                              *
% c     *                       written by : tjt                       *
% c     *                                                              *
% c     *                   last modified: 01/26/15                    *
% c     *                                                              *
% c     *     Calculate partial gamma(s) w.r.t. tau(s) for each slip   *
% c     *     system, for use in J11 in material integration. AADP     *
% c     *                                                              *
% c     ****************************************************************
% c
function [props, np1, n, stress, tt, dgdt] = mm10_dgdt_aadp(props, np1,...
            n, stress, tt)

      dgdt = zeros(1,12);
      
      for i = 1:12

      % Set some properties
      dt = np1.tinc;
      ms = np1.ms(1:6,i);
      rs = stress*ms; % tau^a
      b = props.burgers;
      
      % Get dislocation density
      rho_e = tt(i); % rho^a_edge
      
      % Compute one dependency
      vbar_e = mm10_vbare_aadp(props, np1, n, rs, tt, i,dt); % (A.11)
      
      % Get dislocation density
      rho_s = tt(i+12); % rho^a_screw
      
      % Compute one dependency
      vbar_s = mm10_vbars_aadp(props, np1, n, rs, tt, i,dt); % (A.12)
              
      alpha = i;
      dvedt = mm10_dvedt_aadp(props,np1,n,rs,tt,vbar_e,alpha,dt);
      dvsdt = mm10_dvsdt_aadp(props,np1,n,rs,tt,vbar_s,alpha,dt);  
      
      % Evaluate the equation
      dgdt(i) = dt * (rho_e * dvedt + rho_s * dvsdt) * b * sign(rs); % (A.6)
      
      end
        
end