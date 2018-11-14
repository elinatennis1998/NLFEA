% c
% c     ****************************************************************
% c     *                                                              *
% c     *                 subroutine mm10_dgdh_mts                     *
% c     *                                                              *
% c     *                       written by : tjt                       *
% c     *                                                              *
% c     *                   last modified: 01/26/15                    *
% c     *                                                              *
% c     *     Calculate partial gamma(s) w.r.t. tt for each slip       *
% c     *     system, for use in J12 in material integration. AADP     *
% c     *                                                              *
% c     ****************************************************************
% c
function [props, np1, n, stress, tt, dgammadtt] = mm10_dgdh_aadp(props, np1,...
            n, stress, tt)

        dgammadtt = zeros(12,24);
        
        for alpha = 1:12
            
            for i = 1:24
                
              if i <= 12
                  
              beta = i;
              es = 1;
              % Set some properties
              dt = np1.tinc;
              ms = np1.ms(1:6,alpha);
              rs = stress*ms; % tau^a
              b = props.burgers;

              % Get dislocation density
              rho_e = tt(alpha); % rho^a_edge

              % Compute one dependency
              vbar_e = mm10_vbare_aadp(props, np1, n, rs, tt, alpha,dt); % (A.11)

              % Get dislocation density
              rho_s = tt(alpha+12); % rho^a_screw

              % Compute one dependency
              vbar_s = mm10_vbars_aadp(props, np1, n, rs, tt, alpha,dt); % (A.12)
              
              dvedr = mm10_dvedr_aadp(props,np1,n,rs,tt,vbar_e,alpha,beta,es,dt);
              dvsdr = mm10_dvsdr_aadp(props,np1,n,rs,tt,vbar_s,alpha,beta,es,dt);  
      
              if alpha == beta
                  delta_ab = 1.d0;
              else
                  delta_ab = 0.d0;
              end
              
              % Evaluate the equation
              dgammadtt(alpha,i) = ...
                  dt * (delta_ab * vbar_e + rho_e * dvedr + rho_s * dvsdr) * b * sign(rs); % (A.6)
              
              else
                  
              beta = i - 12;
              es = 2;
              % Set some properties
              dt = np1.tinc;
              ms = np1.ms(1:6,alpha);
              rs = stress*ms; % tau^a
              b = props.burgers;

              % Get dislocation density
              rho_e = tt(alpha); % rho^a_edge

              % Compute one dependency
              vbar_e = mm10_vbare_aadp(props, np1, n, rs, tt, alpha,dt); % (A.11)

              % Get dislocation density
              rho_s = tt(alpha+12); % rho^a_screw

              % Compute one dependency
              vbar_s = mm10_vbars_aadp(props, np1, n, rs, tt, alpha,dt); % (A.12)
              
              dvedr = mm10_dvedr_aadp(props,np1,n,rs,tt,vbar_e,alpha,beta,es,dt);
              dvsdr = mm10_dvsdr_aadp(props,np1,n,rs,tt,vbar_s,alpha,beta,es,dt);  
      
              if alpha == beta
                  delta_ab = 1.d0;
              else
                  delta_ab = 0.d0;
              end
              
              % Evaluate the equation
              dgammadtt(alpha,i) = ...
                  dt * (delta_ab * vbar_s + rho_e * dvedr + rho_s * dvsdr) * b * sign(rs); % (A.6)
              
              end
              
            end
        
        end
        
end