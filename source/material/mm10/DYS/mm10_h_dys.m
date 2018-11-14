% c
% c
% c     DYS:
% c
% c           Actual Roters hardening function
function [props, np1, n, stress, tt, h] = mm10_h_dys(props, np1,...
    n, stress, tt)

        dt = np1.tinc;
        Krho = props.c2;
        K = Krho*props.c1; % page 4 left column

% c
      sliptotal = 0;
      for alpha = 1:12
          
          slipinc = mm10_slipinc_dys(props, np1, n, stress, tt, alpha);
          gammadot = abs(slipinc/dt);
          sliptotal = sliptotal + gammadot;
          
      end
      
      h = n.tau_tilde(1) + dt*K*sliptotal; %(16) page 4
      
end