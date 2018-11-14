% c
% c *****************************************************************************
% c *                                                                           *
% c *         Small constitutive functions                                      *
% c *                                                                           *
% c *****************************************************************************
% c
% c       from mm10_form.f
% c
% c           Calculate the slip increment along system i  
% c           AADP model
function slipinc = ...
    mm10_slipinc_aadp(props, np1, n, stress, tt, i)

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
      
      % Evaluate the equation
      slipinc = dt * (rho_e * vbar_e + rho_s * vbar_s) * b * sign(rs); % (A.6)
          
      end
