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
% c           mrr model
function slipinc = ...
    mm10_slipinc_DJGM(props, np1, n, stress, tt, alpha)

        
        % Compute the shear modulus using Roter's function
        % Load the interaction matrices for parallel and forest dislocs
        [~, gamma_dot_tilde, ~, ~, ~, m, ~] = mm10_DJGM_GH(props);
        tau = transpose(stress*np1.ms);
        g = tt;
        dt = np1.tinc;
        
        % Equation [4]
        slipinc = dt*gamma_dot_tilde(alpha)*abs(tau(alpha)/(g(alpha)))^(1/m(alpha)).*sign(tau(alpha));
      
end
