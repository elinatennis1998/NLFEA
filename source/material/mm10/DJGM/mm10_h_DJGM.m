function [props, np1, n,vec1,vec2, stress, tt, g] = mm10_h_DJGM(props, np1,...
    n,vec1,vec2, stress, tt)

% Deka, Dhyanjyoti, et al. "Crystal plasticity modeling of deformation and 
% creep in polycrystalline Ti-6242." Metallurgical and materials 
% transactions A 37.5 (2006): 1371-1388.

% evolution of hardening parameter 'g'

        % Load material parameters
        dt = np1.tinc;
        slipinc = vec1(1:props.nslip);
        slipinc = sign(real(slipinc)).*slipinc;
        gamma_dot = transpose(slipinc)/dt;
        g = tt;
        gtol = props.tau_hat_y;
        
        [h_0, gamma_dot_tilde, g_tilde, rr, nn, ~, g_0, q] = mm10_DJGM_GH(props);
                
        g_s = zeros(props.nslip,1);
        h = zeros(props.nslip,1);
        
        % Equation [6], g_s equation has to be modified
        for slip_b = 1:props.nslip
            temp = g_tilde(slip_b)*(gamma_dot(slip_b)/gamma_dot_tilde(slip_b))^nn(slip_b);
            g_s(slip_b) = max(temp, gtol*g_0(slip_b)); % threshold to ensure no slip during elastic response
            h(slip_b) = h_0(slip_b)*abs(1-g(slip_b)/g_s(slip_b))^rr(slip_b)*sign(1-g(slip_b)/g_s(slip_b));
        end
        
        % Equation [5]
        g_dot = zeros(props.nslip,1);
        g_n = transpose(n.tau_tilde(1:props.num_hard));
        for slip_a = 1:props.nslip
            for slip_b = 1:props.nslip
                g_dot(slip_a) = g_dot(slip_a) + q(slip_a, slip_b)*h(slip_b)*abs(gamma_dot(slip_b));
            end
            g(slip_a) = g_n(slip_a) + g_dot(slip_a)*dt;
        end
end