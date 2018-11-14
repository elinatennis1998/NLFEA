function [props, np1, n, stress, g, vec1, vec2] = mm10_v_DJGM(props, np1,...
    n, stress, tt, vec1, vec2)

        dt = np1.tinc;
        % Load harderning parameters g
        [~, gamma_dot_tilde, ~, ~, ~, m, ~] = mm10_DJGM_GH(props);
        tau = transpose(stress*np1.ms);
        g = tt;
        
        % Equation [4], slip rate vector
        for slip_a = 1:props.nslip
            vec1(slip_a) = dt * gamma_dot_tilde(slip_a)*abs(tau(slip_a)/(g(slip_a)))^(1/m(slip_a)).*sign(tau(slip_a));
        end

end