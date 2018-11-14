% c
% c     ****************************************************************
% c     *                                                              *
% c     *                 subroutine mm10_dgdt_DJGM                     *
% c     *                                                              *
% c     *                       written by : tjt                       *
% c     *                                                              *
% c     *                   last modified: 06/09/16                    *
% c     *                                                              *
% c     *     Calculate partial gamma(s) w.r.t. tau(s) for each slip   *
% c     *     system, for use in J11 in material integration. DJGM     *
% c     *                                                              *
% c     ****************************************************************
% c
function [props, np1, n, stress, tt, dgdt] = mm10_dgdt_DJGM(props, np1,...
            n, stress, tt)
% Deka, Dhyanjyoti, et al. "Crystal plasticity modeling of deformation and 
% creep in polycrystalline Ti-6242." Metallurgical and materials 
% transactions A 37.5 (2006): 1371-1388.
      dgdt = zeros(1,props.num_hard);

        % Load some material parameters
      dt = np1.tinc;

        % Load harderning parameters g
        [~, gamma_dot_tilde, ~, ~, ~, m, ~] = mm10_DJGM_GH(props);
        tau = transpose(stress*np1.ms);
        
        % Equation [4], slip rate vector
        for slip_a = 1:props.nslip
            dslipinc = dt * (1/m(slip_a))*gamma_dot_tilde(slip_a)*abs(tau(slip_a)/(tt(slip_a)))^(1/m(slip_a) - 1)*(1/tt(slip_a));
            dgdt(slip_a) = dslipinc;
        end
        
        
end