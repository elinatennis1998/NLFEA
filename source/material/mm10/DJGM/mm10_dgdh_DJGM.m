% c
% c     ****************************************************************
% c     *                                                              *
% c     *                 subroutine mm10_dgdh_DJGM                    *
% c     *                                                              *
% c     *                       written by : tjt                       *
% c     *                                                              *
% c     *                   last modified: 06/09/16                    *
% c     *                                                              *
% c     *     Calculate partial gamma(s) w.r.t. tt for each slip       *
% c     *     system, for use in J12 in material integration. DJGM     *
% c     *                                                              *
% c     ****************************************************************
% c
function [props, np1, n, stress, tt, dgammadtt] = mm10_dgdh_DJGM(props, np1,...
            n, stress, tt)
        
      dgammadtt = zeros(props.num_hard,props.num_hard);

        % Load harderning parameters g
        [~, gamma_dot_tilde, ~, ~, ~, m, ~] = mm10_DJGM_GH(props);
        tau = transpose(stress*np1.ms);

        % Load some material parameters
      dt = np1.tinc;
      for slip_a = 1:props.num_hard
        slip_b = slip_a;
        dslipinc = dt * (1/m(slip_a))*gamma_dot_tilde(slip_a)*abs(tau(slip_a)/(tt(slip_a)))^(1/m(slip_a))*(-1/tt(slip_a))*sign(tau(slip_a));
        
          dgammadtt(slip_a,slip_b) = dslipinc;
      end
        
        
end