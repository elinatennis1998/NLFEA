% c
% c     ****************************************************************
% c     *                                                              *
% c     *                 subroutine mm10_dgdd_mts                     *
% c     *                                                              *
% c     *                       written by : tjt                       *
% c     *                                                              *
% c     *                   last modified: 12/17/14                    *
% c     *                                                              *
% c     *     Calculate partial gamma(s) w.r.t. d for each slip        *
% c     *     system, for use in JA in tangent matrix. Voche           *
% c     *                                                              *
% c     ****************************************************************
% c
function [props, np1, n, stress, tt, D, dgammadd] = mm10_dgdd_voche(props, np1,...
            n, stress, tt, D)

        d_mod = D;
        d_mod(4:6) = 0.5 * d_mod(4:6);
        alpha = 2.0/(3.0*np1.dg^2.0);
        dgammadd = alpha*np1.dg * (((stress*np1.ms/tt).^(props.rate_n).*sign(stress*np1.ms)))'*d_mod;
        
end