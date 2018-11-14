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
% c     *     system, for use in JA in tangent matrix. MTS             *
% c     *                                                              *
% c     ****************************************************************
% c
function [props, np1, n, stress, tt, D, dgammadd] = mm10_dgdd_mts_omar(props, np1,...
            n, stress, tt, D)

        d_mod = D;
        d_mod(4:6) = 0.5 * d_mod(4:6);
        alpha = 2.0/(3.0*np1.dg^2.0);
        teff_omar = stress*np1.ms-np1.backstress_omar;
        dgammadd = alpha*np1.dg * (((teff_omar/tt).^(props.rate_n).*sign(teff_omar)))'*d_mod;
        
end