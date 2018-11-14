% c
% c     ****************************************************************
% c     *                                                              *
% c     *                 subroutine mm10_dgdh_mts                     *
% c     *                                                              *
% c     *                       written by : tjt                       *
% c     *                                                              *
% c     *                   last modified: 12/17/14                    *
% c     *                                                              *
% c     *     Calculate partial gamma(s) w.r.t. tt for each slip       *
% c     *     system, for use in J12 in material integration. MTS      *
% c     *                                                              *
% c     ****************************************************************
% c
function [props, np1, n, stress, tt, dgammadtt] = mm10_dgdh_mts_omar(props, np1,...
            n, stress, tt)

        teff_omar = stress*np1.ms-np1.backstress_omar;
        dgammadtt = -props.rate_n/tt*np1.dg * (((teff_omar/tt).^(props.rate_n).*sign(teff_omar)))';
        
end