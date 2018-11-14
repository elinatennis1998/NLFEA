% c
% c     ****************************************************************
% c     *                                                              *
% c     *                 subroutine mm10_dgdt_mts                     *
% c     *                                                              *
% c     *                       written by : tjt                       *
% c     *                                                              *
% c     *                   last modified: 12/17/14                    *
% c     *                                                              *
% c     *     Calculate partial gamma(s) w.r.t. tau(s) for each slip   *
% c     *     system, for use in J11 in material integration. MTS      *
% c     *                                                              *
% c     ****************************************************************
% c
function [props, np1, n, stress, tt, alpha] = mm10_dgdt_mts_omar(props, np1,...
            n, stress, tt)

        teff_omar = stress*np1.ms-np1.backstress_omar;
        alpha = abs(teff_omar).^(props.rate_n-1.0);
        alpha = np1.dg*props.rate_n/tt^(props.rate_n)*alpha;
        
end