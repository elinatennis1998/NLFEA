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
% c     *     system, for use in J12 in material integration. Voche    *
% c     *                                                              *
% c     ****************************************************************
% c
function [props, np1, n, stress, tt, dgammadtt] = mm10_dgdh_voche(props, np1,...
            n, stress, tt)

        dgammadtt = -props.rate_n/tt*np1.dg * (((stress*np1.ms/tt).^(props.rate_n).*sign(stress*np1.ms)))';
        
end