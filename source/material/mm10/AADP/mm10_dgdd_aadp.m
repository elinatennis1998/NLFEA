% c
% c     ****************************************************************
% c     *                                                              *
% c     *                 subroutine mm10_dgdd_mts                     *
% c     *                                                              *
% c     *                       written by : tjt                       *
% c     *                                                              *
% c     *                   last modified: 01/26/15                    *
% c     *                                                              *
% c     *     Calculate partial gamma(s) w.r.t. d for each slip        *
% c     *     system, for use in JA in tangent matrix. AADP            *
% c     *                                                              *
% c     ****************************************************************
% c
function [props, np1, n, stress, tt, D, dgammadd] = mm10_dgdd_aadp(props, np1,...
            n, stress, tt, D)

        dgammadd = zeros(12,6);
        
end