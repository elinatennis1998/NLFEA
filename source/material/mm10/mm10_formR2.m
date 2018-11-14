% c
% c
% c     ****************************************************************
% c     *                                                              *
% c     *                 subroutine mm10_formR2                       *
% c     *                                                              *
% c     *                       written by : mcm                       *
% c     *                                                              *
% c     *                   last modified: 11/26/13                    *
% c     *                                                              *
% c     *     Form R2                                                  *
% c     *                                                              *
% c     ****************************************************************
% c
function [props, np1, n, vec1, vec2, stress, tt, R2] = mm10_formR2(props, np1,...
    n, vec1, vec2, stress, tt)
%       use mm10_defs
%       implicit none
% c
%       type(crystal_props) :: props
%       type(crystal_state) :: np1, n
%       double precision, dimension(6) :: stress
%       double precision :: R2, h
%       double precision :: tt
% c
      if (props.h_type == 1) % voche
        [props, np1, n, stress, tt, h] = mm10_h_voche(props, np1,...
            n, stress, tt);
      elseif (props.h_type == 2) % MTS
        [props, np1, n, stress, tt, h] = mm10_h_mts(props, np1,...
            n, stress, tt);
      elseif (props.h_type == 3) % User
        [props, np1, n, stress, tt, h] = mm10_h_user(props, np1,...
            n, stress, tt);
      elseif (props.h_type == 4) % ORNL
        [props, np1, n, vec1, vec2, stress, tt, h] = mm10_h_ornl(props, np1,...
            n, vec1, vec2, stress, tt);
      elseif (props.h_type == 5) % AADP
        [props, np1, n, stress, tt, h] = mm10_h_aadp(props, np1,...
            n, stress, tt);
elseif (props.h_type == 6) % MTS_omar copy
    [props, np1, n, stress, tt, h] = mm10_h_mts_omar(props, np1,...
        n, stress, tt);
      elseif (props.h_type == 7) % MRR
        [props, np1, n, vec1, vec2, stress, tt, h] = mm10_h_mrr(props, np1,...
            n, vec1, vec2, stress, tt);
      elseif (props.h_type == 8) % DYS
        [props, np1, n, stress, tt, h] = mm10_h_dys(props, np1,...
            n, stress, tt);
      elseif (props.h_type == 9) % Ran_HCP
        [props, np1, n, vec1, vec2, stress, tt, h] = mm10_h_DJGM(props, np1,...
            n, vec1, vec2, stress, tt);
      elseif (props.h_type == 12) % Plugin library
          h = mm10_plugin_lib(4,props,{props, np1,...
            n, vec1, vec2, stress, tt});
      elseif (props.h_type == 14) % halite
        [props, np1, n, vec1, vec2, stress, tt, h] = mm10_h_halite(props, np1,...
            n, vec1, vec2, stress, tt);
      else
        [props] = mm10_unknown_hard_error(props);
      end
% c
      R2 = tt - h;
      % Compute time rate of change of hardening parameters, store for
      % possible prediction during next time/load step
      np1.tt_rate(1:props.num_hard) = (h(1:props.num_hard) - n.tau_tilde(1:props.num_hard)')/np1.tinc;
% c      
%       return
end