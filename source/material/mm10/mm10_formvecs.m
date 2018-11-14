% c
% c
% c     ****************************************************************
% c     *                                                              *
% c     *                 subroutine mm10_formvecs                     *
% c     *                                                              *
% c     *                       written by : tjt                       *
% c     *                                                              *
% c     *                   last modified: 05/29/15                    *
% c     *                                                              *
% c     *     Form R2                                                  *
% c     *                                                              *
% c     ****************************************************************
% c
function [props, np1, n, stress, tt, vec1, vec2] = mm10_formvecs(props, np1,...
    n, stress, tt, vec1, vec2)
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
      elseif (props.h_type == 2) % MTS
      elseif (props.h_type == 3) % User
      elseif (props.h_type == 4) % ORNL
        [props, np1, n, stress, tt, vec1, vec2] = mm10_v_ornl(props, np1,...
            n, stress, tt, vec1, vec2);
      elseif (props.h_type == 5) % AADP
      elseif (props.h_type == 5) % MTS omar
      elseif (props.h_type == 7) % MRR
        [props, np1, n, stress, tt, vec1, vec2] = mm10_v_mrr(props, np1,...
            n, stress, tt, vec1, vec2);
      elseif (props.h_type == 8) % DYS
      elseif (props.h_type == 9) % Ran_HCP
        [props, np1, n, stress, tt, vec1, vec2] = mm10_v_DJGM(props, np1,...
            n, stress, tt, vec1, vec2);
      elseif (props.h_type == 12) % Plugin Library
          AllOutput = mm10_plugin_lib(15,props,{props, np1,...
            n,stress, tt, vec1,vec2});
          vec1 = AllOutput{1};
          vec2 = AllOutput{2};
      elseif (props.h_type == 14) % halite
        [props, np1, n, stress, tt, vec1, vec2] = mm10_v_halite(props, np1,...
            n, stress, tt, vec1, vec2);
      else
        [props] = mm10_unknown_hard_error(props);
      end
% c
end