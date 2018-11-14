% c
% c
% c     ****************************************************************
% c     *                                                              *
% c     *                 subroutine mm10_formarrs                     *
% c     *                                                              *
% c     *                       written by : tjt                       *
% c     *                                                              *
% c     *                   last modified: 05/29/15                    *
% c     *                                                              *
% c     *     Form R2                                                  *
% c     *                                                              *
% c     ****************************************************************
% c
function [props, np1, n, stress, tt, vec1, vec2, arr1, arr2] = mm10_formarrs(props, np1,...
    n, stress, tt, vec1, vec2, arr1, arr2,both)
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
        [props, np1, n, stress, tt, arr1, arr2] = mm10_a_ornl(props, np1,...
            n, stress, tt, arr1, arr2,both);
      elseif (props.h_type == 5) % AADP
      elseif (props.h_type == 6) % mts_omar
      elseif (props.h_type == 7) % MRR
        [props, np1, n, stress, tt, arr1, arr2] = mm10_a_mrr(props, np1,...
            n, stress, tt, arr1, arr2,both);
      elseif (props.h_type == 8) % DYS
      elseif (props.h_type == 9) % Ran_HCP
        [props, np1, n, stress, tt, arr1, arr2] = mm10_a_DJGM(props, np1,...
            n, stress, tt, arr1, arr2,both);
      elseif (props.h_type == 12) % Plugin Library
          AllOutput = mm10_plugin_lib(16,props,{props, np1,...
              n, vec1,vec2,arr1,arr2, stress, tt, both});
          arr1 = AllOutput{1};
          arr2 = AllOutput{2};
      elseif (props.h_type == 14) % halite
          [props, np1, n, stress, tt, arr1, arr2] = mm10_a_halite(props, np1,...
              n, stress, tt, arr1, arr2,both);
      else
        [props] = mm10_unknown_hard_error(props);
      end
% c
end