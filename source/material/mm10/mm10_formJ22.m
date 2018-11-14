% c
% c     ****************************************************************
% c     *                                                              *
% c     *                 subroutine mm10_formJ22                      *
% c     *                                                              *
% c     *                       written by : mcm                       *
% c     *                                                              *
% c     *                   last modified: 11/27/13                    *
% c     *                                                              *
% c     *     Form the hardening varying with hardening part           *
% c     *                                                              *
% c     ****************************************************************
% c
function [props, np1, n,vec1,vec2,arr1,arr2, stress, tt, J22] = mm10_formJ22(props, np1,...
     n,vec1,vec2,arr1,arr2, stress, tt)
%       use mm10_defs
%       implicit none
% c
%       type(crystal_props) :: props
%       type(crystal_state) :: np1, n
%       double precision, dimension(6) :: stress
%       double precision :: tt, J22
% c
%       double precision :: etau
% c
%       J22 = 0.0;
% c
      if (props.h_type == 1) % voche
        [props, np1, n, stress, tt, etau] = mm10_ehard_voche(props,...
            np1, n, stress, tt);
      elseif (props.h_type == 2) % MTS
        [props, np1, n, stress, tt, etau] = mm10_ehard_mts(props, np1,...
            n, stress, tt);
      elseif (props.h_type == 3) % User
        [props, np1, n, stress, tt, etau] = mm10_ehard_user(props,...
            np1, n, stress, tt, etau);
      elseif (props.h_type == 4) % ORNL
        [props, np1, n,vec1,vec2,arr1,arr2, stress, tt, etau] = mm10_ehard_ornl(props, np1,...
            n,vec1,vec2,arr1,arr2, stress, tt);
      elseif (props.h_type == 5) % AADP
        [props, np1, n, stress, tt, etau] = mm10_ehard_aadp(props, np1,...
            n, stress, tt);
      elseif (props.h_type == 6) % MTS_omar copy
          [props, np1, n, stress, tt, etau] = mm10_ehard_mts_omar(props, np1,...
              n, stress, tt);
      elseif (props.h_type == 7) % MRR
        [props, np1, n,vec1,vec2,arr1,arr2, stress, tt, etau] = mm10_ehard_mrr(props, np1,...
            n,vec1,vec2,arr1,arr2, stress, tt);
      elseif (props.h_type == 9) % Ran_HCP
        [props, np1, n,vec1,vec2,arr1,arr2, stress, tt, etau] = mm10_ehard_DJGM(props, np1,...
            n,vec1,vec2,arr1,arr2, stress, tt);
      elseif (props.h_type == 12)% then % Plugin library
          etau = mm10_plugin_lib(14,props,{props, np1,...
            n, vec1,vec2,arr1,arr2, stress, tt});
      elseif (props.h_type == 14) % halite
        [props, np1, n,vec1,vec2,arr1,arr2, stress, tt, etau] = mm10_ehard_halite(props, np1,...
            n,vec1,vec2,arr1,arr2, stress, tt);
      else
        [props] = mm10_unknown_hard_error(props);
      end
% c
      J22 = etau;
      
%       % Complex numerical derivative, also verified.
%       J22i = zeros(props.num_hard,props.num_hard);
%         h = 1e-7;
% %         for j = 1:props.num_hard
%             for k = 1:props.num_hard
%                 A = tt;
%                 hi = h*tt(k);
%                 A(k) = A(k) + 1i*hi;
%                 [~,~,~,~,~,Ri] = mm10_formR2(props,np1,n,stress,A);
%                 J22i(1:props.num_hard,k) = 1/hi*imag(Ri);
%             end
% %         end
%       J22i;
% c
end