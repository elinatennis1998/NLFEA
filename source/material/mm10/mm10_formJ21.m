% c
% c     ****************************************************************
% c     *                                                              *
% c     *                 subroutine mm10_formJ21                      *
% c     *                                                              *
% c     *                       written by : mcm                       *
% c     *                                                              *
% c     *                   last modified: 11/27/13                    *
% c     *                                                              *
% c     *     Form the hardening varying with stress part              *
% c     *                                                              *
% c     ****************************************************************
% c
function [props, np1, n,vec1,vec2,arr1,arr2, stress, tt, J21] = mm10_formJ21(props, np1,...
    n,vec1,vec2,arr1,arr2, stress, tt)
% use mm10_defs
%       implicit none
% c
%       type(crystal_props) :: props
%       type(crystal_state) :: np1, n
%       double precision, dimension(6) :: stress
%       double precision, dimension(6) :: J21
%       double precision :: tt
% c
%       double precision, dimension(6) :: estr
% c
%        J21 = 0.0;
% c

      if (props.h_type == 1) % voche
        [props, np1, n, stress, tt, estr] = mm10_estress_voche(props,...
            np1, n, stress, tt);
      elseif (props.h_type == 2) % MTS
        [props, np1, n, stress, tt, estr] = mm10_estress_mts(props,...
            np1, n, stress, tt);
      elseif (props.h_type == 3) % User
        [props, np1, n, stress, tt, estr] = mm10_estress_user(props,...
            np1, n, stress, tt, estr);
      elseif (props.h_type == 4) % ORNL
        [props, np1, n,vec1,vec2,arr1,arr2, stress, tt, estr] = mm10_estress_ornl(props,...
            np1, n,vec1,vec2,arr1,arr2, stress, tt);
      elseif (props.h_type == 5) % AADP
        [props, np1, n, stress, tt, estr] = mm10_estress_aadp(props,...
            np1, n, stress, tt);
elseif (props.h_type == 6) % MTS_omar copy
    [props, np1, n, stress, tt, estr] = mm10_estress_mts_omar(props,...
        np1, n, stress, tt);
      elseif (props.h_type == 7) % MMR
        [props, np1, n,vec1,vec2,arr1,arr2, stress, tt, estr] = mm10_estress_mrr(props,...
            np1, n,vec1,vec2,arr1,arr2, stress, tt);
      elseif (props.h_type == 9) % Ran_HCP
        [props, np1, n,vec1,vec2,arr1,arr2, stress, tt, estr] = mm10_estress_DJGM(props,...
            np1, n,vec1,vec2,arr1,arr2, stress, tt);
      elseif (props.h_type == 12)% then % Plugin library
          estr = mm10_plugin_lib(13,props,{props, np1,...
            n, vec1,vec2,arr1,arr2, stress, tt});
      elseif (props.h_type == 14) % halite
        [props, np1, n,vec1,vec2,arr1,arr2, stress, tt, estr] = mm10_estress_halite(props,...
            np1, n,vec1,vec2,arr1,arr2, stress, tt);
      else
        [props] = mm10_unknown_hard_error(props);
      end
% c
      J21 = -estr;
      
%       % Complex numerical calculation of derivative; VERIFIED for MTS and
%       % mrr models
%       J22i = zeros(props.num_hard,6);
%         h = 1e-14;
% %         for j = 1:props.num_hard
%             for k = 1:6
%                 A = stress;
%                 hi = h*abs(stress(k));
%                 A(k) = A(k) + 1i*hi;
%                 [~,~,~,~,~,Ri] = mm10_formR2(props,np1,n,A,tt);
%                 J22i(1:props.num_hard,k) = 1/hi*imag(Ri);
%             end
% %         end
%       J22i;
% c
%       return
end