% c
% c     ****************************************************************
% c     *                                                              *
% c     *                 subroutine mm10_formJ12                      *
% c     *                                                              *
% c     *                       written by : mcm                       *
% c     *                                                              *
% c     *                   last modified: 11/27/13                    *
% c     *                                                              *
% c     *     Form the stress varying with hardening part              *
% c     *                                                              *
% c     ****************************************************************
% c
function [props, np1, n,vec1,vec2,arr1,arr2, stress, tt, J12] = mm10_formJ12(props, np1,...
     n,vec1,vec2,arr1,arr2, stress, tt)
%       use mm10_defs
%       implicit none
% c
%       type(crystal_props) :: props
%       type(crystal_state) :: np1, n
%       double precision, dimension(6) :: stress
%       double precision, dimension(6) :: J12
%       double precision :: tt
% c
%       double precision, dimension(6) :: dbar, symTW
%       double precision, dimension(3) :: wp
% c
%       J12 = 0.0;
% c
      symtq = mm10_symSW(stress,...
            np1.qc);
      % Generalization of CP model implementation for other slip rate
      % equations, requiring other forms of d_gamma/d_hardening
      % Vector dgammadtt should be n_slip x n_hardening, which is the
      % derivative of slip rate alpha with respect to hardening variable
      % beta.
      if (props.h_type == 1) % voche
          [props, np1, n, stress, tt, dgammadtt] = mm10_dgdh_voche(props,...
            np1, n, stress, tt);
      elseif (props.h_type == 2) % MTS
          [props, np1, n, stress, tt, dgammadtt] = mm10_dgdh_mts(props, np1,...
            n, stress, tt);
      elseif (props.h_type == 3) % User
          [props, np1, n, stress, tt, dgammadtt] = mm10_dgdh_user(props,...
            np1, n, stress, tt);
      elseif (props.h_type == 4) % ORNL
%           [props, np1, n, stress, tt, dgammadtt] = mm10_dgdh_ornl(props, np1,...
%             n, stress, tt);
           dgammadtt = arr2(1:props.nslip,1:props.num_hard);
      elseif (props.h_type == 5) % AADP
          [props, np1, n, stress, tt, dgammadtt] = mm10_dgdh_aadp(props, np1,...
            n, stress, tt);
elseif (props.h_type == 6) % MTS_omar copy
    [props, np1, n, stress, tt, dgammadtt] = mm10_dgdh_mts_omar(props, np1,...
        n, stress, tt);
      elseif (props.h_type == 7) % MRR
%           [props, np1, n, stress, tt, dgammadtt] = mm10_dgdh_mrr(props, np1,...
%             n, stress, tt);
           dgammadtt = arr2(1:props.nslip,1:props.num_hard);
      elseif (props.h_type == 9) % Ran_HCP
%           [props, np1, n, stress, tt, dgammadtt] = mm10_dgdh_mrr(props, np1,...
%             n, stress, tt);
           dgammadtt = arr2(1:props.nslip,1:props.num_hard);
      elseif (props.h_type == 12)% then % Plugin library
          dgammadtt = mm10_plugin_lib(12,props,{props, np1,...
            n, vec1,vec2,arr1,arr2, stress, tt});
      elseif (props.h_type == 14) % halite
%           [props, np1, n, stress, tt, dgammadtt] = mm10_dgdh_mrr(props, np1,...
%             n, stress, tt);
           dgammadtt = arr2(1:props.nslip,1:props.num_hard);
      else
          [props] = mm10_unknown_hard_error(props);
      end
      J12 = (props.stiffness*np1.ms+ 2.0*symtq)*dgammadtt;
      
%       % Complex numerical derivative, also verified.
%       J22i = zeros(6,props.num_hard);
%         h = 1e-12;
% %         for j = 1:props.num_hard
%             for k = 1:props.num_hard
%                 A = tt;
%                 hi = h*tt(k);
%                 A(k) = A(k) + 1i*hi;
%                 [~,~,~,~,~,Ri] = mm10_formR1(props,np1,n,stress,A);
%                 J22i(1:6,k) = 1/hi*imag(Ri);
%             end
% %         end
%       J22i;
%       % This one agrees with finite difference and analytical
% c
%       return
end