% c
% c     ****************************************************************
% c     *                                                              *
% c     *                 subroutine mm10_formJ11                      *
% c     *                                                              *
% c     *                       written by : mcm                       *
% c     *                                                              *
% c     *                   last modified: 11/26/13                    *
% c     *                                                              *
% c     *     Form the stress varying with stress part                 *
% c     *                                                              *
% c     ****************************************************************
% c
function [props, np1, n,vec1,vec2,arr1,arr2, stress, tt, J11] = mm10_formJ11(props, np1,...
    n,vec1,vec2,arr1,arr2, stress, tt)
%       use mm10_defs
%       implicit none
% c
%       type(crystal_props) :: props
%       type(crystal_state) :: np1, n
%       double precision, dimension(6) :: stress
%       double precision, dimension(6,6) :: J11
%       double precision :: tt
% c
%       double precision :: mm10_rs
% c
%       double precision, dimension(6) :: symtq
%       double precision, dimension(3) :: wp
%       double precision, dimension(6,6) :: Iw
%       double precision :: rs, alpha
%       integer :: i
% c
% c
      symtq = mm10_symSW(stress,...
            np1.qc);
      % Generalization of CP model implementation for other slip rate
      % equations, requiring other forms of d_gamma/d_tau
      % Vector dgammadtau should be 1 x n_slip
      if (props.h_type == 1) % voche
          [props, np1, n, stress, tt, dgammadtau] = mm10_dgdt_voche(props,...
            np1, n, stress, tt);
      elseif (props.h_type == 2) % MTS
          [props, np1, n, stress, tt, dgammadtau] = mm10_dgdt_mts(props, np1,...
            n, stress, tt);
      elseif (props.h_type == 3) % User
          [props, np1, n, stress, tt, dgammadtau] = mm10_dgdt_user(props,...
            np1, n, stress, tt);
      elseif (props.h_type == 4) % ORNL
%           [props, np1, n, stress, tt, dgammadtau] = mm10_dgdt_ornl(props, np1,...
%             n, stress, tt);
          dgammadtau = arr1(1:props.num_hard,1)';
      elseif (props.h_type == 5) % AADP
          [props, np1, n, stress, tt, dgammadtau] = mm10_dgdt_aadp(props, np1,...
            n, stress, tt);
elseif (props.h_type == 6) % MTS_omar copy
    [props, np1, n, stress, tt, dgammadtau] = mm10_dgdt_mts_omar(props, np1,...
        n, stress, tt);
      elseif (props.h_type == 7) % MRR
%           [props, np1, n, stress, tt, dgammadtau] = mm10_dgdt_mrr(props, np1,...
%             n, stress, tt);
          dgammadtau = arr1(1:props.num_hard,1)';
      elseif (props.h_type == 9) % Ran_HCP
%           [props, np1, n, stress, tt, dgammadtau] = mm10_dgdt_mrr(props, np1,...
%             n, stress, tt);
          dgammadtau = arr1(1:props.num_hard,1)';
      elseif (props.h_type == 12)% then % Plugin library
          dgammadtau = mm10_plugin_lib(11,props,{props, np1,...
            n, vec1,vec2,arr1,arr2, stress, tt});
      elseif (props.h_type == 14) % halite
%           [props, np1, n, stress, tt, dgammadtau] = mm10_dgdt_mrr(props, np1,...
%             n, stress, tt);
          dgammadtau = arr1(1:props.num_hard,1)';
      else
          [props] = mm10_unknown_hard_error(props);
      end
      J11 = (props.stiffness*np1.ms+ 2.0*symtq)*((ones(6,1)*dgammadtau).*np1.ms)';
% c
      [props, np1, n,vec1,vec2, stress, tt, wp] = mm10_form_wp(props,...
          np1, n,vec1,vec2, stress, tt);
      Iw = mm10_IW(wp);
      J11 = J11 + Iw;
% c
%	add J of creep strain due to pressure precipitation (Ran)
      if ( props.h_type==14 && abs(props.cp_031-1)<1.0e-5 )
          J11 = mm10_halite_formJpp( props, J11, np1.tinc );
      end
%
      for i=1:6
        J11(i,i) = J11(i,i) + 1.0;
      end
      
%       % Complex numerical derivative, also verified.
%       J22i = zeros(6,6);
%         h = 1e-12;
% %         for j = 1:props.num_hard
%             for k = 1:6
%                 A = stress;
%                 hi = h*abs(stress(k));
%                 A(k) = A(k) + 1i*hi;
%                 [~,~,~,~,~,Ri] = mm10_formR1(props,np1,n,A,tt);
%                 J22i(1:6,k) = 1/hi*imag(Ri);
%             end
% %         end
%       J22i;
%       % Now this one works too for Mark and for Roters
% c
%       return
end