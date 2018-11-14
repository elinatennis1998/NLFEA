% c
% c     ****************************************************************
% c     *                                                              *
% c     *                 subroutine mm10_tangent                      *
% c     *                                                              *
% c     *                       written by : mcm                       *
% c     *                                                              *
% c     *                   last modified: 11/27/13                    *
% c     *                                                              *
% c     *     Calculate the consistent tangent after a converged       *
% c     *     stress update.                                           *
% c     *                                                              *
% c     ****************************************************************
% c
function [props, np1, n] = mm10_tangent(props, np1, n, J)
%       use mm10_defs
%       implicit none
% c
%       type(crystal_props) :: props
%       type(crystal_state) :: np1, n
% c
%       double precision, dimension(6,6) :: J11, JJ, JR, JA, JB
%       double precision, dimension(6) :: J12, J21, d_mod, ed, d_barp, tw
%       double precision, dimension(3) :: w_p
%       double precision :: J22, alpha

% % This stuff is for trying to improve the J for subincrements from the
% % solve_nsteps version of material integration
% [~,max_uhard] = maxparamsCP;
% % Intermediate vectors and arrays to help with pre-computing
% vec1 = zeros(max_uhard,1);
% vec2 = zeros(max_uhard,1);
% arr1 = zeros(max_uhard,max_uhard);
% arr2 = zeros(max_uhard,max_uhard);
% % c
%  [props, np1, n, np1.stress, np1.tau_tilde, vec1,vec2] = mm10_formvecs(props,np1,n,np1.stress,... 
%  np1.tau_tilde,vec1,vec2);
%         [props, np1, n,vec1,vec2,arr1,arr2, np1.stress, np1.tau_tilde, Jb] = mm10_formJ(props, np1,...
%             n,vec1,vec2,arr1,arr2, np1.stress, np1.tau_tilde);
%       J11 = real(Jb(1:6,1:6));
%       J12 = real(Jb(1:6,7:props.num_hard+6));
%       J21 = real(Jb(7:props.num_hard+6,1:6));
%       J22 = real(Jb(7:props.num_hard+6,7:props.num_hard+6));
      % Don't recompute J anymore
      J11 = real(J(1:6,1:6));
      J12 = real(J(1:6,7:props.num_hard+6));
      J21 = real(J(7:props.num_hard+6,1:6));
      J22 = real(J(7:props.num_hard+6,7:props.num_hard+6));
%       Ja = (J+Jb)/2;
%       J11 = real(Ja(1:6,1:6));
%       J12 = real(Ja(1:6,7:props.num_hard+6));
%       J21 = real(Ja(7:props.num_hard+6,1:6));
%       J22 = real(Ja(7:props.num_hard+6,7:props.num_hard+6));
      
% c
      if (props.h_type == 1) % voche
        [props, np1, n, np1.stress, np1.tau_tilde, ed] = ...
            mm10_ed_voche(props, np1, n, np1.stress,...
            np1.tau_tilde);
      elseif (props.h_type == 2) % MTS
        [props, np1, n, np1.stress, np1.tau_tilde, ed] = ...
            mm10_ed_mts(props, np1, n, np1.stress, np1.tau_tilde);
      elseif (props.h_type == 3) % User
        [props, np1, n, np1.stress, np1.tau_tilde, ed] = ...
            mm10_ed_user(props, np1, n, np1.stress, np1.tau_tilde, ed);
      elseif (props.h_type == 4) % ORNL
        [props, np1, n, np1.stress, np1.tau_tilde, ed] = ...
            mm10_ed_ornl(props, np1, n, np1.stress, np1.tau_tilde);
      elseif (props.h_type == 5) % AADP
        [props, np1, n, np1.stress, np1.tau_tilde, ed] = ...
            mm10_ed_aadp(props, np1, n, np1.stress, np1.tau_tilde);
elseif (props.h_type == 6) % MTS_omar copy
    [props, np1, n, np1.stress, np1.tau_tilde, ed] = ...
        mm10_ed_mts_omar(props, np1, n, np1.stress, np1.tau_tilde);
      elseif (props.h_type == 7) % MRR
        [props, np1, n, np1.stress, np1.tau_tilde, ed] = ...
            mm10_ed_mrr(props, np1, n, np1.stress, np1.tau_tilde);
      elseif (props.h_type == 8) % DYS
        [props, np1, n, np1.stress, np1.tau_tilde, ed] = ...
            mm10_ed_dys(props, np1, n, np1.stress, np1.tau_tilde);
      elseif (props.h_type == 9) % Ran_HCP
        [props, np1, n, np1.stress, np1.tau_tilde, ed] = ...
            mm10_ed_DJGM(props, np1, n, np1.stress, np1.tau_tilde);
      elseif (props.h_type == 12) % Plugin Library
        ed = mm10_plugin_lib(19,props,{props, np1,...
        n, np1.stress, np1.tau_tilde});
      elseif (props.h_type == 14) % halite
        [props, np1, n, np1.stress, np1.tau_tilde, ed] = ...
            mm10_ed_halite(props, np1, n, np1.stress, np1.tau_tilde);
      else
        [props] = mm10_unknown_hard_error(props);
      end
% c
% c
%        np1.tangent = 0.0;
% c
      JJ = J11;
      alpha = -inv(J22);
      JJ = J12*alpha*J21 + JJ;
%       JJ = inv(JJ);
% c
      symtq = mm10_symSW(np1.stress,...
            np1.qc);
      % Generalization of CP model implementation for other slip rate
      % equations, requiring other forms of d_gamma/d_d
      % Vector dgammadtt should be n_slip x 6
      if (props.h_type == 1) % voche
          [props, np1, n, np1.stress, np1.tau_tilde, np1.D, dgammadd] = mm10_dgdh_voche(props,...
            np1, n, np1.stress, np1.tau_tilde, np1.D);
      elseif (props.h_type == 2) % MTS
          [props, np1, n, np1.stress, np1.tau_tilde, np1.D, dgammadd] = mm10_dgdd_mts(props, np1,...
            n, np1.stress, np1.tau_tilde, np1.D);
      elseif (props.h_type == 3) % User
          [props, np1, n, np1.stress, np1.tau_tilde, dgammadd] = mm10_dgdd_user(props,...
            np1, n, np1.stress, np1.tau_tilde, np1.D);
      elseif (props.h_type == 4) % ORNL
          [props, np1, n, np1.stress, np1.tau_tilde, np1.D, dgammadd] = mm10_dgdd_ornl(props, np1,...
            n, np1.stress, np1.tau_tilde, np1.D);
      elseif (props.h_type == 5) % AADP
          [props, np1, n, np1.stress, np1.tau_tilde, np1.D, dgammadd] = mm10_dgdd_aadp(props, np1,...
            n, np1.stress, np1.tau_tilde, np1.D);
elseif (props.h_type == 6) % MTS_omar copy
    [props, np1, n, np1.stress, np1.tau_tilde, np1.D, dgammadd] = mm10_dgdd_mts_omar(props, np1,...
        n, np1.stress, np1.tau_tilde, np1.D);
      elseif (props.h_type == 7) % MRR
          [props, np1, n, np1.stress, np1.tau_tilde, np1.D, dgammadd] = mm10_dgdd_mrr(props, np1,...
            n, np1.stress, np1.tau_tilde, np1.D);
      elseif (props.h_type == 8) % DYS
          [props, np1, n, np1.stress, np1.tau_tilde, np1.D, dgammadd] = mm10_dgdd_dys(props, np1,...
            n, np1.stress, np1.tau_tilde, np1.D);
      elseif (props.h_type == 9) % Ran_HCP
          [props, np1, n, np1.stress, np1.tau_tilde, np1.D, dgammadd] = mm10_dgdd_DJGM(props, np1,...
            n, np1.stress, np1.tau_tilde, np1.D);
      elseif (props.h_type == 12) % Plugin Library
            dgammadd = mm10_plugin_lib(20,props,{props, np1,...
            n, np1.stress, np1.tau_tilde, np1.D});
      elseif (props.h_type == 14) % halite
          [props, np1, n, np1.stress, np1.tau_tilde, np1.D, dgammadd] = mm10_dgdd_halite(props, np1,...
            n, np1.stress, np1.tau_tilde, np1.D);
      else
          [props] = mm10_unknown_hard_error(props);
      end
      JA = (props.stiffness*np1.ms+ 2.0*symtq)*dgammadd;

% c
%         JB = 0.0; % kept in program since JB is used in the DGER function
%       alpha = inv(J22);
%       JB = alpha*J12*ed';
      JB = (J12/J22)*ed';
% c
      JR = props.stiffness - JA - JB;
% c
%       np1.tangent = JJ*JR;
      np1.tangent = JJ\JR;
% c
%       return
end