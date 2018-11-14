% c
% c     ****************************************************************
% c     *                                                              *
% c     *                 subroutine mm10_setup                        *
% c     *                                                              *
% c     *                       written by : mcm                       *
% c     *                                                              *
% c     *                   last modified: 11/26/13                    *
% c     *                                                              *
% c     *     setup hardening for a particular state np1               *
% c     *                                                              *
% c     ****************************************************************
% c
function [props, np1, n] = mm10_setup(props, np1, n,subcycle)
%       use mm10_defs
%       implicit none
% c
%       type(crystal_props) :: props
%       type(crystal_state) :: np1, n
% c
%       integer :: i, j, k, s, t
%       integer, external :: l_c
%       double precision, dimension(6,6) :: RE
%       double precision, dimension(3,3) :: RW, RWC
curv = zeros(3,9);
permut = [     0     0     0
     0     0     1
     0    -1     0
     0     0    -1
     0     0     0
     1     0     0
     0     1     0
    -1     0     0
     0     0     0];
% tm = zeros(3,1);%       double precision, dimension(3) :: cn, 
%       double precision :: alpha
% c           The alpha for geometric hardening
alpha=1.0/3.0;
% c
% c           Calculate effective strain increment
      np1.dg = sqrt(2.0/3.0*(np1.D(1:3)*np1.D(1:3)'...
          +0.5*np1.D(4:6)*np1.D(4:6)'));
% c
% c           Calculate the current m and q tensors
% c           Yes, these are supposed to be transposes.  We actually need the
% c           backwards rotation from the lattice state.
      RE = mm10_RT2RVE(transpose(n.Rp));
      RW = mm10_RT2RVW(transpose(n.Rp));
      RWC = ...
          mm10_RT2RVW(np1.R*transpose(n.Rp));
%       for i=1:props.nslip
        np1.ms = RE*props.ms;
        np1.qs = RW*props.qs;
        np1.qc = RWC*props.qs;
%       end
% c
      if (props.h_type == 1) % voche
        [props, np1, n] = mm10_setup_voche(props, np1, n);
      elseif (props.h_type == 2) % MTS
        [props, np1, n] = mm10_setup_mts(props, np1, n);
      elseif (props.h_type == 3) % User
        [props, np1, n] = mm10_setup_user(props, np1, n);
      elseif (props.h_type == 4) % ORNL
        [props, np1, n] = mm10_setup_ornl(props, np1, n);
      elseif (props.h_type == 5) % AADP
        [props, np1, n] = mm10_setup_aadp(props, np1, n);
      elseif (props.h_type == 6) % MTS_ omar copy
        [props, np1, n] = mm10_setup_mts_omar(props, np1, n);
      elseif (props.h_type == 7) % MRR
        [props, np1, n] = mm10_setup_mrr(props, np1, n);
      elseif (props.h_type == 8) % DYS
        [props, np1, n] = mm10_setup_dys(props, np1, n);
      elseif (props.h_type == 9) % Ran_HCP
        [props, np1, n] = mm10_setup_DJGM(props, np1, n);
      elseif (props.h_type == 12) % Plugin library
          AllOutput = mm10_plugin_lib(3,props,{props, np1, n});
          props = AllOutput{1};
          np1 = AllOutput{2};
          n = AllOutput{3};
      elseif (props.h_type == 14) % halite
        [props, np1, n] = mm10_setup_halite(props, np1, n);
      else
        [props] = mm10_unknown_hard_error(props);
      end

% c Compute quantities related to the back stress

if (props.h_type == 1 || props.h_type == 2 || props.h_type == 3 || props.h_type == 6)

% c           Calculate the tau lambdas for geometric hardening
% c           Lattice curvature
% c     
%       curv = 0.0;
        for i = 1:3
          for k = 1:3
            for s = 1:3
                l = (s-1)*3+k;
              curv(i,l) = n.gradFeinv(i,k,s) - n.gradFeinv(i,s,k);
            end
          end
        end
% c           Use Acharya's large strain definition of lambda
        firstterm = props.k_0*props.burgers*(alpha*np1.mu_harden)^2.0/(2.0*props.theta_0);
        if subcycle > 1
        cn = np1.R*(transpose(n.Rp)*props.ns(1:3,:));
        else
        cn = n.R*(transpose(n.Rp)*props.ns(1:3,:));
        end
        tm = 0.5*curv*permut*cn;
        vals = diag(sqrt(tm'*tm));
        np1.tau_l = firstterm*vals(1:props.nslip);

      elseif (props.h_type == 4) % ORNL

% c           REPLACE this with the back stress calculation
% c     
%       curv = 0.0;
        for i = 1:3
          for k = 1:3
            for s = 1:3
                l = (s-1)*3+k;
              curv(i,l) = n.gradFeinv(i,k,s) - n.gradFeinv(i,s,k);
            end
          end
        end
% c           Use Acharya's large strain definition of lambda
        firstterm = props.k_0*props.burgers*(alpha*np1.mu_harden)^2.0/(2.0*props.theta_0);
        if subcycle > 1
        cn = np1.R*(transpose(n.Rp)*props.ns(1:3,:));
        else
        cn = n.R*(transpose(n.Rp)*props.ns(1:3,:));
        end
        tm = 0.5*curv*permut*cn;
        vals = diag(sqrt(tm'*tm));
        np1.tau_l = firstterm*vals(1:props.nslip);

      else
        [props] = mm10_unknown_hard_error(props);
      end
% 
%       return
end
% 
% c     Helper for the above
% function l_c = l_c(i,j,k)
% %         implicit none
% %         integer :: i,j,k
% %         integer :: l_c
% % 
%         l_c = (i-j)*(j-k)*(k-i)/2;
% end