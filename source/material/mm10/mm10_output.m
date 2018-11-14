% c
% c     ****************************************************************
% c     *                                                              *
% c     *                 subroutine mm10_output                       *
% c     *                                                              *
% c     *                       written by : mcm                       *
% c     *                                                              *
% c     *                   last modified: 11/27/13                    *
% c     *                                                              *
% c     *     Calculate various other user output quantities           *
% c     *                                                              *
% c     ****************************************************************
% c
function [props, np1, n,vec1,vec2] = mm10_output(props, np1, n,vec1,vec2)
%       use mm10_defs
%       implicit none
% c
%       type(crystal_props) :: props
%       type(crystal_state) :: np1, n
% c
%       double precision :: mm10_slipinc
% c
%       integer :: i
%       double precision, dimension(6) :: ee, dbarp, ewwe, ep
%       double precision, dimension(3) :: wp
%       double precision, dimension(6,6) :: S, erot
% c
% c     Store the slip increments
% c
% c
      if (props.h_type == 1) % voche
      for i=1:props.nslip
        np1.slip_incs(i) = mm10_slipinc(props.rate_n, np1.dg, np1.ms(1:6,i), np1.stress, np1.tau_tilde);
      end
      elseif (props.h_type == 2) % MTS
      for i=1:props.nslip
        np1.slip_incs(i) = mm10_slipinc(props.rate_n, np1.dg, np1.ms(1:6,i), np1.stress, np1.tau_tilde);
      end
      elseif (props.h_type == 3) % User
      for i=1:props.nslip
        np1.slip_incs(i) = mm10_slipinc(props.rate_n, np1.dg, np1.ms(1:6,i), np1.stress, np1.tau_tilde);
      end
      elseif (props.h_type == 4) % ORNL
          np1.slip_incs(1:props.nslip) = vec1(1:props.nslip);
      elseif (props.h_type == 5) % AADP
        for i = 1:props.nslip
            np1.slip_incs(i) = mm10_slipinc_aadp(props, np1, n, np1.stress, np1.tau_tilde, i);
        end
elseif (props.h_type == 6) % MTS_omar copy
    for i=1:props.nslip
        np1.slip_incs(i) = mm10_slipinc_omar(props.rate_n, np1.dg, np1.ms(1:6,i), np1.stress, np1.tau_tilde);
    end
      elseif (props.h_type == 7) % MRR
          np1.slip_incs(1:props.nslip) = vec1(1:props.nslip);
      elseif (props.h_type == 8) % DYS
        for i = 1:props.nslip
            np1.slip_incs(i) = mm10_slipinc_dys(props, np1, n, np1.stress, np1.tau_tilde, i);
        end
      elseif (props.h_type == 9) % Ran_HCP
          np1.slip_incs(1:props.nslip) = vec1(1:props.nslip);
      elseif (props.h_type == 12) % Plugin Library
            np1.slip_incs(1:props.nslip) = mm10_plugin_lib(18, props, ...
                {props, np1, n, vec1, vec2, np1.stress, np1.tau_tilde});
      elseif (props.h_type == 14) % halite
          np1.slip_incs(1:props.nslip) = vec1(1:props.nslip);
      else
        [props] = mm10_unknown_hard_error(props);
      end
% c
% c     Call a function to store the Euler angles
% c
      [props,np1,n] = mm10_update_euler_angles(props,np1,n);
% c
% c     Take care of miscellaneous things
% c
      np1.work_inc = dot(np1.stress, np1.D);
% c     Plastic strain and work
      S = props.stiffness;
      S = inv(S);
      erot = mm10_RT2RVE(np1.R);
      ee = (erot*(S*np1.stress'));
      [props, np1, n,vec1,vec2, np1.stress, np1.tau_tilde, dbarp] =...
          mm10_form_dbarp(props, np1, n,vec1,vec2, np1.stress,...
          np1.tau_tilde);
      [props, np1, n,~,~, np1.stress, np1.tau_tilde, wp] =...
          mm10_form_wp(props, np1, n,vec1,vec2, np1.stress, np1.tau_tilde);
      ewwe = mm10_symSW(ee, wp);
      np1.elaststrain = ee; % CONFIRM THIS
% c
      ep = dbarp + ewwe;
      np1.p_strain_inc = sqrt(2.0/3.0*(dot(ep(1:3),...
          ep(1:3))+0.5*dot(ep(4:6),ep(4:6))));
      np1.p_work_inc = dot(np1.stress, ep);
% c
%       return
end