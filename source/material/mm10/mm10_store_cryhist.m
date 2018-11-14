% c
% c     ****************************************************************
% c     *                                                              *
% c     *                 subroutine mm10_store_cryhist                *
% c     *                                                              *
% c     *                       written by : mcm                       *
% c     *                                                              *
% c     *                   last modified: 11/26/13                    *
% c     *                                                              *
% c     *    Copy the state np1 struct to the history                  *
% c     *                                                              *
% c     ****************************************************************
% c
function history = mm10_store_cryhist(props,np1,n,history)
%       use mm10_defs
        [max_slip_sys,max_uhard] = maxparamsCP;
%       implicit integer(a-z)
% c
%       double precision, dimension(25+max_uhard) :: history
%       type(crystal_props) :: props
%       type(crystal_state) :: np1, n
% c
      history(1:6) = np1.stress;
      history(7:9) = np1.euler_angles;
      history(10:18) = reshape(np1.Rp,9,1);
      history(18+1:24) = np1.D;
      history(24+1:30) = np1.elaststrain;
      history(30+1:30+max_slip_sys) = np1.slip_incs(1:max_slip_sys);
      history(30+max_slip_sys+1:30+max_slip_sys+props.num_hard) = ...
              np1.tau_tilde(1:props.num_hard);
      history(30+max_slip_sys+max_uhard+1:30+max_slip_sys+2*max_uhard) = ...
              np1.u(1:max_uhard);
      history(30+max_slip_sys+2*max_uhard+1:30+max_slip_sys+2*max_uhard+props.num_hard) = ...
              np1.tt_rate(1:props.num_hard);
% c
%       return
end