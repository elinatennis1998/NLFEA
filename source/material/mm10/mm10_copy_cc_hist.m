% c
% c     ****************************************************************
% c     *                                                              *
% c     *                 subroutine mm10_copy_cc_hist                 *
% c     *                                                              *
% c     *                       written by : mcm                       *
% c     *                                                              *
% c     *                   last modified: 11/26/13                    *
% c     *                                                              *
% c     *    Initialize the state n structure                          *
% c     *                                                              *
% c     ****************************************************************
% c
function [history, gradfe, R, props, n] = mm10_copy_cc_hist(history,...
    gradfe, R, props,backstress_omar)
[max_slip_sys,max_uhard] = maxparamsCP;
%       use mm10_defs
%       implicit integer(a-z)
% c
%       double precision, dimension(25+max_uhard) :: history
%       double precision, dimension(27) :: gradfe
%       double precision, dimension(9) :: R
%       type(crystal_props) :: props
%       type(crystal_state) :: n
n=crystal_state;
% c           Not provided, could be useful...
      n.R(1:3,1:3) = reshape(R, [3,3]);
      n.gradFeinv(1:3,1:3,1:3) = reshape(gradfe, [3,3,3]);
      n.temp = 0.0;
% c           Only used at n+1
      n.tinc = 0.0;
      n.dg = 0.0;
      n.tau_v = 0.0;
      n.tau_y = 0.0;
      n.mu_harden = 0.0;
% c
      n.stress = history(1:6);
      n.euler_angles(1:3) = history(7:9);
% c     
      n.Rp = reshape(history(10:18), [3,3]);
      n.D(1:6) = history(18+1:24);
      n.elaststrain(1:6) = history(24+1:30);
      n.slip_incs(1:max_slip_sys) = history(30+1:30+max_slip_sys);
      n.tau_tilde(1:props.num_hard) = history(30+max_slip_sys+1:30+max_slip_sys+props.num_hard);
% c
      n.u(1:max_uhard) = history(30+max_slip_sys+max_uhard+1:30+max_slip_sys+2*max_uhard);
      n.tt_rate(1:props.num_hard) = history(30+max_slip_sys+2*max_uhard+1: ...
                                            30+max_slip_sys+2*max_uhard+props.num_hard);
      n.backstress_omar = backstress_omar;                                      
% c      
%       return
end