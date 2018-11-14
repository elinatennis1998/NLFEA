% c
% c     ****************************************************************
% c     *                                                              *
% c     *                 subroutine mm10_init_cc_props                *
% c     *                                                              *
% c     *                       written by : mcm                       *
% c     *                                                              *
% c     *                   last modified: 11/26/13                    *
% c     *                                                              *
% c     *    Copy properties over from local_work into the update      *
% c     *    structure                                                 *
% c     *                                                              *
% c     ****************************************************************
% c
function cc_props = ...
    mm10_init_cc_props(inc_props, atype, aconv)
%       use mm10_defs
[max_slip_sys,max_uhard] = maxparamsCP;
%       implicit integer (a-z)
% $add include_sig_up
%       integer :: atype, aconv
%       type(crystal_properties) :: inc_props
%       type(crystal_props) :: cc_props
cc_props = inc_props;
% c
      cc_props.angle_type = atype;
      cc_props.angle_convention = aconv;
% c
%       return
end