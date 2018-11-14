% c
% c
% c     ****************************************************************
% c     *                                                              *
% c     *                 subroutine mm10_dump_props                   *
% c     *                                                              *
% c     *                       written by : mcm                       *
% c     *                                                              *
% c     *                   last modified: 2/25/14                     *
% c     *                                                              *
% c     *     Dump a props variable to try to figure out what's going  *
% c     *     on with these random errors.                             *
% c     *                                                              *
% c     ****************************************************************
% c
function [props] = mm10_dump_props(props)
%       use mm10_defs
%       implicit none
% c
%       type(crystal_props) :: props
% c
%       write(*,*) "PROPS"   % below all "%" from fortran changed to "."
%       write(*,*) props.rate_n, props.tau_hat_y, props.G_0_y, 
%      &                  props.burgers,
%      &                  props.p_v, props.q_v, props.boltzman, 
%      &                  props.theta_0, props.eps_dot_0_v,
%      &                  props.eps_dot_0_y,
%      &                  props.p_y, props.q_y,
%      &                  props.tau_a, props.tau_hat_v, props.G_0_v,
%      &                  props.k_0, props.mu_0, props.D_0, props.T_0, 
%      &                  props.tau_y, props.tau_v, props.voche_m,
%      &                  props.u1, props.u2, props.u3, props.u4, 
%      &                  props.u5, props.u6
%       write(*,*) transpose(props.g)
%       write(*,*) props.ms
%       write(*,*) props.qs
%       write(*,*) props.ns
%       write(*,*) transpose(props.stiffness)
%       write(*,*) props.angle_type, props.angle_convention, props.nslip,
%      &      props.h_type
% c
end