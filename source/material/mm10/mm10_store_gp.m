% c
% c     ****************************************************************
% c     *                                                              *
% c     *                 subroutine mm10_store_gp                     *
% c     *                                                              *
% c     *                       written by : mcm                       *
% c     *                                                              *
% c     *                   last modified: 11/26/13                    *
% c     *                                                              *
% c     *    Store all required gauss point data                       *
% c     *                                                              *
% c     ****************************************************************
% c
function [stress_out, tang_out, ...
    slip_np1, u_new,...
    R_out] = mm10_store_gp(stress_in, tang_in, slip_inc,...
    slip_n, t_work, p_work, p_strain, u_old,...
    R_in)
%       use mm10_defs
%       implicit integer(a-z)
% c
%       double precision, dimension(6) :: stress_in, stress_out
%       double precision, dimension(9) :: R_in, R_out
%       double precision, dimension(6,6) :: tang_in
%       double precision, dimension(36) :: tang_out
%       double precision, dimension(max_slip_sys) :: slip_inc, slip_n, 
%      &      slip_np1
%       double precision :: t_work, p_work, p_strain
%       double precision, dimension(3) :: u_old, u_new
% c
      stress_out(1:6) = stress_in(1:6);
      tang_out(1:36) = reshape(tang_in, 36, 1);
      slip_np1 = slip_n + slip_inc';
% c
      u_new(1) = u_old(1) + t_work;
      u_new(2) = u_old(2) + p_work;
      u_new(3) = u_old(3) + p_strain;
% c
      R_out(1:9) = R_in(1:9);
%       return
end