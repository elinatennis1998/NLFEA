% c
% c     ****************************************************************
% c     *                                                              *
% c     *                 subroutine mm10_dump_state                   *
% c     *                                                              *
% c     *                       written by : mcm                       *
% c     *                                                              *
% c     *                   last modified: 2/25/14                     *
% c     *                                                              *
% c     *     Dump a state variable to try to figure out what's going  *
% c     *     on with these random errors.                             *
% c     *                                                              *
% c     ****************************************************************
% c
function [state] = mm10_dump_state(state)
%       use mm10_defs
%       implicit none
% c
%       type(crystal_state) :: state
% c
%       write(*,*)
%       write(*,*) "step, element, gp"
%       write(*,*) state%step, state%elem, state%gp
%       write(*,*) "R"
%       write(*,*) transpose(state%R)
%       write(*,*) "Rp"
%       write(*,*) transpose(state%Rp)
%       write(*,*) "stress"
%       write(*,*) state%stress
%       write(*,*) "D"
%       write(*,*) state%D
%       write(*,*) "angles"
%       write(*,*) state%euler_angles
%       write(*,*) "tau_l"
%       write(*,*) state%tau_l
%       write(*,*) "slip_incs"
%       write(*,*) state%slip_incs
%       write(*,*) "gradFeinv"
%       write(*,*) state%gradFeInv
%       write(*,*) "tangent"
%       write(*,*) state%tangent
%       write(*,*) "tau_tilde, temp, tinc, dg, tau_v, tau_y"
%       write(*,*) state%tau_tilde, state%temp, state%tinc, state%dg,
%      &            state%tau_v, state%tau_y
%       write(*,*) "mu_harden, work_inc, p_work_inc, p_strain_inc"
%       write(*,*) state%mu_harden, state%work_inc, state%p_work_inc,
%      &            state%p_strain_inc
%       write(*,*) "ms"
%       write(*,*) transpose(state%ms)
%       write(*,*) "qs"
%       write(*,*) transpose(state%qs)
%       write(*,*) "qc"
%       write(*,*) transpose(state%qc)
%       write(*,*) "u"
%       write(*,*) state%u
%       write(*,*)
% c
end