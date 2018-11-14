% c
% c     ****************************************************************
% c     *                                                              *
% c     *                 subroutine mm10_num_tangent                  *
% c     *                                                              *
% c     *                       written by : mcm                       *
% c     *                                                              *
% c     *                   last modified: 11/27/13                    *
% c     *                                                              *
% c     *     Calculate the consistent tangent after a converged       *
% c     *     stress update with a numerical derivative (for check)    *
% c     *                                                              *
% c     ****************************************************************
% c
function [props, np1, n] = mm10_num_tangent(props, np1, n)
%       use mm10_defs
%       implicit none
% c
%       type(crystal_props) :: props
%       type(crystal_state) :: np1, n
% c
%       type(crystal_state) :: np1star
%       integer :: i
%       double precision, dimension(6) :: dD
%       double precision :: eps
%       logical :: cut
 eps = 1e-8;
% c
%       np1.tangent = 0.0
      for i=1:6
%             dD = 0.0
            dD(i) = eps;
            [reshape(np1.R, 9), np1.D +dD, np1.tinc, np1.temp,...
                np1.step, np1.elem, np1.gp, np1star] = ...
                mm10_setup_np1(reshape(np1.R, 9),...
                np1.D+dD, np1.tinc, np1.temp, np1.step, np1.elem,...
                np1.gp, np1star);
            [props, np1star, n, cut] = mm10_solve_strup(props,...
                np1star, n, cut);
%             if (cut) then
%               write (*,*) "Numerical tangent failed"
%               call die_gracefully
%             end if
            np1.tangent(:,i) = (np1star.stress - np1.stress) / eps;
      end
% 
%       return
end