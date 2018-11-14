% c
% c     ****************************************************************
% c     *                                                              *
% c     *                 subroutine mm10_setup_np1                    *
% c     *                                                              *
% c     *                       written by : mcm                       *
% c     *                                                              *
% c     *                   last modified: 11/26/13                    *
% c     *                                                              *
% c     *    Initialize the state np1 structure with new strain/temp/  *
% c     *     time increment                                           *
% c     *                                                              *
% c     ****************************************************************
% c
function np1 = mm10_setup_np1...
    (Rur, dstrain, dt, T, step, elem, gp)
%       use mm10_defs
%       implicit none
% c
%       double precision, dimension(9) :: Rur
%       double precision, dimension(6) :: dstrain
%       double precision :: dt, T
%       integer :: step, elem, gp
%       type(crystal_state) :: np1
% c
np1 = crystal_state;
      np1.R = reshape(Rur(1:9), [3,3]);
      np1.D = dstrain(1:6);
      np1.temp = T;
      np1.tinc = dt;
      np1.step = step;
      np1.elem = elem;
      np1.gp = gp;
% c