% c
% c     ****************************************************************
% c     *                                                              *
% c     *                 subroutine mm10_init_general_hist            *
% c     *                                                              *
% c     *                       written by : mcm                       *
% c     *                                                              *
% c     *                   last modified: 11/26/13                    *
% c     *                                                              *
% c     *    initialize general GP history (grad Fe, tangent, R)       *
% c     *                                                              *
% c     ****************************************************************
% c
function [history] = mm10_init_general_hist(history)
%        implicit none
%       double precision :: history(72)
%       double precision :: eye(3,3)
% c   
      history(1:63) = 0.0;
% c
      eye = 0.0;
      eye(1,1) = 1.0;
      eye(2,2) = 1.0;
      eye(3,3) = 1.0;
      history(64:72) = reshape(eye,9,1);
% return
end
