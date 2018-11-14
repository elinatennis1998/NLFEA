% c
% c     ****************************************************************
% c     *                                                              *
% c     *                 subroutine mm10_init_slip_hist               *
% c     *                                                              *
% c     *                       written by : mcm                       *
% c     *                                                              *
% c     *                   last modified: 11/26/13                    *
% c     *                                                              *
% c     *    initialize the slip totals (output variable)              *
% c     *                                                              *
% c     ****************************************************************
function [history] = mm10_init_slip_hist(history)
max_slip_sys = 12;
%       implicit integer (a-z)
% $add param_def
%       double precision :: history(max_slip_sys)
% c
history(1:max_slip_sys) = 0.0;
%
%return
end