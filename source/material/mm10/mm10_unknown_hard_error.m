% c
% c     ****************************************************************
% c     *                                                              *
% c     *                 subroutine mm10_unknown_hard_error           *
% c     *                                                              *
% c     *                       written by : mcm                       *
% c     *                                                              *
% c     *                   last modified: 1/27/14                     *
% c     *                                                              *
% c     *     A common error message for the general hardening setup   *
% c     *                                                              *
% c     ****************************************************************
% c
function [props] = mm10_unknown_hard_error(props)
%       use mm10_defs
%       implicit none
% c
%       type(crystal_props) :: props;
% c
%       write(props%out,101) props%h_type
%  101  format(
%      &      10x,'>> Error: unknown hardening type ', 'i6', '.',
%      &    /,10x, 'Aborting...')
%      call die_gracefully
end