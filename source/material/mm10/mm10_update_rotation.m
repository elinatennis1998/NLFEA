% c
% c     ****************************************************************
% c     *                                                              *
% c     *                 subroutine mm10_update_rotation              *
% c     *                                                              *
% c     *                       written by : mcm                       *
% c     *                                                              *
% c     *                   last modified: 1/14/14                     *
% c     *                                                              *
% c     *     Update the plastic rotation                              *
% c     *                                                              *
% c     ****************************************************************
% c
function [props,np1,n, vec1, vec2] = mm10_update_rotation(props,np1,n, vec1, vec2)
%       use mm10_defs
%       implicit none
% c
%        type(crystal_props) :: props
%        type(crystal_state) :: np1, n
% c
% c
[props,np1,n, vec1, vec2,np1.stress,np1.tau_tilde,wbarp] = mm10_form_wbarp(props,np1...
    ,n, vec1, vec2,np1.stress,np1.tau_tilde);
wbarp_full = mm10_WV2WT(wbarp);
% c
expw = mm10_expmw3x3(wbarp_full);
% c
np1.Rp = (expw*n.Rp);
% c
end