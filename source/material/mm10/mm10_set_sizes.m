% c
% c     ****************************************************************
% c     *                                                              *
% c     *                 subroutine mm10_set_sizes                    *
% c     *                                                              *
% c     *                       written by : mcm                       *
% c     *                                                              *
% c     *                   last modified: 3/22/12                     *
% c     *                                                              *
% c     *    called by warp3d for each material model to obtain        *
% c     *    various sizes of data for the model                       *
% c     *                                                              *
% c     ****************************************************************
% c
function [size_data, local_el] = mm10_set_sizes(size_data, local_el)
%       use main_data, only: imatprp
max_uhard = 20;
max_slip_sys = 10;
%       implicit integer (a-z)
% $add common.main
%       dimension size_data(*)
%       integer :: local_el, matnum, ncrystals
% c
% c        size_data(1)  :  no. of words of history data for each 
% c                         integration point
% c
% c
% c        in this case sizeof(__)*number of crystals
% c
% c
       matnum = iprops(38,local_el);
       ncrystals = imatprp(101,matnum);
% c
% c       So total history size is going to be:
% c               
% c
      size_data(1) = 76+12+max_slip_sys+ncrystals*(25+max_uhard);
%       return
end
