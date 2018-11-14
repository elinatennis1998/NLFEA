% c
% c     ****************************************************************
% c     *                                                              *
% c     *                 subroutine mm10_formR                        *
% c     *                                                              *
% c     *                       written by : mcm                       *
% c     *                                                              *
% c     *                   last modified: 11/26/13                    *
% c     *                                                              *
% c     *     Form the residual from lower subroutines                 *
% c     *                                                              *
% c     ****************************************************************
% c
function [props, np1, n, vec1, vec2, stress, tt, R] = mm10_formR(props,np1,n,vec1,vec2,... 
stress,tt)
%       use mm10_defs
%       implicit none
% c
%       type(crystal_props) :: props
%       type(crystal_state) :: np1, n
%       double precision, dimension(6) :: stress
R = zeros(props.num_hard+6,1);
%       double precision :: tt
% c
 [props, np1, n, stress, tt, vec1,vec2] = mm10_formvecs(props,np1,n,stress,... 
 tt,vec1,vec2);
 [props, np1, n,vec1,vec2, stress, tt, R(1:6)] = mm10_formR1(props,np1,n,vec1,vec2,stress,... 
 tt);
 [props,np1,n,vec1,vec2,stress,tt,R(7:props.num_hard+6)] = mm10_formR2(props,np1,n,vec1,vec2,stress,tt);
% c
%       return
end