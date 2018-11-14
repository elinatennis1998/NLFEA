% c
% c     ****************************************************************
% c     *                                                              *
% c     *                 subroutine mm10_formJ                        *
% c     *                                                              *
% c     *                       written by : mcm                       *
% c     *                                                              *
% c     *                   last modified: 11/26/13                    *
% c     *                                                              *
% c     *     Form the jacobian from lower subroutines                 *
% c     *                                                              *
% c     ****************************************************************
% c
function [props, np1, n,vec1,vec2,arr1,arr2, stress, tt, J] = mm10_formJr(props, np1, n,vec1,vec2,arr1,arr2,...
    stress,tt)
%       use mm10_defs
%       implicit none
% c
%       type(crystal_props) :: props
%       type(crystal_state) :: np1, n
%       double precision, dimension(6) :: stress
J = zeros(props.num_hard+6,props.num_hard+6);
%       double precision :: tt
% c
[props, np1, n,vec1,vec2,arr1,arr2, stress, tt, J(1:6,1:6)] = mm10_formJ11r(props, np1, n,vec1,vec2,arr1,arr2,...
    stress, tt);
[props, np1, n,vec1,vec2,arr1,arr2, stress, tt, J(1:6,7:props.num_hard+6)] = mm10_formJ12r(props, np1, n,vec1,vec2,arr1,arr2,...
    stress, tt);
[props, np1, n,vec1,vec2,arr1,arr2, stress, tt, J(7:props.num_hard+6, 1:6)] = mm10_formJ21r(props, np1, n,vec1,vec2,arr1,arr2,...
    stress, tt);
[props, np1, n,vec1,vec2,arr1,arr2, stress, tt, J(7:props.num_hard+6,7:props.num_hard+6)] = mm10_formJ22r(props, np1, n,vec1,vec2,arr1,arr2,...
    stress, tt);
% c
%       return
end