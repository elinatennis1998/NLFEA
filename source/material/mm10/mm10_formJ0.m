% c
% c     ****************************************************************
% c     *                                                              *
% c     *                 subroutine mm10_formJ                        *
% c     *                                                              *
% c     *                       written by : tjt                       *
% c     *                                                              *
% c     *                   last modified: 02/17/15                    *
% c     *                                                              *
% c     *     Form the jacobian from lower subroutines                 *
% c     *                                                              *
% c     ****************************************************************
% c
function [props, np1, n, stress, tt, J] = mm10_formJ0(props, np1, n,...
    stress,tt)
%       use mm10_defs
%       implicit none
% c
%       type(crystal_props) :: props
%       type(crystal_state) :: np1, n
%       double precision, dimension(6) :: stress
% J = zeros(props.num_hard,props.num_hard);
%       double precision :: tt
% c
[props, np1, n, stress, tt, J11] = mm10_formJ11(props, np1, n,...
    stress, tt);
[props, np1, n, stress, tt, J12] = mm10_formJ12(props, np1, n,...
    stress, tt);
[props, np1, n, stress, tt, J21] = mm10_formJ21(props, np1, n,...
    stress, tt);
[props, np1, n, stress, tt, J22] = mm10_formJ22(props, np1, n,...
    stress, tt);

J = J22 - J21*(J11\J12);

% c
%       return
end