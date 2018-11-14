% c
% c     ****************************************************************
% c     *                                                              *
% c     *                 subroutine mm10_formR                        *
% c     *                                                              *
% c     *                       written by : mcm                       *
% c     *                                                              *
% c     *                   last modified: 11/26/13                    *
% c     *                                                              *
% c     *     Form residual of dislocation densities                *
% c     *                                                              *
% c     ****************************************************************
% c
function [props, np1, n, stress_n1, tt, R] = mm10_formR0(props,np1,n,... 
stress,tt)
%       use mm10_defs
%       implicit none
% c
%       type(crystal_props) :: props
%       type(crystal_state) :: np1, n
%       double precision, dimension(6) :: stress
R = zeros(props.num_hard,1);
%       double precision :: tt
% c

% Compute current stress prediction using (A.6) and (20 - 23)
 [props, np1, n, stress, tt, R1] = mm10_formR1(props,np1,n,stress,... 
 tt);

 stress_n1 = stress - R1';

% Compute dislocation density residual using (A.11), (A.7), and (A.9)
 [props,np1,n,stress_n1,tt,R(1:props.num_hard)] = mm10_formR2(props,np1,n,stress_n1,tt);
% c
%       return
end