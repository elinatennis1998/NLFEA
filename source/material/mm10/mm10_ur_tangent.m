% c
% c     ****************************************************************
% c     *                                                              *
% c     *                 subroutine mm10_ur_tangent                   *
% c     *                                                              *
% c     *                       written by : mcm                       *
% c     *                                                              *
% c     *                   last modified: 1/16/14                     *
% c     *                                                              *
% c     *     Calculate the consistent tangent after a converged       *
% c     *     stress update.                                           *
% c     *                                                              *
% c     *     This routine uses the old, obsolete unrolling method     *
% c     *     for computing the tangent.  It's just for debug/         *
% c     *     comparison purposes.                                     *
% c     *                                                              *
% c     ****************************************************************
% c
function [props, np1, n] = mm10_ur_tangent(props, np1, n)
%       use mm10_defs
%       implicit none
% c
%       type(crystal_props) :: props
%       type(crystal_state) :: np1, n
% c
%       double precision, dimension(6,6) :: A
%       double precision, dimension(6) :: d_mod, d_barp, tw
%       double precision, dimension(6) :: b
%       double precision, dimension(3) :: w_p
%       double precision :: alpha
%       double precision, dimension(42,42) :: sys
%       double precision, dimension(42) :: rhs
%       double precision, dimension(7,6) :: AA
%       double precision, dimension(7,7) :: Jac
%       integer :: i, j, k, r_ind, c_ind, nes, nrhs, lda, ldb, info
%       integer, dimension(42) :: ipiv
% c
[props, np1, n, np1.stress, np1.tau_tilde, Jac] = ...
    mm10_formJ(props, np1, n, np1.stress, np1.tau_tilde, Jac);
[props, np1, n, np1.stress, np1.tau_tilde, d_barp] = ...
    mm10_form_dbarp(props, np1, n, np1.stress, ...
    np1.tau_tilde, d_barp);
[props, np1, n, np1.stress, np1.tau_tilde, w_p] = ...
    mm10_form_wp(props, np1, n, np1.stress, np1.tau_tilde, w_p);
[np1.stress, ~, tw] = mm10_symSW(np1.stress, w_p, tw);
% c
if (props.h_type == 1) % voche
    [props, np1, n, np1.stress, np1.tau_tilde, b] = ...
        mm10_ed_voche(props, np1, n, np1.stress, np1.tau_tilde, b);
elseif (props.h_type == 2) % MTS
    [props, np1, n, np1.stress, np1.tau_tilde, b] = ...
        mm10_ed_mts(props, np1, n, np1.stress, np1.tau_tilde, b);
elseif (props.h_type == 3) % User
    [props, np1, n, np1.stress, np1.tau_tilde, b] = ...
        mm10_ed_user(props, np1, n, np1.stress, np1.tau_tilde, b);
elseif (props.h_type == 6) % MTS_omar copy
    [props, np1, n, np1.stress, np1.tau_tilde, b] = ...
        mm10_ed_mts_omar(props, np1, n, np1.stress, np1.tau_tilde, b);
else
    [props] = mm10_unknown_hard_error(props);
end
% c
d_mod = np1.D;
d_mod(4:6) = 0.5 * d_mod(4:6);
% c
b = -b;
%       A = 0.0;
A = -props.stiffness;
alpha = 2.0/(3.0*np1.dg^2.0);
[6, 6, ~, props.stiffness*d_barp + 2.0*tw, 1,...
    ~, 1, A, 6] = DGER(6, 6, alpha, props.stiffness...
    *d_barp + 2.0*tw, 1, d_mod, 1, A, 6); % DGER - perform the rank 1 operation   A := alpha*x*y' + A
% c                                             % http://tinyurl.com/o4jxf6u
%       np1.tangent = 0.0;
% c
% c     Unroll
% c
AA(1:6,1:6) = -A;
AA(7,1:6) = -b;
rhs = 0.0;
sys = 0.0;
% c
for i = 1:7
    for j = 1:6
        r_ind = (i-1)*6+j;
        rhs(r_ind) = AA(i,j);
        for k = 1:7
            c_ind = (k-1)*6+j;
            sys(r_ind, c_ind) = Jac(i,k);
        end
    end
end
%
% c     Solve the equation
nes = 42;
nrhs = 1;
lda = 42;
ldb = 42;
[~, ~, ~, ~, ~, rhs, ~, info] = ...
    DGESV(nes, nrhs, sys, lda, ipiv, rhs, ldb, info); % DGESV description here
% c                                                         % http://tinyurl.com/oqk27tm
%       np1.tangent = 0.0;
%
for i = 1:6
    np1.tangent(i,1:6) = rhs(((i-1)*6+1):i*6);
end
%
%
end