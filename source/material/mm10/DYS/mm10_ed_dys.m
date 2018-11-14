% c
% c
% c     mrr:
% c
% c           Derivative of hardening fn wrt strain
function [props, np1, n, stress, tt, ed] = mm10_ed_dys(props, np1,...
    n, stress, tt, ed)
%       use mm10_defs
%       implicit none
% c
%       type(crystal_props) :: props
%       type(crystal_state) :: np1, n
%       double precision, dimension(6) :: stress
%       double precision :: tt
%       double precision, dimension(6) :: ed
% c
%       double precision :: mm10_slipinc
% c
%       double precision :: dslip, lnv, lny, dgc, ty, tv, mnp0, m0np, sc
%       double precision, dimension(6) :: d_mod, dydd, dvdd
%       integer :: s
% c
% c     Form a bunch of simple constitutive things
% c
% c
      ed = zeros(6,1);
% 
%       return
end