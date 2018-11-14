% c
% c     ****************************************************************
% c     *                                                              *
% c     *                 subroutine mm10_formR1                       *
% c     *                                                              *
% c     *                       written by : tjt                       *
% c     *                                                              *
% c     *                   last modified: 02/27/15                    *
% c     *                                                              *
% c     *     Form residual and Jacobian for trust-region solver       *
% c     *                                                              *
% c     ****************************************************************
% c
function [R,J] = mm10_formRfsolve1(x,props,np1,n,tt)
%       use mm10_defs
%       implicit none
% c
%       type(crystal_props) :: props
%       type(crystal_state) :: np1, n
%       double precision, dimension(6) :: stress
stress = x(1:6);
R = zeros(6,1);
%       double precision :: tt
% c
 [~, ~, ~, ~, ~, R(1:6)] = mm10_formR1(props,np1,n,stress,... 
 tt);
% c
%       return
if nargout > 1 % return Jacobian too
[~,~,~,~,~, J] = mm10_formJ11(props, np1, n,...
    stress,tt);
end

end