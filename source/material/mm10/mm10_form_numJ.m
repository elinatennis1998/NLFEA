% c
% c     ****************************************************************
% c     *                                                              *
% c     *                 subroutine mm10_form_numJ                    *
% c     *                                                              *
% c     *                       written by : mcm                       *
% c     *                                                              *
% c     *                   last modified: 11/26/13                    *
% c     *                                                              *
% c     *     Form the jacobian from numerical differentiation         *
% c     *                                                              *
% c     ****************************************************************
% c
function [props, np1, n, stress, tt, J] = mm10_form_numJ(props,...
     np1, n, stress, tt, J)
%       use mm10_defs
%       implicit none
% c
%       type(crystal_props) :: props
%       type(crystal_state) :: np1, n
%       double precision, dimension(6) :: stress
%       double precision, dimension(7,7) :: J
%       double precision :: tt
% c
%       double precision :: eps, pt
%       double precision, dimension(7) :: R, pR
%       double precision, dimension(6) :: pS
%       integer :: i
%       parameter(eps = 1e-6)
% c
%       J = 0.0
      [props, np1, n, stress, tt, R] = mm10_formR(props, np1,...
          n, stress, tt, R);
      for i = 1:7
        if (i ~= 7)
          pS = stress;
          pS(i) = pS(i) + eps;
          [props, np1, n, ~, tt, pR] = mm10_formR(props, np1, n,...
              pS, tt, pR);
        else
          pt = tt + eps;
         [props, np1, n, stress, ~, pR] = mm10_formR(props, np1,...
             n, stress, pt, pR);
        end
        J(:,i) = (pR-R)/eps;
      end
% c
%       return
end
