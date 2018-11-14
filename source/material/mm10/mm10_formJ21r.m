% c
% c     ****************************************************************
% c     *                                                              *
% c     *                 subroutine mm10_formJ21                      *
% c     *                                                              *
% c     *                       written by : mcm                       *
% c     *                                                              *
% c     *                   last modified: 11/27/13                    *
% c     *                                                              *
% c     *     Form the hardening varying with stress part              *
% c     *                                                              *
% c     ****************************************************************
% c
function [props, np1, n,vec1,vec2,arr1,arr2, stress, tt, J21] = mm10_formJ21r(props, np1,...
    n,vec1,vec2,arr1,arr2, stress, tt)
% use mm10_defs
%       implicit none
% c
%       type(crystal_props) :: props
%       type(crystal_state) :: np1, n
%       double precision, dimension(6) :: stress
%       double precision, dimension(6) :: J21
%       double precision :: tt
% c
%       double precision, dimension(6) :: estr
% c
%        J21 = 0.0;
% c
      
      % Complex numerical calculation of derivative; VERIFIED for MTS and
      % mrr models
      J21 = zeros(props.num_hard,6);
        h = 1e-8;
%         for j = 1:props.num_hard
            for k = 1:6
                A = stress;
 [props, np1, n, A, tt, vec1,vec2] = mm10_formvecs(props,np1,n,A,... 
 tt,vec1,vec2);
                [~,~,~,~,~,~,~,R0] = mm10_formR2(props,np1,n, vec1,vec2,A,tt);
                if stress(k) == 0
                hi = h;
                else
                hi = h*abs(stress(k));
                end
                A(k) = A(k) + hi;
 [props, np1, n, A, tt, vec1,vec2] = mm10_formvecs(props,np1,n,A,... 
 tt,vec1,vec2);
                [~,~,~,~,~,~,~,Ri] = mm10_formR2(props,np1,n, vec1,vec2,A,tt);
                J21(1:props.num_hard,k) = 1/hi*(Ri-R0);
            end
%         end
% c
%       return
end