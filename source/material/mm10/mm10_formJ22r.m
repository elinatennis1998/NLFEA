% c
% c     ****************************************************************
% c     *                                                              *
% c     *                 subroutine mm10_formJ22                      *
% c     *                                                              *
% c     *                       written by : mcm                       *
% c     *                                                              *
% c     *                   last modified: 11/27/13                    *
% c     *                                                              *
% c     *     Form the hardening varying with hardening part           *
% c     *                                                              *
% c     ****************************************************************
% c
function [props, np1, n,vec1,vec2,arr1,arr2, stress, tt, J22] = mm10_formJ22r(props, np1,...
     n,vec1,vec2,arr1,arr2, stress, tt)
%       use mm10_defs
%       implicit none
% c
%       type(crystal_props) :: props
%       type(crystal_state) :: np1, n
%       double precision, dimension(6) :: stress
%       double precision :: tt, J22
% c
%       double precision :: etau
% c
%       J22 = 0.0;
% c
      
      % Complex numerical derivative, also verified.
      J22 = zeros(props.num_hard,props.num_hard);
        h = 1e-8;
%         for j = 1:props.num_hard
            for k = 1:props.num_hard
                A = tt;
 [props, np1, n, stress, A, vec1,vec2] = mm10_formvecs(props,np1,n,stress,... 
 A,vec1,vec2);
                [~,~,~,~,~,~,~,R0] = mm10_formR2(props,np1,n, vec1,vec2,stress,A);
                hi = h*tt(k);
                A(k) = A(k) + hi;
 [props, np1, n, stress, A, vec1,vec2] = mm10_formvecs(props,np1,n,stress,... 
 A,vec1,vec2);
                [~,~,~,~,~,~,~,Ri] = mm10_formR2(props,np1,n, vec1,vec2,stress,A);
                J22(1:props.num_hard,k) = 1/hi*(Ri-R0);
                % test to make sure not singular
                if J22(k,k) == 0
                    J22(k,k) = 1;
                end
            end
%         end
% c
end