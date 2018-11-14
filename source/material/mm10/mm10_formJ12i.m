% c
% c     ****************************************************************
% c     *                                                              *
% c     *                 subroutine mm10_formJ12                      *
% c     *                                                              *
% c     *                       written by : mcm                       *
% c     *                                                              *
% c     *                   last modified: 11/27/13                    *
% c     *                                                              *
% c     *     Form the stress varying with hardening part              *
% c     *                                                              *
% c     ****************************************************************
% c
function [props, np1, n,vec1,vec2,arr1,arr2, stress, tt, J12] = mm10_formJ12i(props, np1,...
     n,vec1,vec2,arr1,arr2, stress, tt)
%       use mm10_defs
%       implicit none
% c
%       type(crystal_props) :: props
%       type(crystal_state) :: np1, n
%       double precision, dimension(6) :: stress
%       double precision, dimension(6) :: J12
%       double precision :: tt
% c
%       double precision, dimension(6) :: dbar, symTW
%       double precision, dimension(3) :: wp
% c
%       J12 = 0.0;
% c
      
      % Complex numerical derivative, also verified.
      J12 = zeros(6,props.num_hard);
        h = 1e-12;
%         for j = 1:props.num_hard
            for k = 1:props.num_hard
                A = tt;
                hi = h*tt(k);
                A(k) = A(k) + 1i*hi;
 [props, np1, n, stress, A, vec1,vec2] = mm10_formvecs(props,np1,n,stress,... 
 A,vec1,vec2);
                [~,~,~,~,~,~,~,Ri] = mm10_formR1(props,np1,n, vec1,vec2,stress,A);
                J12(1:6,k) = 1/hi*imag(Ri);
            end
%         end
% c
%       return
end