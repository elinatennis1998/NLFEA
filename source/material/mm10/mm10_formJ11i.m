% c
% c     ****************************************************************
% c     *                                                              *
% c     *                 subroutine mm10_formJ11                      *
% c     *                                                              *
% c     *                       written by : mcm                       *
% c     *                                                              *
% c     *                   last modified: 11/26/13                    *
% c     *                                                              *
% c     *     Form the stress varying with stress part                 *
% c     *                                                              *
% c     ****************************************************************
% c
function [props, np1, n,vec1,vec2,arr1,arr2, stress, tt, J11] = mm10_formJ11i(props, np1,...
    n,vec1,vec2,arr1,arr2, stress, tt)
%       use mm10_defs
%       implicit none
% c
%       type(crystal_props) :: props
%       type(crystal_state) :: np1, n
%       double precision, dimension(6) :: stress
%       double precision, dimension(6,6) :: J11
%       double precision :: tt
% c
%       double precision :: mm10_rs
% c
%       double precision, dimension(6) :: symtq
%       double precision, dimension(3) :: wp
%       double precision, dimension(6,6) :: Iw
%       double precision :: rs, alpha
%       integer :: i
% c
% c
      
      % Complex numerical derivative, also verified.
      J11 = zeros(6,6);
        h = 1e-12;
%         for j = 1:props.num_hard
            for k = 1:6
                A = stress;
                if stress(k) == 0
                hi = h;
                else
                hi = h*abs(stress(k));
                end
                A(k) = A(k) + 1i*hi;
 [props, np1, n, A, tt, vec1,vec2] = mm10_formvecs(props,np1,n,A,... 
 tt,vec1,vec2);
                [~,~,~,~,~,~,~,Ri] = mm10_formR1(props,np1,n,vec1,vec2,A,tt);
                J11(1:6,k) = 1/hi*imag(Ri);
            end
%         end
% c
%       return
end