function [A_iIjJ,C_IJKL] = CSFtoA(C,S,F,ndm)
%
% Function to convert C_IJKL and S_IJ into A_iIjJ

A_iIjJ = zeros(ndm,ndm,ndm,ndm);
C_IJKL = A_iIjJ;
g = eye(ndm);

C_IJKL(1,1,1,1) = C(1,1);
C_IJKL(2,2,2,2) = C(2,2);
C_IJKL(1,2,1,2) = C(ndm+1,ndm+1);
C_IJKL(1,2,2,1) = C(ndm+1,ndm+1);
C_IJKL(2,1,1,2) = C(ndm+1,ndm+1);
C_IJKL(2,1,2,1) = C(ndm+1,ndm+1);
C_IJKL(2,1,2,2) = C(2,ndm+1);
C_IJKL(1,2,2,2) = C(2,ndm+1);
C_IJKL(2,2,1,2) = C(2,ndm+1);
C_IJKL(2,2,2,1) = C(2,ndm+1);
C_IJKL(2,1,1,1) = C(1,ndm+1);
C_IJKL(1,2,1,1) = C(1,ndm+1);
C_IJKL(1,1,1,2) = C(1,ndm+1);
C_IJKL(1,1,2,1) = C(1,ndm+1);
C_IJKL(1,1,2,2) = C(1,2);
C_IJKL(2,2,1,1) = C(1,2);
if ndm == 3
C_IJKL(3,3,3,3) = C(3,3);
C_IJKL(2,3,2,3) = C(5,5);
C_IJKL(2,3,3,2) = C(5,5);
C_IJKL(3,2,2,3) = C(5,5);
C_IJKL(3,2,3,2) = C(5,5);
C_IJKL(3,2,3,3) = C(3,5);
C_IJKL(2,3,3,3) = C(3,5);
C_IJKL(3,3,2,3) = C(3,5);
C_IJKL(3,3,3,2) = C(3,5);
C_IJKL(3,2,2,2) = C(2,5);
C_IJKL(2,3,2,2) = C(2,5);
C_IJKL(2,2,2,3) = C(2,5);
C_IJKL(2,2,3,2) = C(2,5);
C_IJKL(2,2,3,3) = C(2,3);
C_IJKL(3,3,2,2) = C(2,3);
C_IJKL(3,1,3,1) = C(6,6);
C_IJKL(3,1,1,3) = C(6,6);
C_IJKL(1,3,3,1) = C(6,6);
C_IJKL(1,3,1,3) = C(6,6);
C_IJKL(1,3,1,1) = C(1,6);
C_IJKL(3,1,1,1) = C(1,6);
C_IJKL(1,1,3,1) = C(1,6);
C_IJKL(1,1,1,3) = C(1,6);
C_IJKL(1,3,3,3) = C(3,6);
C_IJKL(3,1,3,3) = C(3,6);
C_IJKL(3,3,3,1) = C(3,6);
C_IJKL(3,3,1,3) = C(3,6);
C_IJKL(3,3,1,1) = C(3,1);
C_IJKL(1,1,3,3) = C(3,1);
C_IJKL(3,2,3,1) = C(5,6);
C_IJKL(3,2,1,3) = C(5,6);
C_IJKL(1,3,3,2) = C(5,6);
C_IJKL(1,3,2,3) = C(5,6);
C_IJKL(2,3,3,1) = C(5,6);
C_IJKL(2,3,1,3) = C(5,6);
C_IJKL(3,1,3,2) = C(5,6);
C_IJKL(3,1,2,3) = C(5,6);
C_IJKL(2,3,1,1) = C(1,5);
C_IJKL(3,2,1,1) = C(1,5);
C_IJKL(1,1,3,2) = C(1,5);
C_IJKL(1,1,2,3) = C(1,5);
C_IJKL(1,3,1,2) = C(6,4);
C_IJKL(1,3,2,1) = C(6,4);
C_IJKL(2,1,1,3) = C(6,4);
C_IJKL(2,1,3,1) = C(6,4);
C_IJKL(3,1,1,2) = C(6,4);
C_IJKL(3,1,2,1) = C(6,4);
C_IJKL(1,2,1,3) = C(6,4);
C_IJKL(1,2,3,1) = C(6,4);
C_IJKL(3,1,2,2) = C(2,6);
C_IJKL(1,3,2,2) = C(2,6);
C_IJKL(2,2,1,3) = C(2,6);
C_IJKL(2,2,3,1) = C(2,6);
C_IJKL(2,1,2,3) = C(4,5);
C_IJKL(2,1,3,2) = C(4,5);
C_IJKL(3,2,2,1) = C(4,5);
C_IJKL(3,2,1,2) = C(4,5);
C_IJKL(1,2,2,3) = C(4,5);
C_IJKL(1,2,3,2) = C(4,5);
C_IJKL(2,3,2,1) = C(4,5);
C_IJKL(2,3,1,2) = C(4,5);
C_IJKL(1,2,3,3) = C(3,4);
C_IJKL(2,1,3,3) = C(3,4);
C_IJKL(3,3,2,1) = C(3,4);
C_IJKL(3,3,1,2) = C(3,4);
end

for i = 1:ndm
    for I = 1:ndm
        for j = 1:ndm
            for J = 1:ndm
                F_iK_C_IJKL_F_jL = 0;
                for K = 1:ndm
                    F_iK = F(i,K);
                    for L = 1:ndm
                        F_iK_C_IJKL_F_jL = F_iK_C_IJKL_F_jL + F_iK*C_IJKL(K,I,L,J)*F(j,L);
                    end
                end
                A_iIjJ(i,I,j,J) = S(I,J)*g(i,j) + F_iK_C_IJKL_F_jL;
end,end,end,end