function [EIGPRJ,EIGX,REPEAT] = SPDEC2(X)
%       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
%       PARAMETER
%      1(   MCOMP=4    ,NDIM=2     )
%       LOGICAL REPEAT
%       DIMENSION
%      1    EIGPRJ(MCOMP,NDIM)        ,EIGX(NDIM)                ,
%      2    X(MCOMP)
%       DIMENSION
%      1    AUXMTX(NDIM,NDIM)         ,EIGVEC(NDIM,NDIM)
%       DATA
%      1    R0   ,RP5  ,R1   ,R4   ,SMALL  /
%      2    0.0D0,0.5D0,1.0D0,4.0D0,1.D-5  /
% C***********************************************************************
% C PERFORMS THE CLOSED FORM SPECTRAL DECOMPOSITION OF A
% C SYMMETRIC 2-D TENSOR STORED IN VECTOR FORM
% C
% C REFERENCE: Box A.2
% C***********************************************************************
MCOMP=4;
NDIM=2;
R0 = 0.0D0;
RP5 = 0.5D0;
R1 = 1.0D0;
R4 = 4.0D0;
SMALL=1.D-5;
EIGPRJ = zeros(MCOMP,NDIM);
EIGX = zeros(NDIM,1);
AUXMTX = zeros(NDIM,NDIM);
      REPEAT=0;
% C Compute eigenvalues of X
% C ------------------------
      TRX=X(1)+X(2);
      B=sqrt((X(1)-X(2))^2+R4*X(3)*X(3));
      EIGX(1)=RP5*(TRX+B);
      EIGX(2)=RP5*(TRX-B);
% C Compute eigenprojection tensors
% C -------------------------------
      DIFFER=abs(EIGX(1)-EIGX(2));
      AMXEIG=max(abs(EIGX(1)),abs(EIGX(2)));
      if(AMXEIG~=R0), DIFFER=DIFFER/AMXEIG; end
      if(DIFFER<SMALL)
        REPEAT=1;
% C for repeated (or nearly repeated) eigenvalues, re-compute eigenvalues
% C and compute eigenvectors using the iterative procedure. In such cases,
% C the closed formula for the eigenvectors is singular (or dominated by
% C round-off errors)
        AUXMTX(1,1)=X(1);
        AUXMTX(2,2)=X(2);
        AUXMTX(1,2)=X(3);
        AUXMTX(2,1)=AUXMTX(1,2);
        [EIGX,EIGVEC] = JACOB(AUXMTX,2);
        for IDIR=1:2
          EIGPRJ(1,IDIR)=EIGVEC(1,IDIR)*EIGVEC(1,IDIR);
          EIGPRJ(2,IDIR)=EIGVEC(2,IDIR)*EIGVEC(2,IDIR);
          EIGPRJ(3,IDIR)=EIGVEC(1,IDIR)*EIGVEC(2,IDIR);
          EIGPRJ(4,IDIR)=R0;
        end
      else
% C Use closed formula to compute eigenprojection tensors
        for IDIR=1:2
          B=EIGX(IDIR)-TRX;
          C=R1/(EIGX(IDIR)+B);
          EIGPRJ(1,IDIR)=C*(X(1)+B);
          EIGPRJ(2,IDIR)=C*(X(2)+B);
          EIGPRJ(3,IDIR)=C*X(3);
          EIGPRJ(4,IDIR)=R0;
        end
      end %IF

      