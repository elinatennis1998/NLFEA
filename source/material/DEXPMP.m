function [DEXPX      ,NOCONV] = DEXPMP(X)
%       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
%       PARAMETER
NDIM=3;
NDIM2=9;
% NDIM4=81;
MAXN=100;
MAXN0=MAXN;
% C Arguments
%       LOGICAL  NOCONV
%       DIMENSION
%      1    DEXPX(NDIM,NDIM,NDIM,NDIM),X(NDIM,NDIM)
DEXPX = zeros(NDIM,NDIM,NDIM,NDIM);
% C Local arrays and variables
% C...matrix of powers of X
%       DIMENSION
     R1DFAC = zeros(MAXN,1);
     XMATX = zeros(NDIM,NDIM,MAXN+1);
% C...initialise identity matrix: X to the power 0
     XMATX(:,:,1) = eye(3);

     R0=0.0D0;
     RP5=0.5D0;
     R1=1.0D0;
     TOL=1.0D-10;
     OVER=1.0D+100;
     UNDER=1.0D-100;
% C***********************************************************************
% C COMPUTES THE DERIVATIVE OF THE EXPONENTIAL OF A (GENERALLY
% C UNSYMMETRIC) 3-D TENSOR. USES THE SERIES REPRESENTATION OF THE TENSOR
% C EXPONENTIAL.
% C
% C REFERENCE: Section B.2
% C            Box B.2
% C***********************************************************************
% C Initialise convergence flag
      NOCONV=0;%.FALSE.
% C X to the power 1
      for I=1:NDIM
        for J=1:NDIM
          XMATX(I,J,1+1)=X(I,J);
        end
      end
% C Zero remaining powers of X
%       CALL RVZERO(XMATX(1,1,2),NDIM*NDIM*(MAXN-1))
% C Compute X square
      for I=1:NDIM
        for J=1:NDIM
          for K=1:NDIM
            XMATX(I,J,2+1)=XMATX(I,J,2+1)+X(I,K)*X(K,J);
          end
        end
      end
% C Compute principal invariants of X
      C1=X(1,1)+X(2,2)+X(3,3);
      C2=RP5*(C1*C1-(XMATX(1,1,2+1)+XMATX(2,2,2+1)+XMATX(3,3,2+1)));
      C3=X(1,1)*X(2,2)*X(3,3)+X(1,2)*X(2,3)*X(3,1)+ ...
         X(1,3)*X(2,1)*X(3,2)-X(1,2)*X(2,1)*X(3,3)- ...
         X(1,1)*X(2,3)*X(3,2)-X(1,3)*X(2,2)*X(3,1);
% C Compute X to the powers 3,4,...,NMAX using recursive formula
      R1DFAC(1)=R1;
      R1DFAC(2)=RP5;
      for N=3:MAXN
        R1DFAC(N)=R1DFAC(N-1)/N;
        for I=1:NDIM
          for J=1:NDIM
            XMATX(I,J,N+1)=C1*XMATX(I,J,N-1+1)-C2*XMATX(I,J,N-2+1)+ ...
                         C3*XMATX(I,J,N-3+1);
          end
        end
%         XNNORM=SQRT(SCAPRD(XMATX(:,:,N+1),XMATX(:,:,N+1),NDIM2));
        XNNORM=norm(XMATX(:,:,N+1));
% C...check number of terms required for series convergence
        if(XNNORM>OVER||(XNNORM<UNDER&&XNNORM>R0)...
                                        ||R1DFAC(N)<UNDER) %THEN
% C...numbers are to small or too big: Exit without computing derivative
          NOCONV=1;%.TRUE.
          return
        elseif(XNNORM*R1DFAC(N)<TOL) %THEN
% C...series will converge with NMAX terms:
% C   Carry on to derivative computation
          NMAX=N;
          break
%           GOTO 90
        end %IF
      end
% C...series will not converge for the currently prescribed tolerance
% C   with the currently prescribed maximum number of terms MAXN:
% C   Exit without computing derivative
if N == MAXN0
      NOCONV=1;%.TRUE.
      return
end
%    90 CONTINUE
% C Compute the derivative of exponential map
%       CALL RVZERO(DEXPX,NDIM4)
      for I=1:NDIM
        for J=1:NDIM
          for K=1:NDIM
            for L=1:NDIM
              for N=1:NMAX
                for M=1:N
                  DEXPX(I,J,K,L)=DEXPX(I,J,K,L)+ ...
                                 R1DFAC(N)*XMATX(I,K,M-1+1)*XMATX(L,J,N-M+1);
                end
              end
            end
          end
        end
      end
% C

      
%       function out = SCAPRD(U  ,V  ,N  ) 
%       
% R0 = 0.0D0;
% % C***********************************************************************
% % C SCALAR PRODUCT OF DOUBLE PRECISION VECTORS U AND V OF DIMENSION N
% % C***********************************************************************
%       out=R0;
%       for I=1:N
%         out=out+U(I)*V(I);
%       end