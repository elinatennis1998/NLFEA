function DYDX = DGISO2(DEIGY,EIGPRJ,EIGX,EIGY,OUTOFP,REPEAT)
% 
%       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
%       PARAMETER
%      1(   MCOMP=4    ,MDIM=3     ,NDIM=2     )
% C Arguments
%       LOGICAL OUTOFP, REPEAT
%       DIMENSION
%      1    DEIGY(MDIM,MDIM)  ,DYDX(MCOMP,MCOMP) ,EIGPRJ(MCOMP,NDIM),
%      2    EIGX(NDIM)        ,EIGY(NDIM)
% C Local arrays
%       DIMENSION
%      1    EIGPR3(MCOMP)     ,FOID(MCOMP,MCOMP) ,SOPID(MCOMP)
%       DATA
%      1    FOID(1,1)     ,FOID(1,2)     ,FOID(1,3)     /
%      2    1.0D0         ,0.0D0         ,0.0D0         /
%      3    FOID(2,1)     ,FOID(2,2)     ,FOID(2,3)     /
%      4    0.0D0         ,1.0D0         ,0.0D0         /
%      5    FOID(3,1)     ,FOID(3,2)     ,FOID(3,3)     /
%      6    0.0D0         ,0.0D0         ,0.5D0         /
%       DATA
%      1    SOPID(1)      ,SOPID(2)      ,SOPID(3)      ,SOPID(4)        /
%      2    1.0D0         ,1.0D0         ,0.0D0         ,0.0D0           /
%       DATA
%      1    EIGPR3(1)     ,EIGPR3(2)     ,EIGPR3(3)     ,EIGPR3(4)       /
%      2    0.0D0         ,0.0D0         ,0.0D0         ,1.0D0           /
% C***********************************************************************
% C COMPUTE THE DERIVATIVE OF A GENERAL ISOTROPIC TENSOR FUNCTION OF ONE
% C TENSOR IN 2-D (WITH ONE POSSIBLE OUT-OF-PLANE COMPONENT)
% C
% C REFERENCE: Sections A.3.1-2
% C            Box A.2
% C***********************************************************************
MCOMP=4;
% MDIM=3;
% NDIM=2;
FOID = [1.0D0    ,0.0D0    ,0.0D0
        0.0D0    ,1.0D0    ,0.0D0
        0.0D0    ,0.0D0    ,0.5D0];
SOPID = [1.0D0    ,1.0D0    ,0.0D0    ,0.0D0]';
EIGPR3 = [0.0D0    ,0.0D0    ,0.0D0    ,1.0D0]';

      DYDX = zeros(MCOMP,MCOMP);
      if(REPEAT)
% C Derivative dY/dX for repeated in-plane eigenvalues of X
% C -------------------------------------------------------
% C In-plane component
        for I=1:3
          for J=1:3
            DYDX(I,J)=(DEIGY(1,1)-DEIGY(1,2))*FOID(I,J)+...
                       DEIGY(1,2)*SOPID(I)*SOPID(J);
          end
        end
        if(OUTOFP)
% C out-of-plane components required
          for I=1:4
            for J=1:4
              if(I==4||J==4)
                  DYDX(I,J)=...
                          DEIGY(1,3)*SOPID(I)*EIGPR3(J)+...
                          DEIGY(3,1)*EIGPR3(I)*SOPID(J)+...
                          DEIGY(3,3)*EIGPR3(I)*EIGPR3(J);
              end
            end
          end
        end %IF
      else
% C Derivative dY/dX for distinct in-plane eigenvalues of X
% C -------------------------------------------------------
% C Assemble in-plane DYDX
        A1=(EIGY(1)-EIGY(2))/(EIGX(1)-EIGX(2));
        for I=1:3
          for J=1:3
            DYDX(I,J)=A1*(FOID(I,J)-EIGPRJ(I,1)*EIGPRJ(J,1)-...
                      EIGPRJ(I,2)*EIGPRJ(J,2))+...
                      DEIGY(1,1)*EIGPRJ(I,1)*EIGPRJ(J,1)+...
                      DEIGY(1,2)*EIGPRJ(I,1)*EIGPRJ(J,2)+...
                      DEIGY(2,1)*EIGPRJ(I,2)*EIGPRJ(J,1)+...
                      DEIGY(2,2)*EIGPRJ(I,2)*EIGPRJ(J,2);
          end
        end
        if(OUTOFP)
% C out-of-plane components required
          for I=1:4
            for J=1:4
              if(I==4||J==4)
                  DYDX(I,J)=...
                      DEIGY(1,3)*EIGPRJ(I,1)*EIGPR3(J)+...
                      DEIGY(2,3)*EIGPRJ(I,2)*EIGPR3(J)+...
                      DEIGY(3,1)*EIGPR3(I)*EIGPRJ(J,1)+...
                      DEIGY(3,2)*EIGPR3(I)*EIGPRJ(J,2)+...
                      DEIGY(3,3)*EIGPR3(I)*EIGPR3(J);
              end
            end
          end
        end %IF
      end %IF

      