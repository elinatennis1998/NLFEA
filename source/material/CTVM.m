function DMATX = CTVM(DGAMA,EPFLAG,IPROPS,NTYPE,EPROPS,RPROPS,RSTAVA,STRES)
% 
%       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
%       PARAMETER(IPHARD=4  ,MSTRE=4)
%       LOGICAL EPFLAG
%       DIMENSION
%      1    DMATX(MSTRE,MSTRE),IPROPS(*)           ,RPROPS(*)          ,
%      2    RSTAVA(MSTRE+1)   ,STRES(MSTRE)
%       DIMENSION
%      1    DEVPRJ(MSTRE,MSTRE),FOID(MSTRE,MSTRE)  ,S(MSTRE)           ,
%      2    SOID(MSTRE)
%       DATA
%      1    FOID(1,1),FOID(1,2),FOID(1,3),FOID(1,4)/
%      2    1.0D0    ,0.0D0    ,0.0D0    ,0.0D0    /
%      3    FOID(2,1),FOID(2,2),FOID(2,3),FOID(2,4)/
%      4    0.0D0    ,1.0D0    ,0.0D0    ,0.0D0    /
%      5    FOID(3,1),FOID(3,2),FOID(3,3),FOID(3,4)/
%      6    0.0D0    ,0.0D0    ,0.5D0    ,0.0D0    /
%      7    FOID(4,1),FOID(4,2),FOID(4,3),FOID(4,4)/
%      8    0.0D0    ,0.0D0    ,0.0D0    ,1.0D0    /
%       DATA
%      1    SOID(1)  ,SOID(2)  ,SOID(3)  ,SOID(4)  /
%      2    1.0D0    ,1.0D0    ,0.0D0    ,1.0D0    /
%       DATA
%      1    R1   ,R2   ,R3   ,R6   /
%      2    1.0D0,2.0D0,3.0D0,6.0D0/
% C***********************************************************************
% C COMPUTATION OF THE CONSISTENT TANGENT MODULUS FOR VON MISES TYPE
% C ELASTO-PLASTIC MATERIAL WITH PIECE-WISE LINEAR ISOTROPIC HARDENING.
% C PLANE STRAIN AND AXISYMMETRIC IMPLEMENTATIONS.
% C
% C REFERENCE: Section 7.4.3
% C***********************************************************************
MSTRE=4;
DMATX = zeros(MSTRE,MSTRE);
DEVPRJ = zeros(MSTRE,MSTRE);
S = zeros(MSTRE,1);
FOID = [1.0D0    ,0.0D0    ,0.0D0    ,0.0D0
        0.0D0    ,1.0D0    ,0.0D0    ,0.0D0
        0.0D0    ,0.0D0    ,0.5D0    ,0.0D0
        0.0D0    ,0.0D0    ,0.0D0    ,1.0D0];
SOID = [1.0D0    ,1.0D0    ,0.0D0    ,1.0D0]';
R1 = 1.0D0;
R2 = 2.0D0;
R3 = 3.0D0;
R6 = 6.0D0;
% C Stops program if neither plane strain nor axisymmetric state
      if (NTYPE~=2&&NTYPE~=3), return, end
% C Current accumulated plastic strain
      EPBAR=RSTAVA(MSTRE+1);
% C Set material properties
      YOUNG=EPROPS(1);
      POISS=EPROPS(2);
      NHARD=IPROPS(1);
% C Shear and bulk moduli
      GMODU=YOUNG/(R2*(R1+POISS));
      BULK=YOUNG/(R3*(R1-R2*POISS));
      R2G=R2*GMODU;
      R1D3=R1/R3;
% C Set deviatoric projection tensor
      if(NTYPE==2)
        NSTRE=3;
      elseif(NTYPE==3)
        NSTRE=4;
      end %IF
        for J=1:NSTRE
      for I=1:NSTRE
          DEVPRJ(I,J)=FOID(I,J)-SOID(I)*SOID(J)*R1D3;
      end
        end
      if(EPFLAG)
% C Compute elastoplastic consistent tangent
% C ----------------------------------------
        R3G=R3*GMODU;
        ROO3D2=sqrt(R3/R2);
% C Hydrostatic pressure
        P=(STRES(1)+STRES(2)+STRES(4))*R1D3;
% C Deviatoric stress components
        S(1)=STRES(1)-P;
        S(2)=STRES(2)-P;
        S(3)=STRES(3);
        S(4)=STRES(4)-P;
% C Recover last elastic trial von Mises effective stress
        SNORM=sqrt(S(1)*S(1)+S(2)*S(2)+R2*S(3)*S(3)+S(4)*S(4));
        Q=ROO3D2*SNORM;
        QTRIAL=Q+R3G*DGAMA;
% C Assemble elastoplastic tangent (upper triangle only)
        AFACT=R2G*(R1-R3G*DGAMA/QTRIAL);
        BFACT=R6*GMODU*GMODU*(DGAMA/QTRIAL- ...
              R1/(R3G+DPLFUN(EPBAR,NHARD,RPROPS)))/ ...
              (SNORM*SNORM);
          for J=1:NSTRE
        for I=1:J
            DMATX(I,J)=AFACT*DEVPRJ(I,J)+BFACT*S(I)*S(J)+ ...
                       BULK*SOID(I)*SOID(J);
        end
          end
      else
% C Compute elasticity matrix (upper triangle only)
% C -----------------------------------------------
          for J=1:NSTRE
        for I=1:J
            DMATX(I,J)=R2G*DEVPRJ(I,J)+BULK*SOID(I)*SOID(J);
        end
          end
      end %IF
% C Assemble lower triangle
% C -----------------------
        for I=2:NSTRE
      for J=1:I
          DMATX(I,J)=DMATX(J,I);
      end
        end
%       RETURN
%       END
