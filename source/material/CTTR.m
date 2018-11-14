function DMATX = CTTR(EPFLAG,IPROPS,LALGVA,NTYPE,EPROPS,RPROPS,RSTAVA,STRAT,STRES)
% 
%       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
%       PARAMETER(IPHARD=4  ,MDIM=3,  MSTRE=4)
%       LOGICAL EPFLAG, LALGVA(4), OUTOFP, RIGHT, REPEAT, TWOVEC
%       DIMENSION
%      1    DMATX(MSTRE,MSTRE) ,IPROPS(*)          ,RPROPS(*)          ,
%      2    RSTAVA(MSTRE+1)    ,STRAT(*)           ,STRES(*)
%       DIMENSION
%      1    DPSTRS(MDIM,MDIM)  ,DPSTRE(MDIM,MDIM)  ,EIGPRJ(MSTRE,2)    ,
%      2    FOID(MSTRE,MSTRE)  ,PSTRS(MDIM)        ,PSTRA(MDIM)        ,
%      3    SOID(MSTRE)        ,STRAC(MSTRE)
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
%      1    R0   ,RP5  ,R1   ,R2   ,R3   ,R4   / 
%      2    0.0D0,0.5D0,1.0D0,2.0D0,3.0D0,4.0D0/
% C***********************************************************************
% C COMPUTATION OF CONSISTENT TANGENT MODULUS FOR TRESCA TYPE
% C ELASTO-PLASTIC MATERIAL WITH PIECE-WISE LINEAR ISOTROPIC HARDENING.
% C PLANE STRAIN AND AXISYMMETRIC IMPLEMENTATIONS.
% C
% C REFERENCE: Section 8.1.5
% C***********************************************************************
MSTRE=4;
MDIM=3;
DMATX = zeros(MSTRE,MSTRE);
STRAC = zeros(MSTRE,1);
PSTRS = zeros(MDIM,1);
DPSTRS = zeros(MDIM,MDIM);
DPSTRE = zeros(MDIM,MDIM);
FOID = [1.0D0    ,0.0D0    ,0.0D0    ,0.0D0
        0.0D0    ,1.0D0    ,0.0D0    ,0.0D0
        0.0D0    ,0.0D0    ,0.5D0    ,0.0D0
        0.0D0    ,0.0D0    ,0.0D0    ,1.0D0];
SOID = [1.0D0    ,1.0D0    ,0.0D0    ,1.0D0]';
R0 = 0.0D0;
RP5 = 0.5D0;
R1 = 1.0D0;
R2 = 2.0D0;
R3 = 3.0D0;
R4 = 4.0D0;
% C Stops program if neither plane strain nor axisymmetric state
      if (NTYPE~=2&&NTYPE~=3), return, end
% C Current accumulated plastic strain
      EPBAR=RSTAVA(MSTRE+1);
% C Set material properties
      YOUNG=EPROPS(1);
      POISS=EPROPS(2);
      NHARD=IPROPS(1);
% C Set needed algorithmic variables
      TWOVEC=LALGVA(3);
      RIGHT=LALGVA(4);
% C Set some constants
      GMODU=YOUNG/(R2*(R1+POISS));
      BULK=YOUNG/(R3*(R1-R2*POISS));
      R2G=R2*GMODU;
      R4G=R4*GMODU;
      R1D3=R1/R3;
      R2D3=R2*R1D3;
      if(EPFLAG)
% C Compute elastoplastic consistent tangent
% C ----------------------------------------
% C Spectral decomposition of the elastic trial strain
        STRAC(1)=STRAT(1);
        STRAC(2)=STRAT(2);
        STRAC(3)=STRAT(3)*RP5;
        [EIGPRJ,PSTRA,REPEAT] = SPDEC2(STRAC);
        PSTRA(3)=STRAT(4);
% C and current total stress
        PSTRS(1)=STRES(1)*EIGPRJ(1,1)+STRES(2)*EIGPRJ(2,1)+...
                                   R2*STRES(3)*EIGPRJ(3,1);
        PSTRS(2)=STRES(1)*EIGPRJ(1,2)+STRES(2)*EIGPRJ(2,2)+...
                                   R2*STRES(3)*EIGPRJ(3,2);
        PSTRS(3)=STRES(4);
% C Identify directions of maximum and minimum principal trial stresses
        II=1;
        JJ=1;
        PSTMAX=PSTRA(II);
        PSTMIN=PSTRA(JJ);
        for I=2:3
          if(PSTRA(I)>=PSTMAX)
            II=I;
            PSTMAX=PSTRA(II);
          end %IF
          if(PSTRA(I)<PSTMIN)
            JJ=I;
            PSTMIN=PSTRA(JJ);
          end %IF
        end
        if(II~=1&&JJ~=1), MM=1; end
        if(II~=2&&JJ~=2), MM=2; end
        if(II~=3&&JJ~=3), MM=3; end
        if(TWOVEC)
% C Tangent consistent with two-vector return algorithm
          HSLOPE=DPLFUN(EPBAR,NHARD,RPROPS);
          DAA=R4G+HSLOPE;
          DAB=R2G+HSLOPE;
          DBA=R2G+HSLOPE;
          DBB=R4G+HSLOPE;
          DET=DAA*DBB-DAB*DBA;
          R2GDD=R2G/DET;
          R4G2DD=R2G*R2GDD;
          if(RIGHT)
% C ...returned to right corner
            DPSTRS(II,II)=R2G*(R1-R2GDD*R4G);
            DPSTRS(II,MM)=R4G2DD*(DAA-DAB);
            DPSTRS(II,JJ)=R4G2DD*(DBB-DBA);
            DPSTRS(MM,II)=R4G2DD*R2G;
            DPSTRS(MM,MM)=R2G*(R1-R2GDD*DAA);
            DPSTRS(MM,JJ)=R4G2DD*DBA;
            DPSTRS(JJ,II)=R4G2DD*R2G;
            DPSTRS(JJ,MM)=R4G2DD*DAB;
            DPSTRS(JJ,JJ)=R2G*(R1-R2GDD*DBB);
          else
% C ...returned to left corner
            DPSTRS(II,II)=R2G*(R1-R2GDD*DBB);
            DPSTRS(II,MM)=R4G2DD*DAB;
            DPSTRS(II,JJ)=R4G2DD*R2G;
            DPSTRS(MM,II)=R4G2DD*DBA;
            DPSTRS(MM,MM)=R2G*(R1-R2GDD*DAA);
            DPSTRS(MM,JJ)=R4G2DD*R2G;
            DPSTRS(JJ,II)=R4G2DD*(DBB-DBA);
            DPSTRS(JJ,MM)=R4G2DD*(DAA-DAB);
            DPSTRS(JJ,JJ)=R2G*(R1-R2GDD*R4G);
          end %IF
        else
% C Tangent consistent with one-vector return algorithm
          FACTOR=R2G/(R4G+DPLFUN(EPBAR,NHARD,RPROPS));
          DPSTRS(II,II)=R2G*(R1-FACTOR);
          DPSTRS(II,MM)=R0;
          DPSTRS(II,JJ)=R2G*FACTOR;
          DPSTRS(MM,II)=DPSTRS(II,MM);
          DPSTRS(MM,MM)=R2G;
          DPSTRS(MM,JJ)=R0;
          DPSTRS(JJ,II)=DPSTRS(II,JJ);
          DPSTRS(JJ,MM)=DPSTRS(MM,JJ);
          DPSTRS(JJ,JJ)=DPSTRS(II,II);
        end %IF
        DPSTRE(1,1)=+DPSTRS(1,1)*R2D3-DPSTRS(1,2)*R1D3-DPSTRS(1,3)*R1D3+...
                                                                   BULK;
        DPSTRE(2,1)=+DPSTRS(2,1)*R2D3-DPSTRS(2,2)*R1D3-DPSTRS(2,3)*R1D3+...
                                                                   BULK;
        DPSTRE(3,1)=+DPSTRS(3,1)*R2D3-DPSTRS(3,2)*R1D3-DPSTRS(3,3)*R1D3+...
                                                                   BULK;
        DPSTRE(1,2)=-DPSTRS(1,1)*R1D3+DPSTRS(1,2)*R2D3-DPSTRS(1,3)*R1D3+...
                                                                   BULK;
        DPSTRE(2,2)=-DPSTRS(2,1)*R1D3+DPSTRS(2,2)*R2D3-DPSTRS(2,3)*R1D3+...
                                                                   BULK;
        DPSTRE(3,2)=-DPSTRS(3,1)*R1D3+DPSTRS(3,2)*R2D3-DPSTRS(3,3)*R1D3+...
                                                                   BULK;
        DPSTRE(1,3)=-DPSTRS(1,1)*R1D3-DPSTRS(1,2)*R1D3+DPSTRS(1,3)*R2D3+...
                                                                   BULK;
        DPSTRE(2,3)=-DPSTRS(2,1)*R1D3-DPSTRS(2,2)*R1D3+DPSTRS(2,3)*R2D3+...
                                                                   BULK;
        DPSTRE(3,3)=-DPSTRS(3,1)*R1D3-DPSTRS(3,2)*R1D3+DPSTRS(3,3)*R2D3+...
                                                                   BULK;
        if(NTYPE==2)
          OUTOFP=0;
        elseif(NTYPE==3)
          OUTOFP=1;
        end %IF
        DMATX = DGISO2(DPSTRE,EIGPRJ,PSTRA,PSTRS,OUTOFP,REPEAT);
      else
% C Compute elasticity matrix
% C -------------------------
        if(NTYPE==2)
          NSTRE=3;
        elseif(NTYPE==3)
          NSTRE=4;
        end %IF
% C
        FACTOR=BULK-R2G*R1D3;
        for I=1:NSTRE
          for J=I:NSTRE
            DMATX(I,J)=R2G*FOID(I,J)+FACTOR*SOID(I)*SOID(J);
          end
        end
        for J=1:NSTRE-1
          for I=J+1:NSTRE
            DMATX(I,J)=DMATX(J,I);
          end
        end
      end %IF

