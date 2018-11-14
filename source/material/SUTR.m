function [DGAM,LALGVA,RSTAVA,STRES] = SUTR(IPROPS,NTYPE,EPROPS,RPROPS,RSTAVA,STRAT)
%
%       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
%       PARAMETER(IPHARD=4  ,MSTRE=4)
% C Arguments
%       LOGICAL
%      1    LALGVA(4)
%       DIMENSION
%      1    DGAM(2)            ,IPROPS(*)          ,RPROPS(*)          ,
%      2    RSTAVA(MSTRE+1)    ,STRAT(MSTRE)       ,STRES(MSTRE)
% C Local arrays and variables
%       LOGICAL
%      1    DUMMY, IFPLAS, RIGHT, SUFAIL, TWOVEC
%       DIMENSION
%      1    EIGPRJ(MSTRE,2)    ,PSTRS(3)           ,STREST(3)
% 
%       DATA
%      1    R0   ,R1   ,R2   ,R3   ,R4   ,SMALL ,TOL   / 
%      2    0.0D0,1.0D0,2.0D0,3.0D0,4.0D0,1.D-10,1.D-10/
%       DATA MXITER / 50 /
% C***********************************************************************
% C STRESS UPDATE PROCEDURE FOR TRESCA TYPE ELASTO-PLASTIC MATERIAL WITH
% C PIECE-WISE LINEAR ISOTROPIC HARDENING:
% C IMPLICIT ELASTIC PREDICTOR/RETURN MAPPING ALGORITHM (Boxes 8.1-3).
% C PLANE STRAIN AND AXISYMMETRIC IMPLEMENTATIONS.
% C
% C REFERENCE: Boxes 8.1-3
% C            Section 8.1.2
% C***********************************************************************
MSTRE=4;
STRES = zeros(MSTRE,1);
R0 = 0.0d0;
R1 = 1.0d0;
R2 = 2.0d0;
R3 = 3.0d0;
R4 = 4.0d0;
SMALL = 1.D-10;
TOL = 1.D-10;
MXITER = 50;
% C Stops program if neither plane strain nor axisymmetric state
      if (NTYPE~=2&&NTYPE~=3), return, end
% C Initialize some algorithmic and internal variables
      DGAMA=R0;
      DGAMB=R0;
      IFPLAS=0;
      TWOVEC=0;
      RIGHT=0;
      SUFAIL=0;
      EPBARN=RSTAVA(MSTRE+1);
      EPBAR=EPBARN;
% C Set some material properties
      YOUNG=EPROPS(1);
      POISS=EPROPS(2);
      NHARD=IPROPS(1);
% C Set some constants
      GMODU=YOUNG/(R2*(R1+POISS));
      BULK=YOUNG/(R3*(R1-R2*POISS));
      R2G=R2*GMODU;
      R4G=R4*GMODU;
      R1D3=R1/R3;
% C Compute elastic trial state
% C ---------------------------
% C Volumetric strain and pressure stress
      EEV=STRAT(1)+STRAT(2)+STRAT(4);
      P=BULK*EEV;
% C Spectral decomposition of the elastic trial deviatoric stress
      EEVD3=EEV*R1D3;
      STREST(1)=R2G*(STRAT(1)-EEVD3);
      STREST(2)=R2G*(STRAT(2)-EEVD3);
      STREST(3)=GMODU*STRAT(3);
      [EIGPRJ,PSTRS] = SPDEC2(STREST);
      PSTRS(3)=R2G*(STRAT(4)-EEVD3);
% C Identify maximum (PSTRS1) and minimum (PSTRS3) principal stresses
      II=1;
      JJ=1;
      PSTRS1=PSTRS(II);
      PSTRS3=PSTRS(JJ);
      for I=2:3
        if(PSTRS(I)>=PSTRS1)
          II=I;
          PSTRS1=PSTRS(II);
        end %IF
        if(PSTRS(I)<PSTRS3)
          JJ=I;
          PSTRS3=PSTRS(JJ);
        end %IF
      end
      if(II~=1&&JJ~=1), MM=1; end
      if(II~=2&&JJ~=2), MM=2; end
      if(II~=3&&JJ~=3), MM=3; end
      PSTRS2=PSTRS(MM);
% C Compute trial yield function and check for plastic consistency
% C --------------------------------------------------------------
      SHMAXT=PSTRS1-PSTRS3;
      SIGMAY=PLFUN(EPBARN,NHARD,RPROPS);
      PHIA=SHMAXT-SIGMAY;
      if(PHIA/SIGMAY>TOL)
% C Plastic step: Apply return mapping
% C ==================================
        IFPLAS=1;
% C identify possible two-vector return: right or left of main plane
        SCAPRD=PSTRS1+PSTRS3-PSTRS2*R2;
        if(SCAPRD>=R0)
          RIGHT=1;
        else
          RIGHT=0;
        end %IF
% C Apply one-vector return mapping first (return to main plane)
% C ------------------------------------------------------------
        TWOVEC=0;
% C Start Newton-Raphson iterations
        for NRITER=1:MXITER
% C Compute residual derivative
          DENOM=-R4G-DPLFUN(EPBAR,NHARD,RPROPS);
% C Compute Newton-Raphson increment and update variable DGAMA
          DDGAMA=-PHIA/DENOM;
          DGAMA=DGAMA+DDGAMA;
% C Compute new residual
          EPBAR=EPBARN+DGAMA;
          SIGMAY=PLFUN(EPBAR,NHARD,RPROPS);
          SHMAX=SHMAXT-R4G*DGAMA;
          PHIA=SHMAX-SIGMAY;
% C Check convergence
          RESNOR=abs(PHIA/SIGMAY);
          if(RESNOR<=TOL)
% C Check validity of one-vector return
            S1=PSTRS1-R2G*DGAMA;
            S2=PSTRS2;
            S3=PSTRS3+R2G*DGAMA;
            DELTA=max([abs(S1),abs(S2),abs(S3)])*SMALL;
            if(S1+DELTA>=S2&&S2+DELTA>=S3)
% C converged stress is in the same sextant as trial stress -> 1-vector
% C return is valid. Update EPBAR and principal deviatoric stresses
              RSTAVA(MSTRE+1)=EPBAR;
              PSTRS1=S1;
              PSTRS3=S3;
              goto50=1;
              break
            else
% C 1-vector return is not valid - go to two-vector procedure
              goto50=0;
              break
            end %IF
          end %IF
        end
% C failure of stress update procedure
%         SUFAIL=1
% C Apply two-vector return mapping (return to corner - right or left)
% C ------------------------------------------------------------------
if(goto50==0)
        TWOVEC=1;
        DGAMA=R0;
        DGABAR=R1;
        EPBAR=EPBARN;
        SIGMAY=PLFUN(EPBARN,NHARD,RPROPS);
        SHMXTA=PSTRS1-PSTRS3;
        if(RIGHT)
          SHMXTB=PSTRS1-PSTRS2;
        else
          SHMXTB=PSTRS2-PSTRS3;
        end %IF
        PHIA=SHMXTA-SIGMAY;
        PHIB=SHMXTB-SIGMAY;
% C Start Newton-Raphson iterations
        for NRITER=1:MXITER
% C Compute residual derivative
          HSLOPE=DPLFUN(EPBAR,NHARD,RPROPS);
          DRVAA=-R4G-HSLOPE;
          DRVAB=-R2G-HSLOPE;
          DRVBA=-R2G-HSLOPE;
          DRVBB=-R4G-HSLOPE;
% C Compute Newton-Raphson increment and update variables DGAMA and DGAMB
          R1DDET=R1/(DRVAA*DRVBB-DRVAB*DRVBA);
          DDGAMA=(-DRVBB*PHIA+DRVAB*PHIB)*R1DDET;
          DDGAMB=(DRVBA*PHIA-DRVAA*PHIB)*R1DDET;
          DGAMA=DGAMA+DDGAMA;
          DGAMB=DGAMB+DDGAMB;
% C Compute new residual
          DGABAR=DGAMA+DGAMB;
          EPBAR=EPBARN+DGABAR;
          SIGMAY=PLFUN(EPBAR,NHARD,RPROPS);
          PHIA=SHMXTA-R2G*(R2*DGAMA+DGAMB)-SIGMAY;
          PHIB=SHMXTB-R2G*(DGAMA+R2*DGAMB)-SIGMAY;
% C Check convergence
          RESNOR=(abs(PHIA)+abs(PHIB))/SIGMAY;
          if(RESNOR<=TOL)
% C Update EPBAR and principal deviatoric stresses
            RSTAVA(MSTRE+1)=EPBAR;
            if(RIGHT)
              PSTRS1=PSTRS1-R2G*(DGAMA+DGAMB);
              PSTRS3=PSTRS3+R2G*DGAMA;
              PSTRS2=PSTRS2+R2G*DGAMB;
            else
              PSTRS1=PSTRS1-R2G*DGAMA;
              PSTRS3=PSTRS3+R2G*(DGAMA+DGAMB);
              PSTRS2=PSTRS2-R2G*DGAMB;
            end %IF
            break
          end %IF
        end
% C failure of stress update procedure
%         SUFAIL=1
end %goto50
% C update stress components
% C ------------------------
        PSTRS(II)=PSTRS1;
        PSTRS(JJ)=PSTRS3;
        PSTRS(MM)=PSTRS2;
        STRES(1)=PSTRS(1)*EIGPRJ(1,1)+PSTRS(2)*EIGPRJ(1,2)+P;
        STRES(2)=PSTRS(1)*EIGPRJ(2,1)+PSTRS(2)*EIGPRJ(2,2)+P;
        STRES(3)=PSTRS(1)*EIGPRJ(3,1)+PSTRS(2)*EIGPRJ(3,2);
        STRES(4)=PSTRS(3)+P;
% C and elastic engineering strain
        RSTAVA(1)=(STRES(1)-P)/R2G+EEVD3;
        RSTAVA(2)=(STRES(2)-P)/R2G+EEVD3;
        RSTAVA(3)=STRES(3)/GMODU;
        RSTAVA(4)=PSTRS(3)/R2G+EEVD3;
      else
% C Elastic step: update stress using linear elastic law
% C ====================================================
        STRES(1)=STREST(1)+P;
        STRES(2)=STREST(2)+P;
        STRES(3)=STREST(3);
        STRES(4)=PSTRS(3)+P;
% C elastic engineering strain
        RSTAVA(1)=STRAT(1);
        RSTAVA(2)=STRAT(2);
        RSTAVA(3)=STRAT(3);
        RSTAVA(4)=STRAT(4);
      end %IF
% C Update algorithmic variables before exit
% C ========================================
      DGAM(1)=DGAMA;
      DGAM(2)=DGAMB;
      LALGVA(1)=IFPLAS;
      LALGVA(2)=SUFAIL;
      LALGVA(3)=TWOVEC;
      LALGVA(4)=RIGHT;

      