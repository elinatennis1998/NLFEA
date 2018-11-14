function [DGAMA,LALGVA,RSTAVA,STRES] = SUVM(IPROPS,NTYPE,EPROPS,RPROPS,RSTAVA,STRAT)
% C***********************************************************************
% C STATE UPDATE PROCEDURE FOR THE VON MISES ELASTO-PLASTIC MATERIAL MODEL
% C WITH NON-LINEAR (PIECEWISE LINEAR) ISOTROPIC HARDENING:
% C IMPLICIT ELASTIC PREDICTOR/RETURN MAPPING ALGORITHM (BOXES 7.3-4).
% C PLANE STRAIN AND AXISYMMETRIC IMPLEMENTATIONS.
% C
% C REFERENCE: Section 7.3.5
% C***********************************************************************
%       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
%       PARAMETER(IPHARD=4  ,MSTRE=4)
MSTRE=4;
%       LOGICAL IFPLAS, LALGVA(2), SUFAIL
%       DIMENSION
%      1    IPROPS(*)          ,RPROPS(*)          ,RSTAVA(MSTRE+1)    ,
%      2    STRAT(MSTRE)       ,STRES(MSTRE)
%       DIMENSION
%      1    EET(MSTRE)
%       DATA
%      1    R0   ,RP5  ,R1   ,R2   ,R3   ,TOL   / 
%      2    0.0D0,0.5D0,1.0D0,2.0D0,3.0D0,1.D-06/
%       DATA MXITER / 50 /
STRES = zeros(MSTRE,1);
R0 = 0.0d0;
RP5 = 0.5d0;
R1 = 1.0d0;
R2 = 2.0d0;
R3 = 3.0d0;
TOL = 1.D-06;
MXITER = 50;

% C Stop program if neither plane strain nor axisymmetric state
      if (NTYPE~=2&&NTYPE~=3), return, end
% C Initialise some algorithmic and internal variables
      DGAMA=R0;
      IFPLAS=0;
      SUFAIL=0;
      EPBARN=RSTAVA(MSTRE+1);
% C Set some material properties
      YOUNG=EPROPS(1);
      POISS=EPROPS(2);
      NHARD=IPROPS(1);
% C Shear and bulk moduli and other necessary constants
      GMODU=YOUNG/(R2*(R1+POISS));
      BULK=YOUNG/(R3*(R1-R2*POISS));
      R2G=R2*GMODU;
      R3G=R3*GMODU;
% C Elastic predictor: Compute elastic trial state
% C ----------------------------------------------
% C Volumetric strain and pressure stress
      EEV=STRAT(1)+STRAT(2)+STRAT(4);
      P=BULK*EEV;
% C Elastic trial deviatoric strain
      EEVD3=EEV/R3;
      EET(1)=STRAT(1)-EEVD3;
      EET(2)=STRAT(2)-EEVD3;
      EET(4)=STRAT(4)-EEVD3;
% C Convert engineering shear component into physical component
      EET(3)=STRAT(3)/R2;
% C Compute trial effective stress and uniaxial yield stress
      VARJ2T=R2G*R2G*(EET(3)*EET(3)+RP5*(EET(1)*EET(1)+ ...
                           EET(2)*EET(2)+EET(4)*EET(4)));
      QTRIAL=sqrt(R3*VARJ2T);
      SIGMAY=PLFUN(EPBARN,NHARD,RPROPS);
% C Check for plastic admissibility
% C -------------------------------
      PHI=QTRIAL-SIGMAY;
      if(PHI/SIGMAY>TOL)
% C Plastic step: Apply return mapping - use Newton-Raphson algorithm
% C               to solve the return mapping equation (Box 7.4)
% C -------------------------------------------------------------------
        IFPLAS=1;
        EPBAR=EPBARN;
        for NRITER=1:MXITER
% C Compute residual derivative
          DENOM=-R3G-DPLFUN(EPBAR,NHARD,RPROPS);
% C Compute Newton-Raphson increment and update variable DGAMA
          DDGAMA=-PHI/DENOM;
          DGAMA=DGAMA+DDGAMA;
% C Compute new residual
          EPBAR=EPBAR+DDGAMA;
          SIGMAY=PLFUN(EPBAR,NHARD,RPROPS);
          PHI=QTRIAL-R3G*DGAMA-SIGMAY;
% C Check convergence
          RESNOR=abs(PHI/SIGMAY);
          if(RESNOR<=TOL)
% C update accumulated plastic strain
            RSTAVA(MSTRE+1)=EPBAR;
% C update stress components
            FACTOR=R2G*(R1-R3G*DGAMA/QTRIAL);
            STRES(1)=FACTOR*EET(1)+P;
            STRES(2)=FACTOR*EET(2)+P;
            STRES(3)=FACTOR*EET(3);
            STRES(4)=FACTOR*EET(4)+P;
% C compute converged elastic (engineering) strain components
            FACTOR=FACTOR/R2G;
            RSTAVA(1)=FACTOR*EET(1)+EEVD3;
            RSTAVA(2)=FACTOR*EET(2)+EEVD3;
            RSTAVA(3)=FACTOR*EET(3)*R2;
            RSTAVA(4)=FACTOR*EET(4)+EEVD3;
            break
          end %IF
        end
% C reset failure flag and print warning message if the algorithm fails
        if(NRITER==MXITER)
        SUFAIL=1;
        end
      else
% C Elastic step: Update stress using linear elastic law
% C ----------------------------------------------------
        STRES(1)=R2G*EET(1)+P;
        STRES(2)=R2G*EET(2)+P;
        STRES(3)=R2G*EET(3);
        STRES(4)=R2G*EET(4)+P;
% C elastic engineering strain
        RSTAVA(1)=STRAT(1);
        RSTAVA(2)=STRAT(2);
        RSTAVA(3)=STRAT(3);
        RSTAVA(4)=STRAT(4);
      end %IF
%   999 CONTINUE
% C Update some algorithmic variables before exit
      LALGVA(1)=IFPLAS;
      LALGVA(2)=SUFAIL;
%       RETURN
%       END
