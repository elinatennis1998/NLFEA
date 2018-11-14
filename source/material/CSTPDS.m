function AMATX = CSTPDS(DGAM,EPFLAG,FINCR,IPROPS,LALGVA,NTYPE,RPROPS,RSTAVA,RSTAVN,STRES)
%       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
%       PARAMETER
IPHARD=6; IPHVAR=5; MADIM=4;5; NDIM=2; NGDIM=4;
% NIPROP=3; NLALGV=6; NRALGV=4; NRSTAV=5; NSTRE=4;
NSYST=4;
% Arguments
AMATX = zeros(MADIM,MADIM);
%       LOGICAL
%      1    EPFLAG     ,LALGVA(NLALGV)
%       DIMENSION
%      1    AMATX(MADIM,MADIM) ,DGAM(NRALGV)       ,FINCR(3,3)         ,
%      2    IPROPS(NIPROP)     ,RPROPS(*)          ,RSTAVA(NRSTAV)     ,
%      3    RSTAVN(NRSTAV)     ,STRES(NSTRE)
% Local arrays and variables
%       LOGICAL
%      1    S1ACT  ,S2ACT  ,S3ACT  ,S4ACT
%       DIMENSION
%      1    APMATX(NGDIM,NGDIM),AUXMTX(NGDIM,NGDIM),AUX2ND(NDIM,NDIM)  ,
%      2    AUX4TH(NDIM,NDIM,NDIM,NDIM)            ,BEISO(3,3)         ,
%      3    BMATX(NDIM,NDIM,NSYST)                 ,DELKRO(NDIM,NDIM)  ,
%      4    DEREXP(3,3,3,3)    ,DEVPRJ(NGDIM,NGDIM),DEVSTR(NGDIM)      ,
%      5    DS0M0(NDIM,NDIM)   ,DUMATX(NDIM,NDIM,NDIM,NDIM)            ,
%      6    FE(2,2)            ,FEFETR(NGDIM,NGDIM),FEISO(2,2)         ,
%      7    FEN(2,2)           ,FETISO(2,2)        ,FETRL(2,2)         ,
%      8    FOIDS(NGDIM,NGDIM) ,FPILOG(3,3)        ,GINV(NSYST,NSYST)  ,
%      9    GMATX(NSYST,NSYST) ,IACSET(NSYST)      ,
%      O    SM0MS0(NDIM,NDIM,NSYST)                ,SOID(NGDIM)        ,
%      1    S0M0(NDIM,NDIM,NSYST)                  ,
%      2    UMATX(NDIM,NDIM,NDIM,NDIM)             ,VECM(NDIM,NSYST)   ,
%      3    VECM0(NDIM,NSYST)  ,VECS(NDIM,NSYST)   ,VECS0(NDIM,NSYST)  ,
%      4    VMATX(NGDIM,NGDIM)
AUX4TH = zeros(NDIM,NDIM,NDIM,NDIM);
S0M0 = zeros(NDIM,NDIM,NSYST);
SM0MS0 = zeros(NDIM,NDIM,NSYST);
GMATX = zeros(NSYST,NSYST);
DEVPRJ = zeros(NGDIM,NGDIM);
% C... Kroenecker delta
%       DATA
%      1    DELKRO(1,1),DELKRO(1,2)/
%      2    1.0D0      ,0.D0       /
%      3    DELKRO(2,1),DELKRO(2,2)/
%      4    0.0D0      ,1.D0       /
DELKRO = eye(NDIM);
%... fourth order (symmetric subspace) identity components stored in
%    matrix form using G matrix ordering (11,21,12,22)
%       DATA
%      1    FOIDS(1,1),FOIDS(1,2),FOIDS(1,3),FOIDS(1,4)/
%      2    1.0D0     ,0.0D0     ,0.0D0     ,0.0D0     /
%      3    FOIDS(2,1),FOIDS(2,2),FOIDS(2,3),FOIDS(2,4)/
%      4    0.0D0     ,0.5D0     ,0.5D0     ,0.0D0     /
%      5    FOIDS(3,1),FOIDS(3,2),FOIDS(3,3),FOIDS(3,4)/
%      6    0.0D0     ,0.5D0     ,0.5D0     ,0.0D0     /
%      7    FOIDS(4,1),FOIDS(4,2),FOIDS(4,3),FOIDS(4,4)/
%      8    0.0D0     ,0.0D0     ,0.0D0     ,1.0D0     /
     FOIDS = [1 0 0 0
              0 .5 .5 0
              0 .5 .5 0
              0 0 0 1];
%... second order identity components in stored in vector form using G
%    matrix ordering
%       DATA
%      1    SOID(1)  ,SOID(2)  ,SOID(3)  ,SOID(4)  /
%      2    1.0D0    ,0.0D0    ,0.0D0    ,1.D0     /
     SOID = [1 0 0 1];
%       DATA
%      1    R0   ,R1   ,R2   ,R3   /
%      2    0.0D0,1.0D0,2.0D0,3.0D0/
     R1 = 1; R2 = 2; R3 = 3;
% C***********************************************************************
% C COMPUTATION OF THE CONSISTENT SPATIAL TANGENT MODULUS 'a' FOR
% C THE PLANAR DOULBLE SLIP SINGLE CRYSTAL ELASTO-PLASTIC MODEL.
% C MODEL VALID FOR PLANE STRAIN ONLY.
% C
% C REFERENCE: Expression (16.85) 
% C***********************************************************************
% C Stop program if not plane strain
%       IF(NTYPE.NE.2)CALL ERRPRT('EI0035')
% C Retrieve some state and algorithmic variables
% C ---------------------------------------------
% C... current hardening internal variable
      HRVAR=RSTAVA(IPHVAR);
% C... current elastic deformation gradient
      FE(1,1)=RSTAVA(1);
      FE(2,1)=RSTAVA(2);
      FE(1,2)=RSTAVA(3);
      FE(2,2)=RSTAVA(4);
% C... current active slip-systems logical flags
      S1ACT=LALGVA(3);
      S2ACT=LALGVA(4);
      S3ACT=LALGVA(5);
      S4ACT=LALGVA(6);
% C... elastic deformation gradient at the beginning of the current load
% C    increment
      FEN(1,1)=RSTAVN(1);
      FEN(2,1)=RSTAVN(2);
      FEN(1,2)=RSTAVN(3);
      FEN(2,2)=RSTAVN(4);
% C Retrieve material properties
% C ----------------------------
% C... neo-Hookean constants
      GMODU=RPROPS(2);
      BULK=RPROPS(3);
% C... initial system orientation
      THETA=RPROPS(4);
% C... relative angle between systems
      BETA=RPROPS(5);
% C... number of sampling points on hardening curve
      NHARD=IPROPS(3);
      RPROPS2 = reshape(RPROPS(IPHARD:IPHARD-1+2*NHARD),2,NHARD);
% C Set some constants
% C ------------------
      R1D3=R1/R3;
      R2D3=R2*R1D3;
% C Assemble deviatoric projection tensor (use G matrix ordering)
      for I=1:NGDIM
        for J=1:NGDIM
          DEVPRJ(I,J)=FOIDS(I,J)-R1D3*SOID(I)*SOID(J);
        end
      end
% C Get current Cauchy deviatoric stress and Cauchy hydrostatic pressure
      P=R1D3*(STRES(1)+STRES(2)+STRES(4));
% C... use G matrix component ordering to store in-plane deviatoric
% C    Cauchy stress components
      DEVSTR(1)=STRES(1)-P;
      DEVSTR(2)=STRES(3);
      DEVSTR(3)=STRES(3);
      DEVSTR(4)=STRES(2)-P;
% C Compute isochoric component of Fe and Be
      DETFE=FE(1,1)*FE(2,2)-FE(1,2)*FE(2,1);
      FACTOR=DETFE^(-R1D3);
      FEISO(1,1)=FACTOR*FE(1,1);
      FEISO(1,2)=FACTOR*FE(1,2);
      FEISO(2,1)=FACTOR*FE(2,1);
      FEISO(2,2)=FACTOR*FE(2,2);
%       CALL RVZERO(BEISO,9)
      BEISO = zeros(3,3);
      for I=1:2
        for J=1:2
          for K=1:2
            BEISO(I,J)=BEISO(I,J)+FEISO(I,K)*FEISO(J,K);
          end
        end
      end
      BEISO(3,3)=FACTOR*FACTOR;
% C Trace of isochoric component of Be
      TRBISO=BEISO(1,1)+BEISO(2,2)+BEISO(3,3);
% C
% C
% C Compute ELASTIC tangent modulus
% C ===============================
% C
      GFAC=R2D3*GMODU*TRBISO/DETFE;
      BULFAC=BULK/DETFE;
      R2P=R2*P;
% C... assemble tensorially compact part
      for I=1:NGDIM
        for J=1:NGDIM
          AMATX(I,J)=BULFAC*SOID(I)*SOID(J)-R2P*FOIDS(I,J)+ ...
                     GFAC*DEVPRJ(I,J)- ...
                     R2D3*(DEVSTR(I)*SOID(J)+SOID(I)*DEVSTR(J));
        end
      end
% C... add non-compact part: delta_ik sigma_jl
      AMATX(1,1)=AMATX(1,1)+STRES(1);
      AMATX(3,1)=AMATX(3,1)+STRES(3);
      AMATX(2,2)=AMATX(2,2)+STRES(1);
      AMATX(4,2)=AMATX(4,2)+STRES(3);
      AMATX(1,3)=AMATX(1,3)+STRES(3);
      AMATX(3,3)=AMATX(3,3)+STRES(2);
      AMATX(2,4)=AMATX(2,4)+STRES(3);
      AMATX(4,4)=AMATX(4,4)+STRES(2);
% C
% C
      if(EPFLAG)
% C
% C Compute and add algorithm-consistent PLASTIC contribution to spatial
% C tangent modulus
% C ====================================================================
% C
% C Compute individual terms needed to assemble the plastic contribution
% C --------------------------------------------------------------------
% C
% C Last elastic trial deformation gradient
%         CALL RVZERO(FETRL,4)
        FETRL = zeros(2,2);
        for I=1:2
          for J=1:2
            for K=1:2
              FETRL(I,J)=FETRL(I,J)+FINCR(I,K)*FEN(K,J);
            end
          end
        end
% C... isochoric component
        DETFET=FETRL(1,1)*FETRL(2,2)-FETRL(1,2)*FETRL(2,1);
        VOLFAC=DETFET^(-R1D3);
        FETISO(1,1)=VOLFAC*FETRL(1,1);
        FETISO(2,1)=VOLFAC*FETRL(2,1);
        FETISO(1,2)=VOLFAC*FETRL(1,2);
        FETISO(2,2)=VOLFAC*FETRL(2,2);
% C Assemble relevant fourth order tensor
        for I=1:NDIM
          for J=1:NDIM
            for K=1:NDIM
              for L=1:NDIM
                AUX4TH(I,J,K,L)=FEISO(I,L)*FETISO(J,K);
              end
            end
          end
        end
% C... rearrange in matrix form with G matrix component ordering
%         CALL ARRGO2(AUX4TH,FEFETR)
        FEFETR = ARRGO2(AUX4TH);
% C
% C Compute exact exponential map derivative
% C ----------------------------------------
% C
% C retrieve information on current set of active slip-systems and set up
% C corresponding initial slip-system vectors
        IACSYS=0;
        if(S1ACT)
% C... system 1:
          IACSYS=IACSYS+1;
          IACSET(IACSYS)=1;
          VECS0(1,1)=cos(THETA);
          VECS0(2,1)=sin(THETA);
          VECM0(1,1)=-sin(THETA);
          VECM0(2,1)=cos(THETA);
        end
        if(S2ACT)
% C... system 2:
          IACSYS=IACSYS+1;
          IACSET(IACSYS)=2;
          VECS0(1,2)=cos(THETA+BETA);
          VECS0(2,2)=sin(THETA+BETA);
          VECM0(1,2)=-sin(THETA+BETA);
          VECM0(2,2)=cos(THETA+BETA);
        end
        if(S3ACT)
          IACSYS=IACSYS+1;
          IACSET(IACSYS)=3;
% C... system 3:
          VECS0(1,3)=-cos(THETA);
          VECS0(2,3)=-sin(THETA);
          VECM0(1,3)=-sin(THETA);
          VECM0(2,3)=cos(THETA);
        end
        if(S4ACT)
% C... system 4:
          IACSYS=IACSYS+1;
          IACSET(IACSYS)=4;
          VECS0(1,4)=-cos(THETA+BETA);
          VECS0(2,4)=-sin(THETA+BETA);
          VECM0(1,4)=-sin(THETA+BETA);
          VECM0(2,4)=cos(THETA+BETA);
        end
% C... number of currently active systems
        NACSYS=IACSYS;
% C Compute current elastic push forward of the slip-system vectors
% C of the active systems
%         CALL RVZERO(VECS,NDIM*NSYST)
%         CALL RVZERO(VECM,NDIM*NSYST)
        VECS = zeros(NDIM,NSYST);
        VECM = zeros(NDIM,NSYST);
        for I=1:NDIM
          for J=1:NDIM
            for II=1:NACSYS
              ISYST=IACSET(II);
              VECS(I,ISYST)=VECS(I,ISYST)+FEISO(I,J)*VECS0(J,ISYST);
              VECM(I,ISYST)=VECM(I,ISYST)+FEISO(I,J)*VECM0(J,ISYST);
            end
          end
        end
% C Current slope of Taylor hardening curve
%         HSLOPE=DPLFUN(HRVAR,NHARD,RPROPS(IPHARD));
        HSLOPE=DPLFUN(HRVAR,NHARD,RPROPS2);
% C Compute logarithm of inverse of incremental plastic deformation
% C gradient by summing up contributions from each active slip system
%         CALL RVZERO(FPILOG,9)
        FPILOG = zeros(3,3);
        for II=1:NACSYS
          ISYST=IACSET(II);
          for I=1:2
            for J=1:2
              FPILOG(I,J)=FPILOG(I,J)- ...
                          DGAM(ISYST)*VECS0(I,ISYST)*VECM0(J,ISYST);
            end
          end
        end
% C... and the corresponding exact derivative of the exponential map
%         CALL DEXPMP(   DEREXP      ,NOCONV     ,FPILOG     )
        DEREXP = DEXPMP(FPILOG); %Box 16.3 (iv), (16.70), tensor E
% C
% C Compute jacobian of non-linear system of return mapping equations
% C -----------------------------------------------------------------
% C compute some preliminary matrices
%         CALL RVZERO(BMATX,NDIM*NDIM*NSYST)
        BMATX = zeros(NDIM,NDIM,NSYST);
        for II=1:NACSYS
          ISYST=IACSET(II);
          for I=1:2
            for J=1:2
              S0M0(I,J,ISYST)=VECS0(I,ISYST)*VECM0(J,ISYST);
              SM0MS0(I,J,ISYST)=VECS(I,ISYST)*VECM0(J,ISYST)+ ...
                                VECM(I,ISYST)*VECS0(J,ISYST);
            end
          end
%           CALL RVZERO(DS0M0,NDIM*NDIM)
          DS0M0 = zeros(NDIM,NDIM);
          for I=1:2
            for J=1:2
              for K=1:2
                for L=1:2
                  DS0M0(I,J)=DS0M0(I,J)+ ...
                             DEREXP(I,J,K,L)*S0M0(K,L,ISYST);
                end
              end
            end
          end
          for I=1:2
            for J=1:2
              for K=1:2
                BMATX(I,J,ISYST)=BMATX(I,J,ISYST)+ ...
                                 FETISO(I,K)*DS0M0(K,J);
              end
            end
          end
        end
% C Assemble exact jacobian of non-linear system
        SM0MS0v = reshape(SM0MS0,NDIM*NDIM,NSYST);
        BMATXv = reshape(BMATX,NDIM*NDIM,NSYST);
        for II=1:NACSYS
          ISYST=IACSET(II);
          for JJ=1:NACSYS
            JSYST=IACSET(JJ);
%             GMATX(II,JJ)=GMODU* ...
%                SCAPRD(SM0MS0(1,1,ISYST),BMATX(1,1,JSYST),NDIM*NDIM)+ ...
%                HSLOPE;
            GMATX(II,JJ)=GMODU* ...
               SM0MS0v(:,ISYST)'*BMATXv(:,JSYST)+ ...
               HSLOPE; %Box 16.3 (iv), (16.70)
          end
        end
% C Invert jacobian: Note that for the double slip model only one or two
% C systems may be active
        if(NACSYS == 1)
          if(GMATX(1,1) == 0)
            error('gmat singular');
          end
          GINV(1,1)=R1/GMATX(1,1);
        elseif(NACSYS == 2)
          DETG=GMATX(1,1)*GMATX(2,2)-GMATX(1,2)*GMATX(2,1);
          if(DETG == 0)
            error('gmat singular');
          end
          DETGIN=R1/DETG;
          GINV(1,1)=GMATX(2,2)*DETGIN;
          GINV(2,2)=GMATX(1,1)*DETGIN;
          GINV(1,2)=-GMATX(1,2)*DETGIN;
          GINV(2,1)=-GMATX(2,1)*DETGIN;
        end
% C Compute U matrix
% C ----------------
%         CALL RVZERO(UMATX,NDIM*NDIM*NDIM*NDIM)
        UMATX = zeros(NDIM,NDIM,NDIM,NDIM);
        for I=1:NDIM
          for J=1:NDIM
            for K=1:NDIM
              for L=1:NDIM
                for II=1:NACSYS
                  ISYST=IACSET(II);
                  for JJ=1:NACSYS
                    JSYST=IACSET(JJ);
                    UMATX(I,J,K,L)=UMATX(I,J,K,L)+S0M0(I,J,ISYST)* ...
                                   GINV(II,JJ)*SM0MS0(K,L,JSYST);
                  end
                end
              end
            end
          end
        end
% C Compute product [D:U]
%         CALL RVZERO(DUMATX,NDIM*NDIM*NDIM*NDIM)
        DUMATX = zeros(NDIM,NDIM,NDIM,NDIM);
        for I=1:NDIM
          for J=1:NDIM
            for K=1:NDIM
              for L=1:NDIM
                for M=1:NDIM
                  for N=1:NDIM
                    DUMATX(I,J,K,L)=DUMATX(I,J,K,L)+ ...
                                    DEREXP(I,J,M,N)*UMATX(M,N,K,L);
                  end
                end
              end
            end
          end
        end
% C... and the contribution to a^p involving the product D:U
%         CALL RVZERO(AUX4TH,NDIM*NDIM*NDIM*NDIM)
        AUX4TH = zeros(NDIM,NDIM,NDIM,NDIM);
        for I=1:NDIM
          for J=1:NDIM
            for K=1:NDIM
              for L=1:NDIM
                for M=1:NDIM
                  AUX4TH(I,J,K,L)=AUX4TH(I,J,K,L)+ ...
                                  DUMATX(I,J,K,M)*FEISO(L,M);
                end
              end
            end
          end
        end
%         CALL RVZERO(AUX2ND,NDIM*NDIM)
        AUX2ND = zeros(NDIM,NDIM);
        for I=1:NDIM
          for J=1:NDIM
            for K=1:NDIM
              for L=1:NDIM
                AUX2ND(I,J)=AUX2ND(I,J)+DUMATX(I,J,K,L)*FEISO(K,L);
              end
            end
          end
        end
        for I=1:NDIM
          for J=1:NDIM
            for K=1:NDIM
              for L=1:NDIM
                AUX4TH(I,J,K,L)=AUX4TH(I,J,K,L)- ...
                                R1D3*AUX2ND(I,J)*DELKRO(K,L);
              end
            end
          end
        end
% C... rearrange in matrix form
%         CALL ARRGO2(AUX4TH,VMATX)
        VMATX = ARRGO2(AUX4TH);
%         CALL RVZERO(AUXMTX,NGDIM*NGDIM)
        AUXMTX = zeros(NGDIM,NGDIM);
        for I=1:NGDIM
          for J=1:NGDIM
            for K=1:NGDIM
              AUXMTX(I,J)=AUXMTX(I,J)+FEFETR(I,K)*VMATX(K,J);
            end
          end
        end
% C Compute plastic contribution
%         CALL RVZERO(APMATX,NGDIM*NGDIM)
        APMATX = zeros(NGDIM,NGDIM);
        AUX=R2*GMODU*GMODU/DETFE;
        for I=1:NGDIM
          for J=1:NGDIM
            for K=1:NGDIM
              APMATX(I,J)=APMATX(I,J)-AUX*DEVPRJ(I,K)*AUXMTX(K,J);
            end
          end
        end
% C Add plastic contribution to spatial tangent modulus
        for I=1:NGDIM
          for J=1:NGDIM
            AMATX(I,J)=AMATX(I,J)+APMATX(I,J);
          end
        end
% C
      end
% C
%       RETURN
%       END
