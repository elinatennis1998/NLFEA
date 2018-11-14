function [DGAM,LALGVA,RSTAVA,STRES] = SUPDSC(FINCR,IPROPS,NTYPE,RPROPS,RSTAVA)
%
% The for loops are equivalent to the matrix multiplications
%
%       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
%       PARAMETER
IPHARD=6; IPHVAR=5; MTRIAL=5; NDIM=2; NSYST=4;
% NIPROP=3; NLALGV=6; NRALGV=4; NRSTAV=5; NSTRE=4;
% C Arguments
%       LOGICAL
%      1    LALGVA(NLALGV)
%       DIMENSION
%      1    DGAM(NRALGV)       ,FINCR(3,3)         ,IPROPS(NIPROP)     ,
%      2    RPROPS(*)          ,RSTAVA(NRSTAV)     ,STRES(NSTRE)
% C Local arrays and variables
%       LOGICAL
%      1    IFPLAS ,NOCONV ,SUFAIL ,S1ACT  ,S2ACT  ,S3ACT  ,S4ACT
%       DIMENSION
%      1    BEDEV(4)           ,BEISO(3,3)      ,BMATX(NDIM,NDIM,NSYST),
%      2    DDGAM(NSYST)       ,DEREXP(3,3,3,3)    ,DS0M0(NDIM,NDIM)   ,
%      3    FEISO(2,2)         ,FEN(2,2)           ,FETISO(2,2)        ,
%      4    FETRL(2,2)         ,FPILOG(3,3)        ,FPINCI(3,3)        ,
%      5    GINV(NSYST,NSYST)  ,GMATX(NSYST,NSYST),IACSET(NSYST,MTRIAL),
%      6    IPACT(0:5)         ,NACSYS(MTRIAL)     ,PHI(NSYST)         ,
%      7    SCHMID(NSYST)      ,SM0MS0(NDIM,NDIM,NSYST)                ,
%      8    S0M0(NDIM,NDIM,NSYST),VECM(NDIM,NSYST) ,VECM0(NDIM,NSYST)  ,
%      9    VECS(NDIM,NSYST)   ,VECS0(NDIM,NSYST)
BEDEV = zeros(4,1);
BEISO = zeros(3,3);
% BMATX = zeros(NDIM,NDIM,NSYST);
% DDGAM = zeros(NSYST);
% DEREXP = zeros(3,3,3,3);
% DS0M0 = zeros(NDIM,NDIM);
FEISO = zeros(2,2);
FEN = zeros(2,2);
FETISO = zeros(2,2);
% FETRL = zeros(2,2);
% FPILOG = zeros(3,3);
FPINCI = zeros(3,3);
GINV = zeros(NSYST,NSYST);
GMATX = zeros(NSYST,NSYST);
IACSET = zeros(NSYST,MTRIAL);
IPACT = [4        ,1        ,2        ,3        ,4        ,1];
NACSYS = zeros(MTRIAL,1);
PHI = zeros(NSYST,1);
SCHMID = zeros(NSYST,1);
SM0MS0 = zeros(NDIM,NDIM,NSYST);
S0M0 = zeros(NDIM,NDIM,NSYST);
% VECM = zeros(NDIM,NSYST);
VECM0 = zeros(NDIM,NSYST);
% VECS = zeros(NDIM,NSYST);
VECS0 = zeros(NDIM,NSYST);
%       DATA
R0 = 0;
R1 = 1;
R3 = 3;
SMALL = 1.D-10;
TOL = 1.D-08;
MXITER = 50;
% C***********************************************************************
% C STRESS UPDATE PROCEDURE FOR THE ANISOTROPIC PLANAR DOUBLE-SLIP SINGLE
% C CRYSTAL ELASTO-PLASTIC MODEL WITH PIECE-WISE LINEAR TAYLOR ISOTROPIC
% C HARDENING:
% C MULTI-SURFACE TYPE IMPLICIT ELASTIC PREDICTOR/RETURN MAPPING ALGORITHM
% C BASED ON THE EXPONENTIAL MAP APPROXIMATION OF THE PLASTIC FLOW RULE.
% C
% C REFERENCE: Section 16.6.2
% C            Boxes 16.2-3
% C***********************************************************************
% C Stops program if not plane strain
      if(NTYPE ~= 2)
          error('EI0036')
      end
% C Initialize some algorithmic and internal variables
%       CALL RVZERO(DGAM,NSYST)
      DGAM(1:NSYST) = 0;
      IFPLAS=0;%.FALSE.
      SUFAIL=0;%.FALSE.
      S1ACT=0;%.FALSE.
      S2ACT=0;%.FALSE.
      S3ACT=0;%.FALSE.
      S4ACT=0;%.FALSE.
% % C... hardening internal variable
      HRVARN=RSTAVA(IPHVAR);
      HRVAR=HRVARN;
% C... elastic deformation gradient
      FEN(1,1)=RSTAVA(1);
      FEN(2,1)=RSTAVA(2);
      FEN(1,2)=RSTAVA(3);
      FEN(2,2)=RSTAVA(4);
% C Retrieve material properties
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
% C Set up initial slip systems vectors
% C... system 1:
      VECS0(1,1)=cos(THETA);
      VECS0(2,1)=sin(THETA);
      VECM0(1,1)=-sin(THETA);
      VECM0(2,1)=cos(THETA);
% C... system 2:
      VECS0(1,2)=cos(THETA+BETA);
      VECS0(2,2)=sin(THETA+BETA);
      VECM0(1,2)=-sin(THETA+BETA);
      VECM0(2,2)=cos(THETA+BETA);
% C... system 3:
      VECS0(1,3)=-VECS0(1,1);
      VECS0(2,3)=-VECS0(2,1);
      VECM0(1,3)=VECM0(1,1);
      VECM0(2,3)=VECM0(2,1);
% C... system 4:
      VECS0(1,4)=-VECS0(1,2);
      VECS0(2,4)=-VECS0(2,2);
      VECM0(1,4)=VECM0(1,2);
      VECM0(2,4)=VECM0(2,2);
% C Set some constants
      R1D3=R1/R3;
% C Compute elastic trial state
% C ---------------------------
% C Elastic trial deformation gradient
%       CALL RVZERO(FETRL,NDIM*NDIM)
      FETRL(1:2,1:2) = FINCR(1:2,1:2)*FEN(1:2,1:2); %Box 16.2 (i), (16.27)
%       FETRL = zeros(NDIM,NDIM);%(1:NDIM,1:NDIM) = 0;
%       for I=1:NDIM %30
%         for J=1:NDIM %20
%           for K=1:NDIM %10
%             FETRL(I,J)=FETRL(I,J)+FINCR(I,K)*FEN(K,J);
%           end
%         end
%       end
%    10     CONTINUE
%    20   CONTINUE
%    30 CONTINUE
% C Perform isochoric/volumetric split of elastic trial def. grad.
      DETFET=FETRL(1,1)*FETRL(2,2)-FETRL(1,2)*FETRL(2,1);
      VOLFAC=DETFET^(-R1D3);
      FETISO=VOLFAC*FETRL; %Box 16.2 (ii), (16.58)
%       FETISO(1,1)=VOLFAC*FETRL(1,1);
%       FETISO(2,1)=VOLFAC*FETRL(2,1);
%       FETISO(1,2)=VOLFAC*FETRL(1,2);
%       FETISO(2,2)=VOLFAC*FETRL(2,2);
% C Check plastic consistency
% C -------------------------
% C Compute yield functions values
% C... elastic push forward of slip-systems vectors
      VECS = FETISO*VECS0; %Box 16.2 (iii), (16.62)
      VECM = FETISO*VECM0;
% %       CALL RVZERO(VECS,NDIM*NSYST);
%       VECS = zeros(NDIM,NSYST);
% %       CALL RVZERO(VECM,NDIM*NSYST);
%       VECM = zeros(NDIM,NSYST);
%       for I=1:NDIM %60
%         for J=1:NDIM %50
%           for ISYST=1:NSYST %40
%             VECS(I,ISYST)=VECS(I,ISYST)+FETISO(I,J)*VECS0(J,ISYST);
%             VECM(I,ISYST)=VECM(I,ISYST)+FETISO(I,J)*VECM0(J,ISYST);
%           end
%         end
%       end
%    40     CONTINUE
%    50   CONTINUE
%    60 CONTINUE
% C... current resolved yield stress
      RYIELD=PLFUN(HRVAR,NHARD,RPROPS2);
% C... elastic trial Schmid resolved shear stresses
      SCHMID=GMODU*diag(VECS'*VECM); %Box 16.2 (iii), (16.61)
      for ISYST=1:NSYST %70
        SCHMID(ISYST)=GMODU*VECS(:,ISYST)'*VECM(:,ISYST);
%         SCHMID(ISYST)=GMODU*SCAPRD(VECS(:,ISYST),VECM(:,ISYST),2);
      end
%    70 CONTINUE
% C Check consistency
% C... compute yield functions values and determine set of active systems
% C    at elastic trial state
      IACSYS=0;
      for ISYST=1:NSYST %80 %Box 16.2 (iv), (16.41)
        PHI(ISYST)=SCHMID(ISYST)-RYIELD;
        if(PHI(ISYST)/RYIELD > TOL) %THEN
          IFPLAS=1;%.TRUE.
          IACSYS=IACSYS+1;
          IACSET(IACSYS,1)=ISYST;
        end %IF
      end
%    80 CONTINUE
      NACSYS(1)=IACSYS;
% C... define the other possible tentative sets of active systems
      if(NACSYS(1) == 1) %THEN
        NTENT=3;
        NACSYS(2)=2;
        IACSET(1,2)=IACSET(1,1);
        IACSET(2,2)=IPACT(IACSET(1,1)-1+1);
        NACSYS(3)=2;
        IACSET(1,3)=IACSET(1,1);
        IACSET(2,3)=IPACT(IACSET(1,1)+1+1);
      elseif(NACSYS(1) == 2) %THEN
        NTENT=5;
        NACSYS(2)=1;
        IACSET(1,2)=IACSET(1,1);
        NACSYS(3)=1;
        IACSET(1,3)=IACSET(2,1);
        if(IACSET(1,1) == 1 && IACSET(2,1) == 4) %THEN
          NACSYS(4)=2;
          IACSET(1,4)=1;
          IACSET(2,4)=2;
          NACSYS(5)=2;
          IACSET(1,5)=3;
          IACSET(2,5)=4;
        else
          NACSYS(4)=2;
          IACSET(1,4)=IACSET(1,1);
          IACSET(2,4)=IPACT(IACSET(1,1)-1+1);
          NACSYS(5)=2;
          IACSET(1,5)=IACSET(2,1);
          IACSET(2,5)=IPACT(IACSET(2,1)+1+1);
        end %IF
      end %IF
      if(IFPLAS) %THEN
% C Plastic step: Apply return mapping
% C ==================================
% C Loop over the tentative sets of active systems
% C ----------------------------------------------
        ALLDONE = 0;
        for ITENT=1:NTENT %420
% C re-set elastic push-forward of slip systems vectors
%           CALL RVZERO(VECS,NDIM*NSYST)
          VECS=FETISO*VECS0; %??????????
          VECM=FETISO*VECM0;
%           VECS = zeros(NDIM,NSYST);
% %           CALL RVZERO(VECM,NDIM*NSYST);
%           VECM = zeros(NDIM,NSYST);
%           for I=1:NDIM %110
%             for J=1:NDIM %100
%               for ISYST=1:NSYST %90
%                 VECS(I,ISYST)=VECS(I,ISYST)+FETISO(I,J)*VECS0(J,ISYST);
%                 VECM(I,ISYST)=VECM(I,ISYST)+FETISO(I,J)*VECM0(J,ISYST);
%               end
%             end
%           end
%    90         CONTINUE
%   100       CONTINUE
%   110     CONTINUE
% C re-set hardening variable
          HRVAR=HRVARN;
% C re-set yield function values at trial state
          RYIELD=PLFUN(HRVAR,NHARD,RPROPS2);
% C... elastic trial Schmid resolved shear stresses
      SCHMID=GMODU*diag(VECS'*VECM); %
      PHI=SCHMID-RYIELD;
%           for ISYST=1:NSYST %120
%             SCHMID(ISYST)=GMODU*VECS(:,ISYST)'*VECM(:,ISYST);
%             SCHMID(ISYST)=GMODU*SCAPRD(VECS(1,ISYST),VECM(1,ISYST),NDIM);
%             PHI(ISYST)=SCHMID(ISYST)-RYIELD;
%           end
%   120     CONTINUE
% C Start Newton-Raphson iterations for plastic multipliers
%           CALL RVZERO(DGAM,NSYST)
          DGAM(1:NSYST) = 0;
          for NRITER=1:MXITER %410
            HSLOPE=DPLFUN(HRVAR,NHARD,RPROPS2);
            if (NRITER == 1)
%                 CALL RVZERO(FPILOG,9)
                FPILOG = zeros(3,3);
            end
%             [DEREXP,NOCONV] = DEXPMP(FPILOG); %Box 16.3 (iv), (16.70), tensor E
%             CALL RVZERO(BMATX,NDIM*NDIM*NSYST)
            BMATX = zeros(NDIM,NDIM,NSYST);
            for II=1:NACSYS(ITENT) %220
              ISYST=IACSET(II,ITENT);
              S0M0(:,:,ISYST)=VECS0(:,ISYST)*VECM0(:,ISYST)'; %Box 16.3 (iv), (16.70)
              SM0MS0(:,:,ISYST)=VECS(:,ISYST)*VECM0(:,ISYST)'+ ...
                                VECM(:,ISYST)*VECS0(:,ISYST)'; %Box 16.3 (iv), (16.70), first half of : term
%               for I=1:NDIM %140
%                 for J=1:NDIM %130
%                   S0M0(I,J,ISYST)=VECS0(I,ISYST)*VECM0(J,ISYST);
%                   SM0MS0(I,J,ISYST)=VECS(I,ISYST)*VECM0(J,ISYST)+ ...
%                                     VECM(I,ISYST)*VECS0(J,ISYST);
%                 end
%               end
%   130           CONTINUE
%   140         CONTINUE
%               CALL RVZERO(DS0M0,NDIM*NDIM)
%               DS0M0 = zeros(NDIM,NDIM);
%               for I=1:NDIM %180
%                 for J=1:NDIM %170
%                   for K=1:NDIM %160
%                     for L=1:NDIM %150
%                       DS0M0(I,J)=DS0M0(I,J)+ ...
%                                  DEREXP(I,J,K,L)*S0M0(K,L,ISYST);
%                     end
%                   end
%                 end
%               end
              S0M0a = [S0M0(:,:,ISYST) zeros(2,1); zeros(1,3)];
              [~,DS0M0] = dexpmt(FPILOG,S0M0a); %%%%% CHECK TO SEE IF THIS AGREES
%   150               CONTINUE
%   160             CONTINUE
%   170           CONTINUE
%   180         CONTINUE
              BMATX(:,:,ISYST)=FETISO(1:2,1:2)*DS0M0(1:2,1:2); %Box 16.3 (iv), (16.70), second half of : term
%               for I=1:NDIM %210
%                 for J=1:NDIM %200
%                   for K=1:NDIM %190
%                     BMATX(I,J,ISYST)=BMATX(I,J,ISYST)+ ...
%                                      FETISO(I,K)*DS0M0(K,J);
%                   end
%                 end
%               end
%   190             CONTINUE
%   200           CONTINUE
%   210         CONTINUE
            end
%   220       CONTINUE
            SM0MS0v = reshape(SM0MS0,NDIM*NDIM,NSYST);
            BMATX = reshape(BMATX,NDIM*NDIM,NSYST);
            for II=1:NACSYS(ITENT) %240
              ISYST=IACSET(II,ITENT);
% C Compute exact jacobian of non-linear system of equations
              for JJ=1:NACSYS(ITENT) %230
                JSYST=IACSET(JJ,ITENT);

                GMATX(II,JJ)=GMODU* ...
                   SM0MS0v(:,ISYST)'*BMATX(:,JSYST)+ ...
                   HSLOPE; %Box 16.3 (iv), (16.70)
%                 GMATX(II,JJ)=GMODU* ...
%                    SCAPRD(SM0MS0(:,:,ISYST),BMATX(:,:,JSYST),NDIM*NDIM)+ ...
%                    HSLOPE;

              end
            end
%   230         CONTINUE
%   240       CONTINUE
% C Invert jacobian: Note that for the double slip model only one or two
% C systems may be active
            if(NACSYS(ITENT) == 1) %THEN %Box 16.3 (v), (16.70)
              JACSIG = 0;
              if(GMATX(1,1) < SMALL) %THEN
% C... jacobian is singular: Try another active set or exit if
% C                          all possible sets have already been tried
%                 GOTO 420
                JACSIG = 1;
                break
              end %IF
              GINV(1,1)=R1/GMATX(1,1);
            elseif(NACSYS(ITENT) == 2) %THEN
              DETG=GMATX(1,1)*GMATX(2,2)-GMATX(1,2)*GMATX(2,1);
              JACSIG = 0;
              if (DETG < SMALL) %THEN
% C... jacobian is singular: Try another active set or exit if
% C                          all possible sets have already been tried
%                 GOTO 420
                JACSIG = 1;
                break
              end %IF
              DETGIN=R1/DETG;
              GINV(1,1)=GMATX(2,2)*DETGIN;
              GINV(2,2)=GMATX(1,1)*DETGIN;
              GINV(1,2)=-GMATX(1,2)*DETGIN;
              GINV(2,1)=-GMATX(2,1)*DETGIN;
            end %IF
% C Apply Newton-Raphson correction to plastic multipliers
%             CALL RVZERO(DDGAM,NSYST)
            DDGAM = zeros(NSYST,1);%(1:NSYST) = 0;
            for II=1:NACSYS(ITENT) %260
              ISYST=IACSET(II,ITENT);
              for JJ=1:NACSYS(ITENT) %250
                JSYST=IACSET(JJ,ITENT);
                DDGAM(ISYST)=DDGAM(ISYST)+GINV(II,JJ)*PHI(JSYST); %THEN %Box 16.3 (v), (16.44)
              end
%   250         CONTINUE
              DGAM(ISYST)=DGAM(ISYST)+DDGAM(ISYST); %Box 16.3 (vi), (16.35)
            end
%   260       CONTINUE
% C Compute inverse of incremental plastic deformation gradient
% C... sum up contributions from each active slip system
%             CALL RVZERO(FPILOG,9)
            FPILOG = zeros(3,3);
            for II=1:NACSYS(ITENT) %290
              ISYST=IACSET(II,ITENT);
              for I=1:NDIM %280
                for J=1:NDIM %270
                  FPILOG(I,J)=FPILOG(I,J)- ...
                              DGAM(ISYST)*VECS0(I,ISYST)*VECM0(J,ISYST); %Box 16.3 (vi), (16.67)
                end
              end
            end
%   270           CONTINUE
%   280         CONTINUE
%   290       CONTINUE
% C... use exponential map to update inverse of incremental Fp
%             CALL EXPMAP(   FPINCI     ,NOCONV     ,FPILOG     )
            FPINCI = expm(FPILOG);
%             if(NOCONV) %THEN
% % C... exponential map algorithm failed: Break loop and exit
% %               SUFAIL=1;%.TRUE.
%               error('WE0014')
%             end %IF
% C Update isochoric component of elastic deformation gradient
%             CALL RVZERO(FEISO,4)
            FEISO(1:2,1:2)=FETISO(1:2,1:2)*FPINCI(1:2,1:2); %Box 16.3 (vi), (16.67)
%             FEISO = zeros(4,1);
%             for I=1:NDIM %320
%               for J=1:NDIM %310
%                 for K=1:NDIM %300
%                   FEISO(I,J)=FEISO(I,J)+FETISO(I,K)*FPINCI(K,J);
%                 end
%               end
%             end
%   300           CONTINUE
%   310         CONTINUE
%   320       CONTINUE
% C Update hardening internal variable and yield resolved shear stress
            HRVAR=HRVARN;
            for II=1:NACSYS(ITENT) %330
              ISYST=IACSET(II,ITENT);
              HRVAR=HRVAR+DGAM(ISYST); %Box 16.3 (vi), (16.35)
            end
%   330       CONTINUE
% C Compute yield functions values and check for convergence
% C... elastic push forward of all slip-systems vectors
            VECS=FEISO*VECS0;
            VECM=FEISO*VECM0; %Box 16.3 (vii)
% %             CALL RVZERO(VECS,NDIM*NSYST)
%             VECS = zeros(NDIM,NSYST);
% %             CALL RVZERO(VECM,NDIM*NSYST);
%             VECM = zeros(NDIM,NSYST);
%             for I=1:NDIM %360
%               for J=1:NDIM %350
%                 for ISYST=1:NSYST %340
%                   VECS(I,ISYST)=VECS(I,ISYST)+FEISO(I,J)*VECS0(J,ISYST);
%                   VECM(I,ISYST)=VECM(I,ISYST)+FEISO(I,J)*VECM0(J,ISYST);
%                 end
%               end
%             end
%   340           CONTINUE
%   350         CONTINUE
%   360       CONTINUE
% C... update resolved yield stress
            RYIELD=PLFUN(HRVAR,NHARD,RPROPS2);
% C... Schmid resolved shear stresses for all systems and corresponding
% C    yield function values
      SCHMID=GMODU*diag(VECS'*VECM); %
      PHI=SCHMID-RYIELD;
%             for ISYST=1:NSYST %370
%               SCHMID(ISYST)=GMODU*VECS(:,ISYST)'*VECM(:,ISYST);
%               SCHMID(ISYST)=GMODU* ...
%                             SCAPRD(VECS(1,ISYST),VECM(1,ISYST),NDIM);
%               PHI(ISYST)=SCHMID(ISYST)-RYIELD;
%             end
%   370       CONTINUE
% C... check for convergence
            RESNOR=R0;
            for II=1:NACSYS(ITENT) %380
              ISYST=IACSET(II,ITENT);
              RESNOR=RESNOR+abs(PHI(ISYST));
            end
%   380       CONTINUE
            RESNOR=RESNOR/RYIELD;
            if(RESNOR <= TOL) %THEN
% C... N-R loop converged: check validity of current solution
              CORRECT = 1;
              for ISYST=1:NSYST %390
                if(DGAM(ISYST)<R0|| ...
                   PHI(ISYST)/RYIELD-TOL>R0) %THEN
% C... current solution is not valid: Try another active set or exit if
% C                                   all possible sets have already been
% C                                   tried
%                   GOTO 420
                  CORRECT = 0;
                  break
                end %IF
              end
%   390         CONTINUE
% C... Stress updated converged: Break loop to update necessary variables
% C                              and exit
              if CORRECT
              for II=1:NACSYS(ITENT) %400
                ISYST=IACSET(II,ITENT);
                if(ISYST == 1)
                    S1ACT=1;%.TRUE.
                end
                if(ISYST == 2)
                    S2ACT=1;%.TRUE.
                end
                if(ISYST == 3)
                    S3ACT=1;%.TRUE.
                end
                if(ISYST == 4)
                    S4ACT=1;%.TRUE.
                end
              end
%   400         CONTINUE
%               GOTO 450
              ALLDONE = 1;
              end
              break %%%%%%%%%%% MAY NEED TO BREAK OUT FURTHER
            end %IF
          end
          if ALLDONE
              break
          end
%   410     CONTINUE
        end
%   420   CONTINUE
% C Stress update procedure failed to converge: Break loop and exit
%         SUFAIL=1;%.TRUE.
        if ~ALLDONE
        error('WE0015')
        end
%         GOTO 999
      else
% C Elastic step: Trial state is the actual one
% C ===========================================
        FEISO=FETISO;
%         for I=1:NDIM %440
%           for J=1:NDIM %430
%             FEISO(I,J)=FETISO(I,J);
%           end
%         end
%   430     CONTINUE
%   440   CONTINUE
      end %IF
% C Use neo-Hookean law to update stresses
% C ======================================
%   450 CONTINUE
% C Compute elastic left Cauchy-Green tensor
%       CALL RVZERO(BEISO,9)
      BEISO(1:2,1:2) = FEISO*FEISO'; %Box 16.2 (v), (16.57)
%       BEISO = zeros(9,1);
%       for I=1:NDIM %480
%         for J=1:NDIM %470
%           for K=1:NDIM %460
%             BEISO(I,J)=BEISO(I,J)+FEISO(I,K)*FEISO(J,K);
%           end
%         end
%       end
%   460     CONTINUE
%   470   CONTINUE
%   480 CONTINUE
      BEISO(3,3)=VOLFAC*VOLFAC;
% C Hydrostatic pressure
      P=BULK*log(DETFET); %Box 16.2 (v), (16.59)
% C Deviatoric component of isochoric elastic left Cauchy-Green tensor
      TRACE=BEISO(1,1)+BEISO(2,2)+BEISO(3,3); %Box 16.2 (v), (16.59)
      BEDEV(1)=BEISO(1,1)-R1D3*TRACE;
      BEDEV(2)=BEISO(2,2)-R1D3*TRACE;
      BEDEV(3)=BEISO(1,2);
      BEDEV(4)=BEISO(3,3)-R1D3*TRACE;
% C Update Cauchy stress components
      DETINV=R1/DETFET;
      STRES(1)=(GMODU*BEDEV(1)+P)*DETINV; %Box 16.2 (v), (16.59)
      STRES(2)=(GMODU*BEDEV(2)+P)*DETINV;
      STRES(3)=(GMODU*BEDEV(3))*DETINV;
      STRES(4)=(GMODU*BEDEV(4)+P)*DETINV;
% C Update elastic deformation gradient components
      RSTAVA(1)=FEISO(1,1)/VOLFAC;
      RSTAVA(2)=FEISO(2,1)/VOLFAC;
      RSTAVA(3)=FEISO(1,2)/VOLFAC;
      RSTAVA(4)=FEISO(2,2)/VOLFAC;
% C Store updated hardening variable
      RSTAVA(5)=HRVAR;
%   999 CONTINUE
% C Update some algorithmic variables before exit
% C =============================================
      LALGVA(1)=IFPLAS;
      LALGVA(2)=SUFAIL;
      if(~SUFAIL) %THEN
% C Update active system flags if state update was successful
        LALGVA(3)=S1ACT;
        LALGVA(4)=S2ACT;
        LALGVA(5)=S3ACT;
        LALGVA(6)=S4ACT;
      end %IF
%       RETURN
end
