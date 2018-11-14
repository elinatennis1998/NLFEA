% function ElemK = NL_Elem2_2d(ElemXYZ, nel, ndf, R, S, T, pr, ps, PatchE, Patchv)
%
% Tim Truster
% 01/04/2012
% UIUC

% Master symmetric mixed CST element, 3D
% Verified 01/04/2012 against NL_Elem5_3dCST.m for BC5U4M1.m
% Stress post-processing functions and gives meaningful results, not yet
% verified.
% Verified against files in NL Mixed Elasticity 06/30/2013 using BC5U4M1.m

symmns = 0; % symmetric

switch isw %Task Switch
%%    
    case 1
        
        if ndf > 4
            
            for i = 5:ndf
                lie(i,1) = 0;
            end
            
        end
        
        nh3 = 9;
        
%%
    case 3
        
        ElemK = zeros(nst);
        ElemF = zeros(nst,1);
%         ElemP = zeros(2,1);
%         ElemR = zeros(2,1);
%         Nmat = zeros(3,3*nel);
        Bmat = zeros(10,4*nel);
%         BBmat = zeros(6,3*nel);

        lam = getlam(mateprop);
        
        I1 = [1; 1; 1; 0; 0; 0; 0; 0; 0];
        I2 = eye(3);
        I6 = eye(9);
        P2 = [1 0 0 0 0 0 0 0 0
              0 0 0 1/2 0 0 -1/2 0 0
              0 0 0 0 0 1/2 0 0 1/2
              0 0 0 1/2 0 0 1/2 0 0
              0 1 0 0 0 0 0 0 0
              0 0 0 0 1/2 0 0 -1/2 0
              0 0 0 0 0 1/2 0 0 -1/2
              0 0 0 0 1/2 0 0 1/2 0
              0 0 1 0 0 0 0 0 0];
        P3 = [1     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0
     0     0     0     0     0     0     0     0     1
     0     1     0     1     0     0     0     0     0
     0     0     0     0     0     1     0     1     0
     0     0     1     0     0     0     1     0     0
     0    -1     0     1     0     0     0     0     0
     0     0     0     0     0    -1     0     1     0
     0     0     1     0     0     0    -1     0     0];
        P4 = [1 0 0 0 0 0 0 0 0
              0 0 0 1/2 0 0 1/2 0 0
              0 0 0 0 0 1/2 0 0 -1/2
              0 0 0 1/2 0 0 -1/2 0 0
              0 1 0 0 0 0 0 0 0
              0 0 0 0 1/2 0 0 1/2 0
              0 0 0 0 0 1/2 0 0 1/2
              0 0 0 0 1/2 0 0 -1/2 0
              0 0 1 0 0 0 0 0 0];
        Z2 = zeros(3);
        spvec0 = I1;
        spmat0 = I2;
        cpmat1 = I1*I1';
        cpmat2 = diag([-2,-2,-2,-1,-1,-1,0,0,0]);
        cpmat0 = cpmat1 + cpmat2;
        dpmat1 = [cpmat1; cpmat1; cpmat1; zeros(27,9)];
        dpmat2 =-[6 2 2 0 0 0 0 0 0
                  2 2 0 0 0 0 0 0 0
                  2 0 2 0 0 0 0 0 0
                  0 0 0 1 0 0 0 0 0
                  0 0 0 0 1 0 0 0 0
                  0 0 0 0 0 1 0 0 0
                  zeros(3,9)
                  2 2 0 0 0 0 0 0 0
                  2 6 2 0 0 0 0 0 0
                  0 2 2 0 0 0 0 0 0
                  0 0 0 1 0 0 0 0 0
                  0 0 0 0 1 0 0 0 0
                  0 0 0 0 0 1 0 0 0
                  zeros(3,9)
                  2 0 2 0 0 0 0 0 0
                  0 2 2 0 0 0 0 0 0
                  2 2 6 0 0 0 0 0 0
                  0 0 0 1 0 0 0 0 0
                  0 0 0 0 1 0 0 0 0
                  0 0 0 0 0 1 0 0 0
                  zeros(3,9)
                  0 0 0 1 0 0 0 0 0
                  0 0 0 1 0 0 0 0 0
                  0 0 0 1 0 0 0 0 0
                  1 1 1 0 0 0 0 0 0
                  0 0 0 0 0 0 0 0 0
                  0 0 0 0 0 0 0 0 0
                  zeros(3,9)
                  0 0 0 0 1 0 0 0 0
                  0 0 0 0 1 0 0 0 0
                  0 0 0 0 1 0 0 0 0
                  0 0 0 0 0 0 0 0 0
                  1 1 1 0 0 0 0 0 0
                  0 0 0 0 0 0 0 0 0
                  zeros(3,9)
                  0 0 0 0 0 1 0 0 0
                  0 0 0 0 0 1 0 0 0
                  0 0 0 0 0 1 0 0 0
                  0 0 0 0 0 0 0 0 0
                  0 0 0 0 0 0 0 0 0
                  1 1 1 0 0 0 0 0 0
                  zeros(3,9)];
        dpmat3 = [8 0 0 0 0 0 0 0 0
                  0 0 0 0 0 0 0 0 0
                  0 0 0 0 0 0 0 0 0
                  0 0 0 2 0 0 0 0 0
                  0 0 0 0 0 0 0 0 0
                  0 0 0 0 0 2 0 0 0
                  zeros(3,9)
                  0 0 0 0 0 0 0 0 0
                  0 8 0 0 0 0 0 0 0
                  0 0 0 0 0 0 0 0 0
                  0 0 0 2 0 0 0 0 0
                  0 0 0 0 2 0 0 0 0
                  0 0 0 0 0 0 0 0 0
                  zeros(3,9)
                  0 0 0 0 0 0 0 0 0
                  0 0 0 0 0 0 0 0 0
                  0 0 8 0 0 0 0 0 0
                  0 0 0 0 0 0 0 0 0
                  0 0 0 0 2 0 0 0 0
                  0 0 0 0 0 2 0 0 0
                  zeros(3,9)
                  0 0 0 2 0 0 0 0 0
                  0 0 0 2 0 0 0 0 0
                  0 0 0 0 0 0 0 0 0
                  2 2 0 0 0 0 0 0 0
                  0 0 0 0 0 1 0 0 0
                  0 0 0 0 1 0 0 0 0
                  zeros(3,9)
                  0 0 0 0 0 0 0 0 0
                  0 0 0 0 2 0 0 0 0
                  0 0 0 0 2 0 0 0 0
                  0 0 0 0 0 1 0 0 0
                  0 2 2 0 0 0 0 0 0
                  0 0 0 1 0 0 0 0 0
                  zeros(3,9)
                  0 0 0 0 0 2 0 0 0
                  0 0 0 0 0 0 0 0 0
                  0 0 0 0 0 2 0 0 0
                  0 0 0 0 1 0 0 0 0
                  0 0 0 1 0 0 0 0 0
                  2 0 2 0 0 0 0 0 0
                  zeros(3,9)];
        dpmat0 = dpmat1 + dpmat2 + dpmat3;
        
        % Load Guass Integration Points

        lint = 4;14;36;24;15;11;5; %36 is minimum for full integration of tau, 14 for K and F
        der = 0;
        bf = 1;
        ib = 0;

        %Form tau matrix: compute on first iteration and store; this gives
        %quadratic convergence when tau isn't recalculated (no
        %linearization of tau required).
        if iter <=1 % == 0 % 
        [t11,t12,t13,t21,t22,t23,t31,t32,t33] = TauSCST3(xl,ul,mateprop,nel,nen,nelP,lint);
%         TauList(elem,:) = [t11,t12,t13,t21,t22,t23,t31,t32,t33];
        hr(nhc:nhc+8) = [t11,t12,t13,t21,t22,t23,t31,t32,t33];
        else
%         t11 = TauList(elem,1);
%         t12 = TauList(elem,2);
%         t21 = TauList(elem,4);
%         t22 = TauList(elem,5);
%         t23 = TauList(elem,6);
%         t31 = TauList(elem,7);
%         t32 = TauList(elem,8);
%         t33 = TauList(elem,9);
        t11 = hr(nhc+0);
        t12 = hr(nhc+1);
        t13 = hr(nhc+2);
        t21 = hr(nhc+3);
        t22 = hr(nhc+4);
        t23 = hr(nhc+5);
        t31 = hr(nhc+6);
        t32 = hr(nhc+7);
        t33 = hr(nhc+8);
        end
        
        %Evaluate first derivatives of basis functions at int. point
          [Wgt,ss] =  int3d_t(1,1,ib);
          [shl,shld,shls,be] = shltt(ss,nel,nel,der,bf);
          [Qxy,shgs,Jdet,be,sx] = shgtt(xl+ul(1:3,:),nel,shld,shls,nel,bf,der,be);

        [fi,JxX,F] = kine3d(Qxy,-ul,nel,0); %this is equivalent to ikine2d
        JxX = 1/JxX; %this is equivalent to ikine2d
%         [F,JxX,fi] = ikine2d(Qxy,ul,nel); %this is equivalent to above two lines
        xs = sx;
        Jdet = Jdet/JxX;
        [theta,theta1,theta2,theta3] = ThetaS(JxX,mateprop);
        Pvec = (ul(4,:)*Qxy)';
        Pmat = [Pvec(1) 0 0 Pvec(2) 0       Pvec(3) 0 0 0 
                0 Pvec(2) 0 Pvec(1) Pvec(3) 0       0 0 0
                0 0 Pvec(3) 0       Pvec(2) Pvec(1) 0 0 0];
        [sigmai, cmati] = SigmaCmatNSCST3i(F,JxX,mateprop);
            
        spvec = JxX*theta1*spvec0;
        spmat = JxX*theta1*spmat0;
        cpmat = (theta2*JxX^2 + theta1*JxX)*cpmat1 + theta1*JxX*cpmat2;
        A = JxX*theta1 + 3*JxX^2*theta2 + JxX^3*theta3;
        B = JxX*theta1 + JxX^2*theta2;
        C = JxX*theta1;
        dpmat = (A*dpmat1 + B*dpmat2 + C*dpmat3);
                
        Tsmall = [t11 t12 t13; t21 t22 t23; t31 t32 t33]; %Modified 8/23 for stored tau
            
        % Form BB matrix
        BBmat = [zeros(18,16) 
                 0        0        0        Qxy(1,1) 0        0        0        Qxy(2,1) 0        0        0        Qxy(3,1) 0        0        0        Qxy(4,1)
                 0        0        0        Qxy(1,2) 0        0        0        Qxy(2,2) 0        0        0        Qxy(3,2) 0        0        0        Qxy(4,2)
                 0        0        0        Qxy(1,3) 0        0        0        Qxy(2,3) 0        0        0        Qxy(3,3) 0        0        0        Qxy(4,3)];

        %Integration Loop
        for ll = 1:lint

            %Evaluate first derivatives of basis functions at int. point
            [Wgt,ss] =  int3d_t(ll,lint,ib);
            bb = bubbleCSTT(ss,xs);
            shlp = shltCSTT(ss);

            w = Wgt*Jdet;

            press = ul(4,:)*shlp;
            sigmap = press*spvec;
            cmatp = press*cpmat;
            sigma = sigmai + sigmap;
            cmat = cmati + cmatp;
            
            term32 = spmat*Pvec;
            ElemTR = bb*Tsmall*term32;
            term32M = [term32(1)*I2 term32(2)*I2 term32(3)*I2];
%             term33 = [term32(1) 0 0 term32(2)/2 0           term32(3)/2  term32(2)/2  0           -term32(3)/2
%                       0 term32(2) 0 term32(1)/2 term32(3)/2 0           -term32(1)/2  term32(3)/2  0
%                       0 0 term32(3) 0           term32(2)/2 term32(1)/2  0           -term32(2)/2  term32(1)/2];
            term33 = term32M*P2;
            term34 = Pmat*cpmat;
            BT = [term33+term34 zeros(3,1); zeros(3,10)];
            BBT = [zeros(3,21); zeros(3,18) spmat];
            Svec = [sigma; (theta-press/lam)];
            sig19 = cpmat*Pmat'*ElemTR;
            sigma2 = sigma - sig19;
            Smat = ...
[    sigma2(1),        0,        0,           sigma2(4)/2,                 0,           sigma2(6)/2,           sigma2(4)/2,                 0,          -sigma2(6)/2
         0,    sigma2(2),        0,           sigma2(4)/2,           sigma2(5)/2,                 0,          -sigma2(4)/2,           sigma2(5)/2,                 0
         0,        0,    sigma2(3),                 0,           sigma2(5)/2,           sigma2(6)/2,                 0,          -sigma2(5)/2,           sigma2(6)/2
   sigma2(4)/2,  sigma2(4)/2,        0, sigma2(1)/4 + sigma2(2)/4,           sigma2(6)/4,           sigma2(5)/4, sigma2(2)/4 - sigma2(1)/4,           sigma2(6)/4,          -sigma2(5)/4
         0,  sigma2(5)/2,  sigma2(5)/2,           sigma2(6)/4, sigma2(2)/4 + sigma2(3)/4,           sigma2(4)/4,          -sigma2(6)/4, sigma2(3)/4 - sigma2(2)/4,           sigma2(4)/4
   sigma2(6)/2,        0,  sigma2(6)/2,           sigma2(5)/4,           sigma2(4)/4, sigma2(1)/4 + sigma2(3)/4,           sigma2(5)/4,          -sigma2(4)/4, sigma2(1)/4 - sigma2(3)/4
   sigma2(4)/2, -sigma2(4)/2,        0, sigma2(2)/4 - sigma2(1)/4,          -sigma2(6)/4,           sigma2(5)/4, sigma2(1)/4 + sigma2(2)/4,          -sigma2(6)/4,          -sigma2(5)/4
         0,  sigma2(5)/2, -sigma2(5)/2,           sigma2(6)/4, sigma2(3)/4 - sigma2(2)/4,          -sigma2(4)/4,          -sigma2(6)/4, sigma2(2)/4 + sigma2(3)/4,          -sigma2(4)/4
  -sigma2(6)/2,        0,  sigma2(6)/2,          -sigma2(5)/4,           sigma2(4)/4, sigma2(1)/4 - sigma2(3)/4,          -sigma2(5)/4,          -sigma2(4)/4, sigma2(1)/4 + sigma2(3)/4];
            cvec20 = -dpmat*Pmat'*ElemTR;
            cmat20 = [cvec20(1:9) cvec20(10:18) cvec20(19:27) cvec20(28:36) cvec20(37:45) cvec20(46:54) zeros(9,3)];
%             trp = [ElemTR(1)*Pvec(1) 0 0 ElemTR(2)*Pvec(1)/2  0                   ElemTR(3)*Pvec(1)/2 -ElemTR(2)*Pvec(1)/2  0          -ElemTR(2)*Pvec(1)/2
%                    0 ElemTR(2)*Pvec(2) 0 ElemTR(1)*Pvec(2)/2  ElemTR(1)*Pvec(2)/2
%                    0 0 ElemTR(3)*Pvec(3) ElemTR(1)*Pvec(2)/2  ElemTR(1)*Pvec(2)/2
%                    ElemTR(1)*Pvec(2) ElemTR(2)*Pvec(1) (ElemTR(1)*Pvec(1)+ElemTR(2)*Pvec(2))/2 (ElemTR(1)*Pvec(1)-ElemTR(2)*Pvec(2))/2
%                    zeros(3,9)];
            trp = [ElemTR(1)*I6*Pmat' ElemTR(2)*I6*Pmat' ElemTR(3)*I6*Pmat']*P4;
            term21 = -cpmat*trp - trp'*cpmat;
            ITR2 = [ElemTR(1)*I2; ElemTR(2)*I2; ElemTR(3)*I2];
%             term25 = -cpmat*[ElemTR(1) 0 0
%                              0 ElemTR(2) 0
%                              0 0 ElemTR(3)
%                              ElemTR(2) ElemTR(1) 0
%                              0 ElemTR(3) ElemTR(2)
%                              ElemTR(3) 0 ElemTR(1)
%                              zeros(3,3)];
            term25 = -cpmat*P3*ITR2;
%             term26 = -spmat*[ElemTR(1) 0 0 ElemTR(2)/2 0            ElemTR(3)/2 -ElemTR(2)/2  0            ElemTR(3)/2
%                              0 ElemTR(2) 0 ElemTR(1)/2 ElemTR(3)/2  0            ElemTR(1)/2 -ElemTR(3)/2  0
%                              0 0 ElemTR(3) 0           ElemTR(2)/2  ElemTR(1)/2  0            ElemTR(2)/2  -ElemTR(1)/2];
            term26 = -spmat*ITR2'*P4;
            
            D11 = [Smat+cmat+cmat20+term21 spvec
                   spvec' -1/(lam)];
            D12 = [zeros(9,18) term25+term26'
                   zeros(1,21)];
            D21 = D12';
            Tmat = bb*[Tsmall Z2; Z2 Tsmall];
            Tmat2 = bb*[Z2 Tsmall; Tsmall Z2];
            Tvec = [ElemTR; ElemTR];
            
            % Form B matrix
%             Nmat(:,(ie-1)*3+1:3*ie) = [shlp(ie,1) 0         0
%                                        0         shlp(ie,1) 0
%                                        0         0         shlp(ie,1)];
                
            Bmat = [Qxy(1,1) 0        0        0         Qxy(2,1) 0        0        0         Qxy(3,1) 0        0        0         Qxy(4,1) 0        0        0
                    0        Qxy(1,2) 0        0         0        Qxy(2,2) 0        0         0        Qxy(3,2) 0        0         0        Qxy(4,2) 0        0
                    0        0        Qxy(1,3) 0         0        0        Qxy(2,3) 0         0        0        Qxy(3,3) 0         0        0        Qxy(4,3) 0
                    Qxy(1,2) Qxy(1,1) 0        0         Qxy(2,2) Qxy(2,1) 0        0         Qxy(3,2) Qxy(3,1) 0        0         Qxy(4,2) Qxy(4,1) 0        0
                    0        Qxy(1,3) Qxy(1,2) 0         0        Qxy(2,3) Qxy(2,2) 0         0        Qxy(3,3) Qxy(3,2) 0         0        Qxy(4,3) Qxy(4,2) 0
                    Qxy(1,3) 0        Qxy(1,1) 0         Qxy(2,3) 0        Qxy(2,1) 0         Qxy(3,3) 0        Qxy(3,1) 0         Qxy(4,3) 0        Qxy(4,1) 0
                    Qxy(1,2) -Qxy(1,1) 0        0         Qxy(2,2) -Qxy(2,1) 0        0         Qxy(3,2) -Qxy(3,1) 0        0         Qxy(4,2) -Qxy(4,1) 0        0
                    0        Qxy(1,3) -Qxy(1,2) 0         0        Qxy(2,3) -Qxy(2,2) 0         0        Qxy(3,3) -Qxy(3,2) 0         0        Qxy(4,3) -Qxy(4,2) 0
                    -Qxy(1,3) 0        Qxy(1,1) 0         -Qxy(2,3) 0        Qxy(2,1) 0         -Qxy(3,3) 0        Qxy(3,1) 0         -Qxy(4,3) 0        Qxy(4,1) 0
                    0        0        0        shlp(1,1) 0        0        0        shlp(2,1) 0        0        0        shlp(3,1) 0        0        0        shlp(4,1)];
            
            ElemF = ElemF - w*(Bmat'*(Svec - BT'*Tvec) - BBmat'*BBT'*Tvec);
            
            ElemK = ElemK + w*(Bmat'*D11*Bmat ... + BBmat'*D22*BBmat
                    + BBmat'*D21*Bmat + Bmat'*D12*BBmat ...
                    - Bmat'*BT'*Tmat*BT*Bmat - BBmat'*BBT'*Tmat*BBT*BBmat ...
                    - BBmat'*BBT'*Tmat2*BT*Bmat - Bmat'*BT'*Tmat2*BBT*BBmat);

        end %je
ElemK;
%         if elem == 48
%             [ElemK(1,3) ElemK(3,1)]
%         end
%%
    case 6
        
        ElemF = zeros(nst,1);
%         ElemP = zeros(2,1);
%         ElemR = zeros(2,1);
%         Nmat = zeros(3,3*nel);
        Bmat = zeros(10,4*nel);
%         BBmat = zeros(6,3*nel);

        lam = getlam(mateprop);
        
        I1 = [1; 1; 1; 0; 0; 0; 0; 0; 0];
        I2 = eye(3);
        P2 = [1 0 0 0 0 0 0 0 0
              0 0 0 1/2 0 0 -1/2 0 0
              0 0 0 0 0 1/2 0 0 1/2
              0 0 0 1/2 0 0 1/2 0 0
              0 1 0 0 0 0 0 0 0
              0 0 0 0 1/2 0 0 -1/2 0
              0 0 0 0 0 1/2 0 0 -1/2
              0 0 0 0 1/2 0 0 1/2 0
              0 0 1 0 0 0 0 0 0];
        spvec0 = I1;
        spmat0 = I2;
        
        % Load Guass Integration Points

        lint = 14;36;24;15;11;5;4; %36 is minimum for full integration of tau, 14 for K and F
        der = 0;
        bf = 1;
        ib = 0;

        %Form tau matrix: compute on first iteration and store; this gives
        %quadratic convergence when tau isn't recalculated (no
        %linearization of tau required).
        if iter <=1 % == 0 % 
        [t11,t12,t13,t21,t22,t23,t31,t32,t33] = TauSCST3(xl,ul,mateprop,nel,nen,nelP,lint);
%         TauList(elem,:) = [t11,t12,t13,t21,t22,t23,t31,t32,t33];
        hr(nhc:nhc+8) = [t11,t12,t13,t21,t22,t23,t31,t32,t33];
        else
%         t11 = TauList(elem,1);
%         t12 = TauList(elem,2);
%         t21 = TauList(elem,4);
%         t22 = TauList(elem,5);
%         t23 = TauList(elem,6);
%         t31 = TauList(elem,7);
%         t32 = TauList(elem,8);
%         t33 = TauList(elem,9);
        t11 = hr(nhc+0);
        t12 = hr(nhc+1);
        t13 = hr(nhc+2);
        t21 = hr(nhc+3);
        t22 = hr(nhc+4);
        t23 = hr(nhc+5);
        t31 = hr(nhc+6);
        t32 = hr(nhc+7);
        t33 = hr(nhc+8);
        end
        
        %Evaluate first derivatives of basis functions at int. point
          [Wgt,ss] =  int3d_t(1,1,ib);
          [shl,shld,shls,be] = shltt(ss,nel,nel,der,bf);
          [Qxy,shgs,Jdet,be,sx] = shgtt(xl+ul(1:3,:),nel,shld,shls,nel,bf,der,be);

        [fi,JxX,F] = kine3d(Qxy,-ul,nel,0); %this is equivalent to ikine2d
        JxX = 1/JxX; %this is equivalent to ikine2d
%         [F,JxX,fi] = ikine2d(Qxy,ul,nel); %this is equivalent to above two lines
        xs = sx;
        Jdet = Jdet/JxX;
        [theta,theta1] = ThetaS(JxX,mateprop);
        Pvec = (ul(4,:)*Qxy)';
        Pmat = [Pvec(1) 0 0 Pvec(2) 0       Pvec(3) 0 0 0 
                0 Pvec(2) 0 Pvec(1) Pvec(3) 0       0 0 0
                0 0 Pvec(3) 0       Pvec(2) Pvec(1) 0 0 0];
        sigmai = SigmaCmatNSCST3i(F,JxX,mateprop);
            
        spvec = JxX*theta1*spvec0;
        spmat = JxX*theta1*spmat0;
                
        Tsmall = [t11 t12 t13; t21 t22 t23; t31 t32 t33]; %Modified 8/23 for stored tau
            
        % Form BB matrix
        BBmat = [zeros(18,16) 
                 0        0        0        Qxy(1,1) 0        0        0        Qxy(2,1) 0        0        0        Qxy(3,1) 0        0        0        Qxy(4,1)
                 0        0        0        Qxy(1,2) 0        0        0        Qxy(2,2) 0        0        0        Qxy(3,2) 0        0        0        Qxy(4,2)
                 0        0        0        Qxy(1,3) 0        0        0        Qxy(2,3) 0        0        0        Qxy(3,3) 0        0        0        Qxy(4,3)];

        %Integration Loop
        for ll = 1:lint

            %Evaluate first derivatives of basis functions at int. point
            [Wgt,ss] =  int3d_t(ll,lint,ib);
            bb = bubbleCSTT(ss,xs);
            shlp = shltCSTT(ss);

            w = Wgt*Jdet;

            press = ul(4,:)*shlp;
            sigmap = press*spvec;
            sigma = sigmai + sigmap;
            
            term32 = spmat*Pvec;
            ElemTR = bb*Tsmall*term32;
            term32M = [term32(1)*I2 term32(2)*I2 term32(3)*I2];
%             term33 = [term32(1) 0 0 term32(2)/2 0           term32(3)/2  term32(2)/2  0           -term32(3)/2
%                       0 term32(2) 0 term32(1)/2 term32(3)/2 0           -term32(1)/2  term32(3)/2  0
%                       0 0 term32(3) 0           term32(2)/2 term32(1)/2  0           -term32(2)/2  term32(1)/2];
            term33 = term32M*P2;
            term34 = Pmat*cpmat;
            BT = [term33+term34 zeros(3,1); zeros(3,10)];
            BBT = [zeros(3,21); zeros(3,18) spmat];
            Svec = [sigma; (theta-press/lam)];
            Tvec = [ElemTR; ElemTR];
            
            % Form B matrix
%             Nmat(:,(ie-1)*3+1:3*ie) = [shlp(ie,1) 0         0
%                                        0         shlp(ie,1) 0
%                                        0         0         shlp(ie,1)];
                
            Bmat = [Qxy(1,1) 0        0        0         Qxy(2,1) 0        0        0         Qxy(3,1) 0        0        0         Qxy(4,1) 0        0        0
                    0        Qxy(1,2) 0        0         0        Qxy(2,2) 0        0         0        Qxy(3,2) 0        0         0        Qxy(4,2) 0        0
                    0        0        Qxy(1,3) 0         0        0        Qxy(2,3) 0         0        0        Qxy(3,3) 0         0        0        Qxy(4,3) 0
                    Qxy(1,2) Qxy(1,1) 0        0         Qxy(2,2) Qxy(2,1) 0        0         Qxy(3,2) Qxy(3,1) 0        0         Qxy(4,2) Qxy(4,1) 0        0
                    0        Qxy(1,3) Qxy(1,2) 0         0        Qxy(2,3) Qxy(2,2) 0         0        Qxy(3,3) Qxy(3,2) 0         0        Qxy(4,3) Qxy(4,2) 0
                    Qxy(1,3) 0        Qxy(1,1) 0         Qxy(2,3) 0        Qxy(2,1) 0         Qxy(3,3) 0        Qxy(3,1) 0         Qxy(4,3) 0        Qxy(4,1) 0
                    Qxy(1,2) -Qxy(1,1) 0        0         Qxy(2,2) -Qxy(2,1) 0        0         Qxy(3,2) -Qxy(3,1) 0        0         Qxy(4,2) -Qxy(4,1) 0        0
                    0        Qxy(1,3) -Qxy(1,2) 0         0        Qxy(2,3) -Qxy(2,2) 0         0        Qxy(3,3) -Qxy(3,2) 0         0        Qxy(4,3) -Qxy(4,2) 0
                    -Qxy(1,3) 0        Qxy(1,1) 0         -Qxy(2,3) 0        Qxy(2,1) 0         -Qxy(3,3) 0        Qxy(3,1) 0         -Qxy(4,3) 0        Qxy(4,1) 0
                    0        0        0        shlp(1,1) 0        0        0        shlp(2,1) 0        0        0        shlp(3,1) 0        0        0        shlp(4,1)];
            
            ElemF = ElemF - w*(Bmat'*(Svec - BT'*Tvec) - BBmat'*BBT'*Tvec);

        end %je
%%
    case 11
        
        ElemE = zeros(numEn,1);
        %Set integration number
        lint = 4; IntPoint(nel);
        ib = 0;
        bf = 1;
        der = 1;

        lam = getlam(mateprop);
        
        I1 = [1; 1; 1; 0; 0; 0; 0; 0; 0];
        I2 = eye(3);
        I6 = eye(9);
        
        spvec0 = I1;
        spmat0 = I2;
        
        el2fine = zeros(12,1);
        
        %Form tau matrix: compute on first iteration and store; this gives
        %quadratic convergence when tau isn't recalculated (no
        %linearization of tau required).
%         if iter <=2 % == 0 % 
        [t11,t12,t13,t21,t22,t23,t31,t32,t33] = TauSCST3(xl,ul,mateprop,nel,nen,nelP,lint);
% %         TauList(elem,:) = [t11,t12,t13,t21,t22,t23,t31,t32,t33];
%         hr(nhc:nhc+8) = [t11,t12,t13,t21,t22,t23,t31,t32,t33];
%         else
% %         t11 = TauList(elem,1);
% %         t12 = TauList(elem,2);
% %         t21 = TauList(elem,4);
% %         t22 = TauList(elem,5);
% %         t23 = TauList(elem,6);
% %         t31 = TauList(elem,7);
% %         t32 = TauList(elem,8);
% %         t33 = TauList(elem,9);
%         t11 = hr(nhc+0);
%         t12 = hr(nhc+1);
%         t13 = hr(nhc+2);
%         t21 = hr(nhc+3);
%         t22 = hr(nhc+4);
%         t23 = hr(nhc+5);
%         t31 = hr(nhc+6);
%         t32 = hr(nhc+7);
%         t33 = hr(nhc+8);
%         end
        
        %Evaluate first derivatives of basis functions at int. point
          [Wgt,ss] =  int3d_t(1,1,ib);
          [shl,shld,shls,be] = shltt(ss,nel,nel,der,bf);
          [Qxy,shgs,Jdet,be,sx] = shgtt(xl+ul(1:3,:),nel,shld,shls,nel,bf,der,be);
          sx = inv(sx); % Correct, but only needed when correct derivs of bubble are required

        [fi,JxX,F] = kine3d(Qxy,-ul,nel,0); %this is equivalent to ikine2d
        JxX = 1/JxX; %this is equivalent to ikine2d
%         [F,JxX,fi] = ikine2d(Qxy,ul,nel); %this is equivalent to above two lines
        xs = sx;
        Jdet = Jdet/JxX;
        [theta,theta1,theta2,theta3] = ThetaS(JxX,mateprop);
        Pvec = (ul(4,:)*Qxy)';
                
        Tsmall = [t11 t12 t13; t21 t22 t23; t31 t32 t33]; %Modified 8/23 for stored tau
            
        %Integration Loop
        for ll = 1:lint

            %Evaluate first derivatives of basis functions at int. point
            [Wgt,ss] =  int3d_t(ll,lint,ib);
            [bb,bd] = bubbleCSTT(ss,xs);
            shlp = shltCSTT(ss);

            c1 = Wgt*Jdet;

            xint = xl(1,:)*shlp;
            yint = xl(2,:)*shlp;
            zint = xl(3,:)*shlp;

            %Evaluate residual of equilibrium equation
            rvec = JxX*theta1*Pvec;

            %Evaluate explicit fine scale
            ufine = bb*Tsmall*rvec;
            ufinex = bd(1)*Tsmall*rvec;
            ufiney = bd(2)*Tsmall*rvec;
            ufinez = bd(3)*Tsmall*rvec;

            %Add explicit int. point error to element explicit error
            el2fine = el2fine + c1*[ufine.^2; ufinex.^2; ufiney.^2; ufinez.^2];

        end %je

%         for in= 1:3
            ElemE(1:12) = el2fine;
%             ElemE(in+3) = eprixel(in);
%             ElemE(in+6) = epriyel(in);
%             ElemE(in+9) = el2fine(in);
%             ElemE(in+12) = el2fine(in+3);
%         end
        
        H1up = el2fine(4)+el2fine(5)+el2fine(6)+el2fine(7)+el2fine(8)+el2fine(9)+el2fine(10)+el2fine(11)+el2fine(12);
%         H1u = eprixel(1)+eprixel(2)+epriyel(1)+epriyel(2);
%         Ieffvals(elem,:) = [sqrt(H1up/H1u) H1up H1u];
        Ieffvals(elem,:) = [0 H1up 0];
        
        
        
%%
    case -1 % boundary tractions

        ElemF = zeros(nst,1);
        
        % Reorder nodes for corresponding face of integration
        SurfOrientEtoS
        
        t1 = zeros(3,1);
        t2 = t1;
        t3 = t1;
        
        lint = 13;
        der = 0;
        bf = 1;
        ib = 0;
        
        % Integration Loop
        for je = 1:lint

          [Wgt,ss] =  int3d_t(je,lint,edge);
          [shl,shld,shls,be] = shltt(ss,nel,nel,der,bf);
          [QXY,shgs,Jdet,be,sx] = shgtt(xl,nel,shld,shls,nel,bf,der,be);
          
            %Evaluate tangent and normal vectors
            t1 = sx(:,2);
            [tm1, tu1] = VecNormalize(t1);
            t2 = sx(:,1);
            [tm2, tu2] = VecNormalize(t2);
            t3 = VecCrossProd(t1,t2);
            [tm3, tu3] = VecNormalize(t3);

            %Force components are positive in positive coord. direction
            c1 = Wgt*tm3;
            for o=1:nel
                don = shl(o);
                F = don*traction';

                ElemF(ndf*o-3) = ElemF(ndf*o-3) + F(1)*c1;

                ElemF(ndf*o-2) = ElemF(ndf*o-2) + F(2)*c1;

                ElemF(ndf*o-1) = ElemF(ndf*o-1) + F(3)*c1;

            end %o
            
        end %je
        
        % Reorder nodes back to the orientation of the element
        SurfOrientStoE
        ElemF(1:ndf*nel) = ElemF(ilist2);
%%
    case 22
        
      ElemF;
        
        ElemM = zeros(nst);
        ElemF = zeros(nst,1);

        lam = getlam(mateprop);
        
        I1 = [1; 1; 1; 0; 0; 0; 0; 0; 0];
        spvec0 = I1;
        
        % Load Guass Integration Points

        if nel == 4
            lint = 11; %minimum of 7 for all integrals in deformed state
            der = 0;
        elseif nel == 8
%             lint = 4;
            lint = 8;
            der = 1;
        elseif nel == 10
            lint = 14; %minimum of 13 for all integrals in deformed state
            der = 1;
        else
%             lint = 9;
            lint = 27;
            der = 0;
        end

        %Integration Loop
        for ll = 1:lint

            %Evaluate first derivatives of basis functions at int. point
            if nel == 4 || nel == 10
              [Wgt,ss] =  int3d_t(ll,lint,ib);
              [shl,shld,shls,be] = shltt(ss,nel,nel,der,bf);
              [Qxy,shgs,Jdet,be] = shgtt(xl+ul(1:3,:),nel,shld,shls,nel,bf,der,be);
              shlp = shltt(ss,nelP,nel,der,bf);
            else
              [Wgt,ss] =  intpntb(ll,lint,ib);
              [shl,shld,shls,be] = shlb(ss,nel,nel,der,bf);
              [Qxy,shgs,Jdet,be] = shgb(xl+ul(1:3,:),nel,shld,shls,nel,bf,der,be);
              shlp = shlb(ss,nelP,nel,der,bf);
            end

            if nelS == 4 || nelS == 10
              [shlS,shld,shls] = shltt(ss,nelS,nel,0,0);
            else
              [shlS,shld,shls] = shlb(ss,nelS,nel,0,0);
            end
            
            [fi,JxX,F] = kine3d(Qxy,-ul,nel,0); %this is equivalent to ikine2d
            JxX = 1/JxX; %this is equivalent to ikine2d
    %         [F,JxX,fi] = ikine2d(Qxy,ul,nel); %this is equivalent to above two lines
            Jdet = Jdet/JxX;
            [theta,theta1] = ThetaS(JxX,mateprop);
            sigmai = SigmaCmatNSCST3i(F,JxX,mateprop);

            spvec = JxX*theta1*spvec0;

            press = ul(4,:)*shlp;
            sigmap = press*spvec;
            sigma = sigmai + sigmap;
            
            % Form B matrix
            Nmat = shlS';
            
            if stres <= 6 % stress components
                sigmas = sigma(stres);
            elseif stres >= 8
                if stres <= 10 % principal stresses
                sigma = [sigma(1) sigma(4) sigma(6); sigma(4) sigma(2) sigma(5); sigma(6) sigma(5) sigma(3)];
                psig = eig(sigma);
                sigmas = psig(stres-7);
                else % hydrostatic stress
                sigmas = 1/3*sigma'*I1;
                end
            else % von Mises stress
                trs = sigma'*I1;
                dsig = sigma - 1/3*trs*I1;
                sigmas = sqrt(3/2*(dsig'*dsig));
            end
            
            ElemF = ElemF + w*Nmat'*sigmas;
            
            ElemM = ElemM + w*(Nmat'*Nmat);

        end %je
%%        
    case 25 %Stress Projection2

        ElemS = zeros(nel,npstr+1);
        ElemS2 = zeros(nel,npstr);

        lam = getlam(mateprop);
        
        I1 = [1; 1; 1; 0; 0; 0; 0; 0; 0];
        spvec0 = I1;
        
        thick = 1;
        
        % Load Guass Integration Points

            lint = 1;
            nint = 1;
        
        der = 0;
        bf = 0;
        ib = 0;

        %Stress Loop
        for ll = 1:nint

            %Evaluate first derivatives of basis functions at int. point
              [Wgt,ss] =  int3d_t(ll,nint,ib);
              [shl,shld,shls,be] = shltt(ss,nel,nel,der,bf);
              [Qxy,shgs,Jdet,be] = shgtt(xl+ul(1:3,:),nel,shld,shls,nel,bf,der,be);
            
              shlp = shltt(ss,nelP,nel,0,0);
            
            [fi,JxX,F] = kine3d(Qxy,-ul,nel,0); %this is equivalent to ikine2d
            JxX = 1/JxX; %this is equivalent to ikine2d
    %         [F,JxX,fi] = ikine2d(Qxy,ul,nel); %this is equivalent to above two lines
            Jdet = Jdet/JxX;
            [theta,theta1] = ThetaS(JxX,mateprop);
            sigmai = SigmaCmatNSCST3i(F,JxX,mateprop);

            spvec = JxX*theta1*spvec0;
            press = ul(4,:)*shlp;
            sigmap = press*spvec;
            sigma = (sigmai + sigmap)/JxX;
            
            for stres = 1:npstr
            
            if stres <= 6 % stress components
                sigmas = sigma(stres);
            elseif stres >= 8
                if stres <= 10 % principal stresses
                sigma2 = [sigma(1) sigma(4) sigma(6); sigma(4) sigma(2) sigma(5); sigma(6) sigma(5) sigma(3)];
                psig = eig(sigma2);
                sigmas = psig(stres-7);
                else % hydrostatic stress
                sigmas = 1/3*sigma'*I1;
                end
            else % von Mises stress
                trs = sigma'*I1;
                dsig = sigma - 1/3*trs*I1;
                sigmas = sqrt(3/2*(dsig'*dsig));
            end
            
            ElemS2(ll,stres) = sigmas;
            
            end

        end %je
        
        % interpolate stress at nodes
            plist = [1 0 0 0
                     0 1 0 0
                     0 0 0 1];
        
        for ll = 1:nel
            
            r = plist(1,ll);
            s = plist(2,ll);
            t = plist(3,ll);
            shpS = sshp3d(r,s,t,nint);
            
            for stres = 1:npstr
                
                sigmas = ElemS2(1:nint,stres)'*shpS;
                ElemS(ll,stres) = sigmas;
                
            end
            
        end
        
        %Integration Loop
        Vol = 0;
        for ll = 1:lint

            %Evaluate first derivatives of basis functions at int. point
            if nel == 4 || nel == 10
              [Wgt,ss] =  int3d_t(ll,lint,ib);
              [shl,shld,shls,be] = shltt(ss,nel,nel,der,bf);
              [Qxy,shgs,Jdet,be] = shgtt(xl+ul(1:3,:),nel,shld,shls,nel,bf,der,be);
            else
              [Wgt,ss] =  intpntb(ll,lint,ib);
              [shl,shld,shls,be] = shlb(ss,nel,nel,der,bf);
              [Qxy,shgs,Jdet,be] = shgb(xl+ul(1:3,:),nel,shld,shls,nel,bf,der,be);
            end
            
            [fi,JxX,F] = kine3d(Qxy,-ul,nel,0); %this is equivalent to ikine2d
            JxX = 1/JxX; %this is equivalent to ikine2d
    %         [F,JxX,fi] = ikine2d(Qxy,ul,nel); %this is equivalent to above two lines
            Jdet = Jdet/JxX;

            w = Wgt*Jdet*thick;
            
            Vol = Vol + w;

        end %je

        for i = 1:nel
        ElemS(i,npstr+1) = 1; % use this for simple average Vol; % use this for weighted average 
        end
end