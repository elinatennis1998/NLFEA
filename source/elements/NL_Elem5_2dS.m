% function ElemK = NL_Elem2_2d(ElemXYZ, nel, ndf, R, S, T, pr, ps, PatchE, Patchv)
%
% Tim Truster
% 11/03/2011
% UIUC

% Master symmetric mixed element
% Verified agains NL_Elem5_2d2s.m on 11/16/2011
% Reactions, stresses, and error estimation runs properly 11/16/11
% Submesh fine scales don't look good, though

symmns = 0; % symmetric

switch isw %Task Switch
%%    
    case 1
        
        if ndf > 3
            
            for i = 4:ndf
                lie(i,1) = 0;
            end
            
        end
        
        nh3 = 4;
        
%%
    case 3
        
        ElemK = zeros(nst);
        ElemF = zeros(nst,1);
%         ElemP = zeros(2,1);
%         ElemR = zeros(2,1);
%         Nmat = zeros(3,3*nel);
        Bmat = zeros(4,3*nel);
%         BBmat = zeros(6,3*nel);

        lam = getlam(mateprop);
        
        I1 = [1; 1; 0; 0];
        I2 = eye(2);
        I6 = eye(4);
        P2 = [1 0 0 0
              0 0 1/2 -1/2
              0 0 1/2 1/2
              0 1 0 0];
        P3 = [1 0 0 0
              0 0 0 1
              0 1 1 0
              0 -1 1 0];
        P4 = [1 0 0 0
              0 0 1/2 1/2
              0 0 1/2 -1/2
              0 1 0 0];
        Z2 = zeros(2);
        spvec0 = I1;
        spmat0 = I2;
        cpmat1 = I1*I1';
        cpmat2 = diag([-2,-2,-1,0]);
        cpmat0 = cpmat1 + cpmat2;
        dpmat1 = [cpmat1; cpmat1; zeros(4,4)];
        dpmat2 =-[6 2 0 0
                  2 2 0 0
                  0 0 1 0
                  0 0 0 0
                  2 2 0 0
                  2 6 0 0
                  0 0 1 0
                  0 0 0 0
                  0 0 1 0
                  0 0 1 0
                  1 1 0 0
                  0 0 0 0];
        dpmat3 = [8 0 0 0
                  0 0 0 0
                  0 0 2 0
                  0 0 0 0
                  0 0 0 0
                  0 8 0 0
                  0 0 2 0
                  0 0 0 0
                  0 0 2 0
                  0 0 2 0
                  2 2 0 0
                  0 0 0 0];
        dpmat0 = dpmat1 + dpmat2 + dpmat3;
        
        thick = 1;
        
        % Load Guass Integration Points

        if nel == 3 || nel == 6
            if exist('iprob','var') == 1 && iprob == 6
                lint = 25;
            else
                lint = 13;
            end
        elseif nel == 4
            if exist('iprob','var') == 1 && iprob == 6
                lint = 100;4;
            else
                lint = 4;
            end
        else
            if exist('iprob','var') == 1 && iprob == 6
                lint = 100;
            else
                lint = 9;
            end
        end
        der = 1;
        bf = 1;
        ib = 0;

        %Form tau matrix: compute on first iteration and store; this gives
        %quadratic convergence when tau isn't recalculated (no
        %linearization of tau required).
        if iter <=2 % == 0 % 
        [t11,t12,t21,t22] = TauS(xl,ul,mateprop,nel,nen,nelP,lint);
%         TauList(elem,1:4) = [t11,t12,t21,t22];
        hr(nhc:nhc+3) = [t11,t12,t21,t22];
        else
%         t11 = TauList(elem,1);
%         t12 = TauList(elem,2);
%         t21 = TauList(elem,3);
%         t22 = TauList(elem,4);
        t11 = hr(nhc+0);
        t12 = hr(nhc+1);
        t21 = hr(nhc+2);
        t22 = hr(nhc+3);
        end
                
        Tsmall = [t11 t12; t21 t22]; %Modified 8/23 for stored tau
        
if elem == 25
    t11;
end
        %Integration Loop
        for ll = 1:lint

            %Evaluate first derivatives of basis functions at int. point
            if nel == 3 || nel == 6
              [Wgt,litr,lits] =  intpntt(ll,lint,ib);
              [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
              [Qxy, shgs, Jdet,be] = shgt(xl+ul(1:2,:),nel,shld,shls,nel,bf,der,be);
            else
              [Wgt,litr,lits] =  intpntq(ll,lint,ib);
              [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
              [Qxy, shgs, Jdet,be] = shgq(xl+ul(1:2,:),nel,shld,shls,nel,bf,der,be);
            end

            if nelP == 3 || nelP == 6
              [shlp,shld,shls] = shlt(litr,lits,nelP,nel,0,0);
              [Pxy, shgp] = shgt(xl+ul(1:2,:),nelP,shld,shls,nelP,0,0,be);
            else
              [shlp,shld,shls] = shlq(litr,lits,nelP,nel,0,0);
              [Pxy, shgp] = shgq(xl+ul(1:2,:),nelP,shld,shls,nelP,0,0,be);
            end

            if Jdet < 0
                Jdet
            end


            if nel == 3
            % Form N matrix 
            Nmat = [shl(1)  0        0 shl(2)  0        0 shl(3)  0  0     
                    0         shl(1) 0  0         shl(2) 0 0         shl(3) 0
                    0        0         shlp(1,1) 0        0         shlp(2,1) 0        0         shlp(3,1)];
            % Form B matrix    
            Bmat = [Qxy(1,1)  0        0         Qxy(2,1)  0        0         Qxy(3,1)  0        0
                    0         Qxy(1,2) 0         0         Qxy(2,2) 0         0         Qxy(3,2) 0         
                    Qxy(1,2)  Qxy(1,1) 0         Qxy(2,2)  Qxy(2,1) 0         Qxy(3,2)  Qxy(3,1) 0         
                    Qxy(1,2) -Qxy(1,1) 0         Qxy(2,2) -Qxy(2,1) 0         Qxy(3,2) -Qxy(3,1) 0         
                    0        0         shlp(1,1) 0        0         shlp(2,1) 0        0         shlp(3,1)];
            % Form BB matrix
            BBmat = [zeros(6,9)        
                     0         0        Qxy(1,1) 0         0        Qxy(2,1) 0         0        Qxy(3,1) 
                     0         0        Qxy(1,2) 0         0        Qxy(2,2) 0         0        Qxy(3,2)];
            elseif nel == 4
            % Form N matrix
            Nmat = [shl(1)  0        0 shl(2)  0        0 shl(3)  0        0 shl(4)  0        0
                    0         shl(1) 0 0         shl(2) 0 0         shl(3) 0 0         shl(4) 0
                    0        0         shlp(1,1) 0        0         shlp(2,1) 0        0         shlp(3,1) 0        0         shlp(4,1)]; % Added back for body force BF2U4M0.m
            % Form B matrix
            Bmat = [Qxy(1,1)  0        0         Qxy(2,1)  0        0         Qxy(3,1)  0        0         Qxy(4,1)  0        0
                    0         Qxy(1,2) 0         0         Qxy(2,2) 0         0         Qxy(3,2) 0         0         Qxy(4,2) 0         
                    Qxy(1,2)  Qxy(1,1) 0         Qxy(2,2)  Qxy(2,1) 0         Qxy(3,2)  Qxy(3,1) 0         Qxy(4,2)  Qxy(4,1) 0         
                    Qxy(1,2) -Qxy(1,1) 0         Qxy(2,2) -Qxy(2,1) 0         Qxy(3,2) -Qxy(3,1) 0         Qxy(4,2) -Qxy(4,1) 0         
                    0        0         shlp(1,1) 0        0         shlp(2,1) 0        0         shlp(3,1) 0        0         shlp(4,1)];
            % Form BB matrix
            BBmat = [shgs(1,1) 0         0        shgs(2,1) 0         0        shgs(3,1) 0         0        shgs(4,1) 0         0
                     shgs(1,2) 0         0        shgs(2,2) 0         0        shgs(3,2) 0         0        shgs(4,2) 0         0
                     shgs(1,3) 0         0        shgs(2,3) 0         0        shgs(3,3) 0         0        shgs(4,3) 0         0
                     0         shgs(1,1) 0        0         shgs(2,1) 0        0         shgs(3,1) 0        0         shgs(4,1) 0
                     0         shgs(1,2) 0        0         shgs(2,2) 0        0         shgs(3,2) 0        0         shgs(4,2) 0
                     0         shgs(1,3) 0        0         shgs(2,3) 0        0         shgs(3,3) 0        0         shgs(4,3) 0
                     0         0         Pxy(1,1) 0         0        Pxy(2,1) 0         0        Pxy(3,1) 0         0        Pxy(4,1) 
                     0         0         Pxy(1,2) 0         0        Pxy(2,2) 0         0        Pxy(3,2) 0         0        Pxy(4,2)];
            elseif nel == 6
            % Form N matrix
            Nmat = [shl(1)  0        0 shl(2)  0        0 shl(3)  0        0 shl(4)  0        0 shl(5)  0        0 shl(6)  0        0
                    0         shl(1) 0 0         shl(2) 0 0         shl(3) 0 0         shl(4) 0 0         shl(5) 0 0         shl(6) 0
                    0        0         shlp(1,1) 0        0         shlp(2,1) 0        0         shlp(3,1) 0        0         shlp(4,1) 0        0         shlp(5,1) 0        0         shlp(6,1)];
            % Form B matrix
            Bmat = [Qxy(1,1)  0        0         Qxy(2,1)  0        0         Qxy(3,1)  0        0         Qxy(4,1)  0        0         Qxy(5,1)  0        0         Qxy(6,1)  0        0
                    0         Qxy(1,2) 0         0         Qxy(2,2) 0         0         Qxy(3,2) 0         0         Qxy(4,2) 0         0         Qxy(5,2) 0         0         Qxy(6,2) 0
                    Qxy(1,2)  Qxy(1,1) 0         Qxy(2,2)  Qxy(2,1) 0         Qxy(3,2)  Qxy(3,1) 0         Qxy(4,2)  Qxy(4,1) 0         Qxy(5,2)  Qxy(5,1) 0         Qxy(6,2)  Qxy(6,1) 0
                    Qxy(1,2) -Qxy(1,1) 0         Qxy(2,2) -Qxy(2,1) 0         Qxy(3,2) -Qxy(3,1) 0         Qxy(4,2) -Qxy(4,1) 0         Qxy(5,2) -Qxy(5,1) 0         Qxy(6,2) -Qxy(6,1) 0
                    0        0         shlp(1,1) 0        0         shlp(2,1) 0        0         shlp(3,1) 0        0         shlp(4,1) 0        0         shlp(5,1) 0        0         shlp(6,1)];
            % Form BB matrix
            BBmat = [shgs(1,1) 0         0        shgs(2,1) 0         0        shgs(3,1) 0         0        shgs(4,1) 0         0        shgs(5,1) 0         0        shgs(6,1) 0         0
                     shgs(1,2) 0         0        shgs(2,2) 0         0        shgs(3,2) 0         0        shgs(4,2) 0         0        shgs(5,2) 0         0        shgs(6,2) 0         0
                     shgs(1,3) 0         0        shgs(2,3) 0         0        shgs(3,3) 0         0        shgs(4,3) 0         0        shgs(5,3) 0         0        shgs(6,3) 0         0
                     0         shgs(1,1) 0        0         shgs(2,1) 0        0         shgs(3,1) 0        0         shgs(4,1) 0        0         shgs(5,1) 0        0         shgs(6,1) 0
                     0         shgs(1,2) 0        0         shgs(2,2) 0        0         shgs(3,2) 0        0         shgs(4,2) 0        0         shgs(5,2) 0        0         shgs(6,2) 0
                     0         shgs(1,3) 0        0         shgs(2,3) 0        0         shgs(3,3) 0        0         shgs(4,3) 0        0         shgs(5,3) 0        0         shgs(6,3) 0
                     0         0         Pxy(1,1) 0         0         Pxy(2,1) 0         0         Pxy(3,1) 0         0         Pxy(4,1) 0         0         Pxy(5,1) 0         0         Pxy(6,1)
                     0         0         Pxy(1,2) 0         0         Pxy(2,2) 0         0         Pxy(3,2) 0         0         Pxy(4,2) 0         0         Pxy(5,2) 0         0         Pxy(6,2)];
            elseif nel == 9
            % Form N matrix
            Nmat = [shl(1)  0      0 shl(2)  0      0 shl(3)  0      0 shl(4)  0      0 shl(5)  0      0 shl(6)  0      0 shl(7)  0      0 shl(8)  0      0 shl(9)  0      0     
                    0       shl(1) 0 0       shl(2) 0 0       shl(3) 0 0       shl(4) 0 0       shl(5) 0 0       shl(6) 0 0       shl(7) 0 0       shl(8) 0 0       shl(9) 0
                    0       0      shlp(1,1) 0    0 shlp(2,1) 0    0 shlp(3,1) 0    0 shlp(4,1) 0    0 shlp(5,1) 0    0 shlp(6,1) 0    0 shlp(7,1) 0    0 shlp(8,1) 0    0 shlp(9,1)];
            % Form B matrix
            Bmat = [Qxy(1,1)  0        0         Qxy(2,1)  0        0         Qxy(3,1)  0        0         Qxy(4,1)  0        0         Qxy(5,1)  0        0         Qxy(6,1)  0        0         Qxy(7,1)  0        0         Qxy(8,1)  0        0         Qxy(9,1)  0        0
                    0         Qxy(1,2) 0         0         Qxy(2,2) 0         0         Qxy(3,2) 0         0         Qxy(4,2) 0         0         Qxy(5,2) 0         0         Qxy(6,2) 0         0         Qxy(7,2) 0         0         Qxy(8,2) 0         0         Qxy(9,2) 0
                    Qxy(1,2)  Qxy(1,1) 0         Qxy(2,2)  Qxy(2,1) 0         Qxy(3,2)  Qxy(3,1) 0         Qxy(4,2)  Qxy(4,1) 0         Qxy(5,2)  Qxy(5,1) 0         Qxy(6,2)  Qxy(6,1) 0         Qxy(7,2)  Qxy(7,1) 0         Qxy(8,2)  Qxy(8,1) 0         Qxy(9,2)  Qxy(9,1) 0
                    Qxy(1,2) -Qxy(1,1) 0         Qxy(2,2) -Qxy(2,1) 0         Qxy(3,2) -Qxy(3,1) 0         Qxy(4,2) -Qxy(4,1) 0         Qxy(5,2) -Qxy(5,1) 0         Qxy(6,2) -Qxy(6,1) 0         Qxy(7,2) -Qxy(7,1) 0         Qxy(8,2) -Qxy(8,1) 0         Qxy(9,2) -Qxy(9,1) 0
                    0        0         shlp(1,1) 0        0         shlp(2,1) 0        0         shlp(3,1) 0        0         shlp(4,1) 0        0         shlp(5,1) 0        0         shlp(6,1) 0        0         shlp(7,1) 0        0         shlp(8,1) 0        0         shlp(9,1)];
            % Form BB matrix
            BBmat = [shgs(1,1) 0         0        shgs(2,1) 0         0        shgs(3,1) 0         0        shgs(4,1) 0         0        shgs(5,1) 0         0        shgs(6,1) 0         0        shgs(7,1) 0         0        shgs(8,1) 0         0        shgs(9,1) 0         0
                     shgs(1,2) 0         0        shgs(2,2) 0         0        shgs(3,2) 0         0        shgs(4,2) 0         0        shgs(5,2) 0         0        shgs(6,2) 0         0        shgs(7,2) 0         0        shgs(8,2) 0         0        shgs(9,2) 0         0
                     shgs(1,3) 0         0        shgs(2,3) 0         0        shgs(3,3) 0         0        shgs(4,3) 0         0        shgs(5,3) 0         0        shgs(6,3) 0         0        shgs(7,3) 0         0        shgs(8,3) 0         0        shgs(9,3) 0         0
                     0         shgs(1,1) 0        0         shgs(2,1) 0        0         shgs(3,1) 0        0         shgs(4,1) 0        0         shgs(5,1) 0        0         shgs(6,1) 0        0         shgs(7,1) 0        0         shgs(8,1) 0        0         shgs(9,1) 0
                     0         shgs(1,2) 0        0         shgs(2,2) 0        0         shgs(3,2) 0        0         shgs(4,2) 0        0         shgs(5,2) 0        0         shgs(6,2) 0        0         shgs(7,2) 0        0         shgs(8,2) 0        0         shgs(9,2) 0
                     0         shgs(1,3) 0        0         shgs(2,3) 0        0         shgs(3,3) 0        0         shgs(4,3) 0        0         shgs(5,3) 0        0         shgs(6,3) 0        0         shgs(7,3) 0        0         shgs(8,3) 0        0         shgs(9,3) 0
                     0         0         Pxy(1,1) 0         0         Pxy(2,1) 0         0         Pxy(3,1) 0         0         Pxy(4,1) 0         0         Pxy(5,1) 0         0         Pxy(6,1) 0         0         Pxy(7,1) 0         0         Pxy(8,1) 0         0         Pxy(9,1)
                     0         0         Pxy(1,2) 0         0         Pxy(2,2) 0         0         Pxy(3,2) 0         0         Pxy(4,2) 0         0         Pxy(5,2) 0         0         Pxy(6,2) 0         0         Pxy(7,2) 0         0         Pxy(8,2) 0         0         Pxy(9,2)];
            end
            
            bb = be(3);

            [fi,JxX,F] = kine2d(Qxy,-ul,nel,0); %this is equivalent to ikine2d
            JxX = 1/JxX; %this is equivalent to ikine2d
    %         [F,JxX,fi] = ikine2d(Qxy,ul,nel); %this is equivalent to above two lines
    %         xs = sx;
            Jdet = Jdet/JxX;
            [theta,theta1,theta2,theta3] = ThetaS(JxX,mateprop);
            Pvec = (ul(3,:)*Qxy)';
            Pmat = [Pvec(1) 0 Pvec(2) 0 
                    0 Pvec(2) Pvec(1) 0];
            [sigmai, cmati] = SigmaCmatNSCST2i(F,JxX,mateprop);

            w = Wgt*Jdet*thick;

            spvec = JxX*theta1*spvec0;
            spmat = JxX*theta1*spmat0;
            cpmat = (theta2*JxX^2 + theta1*JxX)*cpmat1 + theta1*JxX*cpmat2;
            A = JxX*theta1 + 3*JxX^2*theta2 + JxX^3*theta3;
            B = JxX*theta1 + JxX^2*theta2;
            C = JxX*theta1;
            dpmat = (A*dpmat1 + B*dpmat2 + C*dpmat3);

            press = ul(3,:)*shlp;
            sigmap = press*spvec;
            cmatp = press*cpmat;
            sigma = sigmai + sigmap;
            cmat = cmati + cmatp;
            
            if exist('iprob','var')
                if iprob == 6
                    mu = 40;
                    X = xl(1,:)*shl;
                    Y = xl(2,:)*shl;
                    fb = [ - (101*mu)/(101*X + 100) - (101*lam*(101*X + 100))/10000
                                                      0; 0];
                else
                    fb = zeros(3,1);
                end
            else
                fb = zeros(3,1);
            end
            fb = fb*lamda; % for proportional loading
            
            term32 = spmat*Pvec;
            ElemTR = bb*Tsmall*(term32 + fb(1:2));% % Added back for body force BF2U4M0.m
            term32M = [term32(1)*I2 term32(2)*I2];
%             term33 = [term32(1) 0 term32(2)/2  term32(2)/2
%                       0 term32(2) term32(1)/2 -term32(1)/2];
            term33 = term32M*P2;
            term34 = Pmat*cpmat;
            BT = [term33+term34 zeros(2,1); zeros(2,5)];
            BBT = [zeros(2,8); zeros(2,6) spmat];
            Svec = [sigma; (theta-press/lam)];
            sig19 = cpmat*Pmat'*ElemTR;
            sigma2 = sigma - sig19;
            Smat = [sigma2(1) 0  sigma2(3)/2 sigma2(3)/2
                    0 sigma2(2)  sigma2(3)/2 -sigma2(3)/2
                    sigma2(3)/2  sigma2(3)/2 (sigma2(2)+sigma2(1))/4 (sigma2(2)-sigma2(1))/4
                    sigma2(3)/2 -sigma2(3)/2 (sigma2(2)-sigma2(1))/4 (sigma2(2)+sigma2(1))/4];
            cvec20 = -dpmat*Pmat'*ElemTR;
            cmat20 = [cvec20(1:4) cvec20(5:8) cvec20(9:12) zeros(4,1)];
%             trp = [ElemTR(1)*Pvec(1) 0 ElemTR(2)*Pvec(1)/2 -ElemTR(2)*Pvec(1)/2
%                    0 ElemTR(2)*Pvec(2) ElemTR(1)*Pvec(2)/2  ElemTR(1)*Pvec(2)/2
%                    ElemTR(1)*Pvec(2) ElemTR(2)*Pvec(1) (ElemTR(1)*Pvec(1)+ElemTR(2)*Pvec(2))/2 (ElemTR(1)*Pvec(1)-ElemTR(2)*Pvec(2))/2
%                    zeros(1,4)];
            trp = [ElemTR(1)*I6*Pmat' ElemTR(2)*I6*Pmat']*P4;
            term21 = -cpmat*trp - trp'*cpmat;
            ITR2 = [ElemTR(1)*I2; ElemTR(2)*I2];
%             term25 = -cpmat*[ElemTR(1) 0
%                              0 ElemTR(2)
%                              ElemTR(2) ElemTR(1)
%                              zeros(1,2)];
            term25 = -cpmat*P3*ITR2;
%             term26 = -spmat*[ElemTR(1) 0 ElemTR(2)/2 -ElemTR(2)/2
%                            0 ElemTR(2) ElemTR(1)/2  ElemTR(1)/2];
            term26 = -spmat*ITR2'*P4;
            
            Tmat = bb*[Tsmall Z2; Z2 Tsmall];
            Tmat2 = bb*[Z2 Tsmall; Tsmall Z2];
            Tvec = [ElemTR; ElemTR];
            
            D11 = [Smat+cmat+cmat20+term21 spvec
                   spvec' -1/(lam)] - BT'*Tmat*BT - BT'*Tmat2*BT;
            D12 = [zeros(4,6) term25+term26'
                   zeros(1,8)] - BT'*Tmat*BBT - BT'*Tmat2*BBT;
            D21 = D12';
            D22 = - BBT'*Tmat*BBT - BBT'*Tmat2*BBT;
            
            ElemF = ElemF - w*(Bmat'*(Svec - BT'*Tvec) - BBmat'*BBT'*Tvec);%  - Nmat'*fb % Added back for BF2U4M0.m
            
            ElemK = ElemK + w*(Bmat'*D11*Bmat + BBmat'*D22*BBmat ...
                    + BBmat'*D21*Bmat + Bmat'*D12*BBmat);

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
        Bmat = zeros(4,3*nel);
%         BBmat = zeros(6,3*nel);

        lam = getlam(mateprop);
        
        I1 = [1; 1; 0; 0];
        I2 = eye(2);
        I6 = eye(4);
        P2 = [1 0 0 0
              0 0 1/2 -1/2
              0 0 1/2 1/2
              0 1 0 0];
        P3 = [1 0 0 0
              0 0 0 1
              0 1 1 0
              0 -1 1 0];
        P4 = [1 0 0 0
              0 0 1/2 1/2
              0 0 1/2 -1/2
              0 1 0 0];
        Z2 = zeros(2);
        spvec0 = I1;
        spmat0 = I2;
        
        thick = 1;
        
        % Load Guass Integration Points

        if nel == 3 || nel == 6
            if exist('iprob','var') == 1 && iprob == 6
                lint = 25;
            else
                lint = 13;
            end
        elseif nel == 4
            if exist('iprob','var') == 1 && iprob == 6
                lint = 100;4;
            else
                lint = 4;
            end
        else
            if exist('iprob','var') == 1 && iprob == 6
                lint = 100;
            else
                lint = 9;
            end
        end
        der = 1;
        bf = 1;
        ib = 0;

        %Form tau matrix: compute on first iteration and store; this gives
        %quadratic convergence when tau isn't recalculated (no
        %linearization of tau required).
        if iter <=2 % == 0 % 
        [t11,t12,t21,t22] = TauS(xl,ul,mateprop,nel,nen,nelP,lint);
%         TauList(elem,1:4) = [t11,t12,t21,t22];
        hr(nhc:nhc+3) = [t11,t12,t21,t22];
        else
%         t11 = TauList(elem,1);
%         t12 = TauList(elem,2);
%         t21 = TauList(elem,3);
%         t22 = TauList(elem,4);
        t11 = hr(nhc+0);
        t12 = hr(nhc+1);
        t21 = hr(nhc+2);
        t22 = hr(nhc+3);
        end
                
        Tsmall = [t11 t12; t21 t22]; %Modified 8/23 for stored tau
        
if elem == 25
    t11;
end
        %Integration Loop
        for ll = 1:lint

            %Evaluate first derivatives of basis functions at int. point
            if nel == 3 || nel == 6
              [Wgt,litr,lits] =  intpntt(ll,lint,ib);
              [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
              [Qxy, shgs, Jdet,be] = shgt(xl+ul(1:2,:),nel,shld,shls,nel,bf,der,be);
            else
              [Wgt,litr,lits] =  intpntq(ll,lint,ib);
              [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
              [Qxy, shgs, Jdet,be] = shgq(xl+ul(1:2,:),nel,shld,shls,nel,bf,der,be);
            end

            if nelP == 3 || nelP == 6
              [shlp,shld,shls] = shlt(litr,lits,nelP,nel,0,0);
              [Pxy, shgp] = shgt(xl+ul(1:2,:),nelP,shld,shls,nelP,0,0,be);
            else
              [shlp,shld,shls] = shlq(litr,lits,nelP,nel,0,0);
              [Pxy, shgp] = shgq(xl+ul(1:2,:),nelP,shld,shls,nelP,0,0,be);
            end

            if Jdet < 0
                Jdet
            end

            if nel == 3            
            % Form B matrix    
            Bmat = [Qxy(1,1)  0        0         Qxy(2,1)  0        0         Qxy(3,1)  0        0
                    0         Qxy(1,2) 0         0         Qxy(2,2) 0         0         Qxy(3,2) 0         
                    Qxy(1,2)  Qxy(1,1) 0         Qxy(2,2)  Qxy(2,1) 0         Qxy(3,2)  Qxy(3,1) 0         
                    Qxy(1,2) -Qxy(1,1) 0         Qxy(2,2) -Qxy(2,1) 0         Qxy(3,2) -Qxy(3,1) 0         
                    0        0         shlp(1,1) 0        0         shlp(2,1) 0        0         shlp(3,1)];
            % Form BB matrix
            BBmat = [zeros(6,9)        
                     0         0        Qxy(1,1) 0         0        Qxy(2,1) 0         0        Qxy(3,1) 
                     0         0        Qxy(1,2) 0         0        Qxy(2,2) 0         0        Qxy(3,2)];
            elseif nel == 4            
            % Form B matrix
            Bmat = [Qxy(1,1)  0        0         Qxy(2,1)  0        0         Qxy(3,1)  0        0         Qxy(4,1)  0        0
                    0         Qxy(1,2) 0         0         Qxy(2,2) 0         0         Qxy(3,2) 0         0         Qxy(4,2) 0         
                    Qxy(1,2)  Qxy(1,1) 0         Qxy(2,2)  Qxy(2,1) 0         Qxy(3,2)  Qxy(3,1) 0         Qxy(4,2)  Qxy(4,1) 0         
                    Qxy(1,2) -Qxy(1,1) 0         Qxy(2,2) -Qxy(2,1) 0         Qxy(3,2) -Qxy(3,1) 0         Qxy(4,2) -Qxy(4,1) 0         
                    0        0         shlp(1,1) 0        0         shlp(2,1) 0        0         shlp(3,1) 0        0         shlp(4,1)];
            % Form BB matrix
            BBmat = [shgs(1,1) 0         0        shgs(2,1) 0         0        shgs(3,1) 0         0        shgs(4,1) 0         0
                     shgs(1,2) 0         0        shgs(2,2) 0         0        shgs(3,2) 0         0        shgs(4,2) 0         0
                     shgs(1,3) 0         0        shgs(2,3) 0         0        shgs(3,3) 0         0        shgs(4,3) 0         0
                     0         shgs(1,1) 0        0         shgs(2,1) 0        0         shgs(3,1) 0        0         shgs(4,1) 0
                     0         shgs(1,2) 0        0         shgs(2,2) 0        0         shgs(3,2) 0        0         shgs(4,2) 0
                     0         shgs(1,3) 0        0         shgs(2,3) 0        0         shgs(3,3) 0        0         shgs(4,3) 0
                     0         0         Pxy(1,1) 0         0        Pxy(2,1) 0         0        Pxy(3,1) 0         0        Pxy(4,1) 
                     0         0         Pxy(1,2) 0         0        Pxy(2,2) 0         0        Pxy(3,2) 0         0        Pxy(4,2)];
            elseif nel == 6
            % Form B matrix
            Bmat = [Qxy(1,1)  0        0         Qxy(2,1)  0        0         Qxy(3,1)  0        0         Qxy(4,1)  0        0         Qxy(5,1)  0        0         Qxy(6,1)  0        0
                    0         Qxy(1,2) 0         0         Qxy(2,2) 0         0         Qxy(3,2) 0         0         Qxy(4,2) 0         0         Qxy(5,2) 0         0         Qxy(6,2) 0
                    Qxy(1,2)  Qxy(1,1) 0         Qxy(2,2)  Qxy(2,1) 0         Qxy(3,2)  Qxy(3,1) 0         Qxy(4,2)  Qxy(4,1) 0         Qxy(5,2)  Qxy(5,1) 0         Qxy(6,2)  Qxy(6,1) 0
                    Qxy(1,2) -Qxy(1,1) 0         Qxy(2,2) -Qxy(2,1) 0         Qxy(3,2) -Qxy(3,1) 0         Qxy(4,2) -Qxy(4,1) 0         Qxy(5,2) -Qxy(5,1) 0         Qxy(6,2) -Qxy(6,1) 0
                    0        0         shlp(1,1) 0        0         shlp(2,1) 0        0         shlp(3,1) 0        0         shlp(4,1) 0        0         shlp(5,1) 0        0         shlp(6,1)];
            % Form BB matrix
            BBmat = [shgs(1,1) 0         0        shgs(2,1) 0         0        shgs(3,1) 0         0        shgs(4,1) 0         0        shgs(5,1) 0         0        shgs(6,1) 0         0
                     shgs(1,2) 0         0        shgs(2,2) 0         0        shgs(3,2) 0         0        shgs(4,2) 0         0        shgs(5,2) 0         0        shgs(6,2) 0         0
                     shgs(1,3) 0         0        shgs(2,3) 0         0        shgs(3,3) 0         0        shgs(4,3) 0         0        shgs(5,3) 0         0        shgs(6,3) 0         0
                     0         shgs(1,1) 0        0         shgs(2,1) 0        0         shgs(3,1) 0        0         shgs(4,1) 0        0         shgs(5,1) 0        0         shgs(6,1) 0
                     0         shgs(1,2) 0        0         shgs(2,2) 0        0         shgs(3,2) 0        0         shgs(4,2) 0        0         shgs(5,2) 0        0         shgs(6,2) 0
                     0         shgs(1,3) 0        0         shgs(2,3) 0        0         shgs(3,3) 0        0         shgs(4,3) 0        0         shgs(5,3) 0        0         shgs(6,3) 0
                     0         0         Pxy(1,1) 0         0         Pxy(2,1) 0         0         Pxy(3,1) 0         0         Pxy(4,1) 0         0         Pxy(5,1) 0         0         Pxy(6,1)
                     0         0         Pxy(1,2) 0         0         Pxy(2,2) 0         0         Pxy(3,2) 0         0         Pxy(4,2) 0         0         Pxy(5,2) 0         0         Pxy(6,2)];
            elseif nel == 9
            % Form B matrix
            Bmat = [Qxy(1,1)  0        0         Qxy(2,1)  0        0         Qxy(3,1)  0        0         Qxy(4,1)  0        0         Qxy(5,1)  0        0         Qxy(6,1)  0        0         Qxy(7,1)  0        0         Qxy(8,1)  0        0         Qxy(9,1)  0        0
                    0         Qxy(1,2) 0         0         Qxy(2,2) 0         0         Qxy(3,2) 0         0         Qxy(4,2) 0         0         Qxy(5,2) 0         0         Qxy(6,2) 0         0         Qxy(7,2) 0         0         Qxy(8,2) 0         0         Qxy(9,2) 0
                    Qxy(1,2)  Qxy(1,1) 0         Qxy(2,2)  Qxy(2,1) 0         Qxy(3,2)  Qxy(3,1) 0         Qxy(4,2)  Qxy(4,1) 0         Qxy(5,2)  Qxy(5,1) 0         Qxy(6,2)  Qxy(6,1) 0         Qxy(7,2)  Qxy(7,1) 0         Qxy(8,2)  Qxy(8,1) 0         Qxy(9,2)  Qxy(9,1) 0
                    Qxy(1,2) -Qxy(1,1) 0         Qxy(2,2) -Qxy(2,1) 0         Qxy(3,2) -Qxy(3,1) 0         Qxy(4,2) -Qxy(4,1) 0         Qxy(5,2) -Qxy(5,1) 0         Qxy(6,2) -Qxy(6,1) 0         Qxy(7,2) -Qxy(7,1) 0         Qxy(8,2) -Qxy(8,1) 0         Qxy(9,2) -Qxy(9,1) 0
                    0        0         shlp(1,1) 0        0         shlp(2,1) 0        0         shlp(3,1) 0        0         shlp(4,1) 0        0         shlp(5,1) 0        0         shlp(6,1) 0        0         shlp(7,1) 0        0         shlp(8,1) 0        0         shlp(9,1)];
            % Form BB matrix
            BBmat = [shgs(1,1) 0         0        shgs(2,1) 0         0        shgs(3,1) 0         0        shgs(4,1) 0         0        shgs(5,1) 0         0        shgs(6,1) 0         0        shgs(7,1) 0         0        shgs(8,1) 0         0        shgs(9,1) 0         0
                     shgs(1,2) 0         0        shgs(2,2) 0         0        shgs(3,2) 0         0        shgs(4,2) 0         0        shgs(5,2) 0         0        shgs(6,2) 0         0        shgs(7,2) 0         0        shgs(8,2) 0         0        shgs(9,2) 0         0
                     shgs(1,3) 0         0        shgs(2,3) 0         0        shgs(3,3) 0         0        shgs(4,3) 0         0        shgs(5,3) 0         0        shgs(6,3) 0         0        shgs(7,3) 0         0        shgs(8,3) 0         0        shgs(9,3) 0         0
                     0         shgs(1,1) 0        0         shgs(2,1) 0        0         shgs(3,1) 0        0         shgs(4,1) 0        0         shgs(5,1) 0        0         shgs(6,1) 0        0         shgs(7,1) 0        0         shgs(8,1) 0        0         shgs(9,1) 0
                     0         shgs(1,2) 0        0         shgs(2,2) 0        0         shgs(3,2) 0        0         shgs(4,2) 0        0         shgs(5,2) 0        0         shgs(6,2) 0        0         shgs(7,2) 0        0         shgs(8,2) 0        0         shgs(9,2) 0
                     0         shgs(1,3) 0        0         shgs(2,3) 0        0         shgs(3,3) 0        0         shgs(4,3) 0        0         shgs(5,3) 0        0         shgs(6,3) 0        0         shgs(7,3) 0        0         shgs(8,3) 0        0         shgs(9,3) 0
                     0         0         Pxy(1,1) 0         0         Pxy(2,1) 0         0         Pxy(3,1) 0         0         Pxy(4,1) 0         0         Pxy(5,1) 0         0         Pxy(6,1) 0         0         Pxy(7,1) 0         0         Pxy(8,1) 0         0         Pxy(9,1)
                     0         0         Pxy(1,2) 0         0         Pxy(2,2) 0         0         Pxy(3,2) 0         0         Pxy(4,2) 0         0         Pxy(5,2) 0         0         Pxy(6,2) 0         0         Pxy(7,2) 0         0         Pxy(8,2) 0         0         Pxy(9,2)];
            end
            
            bb = be(3);

            [fi,JxX,F] = kine2d(Qxy,-ul,nel,0); %this is equivalent to ikine2d
            JxX = 1/JxX; %this is equivalent to ikine2d
    %         [F,JxX,fi] = ikine2d(Qxy,ul,nel); %this is equivalent to above two lines
    %         xs = sx;
            Jdet = Jdet/JxX;
            [theta,theta1,theta2,theta3] = ThetaS(JxX,mateprop);
            Pvec = (ul(3,:)*Qxy)';
            Pmat = [Pvec(1) 0 Pvec(2) 0 
                    0 Pvec(2) Pvec(1) 0];
            sigmai = SigmaCmatNSCST2i(F,JxX,mateprop);

            w = Wgt*Jdet*thick;

            spvec = JxX*theta1*spvec0;
            spmat = JxX*theta1*spmat0;

            press = ul(3,:)*shlp;
            sigmap = press*spvec;
            sigma = sigmai + sigmap;
            
            term32 = spmat*Pvec;
            ElemTR = bb*Tsmall*term32;
            term32M = [term32(1)*I2 term32(2)*I2];
%             term33 = [term32(1) 0 term32(2)/2  term32(2)/2
%                       0 term32(2) term32(1)/2 -term32(1)/2];
            term33 = term32M*P2;
            term34 = Pmat*cpmat;
            BT = [term33+term34 zeros(2,1); zeros(2,5)];
            BBT = [zeros(2,8); zeros(2,6) spmat];
            Svec = [sigma; (theta-press/lam)];
            
            Tvec = [ElemTR; ElemTR];
            
            ElemF = ElemF - w*(Bmat'*(Svec - BT'*Tvec) - BBmat'*BBT'*Tvec);

        end %je
%%
    case -1 % boundary tractions
        
        ElemF = zeros(nst,1);
        
        % Reorder nodes for corresponding face of integration
        SurfOrientEtoS

        % Determine bounds of integration
        
        if nel == 4 || nel == 9
            
        dr = 2;
        ro = -1;
        
        % Upper Limit
        if nodeA == ElemFlag(2)
            eR2 = 1;
        else %
            eR2 = 0;
        end
        % Lower Limit
        if nodeB == ElemFlag(1)
            eR1 = -1;
        else %nodeA == ElemFlag(5)
            eR1 = 0;
        end
        
        elseif nel == 3 || nel == 6
            
        dr = 1;
        ro = 0;
            
        % Upper Limit
        if nodeA == ElemFlag(2)
            eR2 = 1;
        else %nodeA == ElemFlagR(5)
            eR2 = 1/2;
        end
        % Lower Limit
        if nodeB == ElemFlag(1)
            eR1 = 0;
        else %nodeA == ElemFlag(5)
            eR1 = 1/2;
        end
        
        end
        
        % Set jacobian for integration space
        drdr = (eR2 - eR1)/dr;
        
        if exist('iprob','var') == 1 && iprob == 6
            lint = 10; % Use 10 for body force BF2U4M0.m
        else
            lint = 4;
        end
%         % Load Gauss Points for quadrature
%         if enrich == 1
%             [rlist, rWgts, rnum] = GaussPoints(pr+2);
%             slist = -1;
%             sWgts = 1;
%             snum = 1;
%         else
%             [rlist, rWgts, rnum] = GaussPoints(pr+1);
%             slist = -1;
%             sWgts = 1;
%             snum = 1;
%         end

%         lamda = Elemv*ElemE/((1+Elemv)*(1-2*Elemv));
%         mu = ElemE/(2*(1+Elemv));
%         thick = 1;
        ideriv = 0;
        der = 0;
        bf = 0;
        thick = 1;
        
        for je = 1:lint

            if nel == 3 || nel == 6
                [Wgt,r,s] = intpntt(je,lint,1);
            elseif nel == 4 || nel == 9
                [Wgt,r,s] = intpntq(je,lint,1);
            end
            r = drdr*(r-ro)+eR1;

            if nel == 3 || nel == 6
              [shl,shld,shls,be] = shlt(r,s,nel,nel,der,bf);
              [QXY, shgs, Jdet, be, sx] = shgt(xl,nel,shld,shls,nel,bf,der,be);
            else
              [shl,shld,shls,be] = shlq(r,s,nel,nel,der,bf);
              [QXY, shgs, Jdet, be, sx] = shgq(xl,nel,shld,shls,nel,bf,der,be);
            end
                    
            %Evaluate tangent and normal vectors
            t1 = [sx(:,1); 0];
            [tm1, tu1] = VecNormalize(t1);
            t2 = [0; 0; 1];
            tm2 = 1;
            tu2 = t2';
            t3 = VecCrossProd(t1,t2);
            [tm3, tu3] = VecNormalize(t3);
            
            if exist('iprob','var')
                if iprob == 6
                    mu = 40;
                    lam = 40;
                    X = xl(1,:)*shl;
                    Y = xl(2,:)*shl;
                    P = [ (101*X*lam*(101*X + 100))/10000, (101*Y*mu)/(101*X + 100) - (10201*X*Y*lam)/10000, 0
                         (101*Y*mu)/100, mu - (100*mu)/(101*X + 100) + X*((101*lam)/100 + (101*mu)/100),     0
                              0,                                            0, (101*X*lam*(101*X + 100))/10000];
                    Traction = P*tu3';
                else
                    Traction = traction;
                end
            else
                Traction = traction;
            end 
            
            c1 = Wgt*tm3*drdr*thick;
            for o=1:nel
                don = shl(o);
                F = don*Traction';

    %                 Fmtl = F'*t; %Magnitudes of F dot tunit(l=1:3)
    % %                 for l = 1:3
    % %                     for m = 1:3
    %                         Ftl = t*Fmtl'; %Sum of Vectors {F.t(l)}t(m)
    % %                     end
    % %                 end  t*t'*F = eye*F = F

                ElemF(ndf*o-2) = ElemF(ndf*o-2) + F(1)*c1;

                ElemF(ndf*o-1)   = ElemF(ndf*o-1)   + F(2)*c1;

            end %o

        end %ie
        
        % Reorder nodes back to the orientation of the element
        SurfOrientStoE
        ElemF(1:ndf*nel) = ElemF(ilist2);
%%
    case 15 % Body Force
        
        ElemF = zeros(nst,1);

        lam = getlam(mateprop);
        
        thick = 1;
        
        % Load Guass Integration Points

        if nel == 3 || nel == 6
            if exist('iprob','var') == 1 && iprob == 6
                lint = 25;
            else
                lint = 13;
            end
        elseif nel == 4
            if exist('iprob','var') == 1 && iprob == 6
                lint = 100;4;
            else
                lint = 4;
            end
        else
            if exist('iprob','var') == 1 && iprob == 6
                lint = 100;
            else
                lint = 9;
            end
        end
        der = 0;
        bf = 0;
        ib = 0;

        %Integration Loop
        for ll = 1:lint

            %Evaluate first derivatives of basis functions at int. point
            if nel == 3 || nel == 6
              [Wgt,litr,lits] =  intpntt(ll,lint,ib);
              [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
              [Qxy, shgs, Jdet] = shgt(xl,nel,shld,shls,nel,bf,der,be);
            else
              [Wgt,litr,lits] =  intpntq(ll,lint,ib);
              [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
              [Qxy, shgs, Jdet] = shgq(xl,nel,shld,shls,nel,bf,der,be);
            end
            shlp = shl;
                
            if nel == 3
            % Form N matrix 
            Nmat = [shl(1)  0        0 shl(2)  0        0 shl(3)  0  0     
                    0         shl(1) 0  0         shl(2) 0 0         shl(3) 0
                    0        0         shlp(1,1) 0        0         shlp(2,1) 0        0         shlp(3,1)];
            elseif nel == 4
            % Form N matrix
            Nmat = [shl(1)  0        0 shl(2)  0        0 shl(3)  0        0 shl(4)  0        0
                    0         shl(1) 0 0         shl(2) 0 0         shl(3) 0 0         shl(4) 0
                    0        0         shlp(1,1) 0        0         shlp(2,1) 0        0         shlp(3,1) 0        0         shlp(4,1)]; % Added back for body force BF2U4M0.m
            elseif nel == 6
            % Form N matrix
            Nmat = [shl(1)  0        0 shl(2)  0        0 shl(3)  0        0 shl(4)  0        0 shl(5)  0        0 shl(6)  0        0
                    0         shl(1) 0 0         shl(2) 0 0         shl(3) 0 0         shl(4) 0 0         shl(5) 0 0         shl(6) 0
                    0        0         shlp(1,1) 0        0         shlp(2,1) 0        0         shlp(3,1) 0        0         shlp(4,1) 0        0         shlp(5,1) 0        0         shlp(6,1)];
            elseif nel == 9
            % Form N matrix
            Nmat = [shl(1)  0      0 shl(2)  0      0 shl(3)  0      0 shl(4)  0      0 shl(5)  0      0 shl(6)  0      0 shl(7)  0      0 shl(8)  0      0 shl(9)  0      0     
                    0       shl(1) 0 0       shl(2) 0 0       shl(3) 0 0       shl(4) 0 0       shl(5) 0 0       shl(6) 0 0       shl(7) 0 0       shl(8) 0 0       shl(9) 0
                    0       0      shlp(1,1) 0    0 shlp(2,1) 0    0 shlp(3,1) 0    0 shlp(4,1) 0    0 shlp(5,1) 0    0 shlp(6,1) 0    0 shlp(7,1) 0    0 shlp(8,1) 0    0 shlp(9,1)];
            end
                
            c1 = Wgt*Jdet*thick;
            
            if exist('iprob','var')
                if iprob == 6
                    mu = 40;
                    X = xl(1,1:nel)*shl;
                    Y = xl(2,1:nel)*shl;
                    fb = [ - (101*mu)/(101*X + 100) - (101*lam*(101*X + 100))/10000
                                                      0
                                                      0];
                else
                    fb = [bodyf(1:2)'; 0];
                end
            else
                fb = [bodyf(1:2)'; 0];
            end
            
            ElemF = ElemF + c1*Nmat'*fb;

        end %je
    ElemF;
%%
    case 11
        
        ElemE = zeros(numEn,1);
        %Set integration number
        if nel == 3 || nel == 6
            if exist('iprob','var') == 1 && iprob == 6
                lint = 25;
            else
                lint = 13;
            end
        elseif nel == 4
            if exist('iprob','var') == 1 && iprob == 6
                lint = 100;4;
            else
                lint = 16;
            end
        else
            if exist('iprob','var') == 1 && iprob == 6
                lint = 100;
            else
                lint = 16;
            end
        end
        ib = 0;
        bf = 1;
        der = 1;

        lam = getlam(mateprop);
        
        I1 = [1; 1; 0; 0];
        I2 = eye(2);
        I3 = [1 0 0 0; 0 1 0 0; 0 0 2 0];
        I4 = [1 0 0 0 0 1
              0 0 1 0 1 0];
        I5 = diag([1 1 2 1 1 2]);
        I6 = eye(4);
        P1 = [1 0 0 0 0 0
              0 0 0 0 0 1
              0 0 1 1 0 0
              0 0 1 0 0 0
              0 0 0 0 1 0
              0 1 0 0 0 1];
        P2 = [1 0 0 0
              0 0 1/2 -1/2
              0 0 1/2 1/2
              0 1 0 0];
        P3 = [1 0 0 0
              0 0 0 1
              0 1 1 0
              0 -1 1 0];
        P4 = [1 0 0 0
              0 0 1/2 1/2
              0 0 1/2 -1/2
              0 1 0 0];
        Z1 = zeros(1,4);
        Z2 = zeros(2);
        
        spvec0 = I1;
        spmat0 = I2;
        cpmat1 = I1*I1';
        cpmat2 = diag([-2,-2,-1,0]);
        cpmat0 = cpmat1 + cpmat2;
        
        thick = 1;
        fbx = 0;
        fby = 0;
        el2el = zeros(3,1);
        eprixel = zeros(3,1);
        epriyel = zeros(3,1);
        el2fine = zeros(9,1);
        
        %Form tau matrix: compute on first iteration and store; this gives
        %quadratic convergence when tau isn't recalculated (no
        %linearization of tau required).
        if iter == 0
        [t11,t12,t21,t22] = TauS(xl,ul,mateprop,nel,nen,nelP,lint);
%         TauList(elem,1:4) = [t11,t12,t21,t22];
        hr(nhc:nhc+3) = [t11,t12,t21,t22];
        else
%         t11 = TauList(elem,1);
%         t12 = TauList(elem,2);
%         t21 = TauList(elem,3);
%         t22 = TauList(elem,4);
        t11 = hr(nhc+0);
        t12 = hr(nhc+1);
        t21 = hr(nhc+2);
        t22 = hr(nhc+3);
        end
        
        %Evaluate first derivatives of basis functions at int. point
        if nel == 3 || nel == 6
          [shl,shld,shls,be] = shlt(1/3,1/3,nel,nel,der,bf);
          [Qxy, shgs, Jdet, be, sx] = shgt(xl+ul(1:2,:),nel,shld,shls,nen,bf,der,be);
        else
          [shl,shld,shls,be] = shlq(0,0,nel,nel,der,bf);
          [Qxy, shgs, Jdet, be, sx] = shgq(xl+ul(1:2,:),nel,shld,shls,nen,bf,der,be);
        end
        bb = be(3);

        [fi,JxX,F,QXY] = kine2d(Qxy,-ul,nel,1); %this is equivalent to ikine2d
        JxX = 1/JxX; %this is equivalent to ikine2d
%         [F,JxX,fi] = ikine2d(Qxy,ul,nel); %this is equivalent to above two lines
%         xs = sx;
        Jdet = Jdet/JxX;
        [theta,theta1,theta2,theta3] = ThetaS(JxX,mateprop);
                
        Tsmall = [t11 t12; t21 t22]; %Modified 8/23 for stored tau

%         press = ul(3,:)*shlp;
%         sigmai = SigmaCmatNSCST2i(F,JxX,mateprop);
%         sigma = sigmai + sigmap;

        dux = ul*QXY(:,1);
        duy = ul*QXY(:,2);
        u = ul*shl;
        px = JxX*theta1*ul(3,:)*Qxy(:,1);
        py = JxX*theta1*ul(3,:)*Qxy(:,2);
        
        if nel == 3
        % Form N matrix
        Nmat = [shl(1)  0        0 shl(2)  0        0 shl(3)  0        0
                0         shl(1) 0 0         shl(2) 0 0         shl(3) 0
                0        0         shl(1,1) 0        0         shl(2,1) 0        0         shl(3,1)]; % Added back for body force BF2U4M0.m
        % Form B matrix    
        Bmat = [Qxy(1,1)  0        0         Qxy(2,1)  0        0         Qxy(3,1)  0        0
                0         Qxy(1,2) 0         0         Qxy(2,2) 0         0         Qxy(3,2) 0         
                Qxy(1,2)  Qxy(1,1) 0         Qxy(2,2)  Qxy(2,1) 0         Qxy(3,2)  Qxy(3,1) 0         
                Qxy(1,2) -Qxy(1,1) 0         Qxy(2,2) -Qxy(2,1) 0         Qxy(3,2) -Qxy(3,1) 0         
                0        0         shl(1,1) 0        0         shl(2,1) 0        0         shl(3,1)];
        % Form BB matrix
        BBmat = [zeros(6,9)        
                 0         0        Qxy(1,1) 0         0        Qxy(2,1) 0         0        Qxy(3,1) 
                 0         0        Qxy(1,2) 0         0        Qxy(2,2) 0         0        Qxy(3,2)];
        elseif nel == 4
        % Form N matrix
        Nmat = [shl(1)  0        0 shl(2)  0        0 shl(3)  0        0 shl(4)  0        0
                0         shl(1) 0 0         shl(2) 0 0         shl(3) 0 0         shl(4) 0
                0        0         shl(1,1) 0        0         shl(2,1) 0        0         shl(3,1) 0        0         shl(4,1)]; % Added back for body force BF2U4M0.m
        % Form B matrix
        Bmat = [Qxy(1,1)  0        0         Qxy(2,1)  0        0         Qxy(3,1)  0        0         Qxy(4,1)  0        0
                0         Qxy(1,2) 0         0         Qxy(2,2) 0         0         Qxy(3,2) 0         0         Qxy(4,2) 0         
                Qxy(1,2)  Qxy(1,1) 0         Qxy(2,2)  Qxy(2,1) 0         Qxy(3,2)  Qxy(3,1) 0         Qxy(4,2)  Qxy(4,1) 0         
                Qxy(1,2) -Qxy(1,1) 0         Qxy(2,2) -Qxy(2,1) 0         Qxy(3,2) -Qxy(3,1) 0         Qxy(4,2) -Qxy(4,1) 0         
                0        0         shl(1,1) 0        0         shl(2,1) 0        0         shl(3,1) 0        0         shl(4,1)];
        % Form BB matrix
        BBmat = [shgs(1,1) 0         0        shgs(2,1) 0         0        shgs(3,1) 0         0        shgs(4,1) 0         0
                 shgs(1,2) 0         0        shgs(2,2) 0         0        shgs(3,2) 0         0        shgs(4,2) 0         0
                 shgs(1,3) 0         0        shgs(2,3) 0         0        shgs(3,3) 0         0        shgs(4,3) 0         0
                 0         shgs(1,1) 0        0         shgs(2,1) 0        0         shgs(3,1) 0        0         shgs(4,1) 0
                 0         shgs(1,2) 0        0         shgs(2,2) 0        0         shgs(3,2) 0        0         shgs(4,2) 0
                 0         shgs(1,3) 0        0         shgs(2,3) 0        0         shgs(3,3) 0        0         shgs(4,3) 0
                 0         0         Pxy(1,1) 0         0        Pxy(2,1) 0         0        Pxy(3,1) 0         0        Pxy(4,1) 
                 0         0         Pxy(1,2) 0         0        Pxy(2,2) 0         0        Pxy(3,2) 0         0        Pxy(4,2)];
        elseif nel == 6
        % Form B matrix
        Bmat = [Qxy(1,1)  0        0         Qxy(2,1)  0        0         Qxy(3,1)  0        0         Qxy(4,1)  0        0         Qxy(5,1)  0        0         Qxy(6,1)  0        0
                0         Qxy(1,2) 0         0         Qxy(2,2) 0         0         Qxy(3,2) 0         0         Qxy(4,2) 0         0         Qxy(5,2) 0         0         Qxy(6,2) 0
                Qxy(1,2)  Qxy(1,1) 0         Qxy(2,2)  Qxy(2,1) 0         Qxy(3,2)  Qxy(3,1) 0         Qxy(4,2)  Qxy(4,1) 0         Qxy(5,2)  Qxy(5,1) 0         Qxy(6,2)  Qxy(6,1) 0
                Qxy(1,2) -Qxy(1,1) 0         Qxy(2,2) -Qxy(2,1) 0         Qxy(3,2) -Qxy(3,1) 0         Qxy(4,2) -Qxy(4,1) 0         Qxy(5,2) -Qxy(5,1) 0         Qxy(6,2) -Qxy(6,1) 0
                0        0         shl(1,1) 0        0         shl(2,1) 0        0         shl(3,1) 0        0         shl(4,1) 0        0         shl(5,1) 0        0         shl(6,1)];
        % Form BB matrix
        BBmat = [shgs(1,1) 0         0        shgs(2,1) 0         0        shgs(3,1) 0         0        shgs(4,1) 0         0        shgs(5,1) 0         0        shgs(6,1) 0         0
                 shgs(1,2) 0         0        shgs(2,2) 0         0        shgs(3,2) 0         0        shgs(4,2) 0         0        shgs(5,2) 0         0        shgs(6,2) 0         0
                 shgs(1,3) 0         0        shgs(2,3) 0         0        shgs(3,3) 0         0        shgs(4,3) 0         0        shgs(5,3) 0         0        shgs(6,3) 0         0
                 0         shgs(1,1) 0        0         shgs(2,1) 0        0         shgs(3,1) 0        0         shgs(4,1) 0        0         shgs(5,1) 0        0         shgs(6,1) 0
                 0         shgs(1,2) 0        0         shgs(2,2) 0        0         shgs(3,2) 0        0         shgs(4,2) 0        0         shgs(5,2) 0        0         shgs(6,2) 0
                 0         shgs(1,3) 0        0         shgs(2,3) 0        0         shgs(3,3) 0        0         shgs(4,3) 0        0         shgs(5,3) 0        0         shgs(6,3) 0
                 0         0         Pxy(1,1) 0         0         Pxy(2,1) 0         0         Pxy(3,1) 0         0         Pxy(4,1) 0         0         Pxy(5,1) 0         0         Pxy(6,1)
                 0         0         Pxy(1,2) 0         0         Pxy(2,2) 0         0         Pxy(3,2) 0         0         Pxy(4,2) 0         0         Pxy(5,2) 0         0         Pxy(6,2)];
        elseif nel == 9
        % Form B matrix
        Bmat = [Qxy(1,1)  0        0         Qxy(2,1)  0        0         Qxy(3,1)  0        0         Qxy(4,1)  0        0         Qxy(5,1)  0        0         Qxy(6,1)  0        0         Qxy(7,1)  0        0         Qxy(8,1)  0        0         Qxy(9,1)  0        0
                0         Qxy(1,2) 0         0         Qxy(2,2) 0         0         Qxy(3,2) 0         0         Qxy(4,2) 0         0         Qxy(5,2) 0         0         Qxy(6,2) 0         0         Qxy(7,2) 0         0         Qxy(8,2) 0         0         Qxy(9,2) 0
                Qxy(1,2)  Qxy(1,1) 0         Qxy(2,2)  Qxy(2,1) 0         Qxy(3,2)  Qxy(3,1) 0         Qxy(4,2)  Qxy(4,1) 0         Qxy(5,2)  Qxy(5,1) 0         Qxy(6,2)  Qxy(6,1) 0         Qxy(7,2)  Qxy(7,1) 0         Qxy(8,2)  Qxy(8,1) 0         Qxy(9,2)  Qxy(9,1) 0
                Qxy(1,2) -Qxy(1,1) 0         Qxy(2,2) -Qxy(2,1) 0         Qxy(3,2) -Qxy(3,1) 0         Qxy(4,2) -Qxy(4,1) 0         Qxy(5,2) -Qxy(5,1) 0         Qxy(6,2) -Qxy(6,1) 0         Qxy(7,2) -Qxy(7,1) 0         Qxy(8,2) -Qxy(8,1) 0         Qxy(9,2) -Qxy(9,1) 0
                0        0         shl(1,1) 0        0         shl(2,1) 0        0         shl(3,1) 0        0         shl(4,1) 0        0         shl(5,1) 0        0         shl(6,1) 0        0         shl(7,1) 0        0         shl(8,1) 0        0         shl(9,1)];
        % Form BB matrix
        BBmat = [shgs(1,1) 0         0        shgs(2,1) 0         0        shgs(3,1) 0         0        shgs(4,1) 0         0        shgs(5,1) 0         0        shgs(6,1) 0         0        shgs(7,1) 0         0        shgs(8,1) 0         0        shgs(9,1) 0         0
                 shgs(1,2) 0         0        shgs(2,2) 0         0        shgs(3,2) 0         0        shgs(4,2) 0         0        shgs(5,2) 0         0        shgs(6,2) 0         0        shgs(7,2) 0         0        shgs(8,2) 0         0        shgs(9,2) 0         0
                 shgs(1,3) 0         0        shgs(2,3) 0         0        shgs(3,3) 0         0        shgs(4,3) 0         0        shgs(5,3) 0         0        shgs(6,3) 0         0        shgs(7,3) 0         0        shgs(8,3) 0         0        shgs(9,3) 0         0
                 0         shgs(1,1) 0        0         shgs(2,1) 0        0         shgs(3,1) 0        0         shgs(4,1) 0        0         shgs(5,1) 0        0         shgs(6,1) 0        0         shgs(7,1) 0        0         shgs(8,1) 0        0         shgs(9,1) 0
                 0         shgs(1,2) 0        0         shgs(2,2) 0        0         shgs(3,2) 0        0         shgs(4,2) 0        0         shgs(5,2) 0        0         shgs(6,2) 0        0         shgs(7,2) 0        0         shgs(8,2) 0        0         shgs(9,2) 0
                 0         shgs(1,3) 0        0         shgs(2,3) 0        0         shgs(3,3) 0        0         shgs(4,3) 0        0         shgs(5,3) 0        0         shgs(6,3) 0        0         shgs(7,3) 0        0         shgs(8,3) 0        0         shgs(9,3) 0
                 0         0         Pxy(1,1) 0         0         Pxy(2,1) 0         0         Pxy(3,1) 0         0         Pxy(4,1) 0         0         Pxy(5,1) 0         0         Pxy(6,1) 0         0         Pxy(7,1) 0         0         Pxy(8,1) 0         0         Pxy(9,1)
                 0         0         Pxy(1,2) 0         0         Pxy(2,2) 0         0         Pxy(3,2) 0         0         Pxy(4,2) 0         0         Pxy(5,2) 0         0         Pxy(6,2) 0         0         Pxy(7,2) 0         0         Pxy(8,2) 0         0         Pxy(9,2)];
        end

        if exist('iprob','var')
            if iprob == 6
                mu = 40;
                X = xl(1,:)*shl;
                Y = xl(2,:)*shl;
                fbx = - (101*mu)/(101*X + 100) - (101*lam*(101*X + 100))/10000;
            else
            end
        else
        end

        spvec = JxX*theta1*spvec0;
        cpmat = (theta2*JxX^2 + theta1*JxX)*cpmat1 + theta1*JxX*cpmat2;

        press = ul(3,:)*shl;
        sigmap = press*spvec;
        cmatp = press*cpmat;

        [sigmai, cmati] = SigmaCmatNSCST2i(F,JxX,mateprop);
        sigma = sigmai + sigmap;
        cmat = cmati + cmatp;
        ddu = BBmat*reshape(ul,ndf*nel,1);
        ddu = ddu(1:6);

        Fddu = F*ul(1:2,:)*shgs;
        Fddu4 = Fddu*I3;
        diagc = [cmat(1:3,1:3) zeros(3); zeros(3) cmat(1:3,1:3)];
        bigF = [F(1,1) 0 0 F(1,2) 0 0
                0 0 F(2,1) 0 0 F(2,2)
                F(2,1) 0 F(1,1) F(2,2) 0 F(1,2)
                0 0 F(1,1) 0 0 F(1,2)
                0 F(2,1) 0 0 F(2,2) 0
                0 F(1,1) F(2,1) 0 F(1,2) F(2,2)];
        term2 = Fddu4*sigma; %page 6a
        term3 = I4*diagc*bigF*ddu; %page 6a

        %Evaluate residual of equilibrium equation
        rx = px + term2(1) + term3(1) + fbx;
        ry = py + term2(2) + term3(2) + fby;

        %Evaluate explicit fine scale
        bubblevals(elem,1) = (t11*rx+t12*ry);
        bubblevals(elem,2) = (t21*rx+t22*ry);
        JxXm = JxX;
        
        divu = (1) + duy(2);
        
        for ll = 1:lint
                    
            %Evaluate first derivatives of basis functions at int. point
            if nel == 3 || nel == 6
                [Wgt,litr,lits] =  intpntt(ll,lint,ib);
              [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
              [Qxy, shgs, Jdet, be, sx] = shgt(xl+ul(1:2,:),nel,shld,shls,nen,bf,der,be);
            else
                [Wgt,litr,lits] =  intpntq(ll,lint,ib);
              [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
              [Qxy, shgs, Jdet, be, sx] = shgq(xl+ul(1:2,:),nel,shld,shls,nen,bf,der,be);
            end
            bb = be(3);

            [fi,JxX,F,QXY] = kine2d(Qxy,-ul,nel,1); %this is equivalent to ikine2d
            JxX = 1/JxX; %this is equivalent to ikine2d
    %         [F,JxX,fi] = ikine2d(Qxy,ul,nel); %this is equivalent to above two lines
    %         xs = sx;
            Jdet = Jdet/JxX;
            [theta,theta1,theta2,theta3] = ThetaS(JxX,mateprop);

            c1 = Wgt*Jdet*thick;

            xint = xl(1,:)*shl;
            yint = xl(2,:)*shl;
            dux = ul*QXY(:,1);
            duy = ul*QXY(:,2);
            u = ul*shl;
            px = JxX*theta1*ul(3,:)*Qxy(:,1);
            py = JxX*theta1*ul(3,:)*Qxy(:,2);

            if nel == 3
            % Form N matrix
            Nmat = [shl(1)  0        0 shl(2)  0        0 shl(3)  0        0
                    0         shl(1) 0 0         shl(2) 0 0         shl(3) 0
                    0        0         shl(1,1) 0        0         shl(2,1) 0        0         shl(3,1)]; % Added back for body force BF2U4M0.m
            % Form B matrix    
            Bmat = [Qxy(1,1)  0        0         Qxy(2,1)  0        0         Qxy(3,1)  0        0
                    0         Qxy(1,2) 0         0         Qxy(2,2) 0         0         Qxy(3,2) 0         
                    Qxy(1,2)  Qxy(1,1) 0         Qxy(2,2)  Qxy(2,1) 0         Qxy(3,2)  Qxy(3,1) 0         
                    Qxy(1,2) -Qxy(1,1) 0         Qxy(2,2) -Qxy(2,1) 0         Qxy(3,2) -Qxy(3,1) 0         
                    0        0         shl(1,1) 0        0         shl(2,1) 0        0         shl(3,1)];
            % Form BB matrix
            BBmat = [zeros(6,9)        
                     0         0        Qxy(1,1) 0         0        Qxy(2,1) 0         0        Qxy(3,1) 
                     0         0        Qxy(1,2) 0         0        Qxy(2,2) 0         0        Qxy(3,2)];
            elseif nel == 4
            % Form N matrix
            Nmat = [shl(1)  0        0 shl(2)  0        0 shl(3)  0        0 shl(4)  0        0
                    0         shl(1) 0 0         shl(2) 0 0         shl(3) 0 0         shl(4) 0
                    0        0         shl(1,1) 0        0         shl(2,1) 0        0         shl(3,1) 0        0         shl(4,1)]; % Added back for body force BF2U4M0.m
            % Form B matrix
            Bmat = [Qxy(1,1)  0        0         Qxy(2,1)  0        0         Qxy(3,1)  0        0         Qxy(4,1)  0        0
                    0         Qxy(1,2) 0         0         Qxy(2,2) 0         0         Qxy(3,2) 0         0         Qxy(4,2) 0         
                    Qxy(1,2)  Qxy(1,1) 0         Qxy(2,2)  Qxy(2,1) 0         Qxy(3,2)  Qxy(3,1) 0         Qxy(4,2)  Qxy(4,1) 0         
                    Qxy(1,2) -Qxy(1,1) 0         Qxy(2,2) -Qxy(2,1) 0         Qxy(3,2) -Qxy(3,1) 0         Qxy(4,2) -Qxy(4,1) 0         
                    0        0         shl(1,1) 0        0         shl(2,1) 0        0         shl(3,1) 0        0         shl(4,1)];
            % Form BB matrix
            BBmat = [shgs(1,1) 0         0        shgs(2,1) 0         0        shgs(3,1) 0         0        shgs(4,1) 0         0
                     shgs(1,2) 0         0        shgs(2,2) 0         0        shgs(3,2) 0         0        shgs(4,2) 0         0
                     shgs(1,3) 0         0        shgs(2,3) 0         0        shgs(3,3) 0         0        shgs(4,3) 0         0
                     0         shgs(1,1) 0        0         shgs(2,1) 0        0         shgs(3,1) 0        0         shgs(4,1) 0
                     0         shgs(1,2) 0        0         shgs(2,2) 0        0         shgs(3,2) 0        0         shgs(4,2) 0
                     0         shgs(1,3) 0        0         shgs(2,3) 0        0         shgs(3,3) 0        0         shgs(4,3) 0
                     0         0         Pxy(1,1) 0         0        Pxy(2,1) 0         0        Pxy(3,1) 0         0        Pxy(4,1) 
                     0         0         Pxy(1,2) 0         0        Pxy(2,2) 0         0        Pxy(3,2) 0         0        Pxy(4,2)];
            elseif nel == 6
            % Form B matrix
            Bmat = [Qxy(1,1)  0        0         Qxy(2,1)  0        0         Qxy(3,1)  0        0         Qxy(4,1)  0        0         Qxy(5,1)  0        0         Qxy(6,1)  0        0
                    0         Qxy(1,2) 0         0         Qxy(2,2) 0         0         Qxy(3,2) 0         0         Qxy(4,2) 0         0         Qxy(5,2) 0         0         Qxy(6,2) 0
                    Qxy(1,2)  Qxy(1,1) 0         Qxy(2,2)  Qxy(2,1) 0         Qxy(3,2)  Qxy(3,1) 0         Qxy(4,2)  Qxy(4,1) 0         Qxy(5,2)  Qxy(5,1) 0         Qxy(6,2)  Qxy(6,1) 0
                    Qxy(1,2) -Qxy(1,1) 0         Qxy(2,2) -Qxy(2,1) 0         Qxy(3,2) -Qxy(3,1) 0         Qxy(4,2) -Qxy(4,1) 0         Qxy(5,2) -Qxy(5,1) 0         Qxy(6,2) -Qxy(6,1) 0
                    0        0         shl(1,1) 0        0         shl(2,1) 0        0         shl(3,1) 0        0         shl(4,1) 0        0         shl(5,1) 0        0         shl(6,1)];
            % Form BB matrix
            BBmat = [shgs(1,1) 0         0        shgs(2,1) 0         0        shgs(3,1) 0         0        shgs(4,1) 0         0        shgs(5,1) 0         0        shgs(6,1) 0         0
                     shgs(1,2) 0         0        shgs(2,2) 0         0        shgs(3,2) 0         0        shgs(4,2) 0         0        shgs(5,2) 0         0        shgs(6,2) 0         0
                     shgs(1,3) 0         0        shgs(2,3) 0         0        shgs(3,3) 0         0        shgs(4,3) 0         0        shgs(5,3) 0         0        shgs(6,3) 0         0
                     0         shgs(1,1) 0        0         shgs(2,1) 0        0         shgs(3,1) 0        0         shgs(4,1) 0        0         shgs(5,1) 0        0         shgs(6,1) 0
                     0         shgs(1,2) 0        0         shgs(2,2) 0        0         shgs(3,2) 0        0         shgs(4,2) 0        0         shgs(5,2) 0        0         shgs(6,2) 0
                     0         shgs(1,3) 0        0         shgs(2,3) 0        0         shgs(3,3) 0        0         shgs(4,3) 0        0         shgs(5,3) 0        0         shgs(6,3) 0
                     0         0         Pxy(1,1) 0         0         Pxy(2,1) 0         0         Pxy(3,1) 0         0         Pxy(4,1) 0         0         Pxy(5,1) 0         0         Pxy(6,1)
                     0         0         Pxy(1,2) 0         0         Pxy(2,2) 0         0         Pxy(3,2) 0         0         Pxy(4,2) 0         0         Pxy(5,2) 0         0         Pxy(6,2)];
            elseif nel == 9
            % Form B matrix
            Bmat = [Qxy(1,1)  0        0         Qxy(2,1)  0        0         Qxy(3,1)  0        0         Qxy(4,1)  0        0         Qxy(5,1)  0        0         Qxy(6,1)  0        0         Qxy(7,1)  0        0         Qxy(8,1)  0        0         Qxy(9,1)  0        0
                    0         Qxy(1,2) 0         0         Qxy(2,2) 0         0         Qxy(3,2) 0         0         Qxy(4,2) 0         0         Qxy(5,2) 0         0         Qxy(6,2) 0         0         Qxy(7,2) 0         0         Qxy(8,2) 0         0         Qxy(9,2) 0
                    Qxy(1,2)  Qxy(1,1) 0         Qxy(2,2)  Qxy(2,1) 0         Qxy(3,2)  Qxy(3,1) 0         Qxy(4,2)  Qxy(4,1) 0         Qxy(5,2)  Qxy(5,1) 0         Qxy(6,2)  Qxy(6,1) 0         Qxy(7,2)  Qxy(7,1) 0         Qxy(8,2)  Qxy(8,1) 0         Qxy(9,2)  Qxy(9,1) 0
                    Qxy(1,2) -Qxy(1,1) 0         Qxy(2,2) -Qxy(2,1) 0         Qxy(3,2) -Qxy(3,1) 0         Qxy(4,2) -Qxy(4,1) 0         Qxy(5,2) -Qxy(5,1) 0         Qxy(6,2) -Qxy(6,1) 0         Qxy(7,2) -Qxy(7,1) 0         Qxy(8,2) -Qxy(8,1) 0         Qxy(9,2) -Qxy(9,1) 0
                    0        0         shl(1,1) 0        0         shl(2,1) 0        0         shl(3,1) 0        0         shl(4,1) 0        0         shl(5,1) 0        0         shl(6,1) 0        0         shl(7,1) 0        0         shl(8,1) 0        0         shl(9,1)];
            % Form BB matrix
            BBmat = [shgs(1,1) 0         0        shgs(2,1) 0         0        shgs(3,1) 0         0        shgs(4,1) 0         0        shgs(5,1) 0         0        shgs(6,1) 0         0        shgs(7,1) 0         0        shgs(8,1) 0         0        shgs(9,1) 0         0
                     shgs(1,2) 0         0        shgs(2,2) 0         0        shgs(3,2) 0         0        shgs(4,2) 0         0        shgs(5,2) 0         0        shgs(6,2) 0         0        shgs(7,2) 0         0        shgs(8,2) 0         0        shgs(9,2) 0         0
                     shgs(1,3) 0         0        shgs(2,3) 0         0        shgs(3,3) 0         0        shgs(4,3) 0         0        shgs(5,3) 0         0        shgs(6,3) 0         0        shgs(7,3) 0         0        shgs(8,3) 0         0        shgs(9,3) 0         0
                     0         shgs(1,1) 0        0         shgs(2,1) 0        0         shgs(3,1) 0        0         shgs(4,1) 0        0         shgs(5,1) 0        0         shgs(6,1) 0        0         shgs(7,1) 0        0         shgs(8,1) 0        0         shgs(9,1) 0
                     0         shgs(1,2) 0        0         shgs(2,2) 0        0         shgs(3,2) 0        0         shgs(4,2) 0        0         shgs(5,2) 0        0         shgs(6,2) 0        0         shgs(7,2) 0        0         shgs(8,2) 0        0         shgs(9,2) 0
                     0         shgs(1,3) 0        0         shgs(2,3) 0        0         shgs(3,3) 0        0         shgs(4,3) 0        0         shgs(5,3) 0        0         shgs(6,3) 0        0         shgs(7,3) 0        0         shgs(8,3) 0        0         shgs(9,3) 0
                     0         0         Pxy(1,1) 0         0         Pxy(2,1) 0         0         Pxy(3,1) 0         0         Pxy(4,1) 0         0         Pxy(5,1) 0         0         Pxy(6,1) 0         0         Pxy(7,1) 0         0         Pxy(8,1) 0         0         Pxy(9,1)
                     0         0         Pxy(1,2) 0         0         Pxy(2,2) 0         0         Pxy(3,2) 0         0         Pxy(4,2) 0         0         Pxy(5,2) 0         0         Pxy(6,2) 0         0         Pxy(7,2) 0         0         Pxy(8,2) 0         0         Pxy(9,2)];
            end

            if exist('iprob','var')
                if iprob == 6
                    mu = 40;
                    X = xl(1,:)*shl;
                    Y = xl(2,:)*shl;
                    fbx = - (101*mu)/(101*X + 100) - (101*lam*(101*X + 100))/10000;
                else
                end
            else
            end

            spvec = JxX*theta1*spvec0;
            cpmat = (theta2*JxX^2 + theta1*JxX)*cpmat1 + theta1*JxX*cpmat2;

            press = ul(3,:)*shl;
            sigmap = press*spvec;
            cmatp = press*cpmat;

            [sigmai, cmati] = SigmaCmatNSCST2i(F,JxX,mateprop);
            sigma = sigmai + sigmap;
            cmat = cmati + cmatp;
            ddu = BBmat*reshape(ul,ndf*nel,1);
            ddu = ddu(1:6);

            Fddu = F*ul(1:2,:)*shgs;
            Fddu4 = Fddu*I3;
            diagc = [cmat(1:3,1:3) zeros(3); zeros(3) cmat(1:3,1:3)];
            bigF = [F(1,1) 0 0 F(1,2) 0 0
                    0 0 F(2,1) 0 0 F(2,2)
                    F(2,1) 0 F(1,1) F(2,2) 0 F(1,2)
                    0 0 F(1,1) 0 0 F(1,2)
                    0 F(2,1) 0 0 F(2,2) 0
                    0 F(1,1) F(2,1) 0 F(1,2) F(2,2)];
            term2 = Fddu4*sigma; %page 6a
            term3 = I4*diagc*bigF*ddu; %page 6a

            %Evaluate residual of equilibrium equation
            rx = px + term2(1) + term3(1) + fbx;
            ry = py + term2(2) + term3(2) + fby;

            %Evaluate explicit fine scale
            ufinex = (t11*rx+t12*ry)*bb;
            ufiney = (t21*rx+t22*ry)*bb;
            ufinex_x = (t11*rx+t12*ry)*be(1);
            ufiney_x = (t21*rx+t22*ry)*be(1);
            ufinex_y = (t11*rx+t12*ry)*be(2);
            ufiney_y = (t21*rx+t22*ry)*be(2);

                %Compute value of exact fields at int. point
                if iprob == 6
                    [ue,duex,duey] = uexact_BF(xint,yint,lam);
                else
                    ue = zeros(3,1);
                    duex = zeros(3,1);
                    duey = zeros(3,1);
                end

                %Add standard int. point error to element standard error
                for in = 1:3
                    un   = c1 * ( (u(in)-ue(in))^2 );
                    upnx   = c1 * ( (dux(in)-duex(in))^2 );
                    upny   = c1 * ( (duy(in)-duey(in))^2 );
                    el2el(in)   = el2el(in)   + un;
                    eprixel(in) = eprixel(in) + upnx;
                    epriyel(in) = epriyel(in) + upny;
                end

            %Add explicit int. point error to element explicit error
            el2fine(1) = el2fine(1) + c1*ufinex^2;
            el2fine(2) = el2fine(2) + c1*ufiney^2;
            el2fine(3) = el2fine(3) + c1*ufinex_x^2;
            el2fine(4) = el2fine(4) + c1*ufiney_x^2;
            el2fine(5) = el2fine(5) + c1*ufinex_y^2;
            el2fine(6) = el2fine(6) + c1*ufiney_y^2;

        end %je

        for in= 1:3
            ElemE(in) = el2el(in);
            ElemE(in+3) = eprixel(in);
            ElemE(in+6) = epriyel(in);
            ElemE(in+9) = el2fine(in);
            ElemE(in+12) = el2fine(in+3);
        end
        
        H1up = el2fine(3)+el2fine(4)+el2fine(5)+el2fine(6);
%         H1u = eprixel(1)+eprixel(2)+epriyel(1)+epriyel(2);
%         Ieffvals(elem,:) = [sqrt(H1up/H1u) H1up H1u];
        Ieffvals(elem,:) = [0 H1up 0 JxXm divu];
%%
    case 9 %Global error
        
        ElemF = zeros(nst,1);
        xc = zeros(ndm,nen);
        uc = zeros(ndfs,nen);
%         Nmat = zeros(3,3*nel);
        Bmat = zeros(4,3*nel);
%         BBmat = zeros(6,3*nel);

        lam = getlam(mateprop);
        
        I1 = [1; 1; 0; 0];
        I2 = eye(2);
        I3 = [1 0 0 0; 0 1 0 0; 0 0 2 0];
        I4 = [1 0 0 0 0 1
              0 0 1 0 1 0];
        I5 = diag([1 1 2 1 1 2]);
        I6 = eye(4);
        P1 = [1 0 0 0 0 0
              0 0 0 0 0 1
              0 0 1 1 0 0
              0 0 1 0 0 0
              0 0 0 0 1 0
              0 1 0 0 0 1];
        P2 = [1 0 0 0
              0 0 1/2 -1/2
              0 0 1/2 1/2
              0 1 0 0];
        P3 = [1 0 0 0
              0 0 0 1
              0 1 1 0
              0 -1 1 0];
        P4 = [1 0 0 0
              0 0 1/2 1/2
              0 0 1/2 -1/2
              0 1 0 0];
        
        thick = 1;
        
        % Load Guass Integration Points

        if nel == 3 || nel == 6
            if exist('iprob','var') == 1 && iprob == 6
                lint = 25;
            else
                lint = 13;
            end
        elseif nel == 4
            if exist('iprob','var') == 1 && iprob == 6
                lint = 100;4;
            else
                lint = 4;
            end
        else
            if exist('iprob','var') == 1 && iprob == 6
                lint = 100;
            else
                lint = 9;
            end
        end
        der = 1;
        bf = 1;
        ib = 0;

%         %Form tau matrix: compute on first iteration and store; this gives
%         %quadratic convergence when tau isn't recalculated (no
%         %linearization of tau required).
%         if iter <=2 % == 0 % 
%         [t11,t12,t21,t22] = TauSCST(xl,ul,mateprop,nel,nen,nelP,lint);
%         TauList(elem,1:4) = [t11,t12,t21,t22];
%         hr(nhc:nhc+3) = [t11,t12,t21,t22];
%         else
%         t11 = TauList(elem,1);
%         t12 = TauList(elem,2);
%         t21 = TauList(elem,3);
%         t22 = TauList(elem,4);
        t11 = hr(nhc+0);
        t12 = hr(nhc+1);
        t21 = hr(nhc+2);
        t22 = hr(nhc+3);
%         end
        
        %Evaluate first derivatives of basis functions at int. point
          [Wgt,litr,lits] =  intpntt(1,1,ib);
          [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
          [Qxy,shgs,Jdet,be,sx] = shgt(xl+ul(1:2,:),nel,shld,shls,nel,bf,der,be);
%           sx = inv(sx); Correct, but only needed when correct derivs of bubble are required

        [fi,JxX,F] = kine2d(Qxy,-ul,nel,0); %this is equivalent to ikine2d
        JxX = 1/JxX; %this is equivalent to ikine2d
%         [F,JxX,fi] = ikine2d(Qxy,ul,nel); %this is equivalent to above two lines
%         xs = sx;
        Jdet = Jdet/JxX;
        [theta,theta1,theta2,theta3] = ThetaS(JxX,mateprop);
        Pvec = (ul(3,:)*Qxy)';
        Pmat = [Pvec(1) 0 Pvec(2) 0 
                0 Pvec(2) Pvec(1) 0];
        [sigmai, cmati] = SigmaCmatNSCST2i(F,JxX,mateprop);
            
        spvec = JxX*theta1*spvec0;
        spmat = JxX*theta1*spmat0;
        cpmat = (theta2*JxX^2 + theta1*JxX)*cpmat1 + theta1*JxX*cpmat2;
                
        Tsmall = [t11 t12; t21 t22]; %Modified 8/23 for stored tau
            
        % Form BB matrix
        BBmat = [zeros(6,9)        
                 0         0        Qxy(1,1) 0         0        Qxy(2,1) 0         0        Qxy(3,1) 
                 0         0        Qxy(1,2) 0         0        Qxy(2,2) 0         0        Qxy(3,2)];

        %Integration Loop
        for ll = 1:lint

            %Evaluate first derivatives of basis functions at int. point
            if nel == 3 || nel == 6
              [Wgt,litr,lits] =  intpntt(ll,lint,ib);
              [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
              [Qxy, shgs, Jdet,be] = shgt(xl+ul(1:2,:),nel,shld,shls,nel,bf,der,be);
            else
              [Wgt,litr,lits] =  intpntq(ll,lint,ib);
              [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
              [Qxy, shgs, Jdet,be] = shgq(xl+ul(1:2,:),nel,shld,shls,nel,bf,der,be);
            end

            if nelP == 3 || nelP == 6
              [shlp,shld,shls] = shlt(litr,lits,nelP,nel,0,0);
              [Pxy, shgp] = shgt(xl+ul(1:2,:),nelP,shld,shls,nelP,0,0,be);
            else
              [shlp,shld,shls] = shlq(litr,lits,nelP,nel,0,0);
              [Pxy, shgp] = shgq(xl+ul(1:2,:),nelP,shld,shls,nelP,0,0,be);
            end

            if Jdet < 0
                Jdet
            end

            if nel == 3            
            % Form B matrix    
            Bmat = [Qxy(1,1)  0        0         Qxy(2,1)  0        0         Qxy(3,1)  0        0
                    0         Qxy(1,2) 0         0         Qxy(2,2) 0         0         Qxy(3,2) 0         
                    Qxy(1,2)  Qxy(1,1) 0         Qxy(2,2)  Qxy(2,1) 0         Qxy(3,2)  Qxy(3,1) 0         
                    Qxy(1,2) -Qxy(1,1) 0         Qxy(2,2) -Qxy(2,1) 0         Qxy(3,2) -Qxy(3,1) 0         
                    0        0         shlp(1,1) 0        0         shlp(2,1) 0        0         shlp(3,1)];
            % Form BB matrix
            BBmat = [zeros(6,9)        
                     0         0        Qxy(1,1) 0         0        Qxy(2,1) 0         0        Qxy(3,1) 
                     0         0        Qxy(1,2) 0         0        Qxy(2,2) 0         0        Qxy(3,2)];
            elseif nel == 4
            % Form N matrix
            Nmat = [shl(1)  0        0 shl(2)  0        0 shl(3)  0        0 shl(4)  0        0
                    0         shl(1) 0 0         shl(2) 0 0         shl(3) 0 0         shl(4) 0
                    0        0         shlp(1,1) 0        0         shlp(2,1) 0        0         shlp(3,1) 0        0         shlp(4,1)]; % Added back for body force BF2U4M0.m
            % Form B matrix
            Bmat = [Qxy(1,1)  0        0         Qxy(2,1)  0        0         Qxy(3,1)  0        0         Qxy(4,1)  0        0
                    0         Qxy(1,2) 0         0         Qxy(2,2) 0         0         Qxy(3,2) 0         0         Qxy(4,2) 0         
                    Qxy(1,2)  Qxy(1,1) 0         Qxy(2,2)  Qxy(2,1) 0         Qxy(3,2)  Qxy(3,1) 0         Qxy(4,2)  Qxy(4,1) 0         
                    Qxy(1,2) -Qxy(1,1) 0         Qxy(2,2) -Qxy(2,1) 0         Qxy(3,2) -Qxy(3,1) 0         Qxy(4,2) -Qxy(4,1) 0         
                    0        0         shlp(1,1) 0        0         shlp(2,1) 0        0         shlp(3,1) 0        0         shlp(4,1)];
            % Form BB matrix
            BBmat = [shgs(1,1) 0         0        shgs(2,1) 0         0        shgs(3,1) 0         0        shgs(4,1) 0         0
                     shgs(1,2) 0         0        shgs(2,2) 0         0        shgs(3,2) 0         0        shgs(4,2) 0         0
                     shgs(1,3) 0         0        shgs(2,3) 0         0        shgs(3,3) 0         0        shgs(4,3) 0         0
                     0         shgs(1,1) 0        0         shgs(2,1) 0        0         shgs(3,1) 0        0         shgs(4,1) 0
                     0         shgs(1,2) 0        0         shgs(2,2) 0        0         shgs(3,2) 0        0         shgs(4,2) 0
                     0         shgs(1,3) 0        0         shgs(2,3) 0        0         shgs(3,3) 0        0         shgs(4,3) 0
                     0         0         Pxy(1,1) 0         0        Pxy(2,1) 0         0        Pxy(3,1) 0         0        Pxy(4,1) 
                     0         0         Pxy(1,2) 0         0        Pxy(2,2) 0         0        Pxy(3,2) 0         0        Pxy(4,2)];
            elseif nel == 6
            % Form B matrix
            Bmat = [Qxy(1,1)  0        0         Qxy(2,1)  0        0         Qxy(3,1)  0        0         Qxy(4,1)  0        0         Qxy(5,1)  0        0         Qxy(6,1)  0        0
                    0         Qxy(1,2) 0         0         Qxy(2,2) 0         0         Qxy(3,2) 0         0         Qxy(4,2) 0         0         Qxy(5,2) 0         0         Qxy(6,2) 0
                    Qxy(1,2)  Qxy(1,1) 0         Qxy(2,2)  Qxy(2,1) 0         Qxy(3,2)  Qxy(3,1) 0         Qxy(4,2)  Qxy(4,1) 0         Qxy(5,2)  Qxy(5,1) 0         Qxy(6,2)  Qxy(6,1) 0
                    Qxy(1,2) -Qxy(1,1) 0         Qxy(2,2) -Qxy(2,1) 0         Qxy(3,2) -Qxy(3,1) 0         Qxy(4,2) -Qxy(4,1) 0         Qxy(5,2) -Qxy(5,1) 0         Qxy(6,2) -Qxy(6,1) 0
                    0        0         shlp(1,1) 0        0         shlp(2,1) 0        0         shlp(3,1) 0        0         shlp(4,1) 0        0         shlp(5,1) 0        0         shlp(6,1)];
            % Form BB matrix
            BBmat = [shgs(1,1) 0         0        shgs(2,1) 0         0        shgs(3,1) 0         0        shgs(4,1) 0         0        shgs(5,1) 0         0        shgs(6,1) 0         0
                     shgs(1,2) 0         0        shgs(2,2) 0         0        shgs(3,2) 0         0        shgs(4,2) 0         0        shgs(5,2) 0         0        shgs(6,2) 0         0
                     shgs(1,3) 0         0        shgs(2,3) 0         0        shgs(3,3) 0         0        shgs(4,3) 0         0        shgs(5,3) 0         0        shgs(6,3) 0         0
                     0         shgs(1,1) 0        0         shgs(2,1) 0        0         shgs(3,1) 0        0         shgs(4,1) 0        0         shgs(5,1) 0        0         shgs(6,1) 0
                     0         shgs(1,2) 0        0         shgs(2,2) 0        0         shgs(3,2) 0        0         shgs(4,2) 0        0         shgs(5,2) 0        0         shgs(6,2) 0
                     0         shgs(1,3) 0        0         shgs(2,3) 0        0         shgs(3,3) 0        0         shgs(4,3) 0        0         shgs(5,3) 0        0         shgs(6,3) 0
                     0         0         Pxy(1,1) 0         0         Pxy(2,1) 0         0         Pxy(3,1) 0         0         Pxy(4,1) 0         0         Pxy(5,1) 0         0         Pxy(6,1)
                     0         0         Pxy(1,2) 0         0         Pxy(2,2) 0         0         Pxy(3,2) 0         0         Pxy(4,2) 0         0         Pxy(5,2) 0         0         Pxy(6,2)];
            elseif nel == 9
            % Form B matrix
            Bmat = [Qxy(1,1)  0        0         Qxy(2,1)  0        0         Qxy(3,1)  0        0         Qxy(4,1)  0        0         Qxy(5,1)  0        0         Qxy(6,1)  0        0         Qxy(7,1)  0        0         Qxy(8,1)  0        0         Qxy(9,1)  0        0
                    0         Qxy(1,2) 0         0         Qxy(2,2) 0         0         Qxy(3,2) 0         0         Qxy(4,2) 0         0         Qxy(5,2) 0         0         Qxy(6,2) 0         0         Qxy(7,2) 0         0         Qxy(8,2) 0         0         Qxy(9,2) 0
                    Qxy(1,2)  Qxy(1,1) 0         Qxy(2,2)  Qxy(2,1) 0         Qxy(3,2)  Qxy(3,1) 0         Qxy(4,2)  Qxy(4,1) 0         Qxy(5,2)  Qxy(5,1) 0         Qxy(6,2)  Qxy(6,1) 0         Qxy(7,2)  Qxy(7,1) 0         Qxy(8,2)  Qxy(8,1) 0         Qxy(9,2)  Qxy(9,1) 0
                    Qxy(1,2) -Qxy(1,1) 0         Qxy(2,2) -Qxy(2,1) 0         Qxy(3,2) -Qxy(3,1) 0         Qxy(4,2) -Qxy(4,1) 0         Qxy(5,2) -Qxy(5,1) 0         Qxy(6,2) -Qxy(6,1) 0         Qxy(7,2) -Qxy(7,1) 0         Qxy(8,2) -Qxy(8,1) 0         Qxy(9,2) -Qxy(9,1) 0
                    0        0         shlp(1,1) 0        0         shlp(2,1) 0        0         shlp(3,1) 0        0         shlp(4,1) 0        0         shlp(5,1) 0        0         shlp(6,1) 0        0         shlp(7,1) 0        0         shlp(8,1) 0        0         shlp(9,1)];
            % Form BB matrix
            BBmat = [shgs(1,1) 0         0        shgs(2,1) 0         0        shgs(3,1) 0         0        shgs(4,1) 0         0        shgs(5,1) 0         0        shgs(6,1) 0         0        shgs(7,1) 0         0        shgs(8,1) 0         0        shgs(9,1) 0         0
                     shgs(1,2) 0         0        shgs(2,2) 0         0        shgs(3,2) 0         0        shgs(4,2) 0         0        shgs(5,2) 0         0        shgs(6,2) 0         0        shgs(7,2) 0         0        shgs(8,2) 0         0        shgs(9,2) 0         0
                     shgs(1,3) 0         0        shgs(2,3) 0         0        shgs(3,3) 0         0        shgs(4,3) 0         0        shgs(5,3) 0         0        shgs(6,3) 0         0        shgs(7,3) 0         0        shgs(8,3) 0         0        shgs(9,3) 0         0
                     0         shgs(1,1) 0        0         shgs(2,1) 0        0         shgs(3,1) 0        0         shgs(4,1) 0        0         shgs(5,1) 0        0         shgs(6,1) 0        0         shgs(7,1) 0        0         shgs(8,1) 0        0         shgs(9,1) 0
                     0         shgs(1,2) 0        0         shgs(2,2) 0        0         shgs(3,2) 0        0         shgs(4,2) 0        0         shgs(5,2) 0        0         shgs(6,2) 0        0         shgs(7,2) 0        0         shgs(8,2) 0        0         shgs(9,2) 0
                     0         shgs(1,3) 0        0         shgs(2,3) 0        0         shgs(3,3) 0        0         shgs(4,3) 0        0         shgs(5,3) 0        0         shgs(6,3) 0        0         shgs(7,3) 0        0         shgs(8,3) 0        0         shgs(9,3) 0
                     0         0         Pxy(1,1) 0         0         Pxy(2,1) 0         0         Pxy(3,1) 0         0         Pxy(4,1) 0         0         Pxy(5,1) 0         0         Pxy(6,1) 0         0         Pxy(7,1) 0         0         Pxy(8,1) 0         0         Pxy(9,1)
                     0         0         Pxy(1,2) 0         0         Pxy(2,2) 0         0         Pxy(3,2) 0         0         Pxy(4,2) 0         0         Pxy(5,2) 0         0         Pxy(6,2) 0         0         Pxy(7,2) 0         0         Pxy(8,2) 0         0         Pxy(9,2)];
            end

            [fi,JxX,F] = kine2d(Qxy,-ul,nel,0); %this is equivalent to ikine2d
            JxX = 1/JxX; %this is equivalent to ikine2d
    %         [F,JxX,fi] = ikine2d(Qxy,ul,nel); %this is equivalent to above two lines
    %         xs = sx;
            Jdet = Jdet/JxX;
            [theta,theta1,theta2,theta3] = ThetaS(JxX,mateprop);
%             [sigmai, cmati] = SigmaCmatNSCST2i(F,JxX,mateprop);

            w = Wgt*Jdet*thick;
            bb = be(3);

            spvec = JxX*theta1*spvec0;
            spmat = JxX*theta1*spmat0;
            cpmat = (theta2*JxX^2 + theta1*JxX)*cpmat1 + theta1*JxX*cpmat2;
            
            press = ul(3,:)*shlp;
            sigmap = press*spvec;
            cmatp = press*cpmat;
            
            [sigmai, cmati] = SigmaCmatNSCST2i(F,JxX,mateprop);
            sigma = sigmai + sigmap;
            cmat = cmati + cmatp;
            ddu = BBmat*reshape(ul,ndf*nel,1);
            Pvec = ddu(7:8);
            ddu = ddu(1:6);
            Pmat = [Pvec(1) 0 Pvec(2) 0 
                    0 Pvec(2) Pvec(1) 0];
            
            if exist('iprob','var')
                if iprob == 6
                    mu = 40;
                    X = xl(1,:)*shl;
                    Y = xl(2,:)*shl;
                    fb = [ - (101*mu)/(101*X + 100) - (101*lam*(101*X + 100))/10000
                                                      0; 0];
                else
                    fb = zeros(3,1);
                end
            else
                fb = zeros(3,1);
            end
            
            % Residual - NO TERMS 10 or 23, they were skipped in the numbering
            Fddu = F*ul(1:2,:)*shgs;
            Fddu4 = Fddu*I3;
            diagc = [cmat(1:3,1:3) zeros(3); zeros(3) cmat(1:3,1:3)];
            bigF = [F(1,1) 0 0 F(1,2) 0 0
                    0 0 F(2,1) 0 0 F(2,2)
                    F(2,1) 0 F(1,1) F(2,2) 0 F(1,2)
                    0 0 F(1,1) 0 0 F(1,2)
                    0 F(2,1) 0 0 F(2,2) 0
                    0 F(1,1) F(2,1) 0 F(1,2) F(2,2)];
            term2 = Fddu4*sigma; %page 6a
            term3 = I4*diagc*bigF*ddu; %page 6a
            term32 = spmat*Pvec; %page 6a
            ElemTR = bb*Tsmall*(term2 + term3 + term32 + fb(1:2));% % Added back for body force BF2U4M0.m
            term32M = [term32(1)*I2 term32(2)*I2];
%             term33 = [term32(1) 0 term32(2)/2  term32(2)/2
%                       0 term32(2) term32(1)/2 -term32(1)/2];
            term33 = term32M*P2;
            term34 = Pmat*cpmat;
            BT = [term33+term34 zeros(2,1); zeros(2,5)];
            BBT = [zeros(2,8); zeros(2,6) spmat];
            
            Tvec = [ElemTR; ElemTR];
            
            ElemF = ElemF - w*(Bmat'*(-BT'*Tvec) - BBmat'*BBT'*Tvec);

        end %je
    
if implicon == 1
        
        ElemF = zeros(nst,1);
        nele = nel;
        if((nele==3)||(nele==6))
            neleB = 3;
        else
            neleB = 4;
        end
        
%    -----------------------------------------------------
%     Loop over cells in element
%    -----------------------------------------------------
%     If Triangular Element
    
            if(neleB==3)

               celle = 0;

               for cellj = 1:minc
                  for celli = 1:2*(minc-cellj)+1
                 
                  celle = celle + 1;
                
%       Load cell coordinates, dof IDs

                  cell = cel*(elem-1) + celle;
                  ixc = zeros(nen,1);
   
                  for j = 1:nele
                     node = ixs(j,cell);
                     ixc(j) = node;
                     for i = 1:ndm
                         xc(i,j) = xs(i,node);
                     end
                     for i = 1:ndfs
                         uc(i,j) = us(i,node);
                     end
                  end

%       Set element-cell local coordinate map
                  evenodd = floor((celli+1)/2);
                  if(celli/2==evenodd) %even celli

	               %adjust order of cell flags for proper assembly
	               %due to reversal of element local coordinates
                   if(nele==3)
	                  i = ixc(1);
	                  ixc(1) = ixc(2);
	                  ixc(2) = ixc(3);
	                  ixc(3) = i;
                      for i = 1:ndm
	                     temp = xc(i,1);
	                     xc(i,1) = xc(i,2);
	                     xc(i,2) = xc(i,3);
	                     xc(i,3) = temp;
                      end
                      for i = 1:ndm
	                     temp = uc(i,1);
	                     uc(i,1) = uc(i,2);
	                     uc(i,2) = uc(i,3);
	                     uc(i,3) = temp;
                      end
	               else
	                  i = ixc(1);
	                  ixc(1) = ixc(2);
	                  ixc(2) = ixc(3);
	                  ixc(3) = i;
	                  i = ixc(4);
	                  ixc(4) = ixc(5);
	                  ixc(5) = ixc(6);
	                  ixc(6) = i;
                      for i = 1:ndm
	                     temp = xc(i,1);
	                     xc(i,1) = xc(i,2);
	                     xc(i,2) = xc(i,3);
	                     xc(i,3) = temp;
	                     temp = xc(i,4);
	                     xc(i,4) = xc(i,5);
	                     xc(i,5) = xc(i,6);
	                     xc(i,6) = temp;
                      end
                      for i = 1:ndfs
	                     temp = uc(i,1);
	                     uc(i,1) = uc(i,2);
	                     uc(i,2) = uc(i,3);
	                     uc(i,3) = temp;
	                     temp = uc(i,4);
	                     uc(i,4) = uc(i,5);
	                     uc(i,5) = uc(i,6);
	                     uc(i,6) = temp;
                      end
                   end
                     MR = -1.d0/minc;
                     BR = celli/(2.d0*minc);
                     MS = -1.d0/minc;
                     BS = cellj/(minc);
                  else
                     MR = 1.d0/minc;
                     BR = ((celli+1)/2.d0-1.d0)/minc;
                     MS = 1.d0/minc;
                     BS = (cellj-1.d0)/minc;
                  end

%      -----------------------------------------------------
%       Compute cell residual
%      -----------------------------------------------------
         eF = getFcellS(xc,uc,mateprop,iprob,xl,ul,ndf,ndfs,ndm,lint,nele,nen,nummat,MR,BR,MS,BS,t11,t12,t21,t22);
         
         ElemF = ElemF + eF(1:nst);

%    -----------------------------------------------------
%     End loop over cells
%    -----------------------------------------------------
                  end
               end
%%
%     Else Quadrilateral Element
            else
              
               celle = 0;
               for cellj = 1:minc
                  for celli = 1:minc
                 
                  celle = celle + 1;

%       Load cell coordinates, dof IDs

                  cell = cel*(elem-1) + celle;
                  ixc = zeros(nen,1);
   
                  for j = 1:nele
                     node = ixs(j,cell);
                     ixc(j) = node;
                     for i = 1:ndm
                         xc(i,j) = xs(i,node);
                     end
                     for i = 1:ndfs
                         uc(i,j) = us(i,node);
                     end
                  end

%       Set element-cell local coordinate map
                  %MR = 1.d0/dble(m) 
                  MR = 1.d0/minc;
                  BR = -1.d0+(2.d0*celli-1.d0)/minc;
                  MS = 1.d0/minc;
                  BS = -1.d0+(2.d0*cellj-1.d0)/minc;

%      -----------------------------------------------------
%       Compute cell residual
%      -----------------------------------------------------
         eF = getFcellS(xc,uc,mateprop,iprob,xl,ul,ndf,ndfs,ndm,lint,nele,nen,nummat,MR,BR,MS,BS,t11,t12,t21,t22);
         
         ElemF = ElemF + eF(1:nst);

%    -----------------------------------------------------
%     End loop over cells
%    -----------------------------------------------------
                  end
               end

            end
%             ElemF

end
%%
    case 21
        
        ElemK = zeros(nst);
%         ElemP = zeros(2,1);
%         ElemR = zeros(2,1);
%         Nmat = zeros(3,3*nel);
        Bmat = zeros(4,3*nel);
%         BBmat = zeros(6,3*nel);

        lam = getlam(mateprop);
        
        I1 = [1; 1; 0; 0];
        I2 = eye(2);
        I6 = eye(4);
        P2 = [1 0 0 0
              0 0 1/2 -1/2
              0 0 1/2 1/2
              0 1 0 0];
        P3 = [1 0 0 0
              0 0 0 1
              0 1 1 0
              0 -1 1 0];
        P4 = [1 0 0 0
              0 0 1/2 1/2
              0 0 1/2 -1/2
              0 1 0 0];
        Z2 = zeros(2);
        spvec0 = I1;
        spmat0 = I2;
        cpmat1 = I1*I1';
        cpmat2 = diag([-2,-2,-1,0]);
        cpmat0 = cpmat1 + cpmat2;
        dpmat1 = [cpmat1; cpmat1; zeros(4,4)];
        dpmat2 =-[6 2 0 0
                  2 2 0 0
                  0 0 1 0
                  0 0 0 0
                  2 2 0 0
                  2 6 0 0
                  0 0 1 0
                  0 0 0 0
                  0 0 1 0
                  0 0 1 0
                  1 1 0 0
                  0 0 0 0];
        dpmat3 = [8 0 0 0
                  0 0 0 0
                  0 0 2 0
                  0 0 0 0
                  0 0 0 0
                  0 8 0 0
                  0 0 2 0
                  0 0 0 0
                  0 0 2 0
                  0 0 2 0
                  2 2 0 0
                  0 0 0 0];
        dpmat0 = dpmat1 + dpmat2 + dpmat3;
        
        thick = 1;
        
        % Load Guass Integration Points

        if nel == 3 || nel == 6
            lint = 13;25;
        elseif nel == 4
            lint = 100;4;
%              lint = 16;
        else
            lint = 100;9;
        end
        der = 1;
        bf = 1;
        ib = 0;

        %Form tau matrix: compute on first iteration and store; this gives
        %quadratic convergence when tau isn't recalculated (no
        %linearization of tau required).
        if iter <=2 % == 0 % 
        [t11,t12,t21,t22] = TauS(xl,ul,mateprop,nel,nen,nelP,lint);
%         TauList(elem,1:4) = [t11,t12,t21,t22];
        hr(nhc:nhc+3) = [t11,t12,t21,t22];
        else
%         t11 = TauList(elem,1);
%         t12 = TauList(elem,2);
%         t21 = TauList(elem,3);
%         t22 = TauList(elem,4);
        t11 = hr(nhc+0);
        t12 = hr(nhc+1);
        t21 = hr(nhc+2);
        t22 = hr(nhc+3);
        end
        
        %Evaluate first derivatives of basis functions at int. point
          [Wgt,litr,lits] =  intpntt(1,1,ib);
          [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
          [Qxy,shgs,Jdet,be,sx] = shgt(xl+ul(1:2,:),nel,shld,shls,nel,bf,der,be);
%           sx = inv(sx); Correct, but only needed when correct derivs of bubble are required
                
        Tsmall = [t11 t12; t21 t22]; %Modified 8/23 for stored tau
        
if elem == 25
    t11;
end
        %Integration Loop
        for ll = 1:lint

            %Evaluate first derivatives of basis functions at int. point
            if nel == 3 || nel == 6
              [Wgt,litr,lits] =  intpntt(ll,lint,ib);
              [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
              [Qxy, shgs, Jdet,be] = shgt(xl+ul(1:2,:),nel,shld,shls,nel,bf,der,be);
            else
              [Wgt,litr,lits] =  intpntq(ll,lint,ib);
              [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
              [Qxy, shgs, Jdet,be] = shgq(xl+ul(1:2,:),nel,shld,shls,nel,bf,der,be);
            end

            if nelP == 3 || nelP == 6
              [shlp,shld,shls] = shlt(litr,lits,nelP,nel,0,0);
              [Pxy, shgp] = shgt(xl+ul(1:2,:),nelP,shld,shls,nelP,0,0,be);
            else
              [shlp,shld,shls] = shlq(litr,lits,nelP,nel,0,0);
              [Pxy, shgp] = shgq(xl+ul(1:2,:),nelP,shld,shls,nelP,0,0,be);
            end

            if Jdet < 0
                Jdet
            end

            if nel == 3            
            % Form B matrix    
            Bmat = [Qxy(1,1)  0        0         Qxy(2,1)  0        0         Qxy(3,1)  0        0
                    0         Qxy(1,2) 0         0         Qxy(2,2) 0         0         Qxy(3,2) 0         
                    Qxy(1,2)  Qxy(1,1) 0         Qxy(2,2)  Qxy(2,1) 0         Qxy(3,2)  Qxy(3,1) 0         
                    Qxy(1,2) -Qxy(1,1) 0         Qxy(2,2) -Qxy(2,1) 0         Qxy(3,2) -Qxy(3,1) 0         
                    0        0         shlp(1,1) 0        0         shlp(2,1) 0        0         shlp(3,1)];
            % Form BB matrix
            BBmat = [zeros(6,9)        
                     0         0        Qxy(1,1) 0         0        Qxy(2,1) 0         0        Qxy(3,1) 
                     0         0        Qxy(1,2) 0         0        Qxy(2,2) 0         0        Qxy(3,2)];
            elseif nel == 4            
            % Form B matrix
            Bmat = [Qxy(1,1)  0        0         Qxy(2,1)  0        0         Qxy(3,1)  0        0         Qxy(4,1)  0        0
                    0         Qxy(1,2) 0         0         Qxy(2,2) 0         0         Qxy(3,2) 0         0         Qxy(4,2) 0         
                    Qxy(1,2)  Qxy(1,1) 0         Qxy(2,2)  Qxy(2,1) 0         Qxy(3,2)  Qxy(3,1) 0         Qxy(4,2)  Qxy(4,1) 0         
                    Qxy(1,2) -Qxy(1,1) 0         Qxy(2,2) -Qxy(2,1) 0         Qxy(3,2) -Qxy(3,1) 0         Qxy(4,2) -Qxy(4,1) 0         
                    0        0         shlp(1,1) 0        0         shlp(2,1) 0        0         shlp(3,1) 0        0         shlp(4,1)];
            % Form BB matrix
            BBmat = [shgs(1,1) 0         0        shgs(2,1) 0         0        shgs(3,1) 0         0        shgs(4,1) 0         0
                     shgs(1,2) 0         0        shgs(2,2) 0         0        shgs(3,2) 0         0        shgs(4,2) 0         0
                     shgs(1,3) 0         0        shgs(2,3) 0         0        shgs(3,3) 0         0        shgs(4,3) 0         0
                     0         shgs(1,1) 0        0         shgs(2,1) 0        0         shgs(3,1) 0        0         shgs(4,1) 0
                     0         shgs(1,2) 0        0         shgs(2,2) 0        0         shgs(3,2) 0        0         shgs(4,2) 0
                     0         shgs(1,3) 0        0         shgs(2,3) 0        0         shgs(3,3) 0        0         shgs(4,3) 0
                     0         0         Pxy(1,1) 0         0        Pxy(2,1) 0         0        Pxy(3,1) 0         0        Pxy(4,1) 
                     0         0         Pxy(1,2) 0         0        Pxy(2,2) 0         0        Pxy(3,2) 0         0        Pxy(4,2)];
            elseif nel == 6
            % Form B matrix
            Bmat = [Qxy(1,1)  0        0         Qxy(2,1)  0        0         Qxy(3,1)  0        0         Qxy(4,1)  0        0         Qxy(5,1)  0        0         Qxy(6,1)  0        0
                    0         Qxy(1,2) 0         0         Qxy(2,2) 0         0         Qxy(3,2) 0         0         Qxy(4,2) 0         0         Qxy(5,2) 0         0         Qxy(6,2) 0
                    Qxy(1,2)  Qxy(1,1) 0         Qxy(2,2)  Qxy(2,1) 0         Qxy(3,2)  Qxy(3,1) 0         Qxy(4,2)  Qxy(4,1) 0         Qxy(5,2)  Qxy(5,1) 0         Qxy(6,2)  Qxy(6,1) 0
                    Qxy(1,2) -Qxy(1,1) 0         Qxy(2,2) -Qxy(2,1) 0         Qxy(3,2) -Qxy(3,1) 0         Qxy(4,2) -Qxy(4,1) 0         Qxy(5,2) -Qxy(5,1) 0         Qxy(6,2) -Qxy(6,1) 0
                    0        0         shlp(1,1) 0        0         shlp(2,1) 0        0         shlp(3,1) 0        0         shlp(4,1) 0        0         shlp(5,1) 0        0         shlp(6,1)];
            % Form BB matrix
            BBmat = [shgs(1,1) 0         0        shgs(2,1) 0         0        shgs(3,1) 0         0        shgs(4,1) 0         0        shgs(5,1) 0         0        shgs(6,1) 0         0
                     shgs(1,2) 0         0        shgs(2,2) 0         0        shgs(3,2) 0         0        shgs(4,2) 0         0        shgs(5,2) 0         0        shgs(6,2) 0         0
                     shgs(1,3) 0         0        shgs(2,3) 0         0        shgs(3,3) 0         0        shgs(4,3) 0         0        shgs(5,3) 0         0        shgs(6,3) 0         0
                     0         shgs(1,1) 0        0         shgs(2,1) 0        0         shgs(3,1) 0        0         shgs(4,1) 0        0         shgs(5,1) 0        0         shgs(6,1) 0
                     0         shgs(1,2) 0        0         shgs(2,2) 0        0         shgs(3,2) 0        0         shgs(4,2) 0        0         shgs(5,2) 0        0         shgs(6,2) 0
                     0         shgs(1,3) 0        0         shgs(2,3) 0        0         shgs(3,3) 0        0         shgs(4,3) 0        0         shgs(5,3) 0        0         shgs(6,3) 0
                     0         0         Pxy(1,1) 0         0         Pxy(2,1) 0         0         Pxy(3,1) 0         0         Pxy(4,1) 0         0         Pxy(5,1) 0         0         Pxy(6,1)
                     0         0         Pxy(1,2) 0         0         Pxy(2,2) 0         0         Pxy(3,2) 0         0         Pxy(4,2) 0         0         Pxy(5,2) 0         0         Pxy(6,2)];
            elseif nel == 9
            % Form B matrix
            Bmat = [Qxy(1,1)  0        0         Qxy(2,1)  0        0         Qxy(3,1)  0        0         Qxy(4,1)  0        0         Qxy(5,1)  0        0         Qxy(6,1)  0        0         Qxy(7,1)  0        0         Qxy(8,1)  0        0         Qxy(9,1)  0        0
                    0         Qxy(1,2) 0         0         Qxy(2,2) 0         0         Qxy(3,2) 0         0         Qxy(4,2) 0         0         Qxy(5,2) 0         0         Qxy(6,2) 0         0         Qxy(7,2) 0         0         Qxy(8,2) 0         0         Qxy(9,2) 0
                    Qxy(1,2)  Qxy(1,1) 0         Qxy(2,2)  Qxy(2,1) 0         Qxy(3,2)  Qxy(3,1) 0         Qxy(4,2)  Qxy(4,1) 0         Qxy(5,2)  Qxy(5,1) 0         Qxy(6,2)  Qxy(6,1) 0         Qxy(7,2)  Qxy(7,1) 0         Qxy(8,2)  Qxy(8,1) 0         Qxy(9,2)  Qxy(9,1) 0
                    Qxy(1,2) -Qxy(1,1) 0         Qxy(2,2) -Qxy(2,1) 0         Qxy(3,2) -Qxy(3,1) 0         Qxy(4,2) -Qxy(4,1) 0         Qxy(5,2) -Qxy(5,1) 0         Qxy(6,2) -Qxy(6,1) 0         Qxy(7,2) -Qxy(7,1) 0         Qxy(8,2) -Qxy(8,1) 0         Qxy(9,2) -Qxy(9,1) 0
                    0        0         shlp(1,1) 0        0         shlp(2,1) 0        0         shlp(3,1) 0        0         shlp(4,1) 0        0         shlp(5,1) 0        0         shlp(6,1) 0        0         shlp(7,1) 0        0         shlp(8,1) 0        0         shlp(9,1)];
            % Form BB matrix
            BBmat = [shgs(1,1) 0         0        shgs(2,1) 0         0        shgs(3,1) 0         0        shgs(4,1) 0         0        shgs(5,1) 0         0        shgs(6,1) 0         0        shgs(7,1) 0         0        shgs(8,1) 0         0        shgs(9,1) 0         0
                     shgs(1,2) 0         0        shgs(2,2) 0         0        shgs(3,2) 0         0        shgs(4,2) 0         0        shgs(5,2) 0         0        shgs(6,2) 0         0        shgs(7,2) 0         0        shgs(8,2) 0         0        shgs(9,2) 0         0
                     shgs(1,3) 0         0        shgs(2,3) 0         0        shgs(3,3) 0         0        shgs(4,3) 0         0        shgs(5,3) 0         0        shgs(6,3) 0         0        shgs(7,3) 0         0        shgs(8,3) 0         0        shgs(9,3) 0         0
                     0         shgs(1,1) 0        0         shgs(2,1) 0        0         shgs(3,1) 0        0         shgs(4,1) 0        0         shgs(5,1) 0        0         shgs(6,1) 0        0         shgs(7,1) 0        0         shgs(8,1) 0        0         shgs(9,1) 0
                     0         shgs(1,2) 0        0         shgs(2,2) 0        0         shgs(3,2) 0        0         shgs(4,2) 0        0         shgs(5,2) 0        0         shgs(6,2) 0        0         shgs(7,2) 0        0         shgs(8,2) 0        0         shgs(9,2) 0
                     0         shgs(1,3) 0        0         shgs(2,3) 0        0         shgs(3,3) 0        0         shgs(4,3) 0        0         shgs(5,3) 0        0         shgs(6,3) 0        0         shgs(7,3) 0        0         shgs(8,3) 0        0         shgs(9,3) 0
                     0         0         Pxy(1,1) 0         0         Pxy(2,1) 0         0         Pxy(3,1) 0         0         Pxy(4,1) 0         0         Pxy(5,1) 0         0         Pxy(6,1) 0         0         Pxy(7,1) 0         0         Pxy(8,1) 0         0         Pxy(9,1)
                     0         0         Pxy(1,2) 0         0         Pxy(2,2) 0         0         Pxy(3,2) 0         0         Pxy(4,2) 0         0         Pxy(5,2) 0         0         Pxy(6,2) 0         0         Pxy(7,2) 0         0         Pxy(8,2) 0         0         Pxy(9,2)];
            end
            
            bb = be(3);

            [fi,JxX,F] = kine2d(Qxy,-ul,nel,0); %this is equivalent to ikine2d
            JxX = 1/JxX; %this is equivalent to ikine2d
    %         [F,JxX,fi] = ikine2d(Qxy,ul,nel); %this is equivalent to above two lines
    %         xs = sx;
            Jdet = Jdet/JxX;
            [theta,theta1,theta2,theta3] = ThetaS(JxX,mateprop);
            Pvec = (ul(3,:)*Qxy)';
            Pmat = [Pvec(1) 0 Pvec(2) 0 
                    0 Pvec(2) Pvec(1) 0];
            [sigmai, cmati] = SigmaCmatNSCST2i(F,JxX,mateprop);

            w = Wgt*Jdet*thick;

            spvec = JxX*theta1*spvec0;
            spmat = JxX*theta1*spmat0;
            cpmat = (theta2*JxX^2 + theta1*JxX)*cpmat1 + theta1*JxX*cpmat2;
            A = JxX*theta1 + 3*JxX^2*theta2 + JxX^3*theta3;
            B = JxX*theta1 + JxX^2*theta2;
            C = JxX*theta1;
            dpmat = (A*dpmat1 + B*dpmat2 + C*dpmat3);

            press = ul(3,:)*shlp;
            sigmap = press*spvec;
            cmatp = press*cpmat;
            sigma = sigmai + sigmap;
            cmat = cmati + cmatp;
            
            term32 = spmat*Pvec;
            ElemTR = bb*Tsmall*term32;
            term32M = [term32(1)*I2 term32(2)*I2];
%             term33 = [term32(1) 0 term32(2)/2  term32(2)/2
%                       0 term32(2) term32(1)/2 -term32(1)/2];
            term33 = term32M*P2;
            term34 = Pmat*cpmat;
            BT = [term33+term34 zeros(2,1); zeros(2,5)];
            BBT = [zeros(2,8); zeros(2,6) spmat];
            Svec = [sigma; (theta-press/lam)];
            sig19 = cpmat*Pmat'*ElemTR;
            sigma2 = sigma - sig19;
            Smat = [sigma2(1) 0  sigma2(3)/2 sigma2(3)/2
                    0 sigma2(2)  sigma2(3)/2 -sigma2(3)/2
                    sigma2(3)/2  sigma2(3)/2 (sigma2(2)+sigma2(1))/4 (sigma2(2)-sigma2(1))/4
                    sigma2(3)/2 -sigma2(3)/2 (sigma2(2)-sigma2(1))/4 (sigma2(2)+sigma2(1))/4];
            cvec20 = -dpmat*Pmat'*ElemTR;
            cmat20 = [cvec20(1:4) cvec20(5:8) cvec20(9:12) zeros(4,1)];
%             trp = [ElemTR(1)*Pvec(1) 0 ElemTR(2)*Pvec(1)/2 -ElemTR(2)*Pvec(1)/2
%                    0 ElemTR(2)*Pvec(2) ElemTR(1)*Pvec(2)/2  ElemTR(1)*Pvec(2)/2
%                    ElemTR(1)*Pvec(2) ElemTR(2)*Pvec(1) (ElemTR(1)*Pvec(1)+ElemTR(2)*Pvec(2))/2 (ElemTR(1)*Pvec(1)-ElemTR(2)*Pvec(2))/2
%                    zeros(1,4)];
            trp = [ElemTR(1)*I6*Pmat' ElemTR(2)*I6*Pmat']*P4;
            term21 = -cpmat*trp - trp'*cpmat;
            ITR2 = [ElemTR(1)*I2; ElemTR(2)*I2];
%             term25 = -cpmat*[ElemTR(1) 0
%                              0 ElemTR(2)
%                              ElemTR(2) ElemTR(1)
%                              zeros(1,2)];
            term25 = -cpmat*P3*ITR2;
%             term26 = -spmat*[ElemTR(1) 0 ElemTR(2)/2 -ElemTR(2)/2
%                            0 ElemTR(2) ElemTR(1)/2  ElemTR(1)/2];
            term26 = -spmat*ITR2'*P4;
            
            Tmat = bb*[Tsmall Z2; Z2 Tsmall];
            Tmat2 = bb*[Z2 Tsmall; Tsmall Z2];
            
            D11 = [Smat+cmat+cmat20+term21 spvec
                   spvec' -1/(lam)] - BT'*Tmat*BT - BT'*Tmat2*BT;
            D12 = [zeros(4,6) term25+term26'
                   zeros(1,8)] - BT'*Tmat*BBT - BT'*Tmat2*BBT;
            D21 = D12';
            D22 = - BBT'*Tmat*BBT - BBT'*Tmat2*BBT;
            
            ElemK = ElemK + w*(Bmat'*D11*Bmat + BBmat'*D22*BBmat ...
                    + BBmat'*D21*Bmat + Bmat'*D12*BBmat);

        end %je
%%        
    case 22
        
        ElemM = zeros(nst);
        ElemF = zeros(nst,1);

        lam = getlam(mateprop);
        
        I1 = [1; 1; 1; 0; 0; 0; 0; 0; 0];
        stresID = [1 2 4 0 1 3];
        spvec0 = I1;
        
        thick = 1;
        
        % Load Guass Integration Points

        if nel == 3 || nel == 6
            lint = 13;25;
        elseif nel == 4
            lint = 100;4;
%              lint = 16;
        else
            lint = 100;9;
        end
        der = 0;
        bf = 1;
        ib = 0;

        %Integration Loop
        for ll = 1:lint

            %Evaluate first derivatives of basis functions at int. point
            if nel == 3 || nel == 6
              [Wgt,litr,lits] =  intpntt(ll,lint,ib);
              [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
              [Qxy, shgs, Jdet,be] = shgt(xl+ul(1:2,:),nel,shld,shls,nel,bf,der,be);
            else
              [Wgt,litr,lits] =  intpntq(ll,lint,ib);
              [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
              [Qxy, shgs, Jdet,be] = shgq(xl+ul(1:2,:),nel,shld,shls,nel,bf,der,be);
            end

            if nelP == 3 || nelP == 6
              [shlp,shld,shls] = shlt(litr,lits,nelP,nel,0,0);
              [Pxy, shgp] = shgt(xl+ul(1:2,:),nelP,shld,shls,nelP,0,0,be);
            else
              [shlp,shld,shls] = shlq(litr,lits,nelP,nel,0,0);
              [Pxy, shgp] = shgq(xl+ul(1:2,:),nelP,shld,shls,nelP,0,0,be);
            end

            if nelS == 3 || nelS == 6
              [shlS,shld,shls] = shlt(litr,lits,nelS,nel,0,0);
            else
              [shlS,shld,shls] = shlq(litr,lits,nelS,nel,0,0);
            end
            
            [fi,JxX,F] = kine2d(Qxy,-ul,nel,0); %this is equivalent to ikine2d
            JxX = 1/JxX; %this is equivalent to ikine2d
    %         [F,JxX,fi] = ikine2d(Qxy,ul,nel); %this is equivalent to above two lines
            Jdet = Jdet/JxX;
            [theta,theta1,theta2,theta3] = ThetaS(JxX,mateprop);
            spvec = JxX*theta1*spvec0;

            sigmai = SigmaCmatNSCST3i([F zeros(2,1); zeros(1,2) 1],JxX,mateprop);
            press = ul(3,:)*shlp;
            sigmap = press*spvec;
            sigma = (sigmai + sigmap)/JxX;
            
            % Form B matrix
            Nmat = shlS';

            w = Wgt*Jdet*thick;
            
            if stres <= 3 % stress components
                sigmas = sigma(stresID(stres));
            elseif stres >= 5
                if stres <= 6 % principal stresses
                sigma = [sigma(1) sigma(4) sigma(6); sigma(4) sigma(2) sigma(5); sigma(6) sigma(5) sigma(3)];
                psig = eig(sigma);
                sigmas = psig(stresID(stres));
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
        stresID = [1 2 4 0 1 3];
        spvec0 = I1;
        
        thick = 1;
        
        % Load Guass Integration Points

        if nel == 3
            lint = 1;
            nint = 1;
        elseif nel == 4
%             lint = 4;
            lint = 4;
            nint = 1;
        elseif nel == 6
            lint = 7;
            nint = 3;
        else
            lint = 9;
            nint = 4;
        end
        
        der = 0;
        bf = 1;
        ib = 0;

        %Stress Loop
        for ll = 1:nint

            %Evaluate first derivatives of basis functions at int. point
            if nel == 3 || nel == 6
              [Wgt,litr,lits] =  intpntt(ll,nint,ib);
              [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
              [Qxy, shgs, Jdet,be] = shgt(xl+ul(1:2,:),nel,shld,shls,nel,bf,der,be);
            else
              [Wgt,litr,lits] =  intpntq(ll,nint,ib);
              [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
              [Qxy, shgs, Jdet,be] = shgq(xl+ul(1:2,:),nel,shld,shls,nel,bf,der,be);
            end

            if nelP == 3 || nelP == 6
              shlp = shlt(litr,lits,nelP,nel,0,0);
            else
              shlp = shlq(litr,lits,nelP,nel,0,0);
            end
            
            [fi,JxX,F] = kine2d(Qxy,-ul,nel,0); %this is equivalent to ikine2d
            JxX = 1/JxX; %this is equivalent to ikine2d
    %         [F,JxX,fi] = ikine2d(Qxy,ul,nel); %this is equivalent to above two lines
            Jdet = Jdet/JxX;
            [theta,theta1,theta2,theta3] = ThetaS(JxX,mateprop);
            spvec = JxX*theta1*spvec0;

            sigmai = SigmaCmatNSCST3i([F zeros(2,1); zeros(1,2) 1],JxX,mateprop);
            press = ul(3,:)*shlp;
            sigmap = press*spvec;
            sigma = (sigmai + sigmap)/JxX;
            
            for stres = 1:npstr
            
            if stres <= 3 % stress components
                sigmas = sigma(stresID(stres));
            elseif stres >= 5
                if stres <= 6 % principal stresses
                sigma2 = [sigma(1) sigma(4) sigma(6); sigma(4) sigma(2) sigma(5); sigma(6) sigma(5) sigma(3)];
                psig = eig(sigma2);
                sigmas = psig(stresID(stres));
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
        if nel == 3
            plist = [0 1 0
                     0 0 1];
        elseif nel == 4
            plist = [-1 1 1 -1
                     -1 -1 1 1];
        elseif nel == 6
            plist = [-1/3 5/3 -1/3 2/3 2/3 -1/3
                     -1/3 -1/3 5/3 -1/3 2/3 2/3];
        else
            sqr3 = sqrt(3);
            plist = [-sqr3 sqr3 sqr3 -sqr3 0 sqr3 0 -sqr3 0
                     -sqr3 -sqr3 sqr3 sqr3 -sqr3 0 sqr3 0 0];
        end
        
        for ll = 1:nelS
            
            r = plist(1,ll);
            s = plist(2,ll);
            shpS = sshp2d(r,s,nint);
            
            for stres = 1:npstr
                
                sigmas = ElemS2(1:nint,stres)'*shpS;
                ElemS(ll,stres) = sigmas;
                
            end
            
        end
        
        %Integration Loop
        Vol = 0;
        for ll = 1:lint

            %Evaluate first derivatives of basis functions at int. point
            if nel == 3 || nel == 6
              [Wgt,litr,lits] =  intpntt(ll,lint,ib);
              [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
              [Qxy, shgs, Jdet,be] = shgt(xl+ul(1:2,:),nel,shld,shls,nel,bf,der,be);
            else
              [Wgt,litr,lits] =  intpntq(ll,lint,ib);
              [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
              [Qxy, shgs, Jdet,be] = shgq(xl+ul(1:2,:),nel,shld,shls,nel,bf,der,be);
            end
            
            [fi,JxX,F] = kine2d(Qxy,-ul,nel,0); %this is equivalent to ikine2d
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