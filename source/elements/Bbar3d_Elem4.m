% Tim Truster
% CEE 577
% 04/04/2012
% Midterm Exam
% 3-D B-bar element for small strain J2-Plasticity with linear
% isotrop/kinematic hardening
% Bilinear quadrilateral element with consistent tangent

% Tested for single element shear and distorted loading (one nodal force
% removed); quadratic NR convergence rate was maintained for non-Bbar and 
% for Bbar. Also, the condition (4.3.29) in Simo-Hughes was verified
% numerically.

% Verified against NL Mixed Inelasticity on 6/30/2013 using IME_Q3_3d.m

% Set Material Properties

ElemYM = mateprop(2);
Elemv = mateprop(3);
thick = mateprop(1);
Khard = mateprop(4);
Hhard = mateprop(5);
sigy = mateprop(6);

switch isw %Task Switch
%%    
    case 1
        
        if ndf > 3
            
            for i = 4:ndf
                lie(i,1) = 0;
            end
            
        end
        
        nh1 = (6+5+1)*8;
        
%%
    case 3
        
        bulk = ElemYM/(3*(1-2*Elemv));%Elemv*ElemYM/((1+Elemv)*(1-2*Elemv));
        mu = ElemYM/(2*(1+Elemv));
        One = [1.0d0; 1.0d0; 1.0d0; 0.0d0; 0.0d0; 0.0d0];
        I4 = diag([1.0d0 1.0d0 1.0d0 0.5d0 0.5d0 0.5d0]);
        OneOne = (One*One');

        % Initialize Matrix and Vector

        nst = nel*ndf;
        ElemK = zeros(nst);
        ElemF = zeros(nst,1);
        
        % Load Guass Integration Points
        lint = 8;

        der = 0;
        bf = 0;
        ib = 0;

        triagNodes = [1 5 7 8
                      2 5 9 6
                      3 7 6 10
                      4 8 10 9
                      5 7 8 9
                      5 7 9 6
                      7 8 9 10
                      7 6 10 9];

%         Bmat = zeros(6,3*nel);
%         bvecT = zeros(1,3*nel);
        bbarvecT = zeros(1,3*nel);
        Gamvec = 1;
        H = 0;
%         Pdev = I4 - 1/3*OneOne;
        
        % Compute bbarvec
        for l = 1:lint

            nodeT = triagNodes(l,1:4);
            xlT = xl(1:3,nodeT);
%                 if nel == 4 || nel == 10
                  [Wgt,ss] =  int3d_t(1,1,ib);
                  [shl,shld,shls,be] = shltt(ss,4,4,der,bf);
                  [Qxy, shgs, Jdet] = shgtt(xlT,4,shld,shls,4,bf,der,be);
%                 else
%                   [Wgt,ss] =  intpntb(l,lint,ib);
%                   [shl,shld,shls,be] = shlb(ss,nel,nel,der,bf);
%                   [Qxy, shgs, Jdet] = shgb(xl,nel,shld,shls,nel,bf,der,be);
%                 end

                % Form b matrix
                bvecT = zeros(1,3*nel);
                for i = 1:4
                ie = nodeT(i);
                bvecT(1,(ie-1)*3+1:3*ie) = [Qxy(i,1) Qxy(i,2) Qxy(i,3)];
                end

                % Update integration weighting factor
                W = Wgt*Jdet;

                bbarvecT = bbarvecT + W*Gamvec*bvecT;
                H = H + W*(Gamvec*Gamvec');

        %     end %ie
        end %j
        
        % Loop over integration points
        for l = 1:lint

            nodeT = triagNodes(l,1:4);
            xlT = xl(1:3,nodeT);
            
%                 if nel == 4 || nel == 10
                  [Wgt,ss] =  int3d_t(1,1,ib);
                  [shl,shld,shls,be] = shltt(ss,4,4,der,bf);
                  [Qxy, shgs, Jdet] = shgtt(xlT,4,shld,shls,4,bf,der,be);
%                 else
%                   [Wgt,ss] =  intpntb(l,lint,ib);
%                   [shl,shld,shls,be] = shlb(ss,nel,nel,der,bf);
%                   [Qxy, shgs, Jdet] = shgb(xl,nel,shld,shls,nel,bf,der,be);
%                 end

                % Form B matrix
                Bmat = zeros(6,3*nel);
                bvecT = zeros(1,3*nel);
                for i = 1:4
                ie = nodeT(i);
                Bmat(:,(ie-1)*3+1:3*ie) = [Qxy(i,1) 0        0
                                           0        Qxy(i,2) 0
                                           0         0        Qxy(i,3)
                                           Qxy(i,2) Qxy(i,1) 0
                                           0         Qxy(i,3) Qxy(i,2)
                                           Qxy(i,3) 0         Qxy(i,1)];
                bvecT(1,(ie-1)*3+1:3*ie) = [Qxy(i,1) Qxy(i,2) Qxy(i,3)];
                end
                
                Bbar = Bmat - 1/3*One*bvecT + 1/3*One*Gamvec'/H*bbarvecT;
                pn_1 = bulk*Gamvec'/H*bbarvecT*reshape(ul,nst,1);
                
                % Compute input for Radial Return
                eps3d = Bbar*reshape(ul,nst,1); %3D enhanced strain
%                 pn_1 = bulk*One'*eps3d;
                ephr = nh1-1+(l-1)*11; %pointer for plastic strain at pt l
                betahr = nh1-1+(l-1)*11+5;
                ahr = nh1-1+l*11;
                ep_n = [hr(ephr+1); hr(ephr+2); 0; hr(ephr+3); hr(ephr+4); hr(ephr+5)];
                ep_n(3) = -(ep_n(1) + ep_n(2));
                beta_n = [hr(betahr+1); hr(betahr+2); 0; hr(betahr+3); hr(betahr+4); hr(betahr+5)];
                beta_n(3) = -(beta_n(1) + beta_n(2));
                a_n = hr(ahr);
                
                [sig_n1,Cmat,ep_n1,beta_n1,a_n1,s_n1,Cdev_n1] = J2RadialReturn0(eps3d,ep_n,beta_n,a_n,mu,bulk,Khard,Hhard,sigy);
                
                % Retrieve deviatoric output
                sdev3 = s_n1;%sig_n1 - 1/3*OneOne*sig_n1;
                
                % Store history variables
                ephr = nh2-1+(l-1)*11; %pointer for plastic strain at pt l
                betahr = nh2-1+(l-1)*11+5;
                ahr = nh2-1+l*11;
                hr(ephr+1) = ep_n1(1);
                hr(ephr+2) = ep_n1(2);
                hr(ephr+3) = ep_n1(4);
                hr(ephr+4) = ep_n1(5);
                hr(ephr+5) = ep_n1(6);
                hr(betahr+1) = beta_n1(1);
                hr(betahr+2) = beta_n1(2);
                hr(betahr+3) = beta_n1(4);
                hr(betahr+4) = beta_n1(5);
                hr(betahr+5) = beta_n1(6);
                hr(ahr) = a_n1;

                % Update integration weighting factor
                W = Wgt*Jdet;

                ElemF = ElemF - W*Bbar'*(pn_1*One + sdev3);
                ElemK = ElemK + W*Bbar'*(Cdev_n1)*Bbar;
%                 ElemK = ElemK + W*Bmat'*(Cmat)*Bmat;

        %     end %ie
        end %je
        
        Gambbar = One*bbarvecT;
        ElemK = ElemK + bulk/ndm*Gambbar'/H*Gambbar;

    case 24
        
        ElemP = zeros(12,nel);
                
        bulk = ElemYM/(3*(1-2*Elemv));%Elemv*ElemYM/((1+Elemv)*(1-2*Elemv));
        mu = ElemYM/(2*(1+Elemv));
        One = [1.0d0; 1.0d0; 1.0d0; 0.0d0; 0.0d0; 0.0d0];
        I4 = diag([1.0d0 1.0d0 1.0d0 0.5d0 0.5d0 0.5d0]);
        OneOne = (One*One');

        % Initialize Matrix and Vector

        nst = nel*ndf;
        ElemK = zeros(nst);
        ElemF = zeros(nst,1);
        
        % Load Guass Integration Points
        lint = 8;

        der = 0;
        bf = 0;
        ib = 0;

        Bmat = zeros(6,3*nel);
        bvecT = zeros(1,3*nel);
        bbarvecT = zeros(1,3*nel);
        Gamvec = 1;
        H = 0;
        Pdev = I4 - 1/3*OneOne;
        
        % Compute bbarvec
        for l = 1:lint

                if nel == 4 || nel == 10
                  [Wgt,ss] =  int3d_t(l,lint,ib);
                  [shl,shld,shls,be] = shltt(ss,nel,nel,der,bf);
                  [Qxy, shgs, Jdet] = shgtt(xl,nel,shld,shls,nel,bf,der,be);
                else
                  [Wgt,ss] =  intpntb(l,lint,ib);
                  [shl,shld,shls,be] = shlb(ss,nel,nel,der,bf);
                  [Qxy, shgs, Jdet] = shgb(xl,nel,shld,shls,nel,bf,der,be);
                end

                % Form b matrix
                for ie = 1:nel
                bvecT(1,(ie-1)*3+1:3*ie) = [Qxy(ie,1) Qxy(ie,2) Qxy(ie,3)];
                end

                % Update integration weighting factor
                W = Wgt*Jdet;

                bbarvecT = bbarvecT + W*Gamvec*bvecT;
                H = H + W*(Gamvec*Gamvec');

        %     end %ie
        end %j
        
        % Loop over integration points
        for l = 1:lint

                if nel == 4 || nel == 10
                  [Wgt,ss] =  int3d_t(l,lint,ib);
                  [shl,shld,shls,be] = shltt(ss,nel,nel,der,bf);
                  [Qxy, shgs, Jdet] = shgtt(xl,nel,shld,shls,nel,bf,der,be);
                else
                  [Wgt,ss] =  intpntb(l,lint,ib);
                  [shl,shld,shls,be] = shlb(ss,nel,nel,der,bf);
                  [Qxy, shgs, Jdet] = shgb(xl,nel,shld,shls,nel,bf,der,be);
                end

                % Form B matrix
                for ie = 1:nel
                Bmat(:,(ie-1)*3+1:3*ie) = [Qxy(ie,1) 0        0
                                           0        Qxy(ie,2) 0
                                           0         0        Qxy(ie,3)
                                           Qxy(ie,2) Qxy(ie,1) 0
                                           0         Qxy(ie,3) Qxy(ie,2)
                                           Qxy(ie,3) 0         Qxy(ie,1)];
                bvec(1,(ie-1)*3+1:3*ie) = [Qxy(ie,1) Qxy(ie,2) Qxy(ie,3)];
                end
                
                Bbar = Bmat - 1/3*One*bvecT + 1/3*One*Gamvec'/H*bbarvecT;
                pn_1 = bulk*Gamvec'/H*bbarvecT*reshape(ul,nst,1);
                
                % Compute input for Radial Return
                eps3d = Bbar*reshape(ul,nst,1); %3D enhanced strain
                ephr = nh1-1+(l-1)*11; %pointer for plastic strain at pt l
                betahr = nh1-1+(l-1)*11+5;
                ahr = nh1-1+l*11;
                ep_n = [hr(ephr+1); hr(ephr+2); 0; hr(ephr+3); hr(ephr+4); hr(ephr+5)];
                ep_n(3) = -(ep_n(1) + ep_n(2));
                beta_n = [hr(betahr+1); hr(betahr+2); 0; hr(betahr+3); hr(betahr+4); hr(betahr+5)];
                beta_n(3) = -(beta_n(1) + beta_n(2));
                a_n = hr(ahr);
                
                [sig_n1,Cmat,ep_n1,beta_n1,a_n1] = J2RadialReturn0(eps3d,ep_n,beta_n,a_n,mu,bulk,Khard,Hhard,sigy);
                
                % Convert output from 3D to 2D
                sdev3 = Pdev*sig_n1;
                
                ElemP(1,l) = pn_1 + sdev3(1);
                ElemP(2,l) = pn_1 + sdev3(2);
                ElemP(3,l) = sdev3(4);
                ElemP(4,l) = pn_1 + sdev3(3);
                ElemP(5,l) = a_n1;
                ElemP(6,l) = beta_n1(1);
                ElemP(7,l) = beta_n1(2);
                ElemP(8,l) = beta_n1(4);
                ElemP(9,l) = beta_n1(3);
                ElemP(10,l) = Cmat(1,1);
                ElemP(11,l) = Cmat(2,2);
                ElemP(12,l) = Cmat(4,4);

        end %je

    case -1
        
        ElemF = zeros(nst,1);
        
        % Reorder nodes for corresponding face of integration
        SurfOrientEtoS
        
        t1 = zeros(3,1);
        t2 = t1;
        t3 = t1;

        triagNodes = [1 5 7 8
                      2 6 5 9
                      3 7 6 10
                      5 6 7 9];
        
            lint = 3;
        der = 0;
        bf = 0;
        ib = 0;
        
        % Integration Loop
        for je = 1:lint
            
            for ke = 1:4 % triangle loop
                
            nodeT = triagNodes(ke,1:4);
            xlT = xl(1:3,nodeT);

%             if nel == 4 || nel == 10
              [Wgt,ss] =  int3d_t(je,lint,edge);
              [shl,shld,shls,be] = shltt(ss,4,4,der,bf);
              [QXY,shgs,Jdet,be,sx] = shgtt(xlT,4,shld,shls,4,bf,der,be);
%             else
%               [Wgt,ss] =  intpntb(je,lint,5);
%               [shl,shld,shls,be] = shlb(ss,nel,nel,der,bf);
%               [QXY,shgs,Jdet,be,sx] = shgb(xl,nel,shld,shls,nel,bf,der,be);
%             end
          
            %Evaluate tangent and normal vectors
            t1 = sx(:,2);
            [tm1, tu1] = VecNormalize(t1);
            t2 = sx(:,1);
            [tm2, tu2] = VecNormalize(t2);
            t3 = VecCrossProd(t1,t2);
            [tm3, tu3] = VecNormalize(t3);

            %Force components are positive in positive coord. direction
            c1 = Wgt*tm3;
            for o=1:4

                ie = nodeT(o);
                don = shl(o);
                F = don*traction';

                ElemF(ndf*(ie-1)+1) = ElemF(ndf*(ie-1)+1) + F(1)*c1;

                ElemF(ndf*(ie-1)+2) = ElemF(ndf*(ie-1)+2) + F(2)*c1;

                ElemF(ndf*(ie-1)+3) = ElemF(ndf*(ie-1)+3) + F(3)*c1;

            end %o
            
            end
            
        end %je
        ElemF;
        
        % Reorder nodes back to the orientation of the element
        SurfOrientStoE
        ElemF(1:ndf*nel) = ElemF(ilist2);
        
    case 26
        
        ElemS = zeros(nestr,1);        
        
        bulk = ElemYM/(3*(1-2*Elemv));%Elemv*ElemYM/((1+Elemv)*(1-2*Elemv));
        mu = ElemYM/(2*(1+Elemv));
        One = [1.0d0; 1.0d0; 1.0d0; 0.0d0; 0.0d0; 0.0d0];
        I4 = diag([1.0d0 1.0d0 1.0d0 0.5d0 0.5d0 0.5d0]);
        OneOne = (One*One');

        % Initialize Matrix and Vector
        ElemS2 = zeros(nestr,8);
        
        % Load Guass Integration Points
        lint = 8;

        der = 0;
        bf = 0;
        ib = 0;

        triagNodes = [1 5 7 8
                      2 5 9 6
                      3 7 6 10
                      4 8 10 9
                      5 7 8 9
                      5 7 9 6
                      7 8 9 10
                      7 6 10 9];

%         Bmat = zeros(6,3*nel);
%         bvecT = zeros(1,3*nel);
        bbarvecT = zeros(1,3*nel);
        Gamvec = 1;
        H = 0;
%         Pdev = I4 - 1/3*OneOne;
        
        % Compute bbarvec
        for l = 1:lint

            nodeT = triagNodes(l,1:4);
            xlT = xl(1:3,nodeT);
%                 if nel == 4 || nel == 10
                  [Wgt,ss] =  int3d_t(1,1,ib);
                  [shl,shld,shls,be] = shltt(ss,4,4,der,bf);
                  [Qxy, shgs, Jdet] = shgtt(xlT,4,shld,shls,4,bf,der,be);
%                 else
%                   [Wgt,ss] =  intpntb(l,lint,ib);
%                   [shl,shld,shls,be] = shlb(ss,nel,nel,der,bf);
%                   [Qxy, shgs, Jdet] = shgb(xl,nel,shld,shls,nel,bf,der,be);
%                 end

                % Form b matrix
                bvecT = zeros(1,3*nel);
                for i = 1:4
                ie = nodeT(i);
                bvecT(1,(ie-1)*3+1:3*ie) = [Qxy(i,1) Qxy(i,2) Qxy(i,3)];
                end

                % Update integration weighting factor
                W = Wgt*Jdet;

                bbarvecT = bbarvecT + W*Gamvec*bvecT;
                H = H + W*(Gamvec*Gamvec');

        %     end %ie
        end %j
        
        % Loop over integration points
        for l = 1:lint

            nodeT = triagNodes(l,1:4);
            xlT = xl(1:3,nodeT);
            
%                 if nel == 4 || nel == 10
                  [Wgt,ss] =  int3d_t(1,1,ib);
                  [shl,shld,shls,be] = shltt(ss,4,4,der,bf);
                  [Qxy, shgs, Jdet] = shgtt(xlT,4,shld,shls,4,bf,der,be);
%                 else
%                   [Wgt,ss] =  intpntb(l,lint,ib);
%                   [shl,shld,shls,be] = shlb(ss,nel,nel,der,bf);
%                   [Qxy, shgs, Jdet] = shgb(xl,nel,shld,shls,nel,bf,der,be);
%                 end

                % Form B matrix
                Bmat = zeros(6,3*nel);
                bvecT = zeros(1,3*nel);
                for i = 1:4
                ie = nodeT(i);
                Bmat(:,(ie-1)*3+1:3*ie) = [Qxy(i,1) 0        0
                                           0        Qxy(i,2) 0
                                           0         0        Qxy(i,3)
                                           Qxy(i,2) Qxy(i,1) 0
                                           0         Qxy(i,3) Qxy(i,2)
                                           Qxy(i,3) 0         Qxy(i,1)];
                bvecT(1,(ie-1)*3+1:3*ie) = [Qxy(i,1) Qxy(i,2) Qxy(i,3)];
                end
                
                Bbar = Bmat - 1/3*One*bvecT + 1/3*One*Gamvec'/H*bbarvecT;
                pn_1 = bulk*Gamvec'/H*bbarvecT*reshape(ul,nst,1);
                
                % Compute input for Radial Return
                eps3d = Bbar*reshape(ul,nst,1); %3D enhanced strain
%                 pn_1 = bulk*One'*eps3d;
                ephr = nh1-1+(l-1)*11; %pointer for plastic strain at pt l
                betahr = nh1-1+(l-1)*11+5;
                ahr = nh1-1+l*11;
                ep_n = [hr(ephr+1); hr(ephr+2); 0; hr(ephr+3); hr(ephr+4); hr(ephr+5)];
                ep_n(3) = -(ep_n(1) + ep_n(2));
                beta_n = [hr(betahr+1); hr(betahr+2); 0; hr(betahr+3); hr(betahr+4); hr(betahr+5)];
                beta_n(3) = -(beta_n(1) + beta_n(2));
                a_n = hr(ahr);
                
                [sig_n1,Cmat,ep_n1,beta_n1,a_n1,s_n1,Cdev_n1] = J2RadialReturn0(eps3d,ep_n,beta_n,a_n,mu,bulk,Khard,Hhard,sigy);
                
                % Retrieve deviatoric output
                sdev3 = s_n1;%sig_n1 - 1/3*OneOne*sig_n1;

                sigma = (pn_1*One + sdev3);
            
            for stres = 1:npstr
            
            if stres <= 6 % stress components
                sigmas = sigma(stres);
            elseif stres >= 8
                if stres <= 10 % principal stresses
                sigma2 = [sigma(1) sigma(4) sigma(6); sigma(4) sigma(2) sigma(5); sigma(6) sigma(5) sigma(3)];
                psig = eig(sigma2);
                sigmas = psig(stres-7);
                else % hydrostatic stress
                sigmas = 1/3*sigma'*One;
                end
            else % von Mises stress
                trs = sigma'*One;
                dsig = sigma - 1/3*trs*One;
                sigmas = sqrt(3/2*(dsig'*dsig));
            end
            
            ElemS(stres) = ElemS(stres) + sigmas;
            
            end
%                 ElemK = ElemK + W*Bmat'*(Cmat)*Bmat;

        %     end %ie
        end %je
        
        ElemS = ElemS/8;
        
end %Task Switch