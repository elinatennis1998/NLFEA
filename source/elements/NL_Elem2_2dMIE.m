
% Tim Truster
% 01/12/2014
% 2-D element for large strain plasticity 
% Bilinear quadrilateral element with consistent tangent
% Uses elastic strain as history variable

% Uses deSouza implementation

% Tested using IFbar2Dtest.m for general loading; viscosity has not been
% tested yet, but isotropic and kinematic hardening have

% Material properties are called as follows:
% mateprop = [mati matv matp E v 0 0 sigy Khard Hhard];
% mati = flag for isochoric elastic material model
% matv = flag for volumetric elastic material model
% matp = flag for plasticity model

switch isw %Task Switch
%%    
    case 1
        
        if ndf > 2
            
            for i = 3:ndf
                lie(i,1) = 0;
            end
            
        end
        
        nh1 = (8)*4;
        
%%
    case 3
        
        ElemK = zeros(nst);
        ElemF = zeros(nst,1);
        Bmat = zeros(4,2*nel);

        lam = getlam(mateprop); % Modified function numbering to agree with plasticity material models
        One = [1; 1; 1; 0; 0; 0; 0; 0; 0];
        mat1 = One*One';
        matE = diag([2,2,2,1,1,1,0,0,0]);
        dt = s_del_a;
        
        ind3to2 = [1 2 4 7];
        
        
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
        cpmat2 = diag([2,2,1,0]);
        cpmat0 = cpmat1 + cpmat2;
        
        % Load Guass Integration Points

            lint = 4;
        der = 0;
        bf = 0;
        ib = 0;


        %Integration Loop
        for l = 1:lint

                  [Wgt,litr,lits] =  intpntq(l,lint,ib);
                  [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
                  [QXY, shgs, Jdet] = shgq(xl,nel,shld,shls,nel,bf,der,be);

            % Push forward kinematics
            [Fn1,JxX,fi,df,Fn,Qxy] = kine2m(QXY,ul,uld,nel,1);
            df = inv(Fn*fi);
            JxXi = det(df);
            
                Fn1 = [Fn1 zeros(2,1); zeros(1,2) 1]; %#ok<AGROW>
                Fn = [Fn zeros(2,1); zeros(1,2) 1]; %#ok<AGROW>
                df = [df zeros(2,1); zeros(1,2) 1]; %#ok<AGROW>
        
            for mm = 1:nel    
 Bmat(:,2*mm-1:2*mm) = [Qxy(mm,1) 0        
                        0         Qxy(mm,2)
                        Qxy(mm,2) Qxy(mm,1)
                        Qxy(mm,2) -Qxy(mm,1)];
            end

            % Compute input for Radial Return
            ephr = nh1-1+(l-1)*8; %pointer for plastic strain at pt l
            betahr = nh1-1+(l-1)*8+4;
            ahr = nh1-1+l*8;
            ee_n = [hr(ephr+1); hr(ephr+2); hr(ephr+3); hr(ephr+4); 0; 0];
            beta_n = [hr(betahr+1); hr(betahr+2); 0; hr(betahr+3); 0; 0];
            beta_n(3) = -(beta_n(1) + beta_n(2));
            a_n = hr(ahr);
                
            logplasfd
            
            % Store history variables
            ephr = nh2-1+(l-1)*8; %pointer for plastic strain at pt l
            betahr = nh2-1+(l-1)*8+4;
            ahr = nh2-1+l*8;
            hr(ephr+1) = ee_n1(1);
            hr(ephr+2) = ee_n1(2);
            hr(ephr+3) = ee_n1(3);
            hr(ephr+4) = ee_n1(4);
            hr(betahr+1) = beta_n1(1);
            hr(betahr+2) = beta_n1(2);
            hr(betahr+3) = beta_n1(4);
            hr(ahr) = a_n1;
            
            sdev = s_n1(ind3to2);
            Cmatdev = cdev_n1(ind3to2,ind3to2);
            
            % Append volumetric stress and material moduli
            [theta1,theta2] = ThetaNS(JxX(1),mateprop); % Modified function numbering to agree with plasticity material models
            sigmap = lam*theta1*JxX(1)*spvec0;
            cpmat = lam*((theta2*JxX(1)^2 + theta1*JxX(1))*cpmat1 - theta1*JxX(1)*cpmat2);
            sigma = sdev + sigmap;
            cmat = Cmatdev + cpmat;

            sigma2 = sigma;
            Smat = [sigma2(1) 0  sigma2(3)/2 sigma2(3)/2
                    0 sigma2(2)  sigma2(3)/2 -sigma2(3)/2
                    sigma2(3)/2  sigma2(3)/2 (sigma2(2)+sigma2(1))/4 (sigma2(2)-sigma2(1))/4
                    sigma2(3)/2 -sigma2(3)/2 (sigma2(2)-sigma2(1))/4 (sigma2(2)+sigma2(1))/4];

            % Update integration weighting factor
            W = Wgt*Jdet;

            ElemF = ElemF - W*Bmat'*(sigma);
            ElemK = ElemK + W*Bmat'*(cmat + Smat)*Bmat;
            
        end %je
   ElemK;     
%%
    case 6
        
        
        ElemK = zeros(nst);
        ElemF = zeros(nst,1);
        Bmat = zeros(4,2*nel);

        lam = getlam(mateprop); % Modified function numbering to agree with plasticity material models
        One = [1; 1; 1; 0; 0; 0; 0; 0; 0];
        mat1 = One*One';
        matE = diag([2,2,2,1,1,1,0,0,0]);
        dt = s_del_a;
        
        ind3to2 = [1 2 4 7];
        
        
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
        
        % Load Guass Integration Points

            lint = 4;
        der = 0;
        bf = 0;
        ib = 0;


        %Integration Loop
        for ll = 1:lint

                  [Wgt,litr,lits] =  intpntq(ll,lint,ib);
                  [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
                  [QXY, shgs, Jdet] = shgq(xl,nel,shld,shls,nel,bf,der,be);

            % Push forward kinematics
            [Fn1,JxX,fi,df,Fn,Qxy] = kine2m(QXY,ul,uld,nel,1);
            df = inv(Fn*fi);
            JxXi = det(df);
            
                Fn1 = [Fn1 zeros(2,1); zeros(1,2) 1]; %#ok<AGROW>
                Fn = [Fn zeros(2,1); zeros(1,2) 1]; %#ok<AGROW>
                df = [df zeros(2,1); zeros(1,2) 1]; %#ok<AGROW>
        
            for mm = 1:nel    
 Bmat(:,2*mm-1:2*mm) = [Qxy(mm,1) 0        
                        0         Qxy(mm,2)
                        Qxy(mm,2) Qxy(mm,1)
                        Qxy(mm,2) -Qxy(mm,1)];
            end

            % Compute input for Radial Return
            ephr = nh1-1+(l-1)*8; %pointer for plastic strain at pt l
            betahr = nh1-1+(l-1)*8+4;
            ahr = nh1-1+l*8;
            ee_n = [hr(ephr+1); hr(ephr+2); hr(ephr+3); hr(ephr+4); 0; 0];
            beta_n = [hr(betahr+1); hr(betahr+2); 0; hr(betahr+3); 0; 0];
            beta_n(3) = -(beta_n(1) + beta_n(2));
            a_n = hr(ahr);
                
            logplasfd2
            
            % Store history variables
            ephr = nh2-1+(l-1)*8; %pointer for plastic strain at pt l
            betahr = nh2-1+(l-1)*8+4;
            ahr = nh2-1+l*8;
            hr(ephr+1) = ee_n1(1);
            hr(ephr+2) = ee_n1(2);
            hr(ephr+3) = ee_n1(3);
            hr(ephr+4) = ee_n1(4);
            hr(betahr+1) = beta_n1(1);
            hr(betahr+2) = beta_n1(2);
            hr(betahr+3) = beta_n1(4);
            hr(ahr) = a_n1;
            
            sdev = s_n1(ind3to2);
            
            % Append volumetric stress and material moduli
            [theta1,theta2] = ThetaNS(JxX(1),mateprop); % Modified function numbering to agree with plasticity material models
            sigmap = lam*theta1*JxX(1)*spvec0;
            sigma = sdev + sigmap;

            % Update integration weighting factor
            W = Wgt*Jdet;

            ElemF = ElemF - W*Bmat'*(sigma);
            
        end %je
   ElemK;     
%%  
    case 24
        
        ElemP = zeros(12,4);
        
        % Initialize Matrix and Vector

        nst = nel*ndf;
        ElemK = zeros(nst);
        ElemF = zeros(nst,1);
        
        % Load Guass Integration Points
        lint = 4;

        der = 0;
        bf = 0;
        ib = 0;
        
        % Loop over integration points
        for l = 1:lint

                  [Wgt,litr,lits] =  intpntq(l,lint,ib);
                  [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
                  [QXY, shgs, Jdet] = shgq(xl,nel,shld,shls,nel,bf,der,be);

            % Push forward kinematics
            [Fn1,JxX,fi,df,Fn,Qxy] = kine2m(QXY,ul,uld,nel,1);
            df = inv(Fn*fi);
            JxXi = det(df);
            
                Fn1 = [Fn1 zeros(2,1); zeros(1,2) 1]; %#ok<AGROW>
                Fn = [Fn zeros(2,1); zeros(1,2) 1]; %#ok<AGROW>
                df = [df zeros(2,1); zeros(1,2) 1]; %#ok<AGROW>
                
            % Compute input for Radial Return
            ephr = nh1-1+(l-1)*8; %pointer for plastic strain at pt l
            betahr = nh1-1+(l-1)*8+4;
            ahr = nh1-1+l*8;
            ee_n = [hr(ephr+1); hr(ephr+2); hr(ephr+3); hr(ephr+4); 0; 0];
            beta_n = [hr(betahr+1); hr(betahr+2); 0; hr(betahr+3); 0; 0];
            beta_n(3) = -(beta_n(1) + beta_n(2));
            a_n = hr(ahr);
                
%             [tau_n1,c_n1,ee_n1,beta_n1,a_n1,s_n1,Cdev_n1] = logplasfd(ee_n,df,beta_n,a_n,mu,bulk,K,H,sigy);
                logplasfd

                % Append volumetric stress and material moduli
            [theta1,theta2] = ThetaNS(JxX(1),mateprop);
                sigmap = lam*theta1*JxX(1)*One;
                cpmat = lam*((theta2*JxX(1)^2 + theta1*JxX(1))*mat1 - theta1*JxX(1)*matE);
                sigma = s_n1 + sigmap;
                cmat = cdev_n1 + cpmat;
                
                % Retrieve deviatoric output
                
                ElemP(1,l) = sigma(1);
                ElemP(2,l) = sigma(2);
                ElemP(3,l) = sigma(4);
                ElemP(4,l) = sigma(3);
                ElemP(5,l) = a_n1;
                ElemP(6,l) = beta_n1(1);
                ElemP(7,l) = beta_n1(2);
                ElemP(8,l) = beta_n1(4);
                ElemP(9,l) = beta_n1(3);
                ElemP(10,l) = cmat(1,1);
                ElemP(11,l) = cmat(2,2);
                ElemP(12,l) = cmat(4,4);

        end %je
                
        
    case 40 % Initialize Intermediate Configuration
        
        lint = 8;
        
        % Loop over integration points
        for l = 1:lint
                
                % Store history variables
                ephr = nh1-1+(l-1)*13; %pointer for plastic strain at pt l
                betahr = nh1-1+(l-1)*13+6;
                ahr = nh1-1+l*13;
%                 hr(ephr+1) = 1.d0;
%                 hr(ephr+2) = 1.d0;
%                 hr(ephr+3) = 1.d0;

        end %je
end