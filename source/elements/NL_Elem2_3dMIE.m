% function ElemK = NL_Elem2_2d(ElemXYZ, nel, ndf, R, S, T, pr, ps, PatchE, Patchv)
%
% Tim Truster
% Tim Truster
% 05/14/2012
% 3-D element for large strain plasticity 
% Bilinear brick element with consistent tangent
% Uses elastic strain as history variable

% Uses deSouza implementation

% Master pure-displacement element, 3D
% Stress post-processing functions and gives meaningful results, not yet
% verified.

% Verified against NL Mixed Inelasticity on 6/30/2013 using PT2Biaxial.m

% Material properties are called as follows:
% mateprop = [mati matv matp E v 0 0 sigy Khard Hhard];
% mati = flag for isochoric elastic material model
% matv = flag for volumetric elastic material model
% matp = flag for plasticity model

switch isw %Task Switch
%%    
    case 1
        
        if ndf > 3
            
            for i = 4:ndf
                lie(i,1) = 0;
            end
            
        end
        
        nh1 = (6+6+1)*8;
        
%%
    case 3
        
        ElemK = zeros(nst);
        ElemF = zeros(nst,1);
        Bmat = zeros(9,3*nel);

        lam = getlam(mateprop); % Modified function numbering to agree with plasticity material models
        One = [1; 1; 1; 0; 0; 0; 0; 0; 0];
        mat1 = One*One';
        matE = diag([2,2,2,1,1,1,0,0,0]);
        dt = s_del_a;
        
        % Load Guass Integration Points

        if nel == 4
            lint = 4;
        elseif nel == 8
            lint = 8;
        elseif nel == 10
            lint = 14;
        else
            lint = 27;
        end
        der = 0;
        bf = 0;
        ib = 0;

        %Integration Loop
        for ll = 1:lint

            if nel == 4 || nel == 10
              [Wgt,ss] =  int3d_t(ll,lint,ib);
              [shl,shld,shls,be] = shltt(ss,nel,nel,der,bf);
              [QXY, shgs, Jdet] = shgtt(xl,nel,shld,shls,nel,bf,der,be);
            else
              [Wgt,ss] =  intpntb(ll,lint,ib);
              [shl,shld,shls,be] = shlb(ss,nel,nel,der,bf);
              [QXY, shgs, Jdet] = shgb(xl,nel,shld,shls,nel,bf,der,be);
            end

            % Push forward kinematics
            [Fn1,JxX,fi,df,Fn,Qxy] = kine3m(QXY,ul,uld,nel,1);
            df = inv(Fn*fi);
                   
            for ie = 1:nel
                Bmat(:,(ie-1)*3+1:3*ie) = [Qxy(ie,1) 0        0
                                           0        Qxy(ie,2) 0
                                           0         0        Qxy(ie,3)
                                           Qxy(ie,2) Qxy(ie,1) 0
                                           0         Qxy(ie,3) Qxy(ie,2)
                                           Qxy(ie,3) 0         Qxy(ie,1)
                                           Qxy(ie,2) -Qxy(ie,1) 0
                                           0         Qxy(ie,3) -Qxy(ie,2)
                                           -Qxy(ie,3) 0         Qxy(ie,1)];
            end

            % Compute input for Radial Return
            ephr = nh1-1+(ll-1)*13; %pointer for plastic strain at pt l
            betahr = nh1-1+(ll-1)*13+6;
            ahr = nh1-1+ll*13;
            ee_n = [hr(ephr+1); hr(ephr+2); hr(ephr+3); hr(ephr+4); hr(ephr+5); hr(ephr+6)];
            beta_n = [hr(betahr+1); hr(betahr+2); hr(betahr+3); hr(betahr+4); hr(betahr+5); hr(betahr+6)];
            a_n = hr(ahr);
                
%             [tau_n1,c_n1,ee_n1,beta_n1,a_n1,s_n1,Cdev_n1] = logplasfd(ee_n,df,beta_n,a_n,mu,bulk,K,H,sigy);
            logplasfd

            % Store history variables
            ephr = nh2-1+(ll-1)*13; %pointer for plastic strain at pt l
            betahr = nh2-1+(ll-1)*13+6;
            ahr = nh2-1+ll*13;
            hr(ephr+1) = ee_n1(1);
            hr(ephr+2) = ee_n1(2);
            hr(ephr+3) = ee_n1(3);
            hr(ephr+4) = ee_n1(4);
            hr(ephr+5) = ee_n1(5);
            hr(ephr+6) = ee_n1(6);
            hr(betahr+1) = beta_n1(1);
            hr(betahr+2) = beta_n1(2);
            hr(betahr+3) = beta_n1(3);
            hr(betahr+4) = beta_n1(4);
            hr(betahr+5) = beta_n1(5);
            hr(betahr+6) = beta_n1(6);
            hr(ahr) = a_n1;
            
            % Append volumetric stress and material moduli
            [theta1,theta2] = ThetaNS(JxX(1),mateprop); % Modified function numbering to agree with plasticity material models
            sigmap = lam*theta1*JxX(1)*One;
            cpmat = lam*((theta2*JxX(1)^2 + theta1*JxX(1))*mat1 - theta1*JxX(1)*matE);
            sigma = s_n1 + sigmap;
            cmat = cdev_n1 + cpmat;

            sigma2 = sigma;
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

            % Update integration weighting factor
            W = Wgt*Jdet;

            ElemF = ElemF - W*Bmat'*(sigma);
            ElemK = ElemK + W*Bmat'*(cmat + Smat)*Bmat;
            
        end %je
   ElemK;     
%%
    case -1

        ElemF = zeros(nst,1);
        
        % Reorder nodes for corresponding face of integration
        SurfOrientEtoS
        
        t1 = zeros(3,1);
        t2 = t1;
        t3 = t1;
        
        if nel == 4 || nel == 10
            lint = 13;
        else
            lint = 16;
        end
        der = 0;
        bf = 0;
        ib = 0;
        
        % Integration Loop
        for je = 1:lint

            if nel == 4 || nel == 10
              [Wgt,ss] =  int3d_t(je,lint,edge);
              [shl,shld,shls,be] = shltt(ss,nel,nen,der,bf);
              [QXY,shgs,Jdet,be,sx] = shgtt(xl,nel,shld,shls,nel,bf,der,be);
            else
              [Wgt,ss] =  intpntb(je,lint,5);
              [shl,shld,shls,be] = shlb(ss,nel,nen,der,bf);
              [QXY,shgs,Jdet,be,sx] = shgb(xl,nel,shld,shls,nel,bf,der,be);
            end
          
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

                ElemF(ndf*o-2) = ElemF(ndf*o-2) + F(1)*c1;

                ElemF(ndf*o-1) = ElemF(ndf*o-1) + F(2)*c1;

                ElemF(ndf*o-0) = ElemF(ndf*o-0) + F(3)*c1;

            end %o
            
        end %je
        
        % Reorder nodes back to the orientation of the element
        SurfOrientStoE
        ElemF(1:ndf*nel) = ElemF(ilist2);
%%
    case 6
        
        ElemF = zeros(nst,1);
        Bmat = zeros(9,3*nel);

        lam = getlam(mateprop);
        dt = s_del_a;
        
        % Load Guass Integration Points

        if nel == 4
            lint = 11;5;16;
        elseif nel == 8
            lint = 8;
        elseif nel == 10
            lint = 14;
%             lint = 27;
        else
            lint = 27;
        end
        der = 0;
        bf = 0;
        ib = 0;

        %Integration Loop
        for ll = 1:lint

            if nel == 4 || nel == 10
              [Wgt,ss] =  int3d_t(ll,lint,ib);
              [shl,shld,shls,be] = shltt(ss,nel,nel,der,bf);
              [QXY, shgs, Jdet] = shgtt(xl,nel,shld,shls,nel,bf,der,be);
            else
              [Wgt,ss] =  intpntb(ll,lint,ib);
              [shl,shld,shls,be] = shlb(ss,nel,nel,der,bf);
              [QXY, shgs, Jdet] = shgb(xl,nel,shld,shls,nel,bf,der,be);
            end

            % Push forward kinematics
            [Fn1,JxX,fi,df,Fn,Qxy] = kine3m(QXY,ul,uld,nel,1);
            df = inv(Fn*fi);
                   
            for ie = 1:nel
                Bmat(:,(ie-1)*3+1:3*ie) = [Qxy(ie,1) 0        0
                                           0        Qxy(ie,2) 0
                                           0         0        Qxy(ie,3)
                                           Qxy(ie,2) Qxy(ie,1) 0
                                           0         Qxy(ie,3) Qxy(ie,2)
                                           Qxy(ie,3) 0         Qxy(ie,1)
                                           Qxy(ie,2) -Qxy(ie,1) 0
                                           0         Qxy(ie,3) -Qxy(ie,2)
                                           -Qxy(ie,3) 0         Qxy(ie,1)];
            end

            % Compute input for Radial Return
            ephr = nh1-1+(ll-1)*13; %pointer for plastic strain at pt l
            betahr = nh1-1+(ll-1)*13+6;
            ahr = nh1-1+ll*13;
            ee_n = [hr(ephr+1); hr(ephr+2); hr(ephr+3); hr(ephr+4); hr(ephr+5); hr(ephr+6)];
            beta_n = [hr(betahr+1); hr(betahr+2); hr(betahr+3); hr(betahr+4); hr(betahr+5); hr(betahr+6)];
            a_n = hr(ahr);
                
%             [tau_n1,c_n1,ee_n1,beta_n1,a_n1,s_n1,Cdev_n1] = logplasfd(ee_n,df,beta_n,a_n,mu,bulk,K,H,sigy);
            logplasfd

            % Store history variables
            ephr = nh2-1+(ll-1)*13; %pointer for plastic strain at pt l
            betahr = nh2-1+(ll-1)*13+6;
            ahr = nh2-1+ll*13;
            hr(ephr+1) = ee_n1(1);
            hr(ephr+2) = ee_n1(2);
            hr(ephr+3) = ee_n1(3);
            hr(ephr+4) = ee_n1(4);
            hr(ephr+5) = ee_n1(5);
            hr(ephr+6) = ee_n1(6);
            hr(betahr+1) = beta_n1(1);
            hr(betahr+2) = beta_n1(2);
            hr(betahr+3) = beta_n1(3);
            hr(betahr+4) = beta_n1(4);
            hr(betahr+5) = beta_n1(5);
            hr(betahr+6) = beta_n1(6);
            hr(ahr) = a_n1;
            
            % Append volumetric stress and material moduli
            [theta1,theta2] = ThetaNS(JxX(1),mateprop);
            sigmap = lam*theta1*JxX(1)*One;
            sigma = s_n1 + sigmap;

            % Update integration weighting factor
            W = Wgt*Jdet;

            ElemF = ElemF - W*Bmat'*(sigma);
            
        end %je
%%  
    case 22
        
        ElemM = zeros(nst);
        ElemF = zeros(nst,1);

        lam = getlam(mateprop);
        dt = s_del_a;
        
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
              [shl,shld,shls,be] = shltt(ss,nel,nen,der,bf);
              [QXY, shgs, Jdet] = shgtt(xl,nel,shld,shls,nel,bf,der,be);
            else
              [Wgt,ss] =  intpntb(ll,lint,ib);
              [shl,shld,shls,be] = shlb(ss,nel,nen,der,bf);
              [QXY, shgs, Jdet] = shgb(xl,nel,shld,shls,nel,bf,der,be);
            end

            if nelP == 4 || nelP == 10
              [shlS,shld,shls] = shltt(ss,nelS,nel,0,0);
            else
              [shlS,shld,shls] = shlb(ss,nelS,nel,0,0);
            end
            
            % Push forward kinematics
            [Fn1,JxX,fi,df,Fn,Qxy] = kine3m(QXY,ul,uld,nel,1);
            df = inv(Fn*fi);

            % Compute input for Radial Return
            ephr = nh1-1+(ll-1)*13; %pointer for plastic strain at pt l
            betahr = nh1-1+(ll-1)*13+6;
            ahr = nh1-1+ll*13;
            be_n = [hr(ephr+1); hr(ephr+2); hr(ephr+3); hr(ephr+4); hr(ephr+5); hr(ephr+6)];
            beta_n = [hr(betahr+1); hr(betahr+2); hr(betahr+3); hr(betahr+4); hr(betahr+5); hr(betahr+6)];
            a_n = hr(ahr);

%             [tau_n1,c_n1,ee_n1,beta_n1,a_n1,s_n1,Cdev_n1] = logplasfd(ee_n,df,beta_n,a_n,mu,bulk,K,H,sigy);
            logplasfd

            % Append volumetric stress and material moduli
            [theta1,theta2] = ThetaNS(JxX(1),mateprop);
            sigmap = lam*theta1*JxX(1)*One;
            sigma = s_n1 + sigmap;
            sigma = sigma/JxX(1);
            
            % Form B matrix
            Nmat = shlS';

            w = Wgt*Jdet;
            
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
        dt = s_del_a;
        
        I1 = [1; 1; 1; 0; 0; 0; 0; 0; 0];
        spvec0 = I1;
        
        thick = 1;
        
        % Load Guass Integration Points

        if nel == 4
            lint = 1;
            nint = 1;
        elseif nel == 8
%             lint = 4;
            lint = 8;
            nint = 1;
        elseif nel == 10
            lint = 11;
            nint = 4;
        else
            lint = 27;
            nint = 8;
        end
        
        der = 0;
        bf = 0;
        ib = 0;

        %Stress Loop
        for ll = 1:nint

            %Evaluate first derivatives of basis functions at int. point
            if nel == 4 || nel == 10
              [Wgt,ss] =  int3d_t(ll,nint,ib);
              [shl,shld,shls,be] = shltt(ss,nel,nel,der,bf);
              [QXY,shgs,Jdet,be] = shgtt(xl,nel,shld,shls,nel,bf,der,be);
            else
              [Wgt,ss] =  intpntb(ll,nint,ib);
              [shl,shld,shls,be] = shlb(ss,nel,nel,der,bf);
              [QXY,shgs,Jdet,be] = shgb(xl,nel,shld,shls,nel,bf,der,be);
            end
            
            if nelP == 4 || nelP == 10
              shlp = shltt(ss,nelP,nel,0,0);
            else
              shlp = shlb(ss,nelP,nel,0,0);
            end

            % Push forward kinematics
            [Fn1,JxX,fi,df,Fn,Qxy] = kine3m(QXY,ul,uld,nel,1);
            df = inv(Fn*fi);
            
            % Compute input for Radial Return
            ephr = nh1-1+(ll-1)*13; %pointer for plastic strain at pt l
            betahr = nh1-1+(ll-1)*13+6;
            ahr = nh1-1+ll*13;
            be_n = [hr(ephr+1); hr(ephr+2); hr(ephr+3); hr(ephr+4); hr(ephr+5); hr(ephr+6)];
            beta_n = [hr(betahr+1); hr(betahr+2); hr(betahr+3); hr(betahr+4); hr(betahr+5); hr(betahr+6)];
            a_n = hr(ahr);

%             [tau_n1,c_n1,ee_n1,beta_n1,a_n1,s_n1,Cdev_n1] = logplasfd(ee_n,df,beta_n,a_n,mu,bulk,K,H,sigy);
            logplasfd

            % Append volumetric stress and material moduli
            [theta1,theta2] = ThetaNS(JxX(1),mateprop);
            sigmap = lam*theta1*JxX(1)*One;
            sigma = s_n1 + sigmap;
            sigma = sigma/JxX(1);
            
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
        if nel == 4
            plist = [1 0 0 0
                     0 1 0 0
                     0 0 0 1];
        elseif nel == 8
            plist = [-1 1 1 -1 -1 1 1 -1
                     -1 -1 1 1 -1 -1 1 1
                     -1 -1 -1 -1 1 1 1 1];
        elseif nel == 10
            plist = [ 1.927050983124842 -0.309016994374947 -0.309016994374947 -0.309016994374947  0.809016994374947 -0.309016994374947  0.809016994374947  0.809016994374947 -0.309016994374947 -0.309016994374947
                     -0.309016994374947  1.927050983124842 -0.309016994374947 -0.309016994374947  0.809016994374947  0.809016994374947 -0.309016994374947 -0.309016994374947  0.809016994374947 -0.309016994374947
                     -0.309016994374947 -0.309016994374947 -0.309016994374947  1.927050983124842 -0.309016994374947 -0.309016994374947 -0.309016994374947  0.809016994374947  0.809016994374947  0.809016994374947];
        else
            sqr3 = sqrt(3);
            plist = [-sqr3 sqr3 sqr3 -sqr3 -sqr3 sqr3 sqr3 -sqr3 0 sqr3 0 -sqr3 0 sqr3 0 -sqr3 -sqr3 sqr3 sqr3 -sqr3 0 0 -sqr3 sqr3 0 0 0
                     -sqr3 -sqr3 sqr3 sqr3 -sqr3 -sqr3 sqr3 sqr3 -sqr3 0 sqr3 0 -sqr3 0 sqr3 0 -sqr3 -sqr3 sqr3 sqr3 0 0 0 0 -sqr3 sqr3 0
                     -sqr3 -sqr3 -sqr3 -sqr3 sqr3 sqr3 sqr3 sqr3 -sqr3 -sqr3 -sqr3 -sqr3 sqr3 sqr3 sqr3 sqr3 0 0 0 0 -sqr3 sqr3 0 0 0 0 0];
        end
        
        for ll = 1:nelS
            
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

    case 24
        
        ElemP = zeros(12,nel);
                
%         One = [1.0d0; 1.0d0; 1.0d0; 0.0d0; 0.0d0; 0.0d0];
%         I4 = diag([1.0d0 1.0d0 1.0d0 0.5d0 0.5d0 0.5d0]);
%         OneOne = (One*One');

        % Initialize Matrix and Vector

        nst = nel*ndf;
        ElemK = zeros(nst);
        ElemF = zeros(nst,1);
        
        % Load Guass Integration Points
        lint = 8;

        der = 0;
        bf = 0;
        ib = 0;

        Bmat = zeros(9,3*nel);
%         bvecT = zeros(1,3*nel);
%         bbarvecT = zeros(1,3*nel);
%         Gamvec = 1;
%         H = 0;
%         Pdev = I4 - 1/3*OneOne;
        
%         % Compute bbarvec
%         for l = 1:lint
% 
%                 if nel == 4 || nel == 10
%                   [Wgt,ss] =  int3d_t(l,lint,ib);
%                   [shl,shld,shls,be] = shltt(ss,nel,der,bf);
%                   [Qxy, shgs, Jdet] = shgtt(xl,nel,shld,shls,nel,bf,der,be);
%                 else
%                   [Wgt,ss] =  intpntb(l,lint,ib);
%                   [shl,shld,shls,be] = shlb(ss,nel,der,bf);
%                   [Qxy, shgs, Jdet] = shgb(xl,nel,shld,shls,nel,bf,der,be);
%                 end
% 
%                 % Form b matrix
%                 for ie = 1:nel
%                 bvecT(1,(ie-1)*3+1:3*ie) = [Qxy(ie,1) Qxy(ie,2) Qxy(ie,3)];
%                 end
% 
%                 % Update integration weighting factor
%                 W = Wgt*Jdet;
% 
%                 bbarvecT = bbarvecT + W*Gamvec*bvecT;
%                 H = H + W*(Gamvec*Gamvec');
% 
%         %     end %ie
%         end %j
        
        % Loop over integration points
        for ll = 1:lint

                if nel == 4 || nel == 10
                  [Wgt,ss] =  int3d_t(ll,lint,ib);
                  [shl,shld,shls,be] = shltt(ss,nel,nel,der,bf);
                  [QXY, shgs, Jdet] = shgtt(xl,nel,shld,shls,nel,bf,der,be);
                else
                  [Wgt,ss] =  intpntb(ll,lint,ib);
                  [shl,shld,shls,be] = shlb(ss,nel,nel,der,bf);
                  [QXY, shgs, Jdet] = shgb(xl,nel,shld,shls,nel,bf,der,be);
                end

                % Push forward kinematics
                [Fn1,JxX,fi,df,Fn,Qxy] = kine3m(QXY,ul,uld,nel,1);
                df = inv(Fn*fi);
                
                % Form B matrix
                for ie = 1:nel
                Bmat(:,(ie-1)*3+1:3*ie) = [Qxy(ie,1) 0        0
                                           0        Qxy(ie,2) 0
                                           0         0        Qxy(ie,3)
                                           Qxy(ie,2) Qxy(ie,1) 0
                                           0         Qxy(ie,3) Qxy(ie,2)
                                           Qxy(ie,3) 0         Qxy(ie,1)
                                           Qxy(ie,2) -Qxy(ie,1) 0
                                           0         Qxy(ie,3) -Qxy(ie,2)
                                           -Qxy(ie,3) 0         Qxy(ie,1)];
%                 bvec(1,(ie-1)*3+1:3*ie) = [Qxy(ie,1) Qxy(ie,2) Qxy(ie,3)];
                end
                
%                 Bbar = Bmat - 1/3*One*bvecT + 1/3*One*Gamvec'/H*bbarvecT;
%                 pn_1 = bulk*Gamvec'/H*bbarvecT*reshape(ul,nst,1);
                
                % Compute input for Radial Return
                ephr = nh1-1+(ll-1)*13; %pointer for plastic strain at pt l
                betahr = nh1-1+(ll-1)*13+6;
                ahr = nh1-1+ll*13;
                be_n = [hr(ephr+1); hr(ephr+2); hr(ephr+3); hr(ephr+4); hr(ephr+5); hr(ephr+6)];
                beta_n = [hr(betahr+1); hr(betahr+2); hr(betahr+3); hr(betahr+4); hr(betahr+5); hr(betahr+6)];
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
                
                ElemP(1,ll) = sigma(1);
                ElemP(2,ll) = sigma(2);
                ElemP(3,ll) = sigma(4);
                ElemP(4,ll) = sigma(3);
                ElemP(5,ll) = a_n1;
                ElemP(6,ll) = beta_n1(1);
                ElemP(7,ll) = beta_n1(2);
                ElemP(8,ll) = beta_n1(4);
                ElemP(9,ll) = beta_n1(3);
                ElemP(10,ll) = cmat(1,1);
                ElemP(11,ll) = cmat(2,2);
                ElemP(12,ll) = cmat(4,4);

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