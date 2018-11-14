
% Tim Truster
% 2/11/2014
% UIUC
% Generalized alpha method dynamic element, 3-D


% Set Material Properties

PatchE = mateprop(1);
Patchv = mateprop(2);
rho0 = mateprop(4);

Bcol1 = [1; 4; 6];
Bcol2 = [2; 4; 5];
Bcol3 = [3; 5; 6];
col1 = [1; 2; 3];
col2 = [2; 1; 3];
col3 = [3; 2; 1];
iemat = ieFEAP(1:ndf,ma);

switch isw %Task Switch
%%    
    case 1
        
        if ndf > 3
            
            for i = 4:ndf
                lie(i,1) = 0;
            end
            
        end
        
%%
    case 3 %Compute Stiffness and Residual
        
        ElemK = zeros(nst);
        ElemF = zeros(nst,1);
        ElemFn = zeros(nst,1);

        %Set integration number
        if nel ~= 6
        lint = IntPoint3(nel);
        else
            linttw = 7;
            lintlw = 2;
            lint = linttw*lintlw;
        end
        ib = 0;
        bf = 0;
        der = 0;

        lam = Patchv*PatchE/((1+Patchv)*(1-2*Patchv));
        mu = PatchE/(2*(1+Patchv));
        fbx = 0;
        fby = 0;
        fbz = 0;
        Dmat = mu*diag([2 2 2 1 1 1]) + lam*[1; 1; 1; 0; 0; 0]*[1 1 1 0 0 0];
        Nmat = zeros(3,nst);
        Bmat = zeros(6,nst);

        % Compute kinematic fields at intermediate time levels
        ul_n_af = (1 - Nalphaf)*ul + Nalphaf*ul_n;
        al_n_am = (1 - Nalpham)*al + Nalpham*al_n;
        ulres = reshape(ul_n_af,ndf*nen,1);
        alres = reshape(al_n_am,ndf*nen,1);
        
        for l = 1:lint

            % Evaluate 1-D basis functions at integration points
            if nel == 4 || nel == 10
              [w,ss] =  int3d_t(l,lint,ib);
              [shl,shld,shls,be] = shltt(ss,nel,nel,der,bf);
              [shg,shgs,Jdet,be] = shgtt(xl,nel,shld,shls,nel,bf,der,be);
            elseif nel == 6 % wedge
              [w,ss] =  intpntw(l,linttw,lintlw,ib);
              [shl,shld,shls,be] = shlw(ss,nel,nel,der,bf);
              [shg,shgs,Jdet,be,xs] = shgw(xl,nel,shld,shls,nel,bf,der,be);
            else
              [w,ss] =  intpntb(l,lint,ib);
              [shl,shld,shls,be] = shlb(ss,nel,nel,der,bf);
              [shg,shgs,Jdet,be] = shgb(xl,nel,shld,shls,nel,bf,der,be);
            end
            
            if Jdet < 0
                elem
                Jdet %#ok<NOPTS>
            end
            c1 = Jdet*w;
            
            Fvec = [fbx; fby; fbz];
            
            % Form B matrix
            for i = 1:nel
                
              Nmat(1,(i-1)*ndf+iemat(1)) = shl(i);
              Nmat(2,(i-1)*ndf+iemat(2)) = shl(i);
              Nmat(3,(i-1)*ndf+iemat(3)) = shl(i);
                
              Bmat(Bcol1,(i-1)*ndf+iemat(1)) = shg(i,col1);
              Bmat(Bcol2,(i-1)*ndf+iemat(2)) = shg(i,col2);
              Bmat(Bcol3,(i-1)*ndf+iemat(3)) = shg(i,col3);
                 
            end

            %Compute rho, accel
            rho = rho0; %rho = f(x,y);
            accel = al_n_am(1:ndf,1:nel)*shl;
            sigma = Dmat*(Bmat*ulres);
            
            ElemF = ElemF + c1*(Nmat'*Fvec - rho*Nmat'*accel - Bmat'*sigma);
            
            ElemK = ElemK + c1*(Bmat'*Dmat*Bmat);

        end %je
        
    case 5 %Compute Mass
        
        ElemM = zeros(nst);
        
        % Load Guass Integration Points
        if nel ~= 6
        lint = IntPoint3(nel);
        else
            linttw = 7;
            lintlw = 2;
            lint = linttw*lintlw;
        end
        der = 0;
        bf = 0;
        ib = 0;
        Nmat = zeros(3,nst);
        
        for l = 1:lint                    

            % Evaluate 1-D basis functions at integration points
            if nel == 4 || nel == 10
              [w,ss] =  int3d_t(l,lint,ib);
              [shl,shld,shls,be] = shltt(ss,nel,nel,der,bf);
              [shg,shgs,Jdet,be] = shgtt(xl,nel,shld,shls,nel,bf,der,be);
            elseif nel == 6 % wedge
              [w,ss] =  intpntw(l,linttw,lintlw,ib);
              [shl,shld,shls,be] = shlw(ss,nel,nel,der,bf);
              [shg,shgs,Jdet,be,xs] = shgw(xl,nel,shld,shls,nel,bf,der,be);
            else
              [w,ss] =  intpntb(l,lint,ib);
              [shl,shld,shls,be] = shlb(ss,nel,nel,der,bf);
              [shg,shgs,Jdet,be] = shgb(xl,nel,shld,shls,nel,bf,der,be);
            end
            
            if Jdet < 0
                elem
                Jdet %#ok<NOPTS>
            end
            c1 = Jdet*w;
            
            % Form B matrix
            for i = 1:nel
                
              Nmat(1,(i-1)*ndf+iemat(1)) = shl(i);
              Nmat(2,(i-1)*ndf+iemat(2)) = shl(i);
              Nmat(3,(i-1)*ndf+iemat(3)) = shl(i);
                 
            end

            %Compute rho, accel
            rho = rho0; %rho = f(x,y);
            
            ElemM = ElemM + c1*rho*(Nmat'*Nmat);

        end %je

    case 6 %Compute Residual
        
        ElemF = zeros(nst,1);
        ElemFn = zeros(nst,1);

        %Set integration number
        if nel ~= 6
        lint = IntPoint3(nel);
        else
            linttw = 7;
            lintlw = 2;
            lint = linttw*lintlw;
        end
        ib = 0;
        bf = 0;
        der = 0;

        lam = Patchv*PatchE/((1+Patchv)*(1-2*Patchv));
        mu = PatchE/(2*(1+Patchv));
        fbx = 0;
        fby = 0;
        fbz = 0;
        Dmat = mu*diag([2 2 2 1 1 1]) + lam*[1; 1; 1; 0; 0; 0]*[1 1 1 0 0 0];
        Nmat = zeros(3,nst);
        Bmat = zeros(6,nst);

        % Compute kinematic fields at intermediate time levels
        ul_n_af = (1 - Nalphaf)*ul + Nalphaf*ul_n;
        al_n_am = (1 - Nalpham)*al + Nalpham*al_n;
        ulres = reshape(ul_n_af,ndf*nen,1);
        alres = reshape(al_n_am,ndf*nen,1);
        
        for l = 1:lint

            % Evaluate 1-D basis functions at integration points
            if nel == 4 || nel == 10
              [w,ss] =  int3d_t(l,lint,ib);
              [shl,shld,shls,be] = shltt(ss,nel,nel,der,bf);
              [shg,shgs,Jdet,be] = shgtt(xl,nel,shld,shls,nel,bf,der,be);
            elseif nel == 6 % wedge
              [w,ss] =  intpntw(l,linttw,lintlw,ib);
              [shl,shld,shls,be] = shlw(ss,nel,nel,der,bf);
              [shg,shgs,Jdet,be,xs] = shgw(xl,nel,shld,shls,nel,bf,der,be);
            else
              [w,ss] =  intpntb(l,lint,ib);
              [shl,shld,shls,be] = shlb(ss,nel,nel,der,bf);
              [shg,shgs,Jdet,be] = shgb(xl,nel,shld,shls,nel,bf,der,be);
            end
            
            if Jdet < 0
                elem
                Jdet %#ok<NOPTS>
            end
            c1 = Jdet*w;
            
            Fvec = [fbx; fby; fbz];
            
            % Form B matrix
            for i = 1:nel
                
              Nmat(1,(i-1)*ndf+iemat(1)) = shl(i);
              Nmat(2,(i-1)*ndf+iemat(2)) = shl(i);
              Nmat(3,(i-1)*ndf+iemat(3)) = shl(i);
                
              Bmat(Bcol1,(i-1)*ndf+iemat(1)) = shg(i,col1);
              Bmat(Bcol2,(i-1)*ndf+iemat(2)) = shg(i,col2);
              Bmat(Bcol3,(i-1)*ndf+iemat(3)) = shg(i,col3);
                 
            end

            %Compute rho, accel
            rho = rho0; %rho = f(x,y);
            accel = al_n_am(1:ndf,1:nel)*shl;
            sigma = Dmat*(Bmat*ulres);
            
            ElemF = ElemF + c1*(Nmat'*Fvec - rho*Nmat'*accel - Bmat'*sigma);

        end %je

    case 21 %Compute Stiffness
        
        ElemK = zeros(nst);

        %Set integration number
        if nel ~= 6
        lint = IntPoint3(nel);
        else
            linttw = 7;
            lintlw = 2;
            lint = linttw*lintlw;
        end
        ib = 0;
        bf = 0;
        der = 0;

        lam = Patchv*PatchE/((1+Patchv)*(1-2*Patchv));
        mu = PatchE/(2*(1+Patchv));
        fbx = 0;
        fby = 0;
        fbz = 0;
        Dmat = mu*diag([2 2 2 1 1 1]) + lam*[1; 1; 1; 0; 0; 0]*[1 1 1 0 0 0];
        Nmat = zeros(3,nst);
        Bmat = zeros(6,nst);
        
        for l = 1:lint

            % Evaluate 1-D basis functions at integration points
            if nel == 4 || nel == 10
              [w,ss] =  int3d_t(l,lint,ib);
              [shl,shld,shls,be] = shltt(ss,nel,nel,der,bf);
              [shg,shgs,Jdet,be] = shgtt(xl,nel,shld,shls,nel,bf,der,be);
            elseif nel == 6 % wedge
              [w,ss] =  intpntw(l,linttw,lintlw,ib);
              [shl,shld,shls,be] = shlw(ss,nel,nel,der,bf);
              [shg,shgs,Jdet,be,xs] = shgw(xl,nel,shld,shls,nel,bf,der,be);
            else
              [w,ss] =  intpntb(l,lint,ib);
              [shl,shld,shls,be] = shlb(ss,nel,nel,der,bf);
              [shg,shgs,Jdet,be] = shgb(xl,nel,shld,shls,nel,bf,der,be);
            end
            
            if Jdet < 0
                elem
                Jdet %#ok<NOPTS>
            end
            c1 = Jdet*w;
            
            % Form B matrix
            for i = 1:nel
                
              Bmat(Bcol1,(i-1)*ndf+iemat(1)) = shg(i,col1);
              Bmat(Bcol2,(i-1)*ndf+iemat(2)) = shg(i,col2);
              Bmat(Bcol3,(i-1)*ndf+iemat(3)) = shg(i,col3);
                 
            end
            
            ElemK = ElemK + c1*(Bmat'*Dmat*Bmat);

        end %je
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
              [shl,shld,shls,be] = shltt(ss,nel,nel,der,bf);
              [QXY,shgs,Jdet,be,sx] = shgtt(xl,nel,shld,shls,nel,bf,der,be);
            else
              [Wgt,ss] =  intpntb(je,lint,5);
              [shl,shld,shls,be] = shlb(ss,nel,nel,der,bf);
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

                ElemF(ndf*(o-1)+iemat(1)) = ElemF(ndf*(o-1)+iemat(1)) + F(1)*c1;

                ElemF(ndf*(o-1)+iemat(2)) = ElemF(ndf*(o-1)+iemat(2)) + F(2)*c1;

                ElemF(ndf*(o-1)+iemat(3)) = ElemF(ndf*(o-1)+iemat(3)) + F(3)*c1;

            end %o
            
        end %je
        ElemF;
        
        % Reorder nodes back to the orientation of the element
        SurfOrientStoE
        ElemF(1:ndf*nel) = ElemF(ilist2);
%%        
    case 25 %Stress Projection2

        ElemS = zeros(nel,npstr+1);
        ElemS2 = zeros(nel,npstr);

        lam = Patchv*PatchE/((1+Patchv)*(1-2*Patchv));
        mu = PatchE/(2*(1+Patchv));
        thick = 1;
        fbx = 0;
        fby = 0;
        fbz = 0;
        Dmat = mu*diag([2 2 2 1 1 1]) + lam*[1; 1; 1; 0; 0; 0]*[1 1 1 0 0 0];
        Nmat = zeros(3,3*nel);
        Bmat = zeros(6,3*nel);
        I1 = [1; 1; 1; 0; 0; 0];
        
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
              [shg,shgs,Jdet,be] = shgtt(xl,nel,shld,shls,nel,bf,der,be);
            else
              [Wgt,ss] =  intpntb(ll,nint,ib);
              [shl,shld,shls,be] = shlb(ss,nel,nel,der,bf);
              [shg,shgs,Jdet,be] = shgb(xl,nel,shld,shls,nel,bf,der,be);
            end
            
            % Form B matrix
            for ie = 1:nel
              
              Bmat(Bcol1,(ie-1)*3+1) = shg(ie,col1);
              Bmat(Bcol2,(ie-1)*3+2) = shg(ie,col2);
              Bmat(Bcol3,3*ie      ) = shg(ie,col3);
                 
            end
            
            epsil = Bmat*reshape(ul(:,1:nel),ndf*nel,1);
            sigma = Dmat*epsil;
            
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
              [Qxy,shgs,Jdet,be] = shgtt(xl,nel,shld,shls,nel,bf,der,be);
            else
              [Wgt,ss] =  intpntb(ll,lint,ib);
              [shl,shld,shls,be] = shlb(ss,nel,nel,der,bf);
              [Qxy,shgs,Jdet,be] = shgb(xl,nel,shld,shls,nel,bf,der,be);
            end

            w = Wgt*Jdet*thick;
            
            Vol = Vol + w;

        end %je

        for i = 1:nel
        ElemS(i,npstr+1) = 1; % use this for simple average Vol; % use this for weighted average 
        end
        
    case 60
        
        numhr = 4;
        ElemI = zeros(10,numhr);
        
    case 12 % energy
        
        ElemE = 0;
        Nmat = zeros(3,nst);
        Bmat = zeros(6,nst);

        lam = Patchv*PatchE/((1+Patchv)*(1-2*Patchv));
        mu = PatchE/(2*(1+Patchv));
        
        % Load Guass Integration Points

        %Set integration number
        if nel ~= 6
        lint = IntPoint3(nel);
        else
            linttw = 7;
            lintlw = 2;
            lint = linttw*lintlw;
        end
        der = 0;
        bf = 0;
        ib = 0;
        Dmat = mu*diag([2 2 2 1 1 1]) + lam*[1; 1; 1; 0; 0; 0]*[1 1 1 0 0 0];

        %Integration Loop
        for l = 1:lint

                % Evaluate 1-D basis functions at integration points
                if nel == 4 || nel == 10
                  [w,ss] =  int3d_t(l,lint,ib);
                  [shl,shld,shls,be] = shltt(ss,nel,nel,der,bf);
                  [shg,shgs,Jdet,be] = shgtt(xl,nel,shld,shls,nel,bf,der,be);
                elseif nel == 6 % wedge
                  [w,ss] =  intpntw(l,linttw,lintlw,ib);
                  [shl,shld,shls,be] = shlw(ss,nel,nel,der,bf);
                  [shg,shgs,Jdet,be,xs] = shgw(xl,nel,shld,shls,nel,bf,der,be);
                else
                  [w,ss] =  intpntb(l,lint,ib);
                  [shl,shld,shls,be] = shlb(ss,nel,nel,der,bf);
                  [shg,shgs,Jdet,be] = shgb(xl,nel,shld,shls,nel,bf,der,be);
                end
                
                for i = 1:nel
                  Bmat(Bcol1,(i-1)*ndf+iemat(1)) = shg(i,col1);
                  Bmat(Bcol2,(i-1)*ndf+iemat(2)) = shg(i,col2);
                  Bmat(Bcol3,(i-1)*ndf+iemat(3)) = shg(i,col3);
                end
                %Compute F, sigma, D
                E = Bmat*reshape(ul,nst,1);
                
                c1 = w*Jdet;
                ElemE = ElemE + c1/2*E'*Dmat*E;

        end %je

end