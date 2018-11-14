
% Tim Truster
% 10/03/2013
% UIUC
% Q4 linear elastothermodynamics
%
% Fint (isw=6), Tractions (isw=-1), Stiffness (isw=3), and body force
% (isw=15) were verified using Itherm.m and varying the initialization of
% the displacement and theta state variables. In all cases, the algorithm
% converges to the correct solution after one iteration. Exact solution is
% derived in temperprob.m.
%
% Naming convention follows Masud's linear thermoelasticity notes.

% Set Material Properties

PatchE = mateprop(2);
Patchv = mateprop(3);
rho0 = mateprop(4);
ct = mateprop(5);
kt = mateprop(6);
alpt = mateprop(7);

Bcol1 = [1; 3];
Bcol2 = [2; 3];
col1 = [1; 2];
col2 = [2; 1];

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
        lint = IntPoint(nel);
        ib = 0;
        bf = 1;
        der = 0;

        lam = Patchv*PatchE/((1+Patchv)*(1-2*Patchv));
        mu = PatchE/(2*(1+Patchv));
        thick = 1;
        Dmat = mu*diag([2 2 1]) + lam*[1; 1; 0]*[1 1 0];
        mvec = 3*alpt*(lam+2/3*mu)*[1; 1; 0]; % thermal coupling tensor in vector form
        mmatT = 3*alpt*(lam+2/3*mu)*eye(2); % thermal coupling tensor in matrix form
        Nmat = zeros(3,nst);
        Bmat = zeros(3,nst);
        Mmat = diag([rho/(Nbeta*tstep^2) rho/(Nbeta*tstep^2) ct/(Nalphat*tstep)]); % mass matrix and thermal evolution
        DmatT = kt*eye(2); % Fourier matrix
        BmatT = zeros(2,nst);

        % Compute kinematic fields at intermediate time levels
        % ONLY VERIFIED FOR alphaf = alpham = 0
        ul_n_af = (1 - Nalphaf)*ul + Nalphaf*ul_n;
        al_n_am = (1 - Nalpham)*al + Nalpham*al_n;
        ulres = reshape(ul_n_af,ndf*nel,1);
        alres = reshape(al_n_am,ndf*nel,1);
        vlres = reshape(vl,ndf*nel,1);
        
        for ll = 1:lint                    

%             Evaluate 1-D basis functions at integration points
            if nel == 3 || nel == 6
              [w,litr,lits] =  intpntt(ll,lint,ib);
              [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
              [shg,shgs,Jdet,be] = shgt(xl,nel,shld,shls,nel,bf,der,be);
            else
              [w,litr,lits] =  intpntq(ll,lint,ib);
              [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
              [shg,shgs,Jdet,be] = shgq(xl,nel,shld,shls,nel,bf,der,be);
            end
            
            if Jdet < 0
                elem
                Jdet %#ok<NOPTS>
            end
            c1 = Jdet*w*thick;
            
            % Form B matrix
            for ie = 1:nel
                
              Nmat(1,(ie-1)*ndf+1) = shl(ie);
              Nmat(2,(ie-1)*ndf+2) = shl(ie);
              Nmat(3,(ie-1)*ndf+3) = shl(ie);
                
              Bmat(Bcol1,(ie-1)*ndf+1) = shg(ie,col1);
              Bmat(Bcol2,(ie-1)*ndf+2) = shg(ie,col2);
              BmatT(1,(ie-1)*ndf+3) = shg(ie,1);
              BmatT(2,(ie-1)*ndf+3) = shg(ie,2);
                 
            end

            %Compute rho, accel
            rho = rho0; %rho = f(x,y);
            accel = Nmat*alres;
            timeder = Nmat*vlres;
            theta = [0 0 1]*Nmat*ulres; % temperature
            veloc = timeder(1:2); % velocity
            thetadot = timeder(3); % temperature rate
            sigma = Dmat*(Bmat*ulres);
            qvec = DmatT*(BmatT*ulres);
            
            ElemF = ElemF - c1*(Nmat'*rho*diag([1 1 0])*accel + Nmat'*ct*[0 0 1]'*thetadot ...
                  + Bmat'*(sigma - mvec*theta) + BmatT'*(qvec + mmatT*veloc));
            
            ElemK = ElemK + c1*(Nmat'*Mmat*Nmat + Bmat'*Dmat*Bmat ...
                  - Bmat'*mvec*[0 0 1]*Nmat + BmatT'*DmatT*BmatT ...
                  + BmatT'*mmatT*Ngamma/(Nbeta*tstep)*[1 0 0; 0 1 0]*Nmat);

        end %je
        ElemF;
    case 5 %Compute Mass
        
        ElemM = zeros(nst);
        
        thick = 1;
        
        % Load Guass Integration Points
        lint = IntPoint(nel);
        der = 0;
        bf = 1;
        ib = 0;
        Nmat = zeros(3,nst);
        Mmat = diag([rho0 rho0 ct]);
        
        for ll = 1:lint                    

%             Evaluate 1-D basis functions at integration points
            if nel == 3 || nel == 6
              [w,litr,lits] =  intpntt(ll,lint,ib);
              [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
              [shg,shgs,Jdet,be] = shgt(xl,nel,shld,shls,nel,bf,der,be);
            else
              [w,litr,lits] =  intpntq(ll,lint,ib);
              [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
              [shg,shgs,Jdet,be] = shgq(xl,nel,shld,shls,nel,bf,der,be);
            end
            
            if Jdet < 0
                elem
                Jdet %#ok<NOPTS>
            end
            c1 = Jdet*w*thick;
            
            % Form B matrix
            for ie = 1:nel
                
              Nmat(1,(ie-1)*ndf+1) = shl(ie);
              Nmat(2,(ie-1)*ndf+2) = shl(ie);
              Nmat(3,(ie-1)*ndf+3) = shl(ie);
                 
            end

            %Compute rho, accel
            rho = rho0; %rho = f(x,y);
            
            ElemM = ElemM + c1*(Nmat'*Mmat*Nmat);

        end %je

    case 6 %Compute Residual
        
        ElemK = zeros(nst);
        ElemF = zeros(nst,1);
        ElemFn = zeros(nst,1);

        %Set integration number
        lint = IntPoint(nel);
        ib = 0;
        bf = 1;
        der = 0;

        lam = Patchv*PatchE/((1+Patchv)*(1-2*Patchv));
        mu = PatchE/(2*(1+Patchv));
        thick = 1;
        Dmat = mu*diag([2 2 1]) + lam*[1; 1; 0]*[1 1 0];
        mvec = 3*alpt*(lam+2/3*mu)*[1; 1; 0];
        mmatT = 3*alpt*(lam+2/3*mu)*eye(2);
        Nmat = zeros(3,nst);
        Bmat = zeros(3,nst);
        Mmat = diag([rho/(Nbeta*tstep^2) rho/(Nbeta*tstep^2) ct]);
        DmatT = kt*eye(2);
        BmatT = zeros(2,nst);

        % Compute kinematic fields at intermediate time levels
        % ONLY VERIFIED FOR alphaf = alpham = 0
        ul_n_af = (1 - Nalphaf)*ul + Nalphaf*ul_n;
        al_n_am = (1 - Nalpham)*al + Nalpham*al_n;
        ulres = reshape(ul_n_af,ndf*nel,1);
        alres = reshape(al_n_am,ndf*nel,1);
        vlres = reshape(vl,ndf*nel,1);
        
        for ll = 1:lint                    

%             Evaluate 1-D basis functions at integration points
            if nel == 3 || nel == 6
              [w,litr,lits] =  intpntt(ll,lint,ib);
              [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
              [shg,shgs,Jdet,be] = shgt(xl,nel,shld,shls,nel,bf,der,be);
            else
              [w,litr,lits] =  intpntq(ll,lint,ib);
              [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
              [shg,shgs,Jdet,be] = shgq(xl,nel,shld,shls,nel,bf,der,be);
            end
            
            if Jdet < 0
                elem
                Jdet %#ok<NOPTS>
            end
            c1 = Jdet*w*thick;
            
            % Form B matrix
            for ie = 1:nel
                
              Nmat(1,(ie-1)*ndf+1) = shl(ie);
              Nmat(2,(ie-1)*ndf+2) = shl(ie);
              Nmat(3,(ie-1)*ndf+3) = shl(ie);
                
              Bmat(Bcol1,(ie-1)*ndf+1) = shg(ie,col1);
              Bmat(Bcol2,(ie-1)*ndf+2) = shg(ie,col2);
              BmatT(1,(ie-1)*ndf+3) = shg(ie,1);
              BmatT(2,(ie-1)*ndf+3) = shg(ie,2);
                 
            end

            %Compute rho, accel
            rho = rho0; %rho = f(x,y);
            accel = Nmat*alres;
            timeder = Nmat*vlres;
            theta = [0 0 1]*Nmat*ulres;
            xint = xl(1,:)*shl;
            yint = xl(2,:)*shl;
            veloc = timeder(1:2);
            thetadot = timeder(3);
            sigma = Dmat*(Bmat*ulres);
            qvec = DmatT*(BmatT*ulres);
            
            ElemF = ElemF - c1*(Nmat'*rho*diag([1 1 0])*accel + Nmat'*ct*[0 0 1]'*thetadot ...
                  + Bmat'*(sigma - mvec*theta) + BmatT'*(qvec + mmatT*veloc));

        end %je
        
        ElemF;
%%
    case 15 % body force for linear element
        
%         ElemF = zeros(nst,1); In FormFE
        
%         if iprob == 1
            
        Nmat = zeros(3,nst);
        
        % Load Gauss Points for quadrature
        if nel == 3
            lint = lintt3;%13;
        elseif nel == 4
            lint = lintq4;
        elseif nel == 6
            lint = lintt6;%13;
        elseif nel == 9
            lint = lintq9;
        end
        
        thick = 1;

        der = 0;
        bf = 0;
        
%         if iprob == 5
%             fbx = 0;
%             grav = 9.81;
%             fby = -rho*grav;
%             fb = [fbx; fby];
%         else
            fb = bodyf';
%         end
        
        for je = 1:lint

            if nel == 3 || nel == 6
                [Wgt,r,s] = intpntt(je,lint,0);
            else %if nel == 4 || nel == 9
                [Wgt,r,s] = intpntq(je,lint,0);
            end
            
            if nel == 3 || nel == 6
                [shl,shld,shls,be] = shlt(r,s,nel,nel,der,bf);
                [shg, shgs, Jdet, be, xs] = shgt(xl,nel,shld,shls,nen,bf,der,be);
            elseif nel == 4 || nel == 9
                [shl,shld,shls,be] = shlq(r,s,nel,nel,der,bf);
                [shg, shgs, Jdet, be, xs] = shgq(xl,nel,shld,shls,nen,bf,der,be);
            end
%             % Evaluate 1-D basis functions at integration points
%             [shp,shp2,be] = shpl_2d(r,s,nel,ideriv,ep,enrich);
%             %Evaluate first derivatives of basis functions at int. point
%             [Qxy, cartd2, Jdet] = shpg_2d(shp,shp2,xl,nel2,ideriv,be);
            c1 = Wgt*Jdet*thick;
            
            % Form B matrix
            for ie = 1:nel
                
              Nmat(1,(ie-1)*ndf+1) = shl(ie);
              Nmat(2,(ie-1)*ndf+2) = shl(ie);
              Nmat(3,(ie-1)*ndf+3) = shl(ie);
                                    
            end
            
            ElemF = ElemF + c1*(Nmat'*fb);
                
        end %je
        
%         end
        ElemF;

    case 21 %Compute Stiffness
        
        ElemK = zeros(nst);

        %Set integration number
        lint = IntPoint(nel);
        ib = 0;
        bf = 1;
        der = 0;

        lam = Patchv*PatchE/((1+Patchv)*(1-2*Patchv));
        mu = PatchE/(2*(1+Patchv));
        thick = 1;
        Dmat = mu*diag([2 2 1]) + lam*[1; 1; 0]*[1 1 0];
        mvec = 3*alpt*(lam+2/3*mu)*[1; 1; 0];
        mmatT = 3*alpt*(lam+2/3*mu)*eye(2);
        Nmat = zeros(3,nst);
        Bmat = zeros(3,nst);
        Mmat = diag([rho/(Nbeta*tstep^2) rho/(Nbeta*tstep^2) ct/(Nalphat*tstep)]);
        DmatT = kt*eye(2);
        BmatT = zeros(2,nst);

        % Compute kinematic fields at intermediate time levels
        % ONLY VERIFIED FOR alphaf = alpham = 0
        ul_n_af = (1 - Nalphaf)*ul + Nalphaf*ul_n;
        al_n_am = (1 - Nalpham)*al + Nalpham*al_n;
        ulres = reshape(ul_n_af,ndf*nel,1);
        alres = reshape(al_n_am,ndf*nel,1);
        vlres = reshape(vl,ndf*nel,1);
        
        for ll = 1:lint                    

%             Evaluate 1-D basis functions at integration points
            if nel == 3 || nel == 6
              [w,litr,lits] =  intpntt(ll,lint,ib);
              [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
              [shg,shgs,Jdet,be] = shgt(xl,nel,shld,shls,nel,bf,der,be);
            else
              [w,litr,lits] =  intpntq(ll,lint,ib);
              [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
              [shg,shgs,Jdet,be] = shgq(xl,nel,shld,shls,nel,bf,der,be);
            end
            
            if Jdet < 0
                elem
                Jdet %#ok<NOPTS>
            end
            c1 = Jdet*w*thick;
            
            % Form B matrix
            for ie = 1:nel
                
              Nmat(1,(ie-1)*ndf+1) = shl(ie);
              Nmat(2,(ie-1)*ndf+2) = shl(ie);
              Nmat(3,(ie-1)*ndf+3) = shl(ie);
                
              Bmat(Bcol1,(ie-1)*ndf+1) = shg(ie,col1);
              Bmat(Bcol2,(ie-1)*ndf+2) = shg(ie,col2);
              BmatT(1,(ie-1)*ndf+3) = shg(ie,1);
              BmatT(2,(ie-1)*ndf+3) = shg(ie,2);
                 
            end

            %Compute rho, accel
            rho = rho0; %rho = f(x,y);
            accel = Nmat*alres;
            timeder = Nmat*vlres;
            theta = [0 0 1]*Nmat*ulres;
            veloc = timeder(1:2);
            thetadot = timeder(3);
            sigma = Dmat*(Bmat*ulres);
            qvec = DmatT*(BmatT*ulres);
            
            ElemK = ElemK + c1*(Nmat'*Mmat*Nmat + Bmat'*Dmat*Bmat ...
                  - Bmat'*mvec*[0 0 1]*Nmat + BmatT'*DmatT*BmatT ...
                  + BmatT'*mmatT*Ngamma/(Nbeta*tstep)*[1 0 0; 0 1 0]*Nmat);

        end %je

    case -1 % boundary tractions
        
        ElemF = zeros(nst,1);
        Nmat = zeros(3,nst);
        
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
        
        lint = 4;
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
        der = 0;
        bf = 0;
        
        for je = 1:lint

            if nel == 3
                [Wgt,r,s] = intpntt(je,lint,1);
            else %if nel == 4
                [Wgt,r,s] = intpntq(je,lint,1);
            end
            r = drdr*(r-ro)+eR1;

            if nel == 3 || nel == 6
              [shl,shld,shls,be] = shlt(r,s,nel,nel,der,bf);
              [Qxy,shgs,Jdet,bubble,xs] = shgt(xl(:,1:nel),nel,shld,shls,nel,bf,der,be);
            else
              [shl,shld,shls,be] = shlq(r,s,nel,nel,der,bf);
              [Qxy,shgs,Jdet,bubble,xs] = shgq(xl(:,1:nel),nel,shld,shls,nel,bf,der,be);
            end
                    
            %Evaluate tangent and normal vectors
            t1 = [xs(:,1); 0];
            [tm1, tu1] = VecNormalize(t1);
            t2 = [0; 0; 1];
            tm2 = 1;
            tu2 = t2';
            t3 = VecCrossProd(t1,t2);
            [tm3, tu3] = VecNormalize(t3);
            
            xint = xl(1,:)*shl;
            yint = xl(2,:)*shl;
            
            if iprob == 1
                if(strcmp(AssemQuant,'AssemLoad'))
                    lam = Patchv*PatchE/((1+Patchv)*(1-2*Patchv));
                    mu = PatchE/(2*(1+Patchv));
                    sigma = [((lam+2*mu)-mu*Patchv/(1-Patchv))*h 0 0]; % Dmat*epsil(t)
                    qvec = zeros(1,2);
                else % non-proportional
                    lam = Patchv*PatchE/((1+Patchv)*(1-2*Patchv));
                    mu = PatchE/(2*(1+Patchv));
                    sigma = -[3*alpt*(lam+2/3*mu)*g*xint 3*alpt*(lam+2/3*mu)*g*xint 0]; %theta*m_ij
                    qvec = [kt*g 0]; % kt_ij*grad(theta)
                    qvec = qvec +3*alpt*(lam+2/3*mu)*h*[xint -yint*Patchv/(1-Patchv)]; %m_ij*velo_i
                end
                traction = [sigma(1) sigma(3); sigma(3) sigma(2)]*tu3(1:2)';
                tflux = qvec*tu3(1:2)';
            end
                    
            c1 = Wgt*tm3*drdr*thick;
            
            for ie = 1:nel
                
              Nmat(1,(ie-1)*ndf+1) = shl(ie);
              Nmat(2,(ie-1)*ndf+2) = shl(ie);
              Nmat(3,(ie-1)*ndf+3) = shl(ie);
                 
            end
            
            ElemF = ElemF + c1*Nmat'*[traction; tflux];

        end %ie
        
        % Reorder nodes back to the orientation of the element
        SurfOrientStoE
        ElemF(1:ndf*nel) = ElemF(ilist2);
        ElemF;
        
    case 51 % Volume stress/strain homogenization
        
        ElemSS = zeros(13,1);

    case 52 % Surface strain homogenization
        
        ElemSS = zeros(13,1);
%%        
    case 25 %Stress Projection2

        ElemS = zeros(nel,npstr+1);
        ElemS2 = zeros(nel,npstr);
        stresID = [1 2 4 0 1 3];

        lam = Patchv*PatchE/((1+Patchv)*(1-2*Patchv));
        mu = PatchE/(2*(1+Patchv));
        thick = 1;
        fbx = 0;
        fby = 0;
        fbz = 0;
        Dmat = mu*diag([2 2 1]) + lam*[1; 1; 0]*[1 1 0];
        Nmat = zeros(2,2*nel);
        Bmat = zeros(3,2*nel);
        I1 = [1; 1; 0];
        
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
        bf = 0;
        ib = 0;
        
        ulres = reshape(ul,ndf*nel,1);

        %Stress Loop
        for ll = 1:nint

            %Evaluate first derivatives of basis functions at int. point
            if nel == 3 || nel == 6
              [w,litr,lits] =  intpntt(ll,nint,ib);
              [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
              [shg,shgs,Jdet,be] = shgt(xl,nel,shld,shls,nel,bf,der,be);
            else
              [w,litr,lits] =  intpntq(ll,nint,ib);
              [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
              [shg,shgs,Jdet,be] = shgq(xl,nel,shld,shls,nel,bf,der,be);
            end
            
            % Form B matrix
            for ie = 1:nel
              
              Bmat(Bcol1,(ie-1)*2+1) = shg(ie,col1);
              Bmat(Bcol2,(ie-1)*2+2) = shg(ie,col2);
                 
            end
            
            epsil = Bmat*ulres;
            sigma = Dmat*epsil;
            
%             for stres = 1:npstr
%             
%             if stres <= 3 % stress components
%                 sigmas = sigma(stresID(stres));
% %             elseif stres >= 5
% %                 if stres <= 6 % principal stresses
% %                 sigma2 = [sigma(1) sigma(4) sigma(6); sigma(4) sigma(2) sigma(5); sigma(6) sigma(5) sigma(3)];
% %                 psig = eig(sigma2);
% %                 sigmas = psig(stresID(stres));
% %                 else % hydrostatic stress
% %                 sigmas = 1/3*sigma'*I1;
% %                 end
% %             else % von Mises stress
% %                 trs = sigma'*I1;
% %                 dsig = sigma - 1/3*trs*I1;
% %                 sigmas = sqrt(3/2*(dsig'*dsig));
%             end
%             
%             ElemS2(ll,stres) = sigmas;
%             
%             end
            ElemS2(ll,1:3) = sigma';

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
            
%             for stres = 1:npstr
%                 
%                 sigmas = ElemS2(1:nint,stres)'*shpS;
%                 ElemS(ll,stres) = sigmas;
%                 
%             end
            ElemS(ll,1:npstr) = (ElemS2(1:nint,:)'*shpS)';
            
        end
        
%         %Integration Loop
%         Vol = 0;
%         for ll = 1:lint
% 
%             %Evaluate first derivatives of basis functions at int. point
%             if nel == 3 || nel == 6
%               [Wgt,litr,lits] =  intpntt(ll,lint,ib);
%               [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
%               [Qxy, shgs, Jdet,be] = shgt(xl+ul(1:2,:),nel,shld,shls,nel,bf,der,be);
%             else
%               [Wgt,litr,lits] =  intpntq(ll,lint,ib);
%               [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
%               [Qxy, shgs, Jdet,be] = shgq(xl+ul(1:2,:),nel,shld,shls,nel,bf,der,be);
%             end
% 
%             w = Wgt*Jdet*thick;
%             
%             Vol = Vol + w;
% 
%         end %je

        for i = 1:nel
        ElemS(i,npstr+1) = 1; % use this for simple average Vol; % use this for weighted average 
        end
        
    case 60
        
        numhr = 4;
        ElemI = zeros(10,numhr);

end