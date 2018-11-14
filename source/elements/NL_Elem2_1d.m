
% Tim Truster
% 6/20/2013
% UIUC
% Third attempt at multiscale dynamics
% Verified 6/23/13 to conserve energy for Newmark
% Use with transient=6, example file IMSDYN1.m


% Set Material Properties

Emod = mateprop(1);
Aelem = mateprop(2);
rho = mateprop(3);
% Et = mateprop(4);
% eps_y = mateprop(5);
% sig_y = Emod*eps_y;

Bcol1 = [1; 3];
Bcol2 = [2; 3];
col1 = [1; 2];
col2 = [2; 1];

switch isw %Task Switch
%%    
    case 1
        
        if ndf > 1
            
            for i = 2:ndf
                lie(i,1) = 0;
            end
            
        end
        
        nh1 = nen*ndm*3;
        
%%
    case 3 %Compute Stiffness and Residual
        
        ElemK = zeros(ndf*nel);
        ElemF = zeros(ndf*nel,1);

        %Set integration number
        lint = nel;
        
        ib = 0;
        bf = 1;
        der = 1;
        
        fb = 0;
        Dmat = Emod;
        Nmat = zeros(1,ndf*nel);
        Bmat = zeros(1,ndf*nel);
        BBmat = zeros(1,nst);
        
        [tau,intb,Mp,Kp,vol] = Tau1_1d(xl,Dmat,rho,Nbeta,tstep,nel,ndm,nen);

        % Compute kinematic fields at intermediate time levels
        ul_n_af = (1 - Nalphaf)*ul + Nalphaf*ul_n;
        al_n_am = (1 - Nalpham)*al + Nalpham*al_n;
        ulres = reshape(ul_n_af,ndf*nel,1);
        alres = reshape(al_n_am,ndf*nel,1);
        
        sw = int1d(lint);
        [shl,shld,shls,be] = shl1d(sw(1,:),lint,nel-1);
        
        for ll = 1:lint                    

            % Evaluate 1-D basis functions at integration points
            Wgt = sw(2,ll);
            [shg, shgs, Jdet, be] = shg1d(xl,ndm,nel,shld(ll,:),shls(ll,:),nen,bf,der,be);
            
            c1 = Wgt*Jdet*Aelem;
            
            Nmat(1:ndf:ndf*(nel-1)+1) = shl(ll,:);
                
            Bmat(1:ndf:ndf*(nel-1)+1) = shg';
            
            BBmat(1:ndf:ndf*(nel-1)+1) = shgs';

            %Compute rho, accel
            accel = Nmat*alres;
            sigma = Dmat*(Bmat*ulres);
            
            Fvec = fb;
            Tmat = tau*intb;
            bave = intb/vol;
          
            % load fine scale
            beta_n1 = hr(nha-1+(ll-1)*3+1);
            betadot_n1 = hr(nha-1+(ll-1)*3+2);
            betaddot_n1 = hr(nha-1+(ll-1)*3+3);
            
            divsig = Dmat*(BBmat*ulres);
            ElemTR = Tmat*(-rho*accel + divsig + Fvec);
            
            % update fine scales
            beta_n = ElemTR + tau*(1/(2*Nbeta*tstep^2)*Mp*(betaddot_n1*tstep^2 + 2*betadot_n1*tstep + 2*beta_n1) - Mp*betaddot_n1);
            betadot_n = tstep*(1-Ngamma/(2*Nbeta))*betaddot_n1 + (1-Ngamma/Nbeta)*betadot_n1 + Ngamma/(Nbeta*tstep)*(beta_n - beta_n1);
            betaddot_n = (1-1/(2*Nbeta))*betaddot_n1 - 1/(Nbeta*tstep)*betadot_n1 + 1/(Nbeta*tstep^2)*(beta_n - beta_n1);
            
            % store fine scale
            hr(nhb-1+(ll-1)*3+1) = beta_n;
            hr(nhb-1+(ll-1)*3+2) = betadot_n;
            hr(nhb-1+(ll-1)*3+3) = betaddot_n;
            
            ElemF = ElemF + c1*(Nmat'*Fvec - rho*Nmat'*accel - Bmat'*sigma ...
                                - Nmat'*rho*bave*betaddot_n);
            
            ElemK = ElemK + c1*(Nmat'*rho/(Nbeta*tstep^2)*Nmat + Bmat'*Dmat*Bmat ...
                                - (-rho/(Nbeta*tstep^2)*Nmat)'*bave*Tmat*(-rho/(Nbeta*tstep^2)*Nmat));

        end %je
        
    case 5 %Compute Mass
        
        ElemM = zeros(nst);
        
        thick = 1;
        
        % Load Guass Integration Points
        lint = nel;
        der = 0;
        bf = 1;
        ib = 0;
        Nmat = zeros(1,nst);
        
        [tau,intb,Mp,Kp,vol] = Tau1_1d(xl,0,rho,1,1,nel,ndm,nen);
        Tmat = intb^2/vol*tau;
        
        sw = int1d(lint);
        [shl,shld,shls,be] = shl1d(sw(1,:),lint,nel-1);
        
        for ll = 1:lint                    

            % Evaluate 1-D basis functions at integration points
            Wgt = sw(2,ll);
            [shg, shgs, Jdet, be] = shg1d(xl,ndm,nel,shld(ll,:),shls(ll,:),nen,bf,der,be);
            
            c1 = Wgt*Jdet*Aelem;
            
            Nmat(1:ndf:ndf*(nel-1)+1) = shl(ll,:);
            
            ElemM = ElemM + c1*rho*(Nmat'*Nmat - rho*(Nmat'*Tmat*Nmat));
%             ElemM = ElemM + c1*rho*(Nmat'*Nmat);

        end %je
        ElemM;

    case 6 %Compute Residual
        
        ElemF = zeros(ndf*nel,1);

        %Set integration number
        lint = nel;
        
        ib = 0;
        bf = 1;
        der = 1;
        
        fb = 0;
        Dmat = Emod;
        Nmat = zeros(1,ndf*nel);
        Bmat = zeros(1,ndf*nel);
        BBmat = zeros(1,nst);
        
        [tau,intb,Mp,Kp,vol] = Tau1_1d(xl,Dmat,rho,Nbeta,tstep,nel,ndm,nen);

        % Compute kinematic fields at intermediate time levels
        ul_n_af = (1 - Nalphaf)*ul + Nalphaf*ul_n;
        al_n_am = (1 - Nalpham)*al + Nalpham*al_n;
        ulres = reshape(ul_n_af,ndf*nel,1);
        alres = reshape(al_n_am,ndf*nel,1);
        
        sw = int1d(lint);
        [shl,shld,shls,be] = shl1d(sw(1,:),lint,nel-1);
        
        for ll = 1:lint                    

            % Evaluate 1-D basis functions at integration points
            Wgt = sw(2,ll);
            [shg, shgs, Jdet, be] = shg1d(xl,ndm,nel,shld(ll,:),shls(ll,:),nen,bf,der,be);
            
            c1 = Wgt*Jdet*Aelem;
            
            Nmat(1:ndf:ndf*(nel-1)+1) = shl(ll,:);
                
            Bmat(1:ndf:ndf*(nel-1)+1) = shg';
            
            BBmat(1:ndf:ndf*(nel-1)+1) = shgs';

            %Compute rho, accel
            accel = Nmat*alres;
            sigma = Dmat*(Bmat*ulres);
            
            Fvec = fb;
            Tmat = tau*intb;
            bave = intb/vol;
          
            % load fine scale
            beta_n1 = hr(nha-1+(ll-1)*3+1);
            betadot_n1 = hr(nha-1+(ll-1)*3+2);
            betaddot_n1 = hr(nha-1+(ll-1)*3+3);
            
            divsig = Dmat*(BBmat*ulres);
            ElemTR = Tmat*(-rho*accel + divsig + Fvec);
            
            % update fine scales
            beta_n = ElemTR + tau*(1/(2*Nbeta*tstep^2)*Mp*(betaddot_n1*tstep^2 + 2*betadot_n1*tstep + 2*beta_n1) - Mp*betaddot_n1);
            betadot_n = tstep*(1-Ngamma/(2*Nbeta))*betaddot_n1 + (1-Ngamma/Nbeta)*betadot_n1 + Ngamma/(Nbeta*tstep)*(beta_n - beta_n1);
            betaddot_n = (1-1/(2*Nbeta))*betaddot_n1 - 1/(Nbeta*tstep)*betadot_n1 + 1/(Nbeta*tstep^2)*(beta_n - beta_n1);
            
            % store fine scale
            hr(nhb-1+(ll-1)*3+1) = beta_n;
            hr(nhb-1+(ll-1)*3+2) = betadot_n;
            hr(nhb-1+(ll-1)*3+3) = betaddot_n;
            
            ElemF = ElemF + c1*(Nmat'*Fvec - rho*Nmat'*accel - Bmat'*sigma ...
                                - Nmat'*rho*bave*betaddot_n);

        end %je

    case 21 %Compute Stiffness
        
        ElemK = zeros(nst);

        %Set integration number
        lint = IntPoint(nel);
        ib = 0;
        bf = 1;
        der = 1;

        lam = Patchv*PatchE/((1+Patchv)*(1-2*Patchv));
        mu = PatchE/(2*(1+Patchv));
        thick = 1;
        Dmat = mu*diag([2 2 1]) + lam*[1; 1; 0]*[1 1 0];
        Nmat = zeros(2,nst);
        Bmat = zeros(3,nst);
        BBmat = zeros(6,nst);
        cmat = Dmat;
        I4 = [1 0 0 0 0 1
              0 0 1 0 1 0];
        P1 = [1 0 0 0 0 0
              0 0 0 0 0 1
              0 0 1 1 0 0
              0 0 1 0 0 0
              0 0 0 0 1 0
              0 1 0 0 0 1];
        
        [tau,intb,Mp,Kp,vol] = Tau9_2d(xl,Dmat,rho,Nbeta,tstep,nel);
        
        for ll = 1:lint                    

            % Evaluate 1-D basis functions at integration points
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
                
              Nmat(1,(ie-1)*2+1) = shl(ie);
              Nmat(2,(ie-1)*2+2) = shl(ie);
                
              Bmat(Bcol1,(ie-1)*2+1) = shg(ie,col1);
              Bmat(Bcol2,(ie-1)*2+2) = shg(ie,col2);
              
              BBmat(:,(ie-1)*2+1:2*ie) = [shgs(ie,1) 0               
                     shgs(ie,2) 0             
                     shgs(ie,3) 0               
                     0         shgs(ie,1)       
                     0         shgs(ie,2)       
                     0         shgs(ie,3)];
                 
            end

            %Compute rho, accel
            rho = rho0; %rho = f(x,y);
            accel = al_n_am*shl;
            sigma = Dmat*(Bmat*ulres);
            
            Fvec = [fbx; fby];
            Tmat = tau*intb;
            bave = intb/vol;
            
            ElemK = ElemK + c1*(Nmat'*rho/(Nbeta*tstep^2)*Nmat + Bmat'*Dmat*Bmat ...
                                - (-rho/(Nbeta*tstep^2)*Nmat + BBT*BBmat)'*bave*Tmat*(-rho/(Nbeta*tstep^2)*Nmat + BBT*BBmat));

        end %je
        
    case 40 % Initialize FS acceleration

        %Set integration number
        lint = nel;
        ib = 0;
        bf = 1;
        der = 1;
        
        fb = 0;
        Dmat = Emod;
        Nmat = zeros(1,ndf*nel);
        Bmat = zeros(1,ndf*nel);
        BBmat = zeros(1,ndf*nel);
        
        [tau,intb,Mp,Kp,vol] = Tau1_1d(xl,0,rho,1,1,nel,ndm,nen);
        Tmat = intb*tau;

        % Compute kinematic fields at intermediate time levels
        ul_n_af = ul_n;
        al_n_am = al_n;
        ulres = reshape(ul_n_af,ndf*nel,1);
        alres = reshape(al_n_am,ndf*nel,1);
        
        sw = int1d(lint);
        [shl,shld,shls,be] = shl1d(sw(1,:),lint,nel-1);
        
        for ll = 1:lint                    

            % Evaluate 1-D basis functions at integration points
            Wgt = sw(2,ll);
            [shg, shgs, Jdet, be] = shg1d(xl,ndm,nel,shld(ll,:),shls(ll,:),nen,bf,der,be);
            
            c1 = Wgt*Jdet*Aelem;
            
            Nmat(1:ndf:ndf*(nel-1)+1) = shl(ll,:);
                
            Bmat(1:ndf:ndf*(nel-1)+1) = shg';
            
            BBmat(1:ndf:ndf*(nel-1)+1) = shgs';

            accel = Nmat*alres;
            
            divsig = Dmat*(BBmat*ulres);
            
            Fvec = fb;
            
            ElemTR = Tmat*(-rho*accel + divsig + Fvec);
            
            hr(nha-1+(ll-1)*3+3) = ElemTR;

        end %je

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
                    
            c1 = Wgt*tm3*drdr*thick;
            for o=1:nel
                don = shl(o);
                F = don*traction';

    %                 Fmtl = F'*t; %Magnitudes of F dot tunit(l=1:3)
    % %                 for l = 1:3
    % %                     for m = 1:3
    %                         Ftl = t*Fmtl'; %Sum of Vectors {F.t(l)}t(m)
    % %                     end
    % %                 end  t*t'*F = eye*F = F

                ElemF(ndf*o-1) = ElemF(ndf*o-1) + F(1)*c1;

                ElemF(ndf*o-0)   = ElemF(ndf*o-0)   + F(2)*c1;

            end %o

        end %ie
        
        % Reorder nodes back to the orientation of the element
        SurfOrientStoE
        ElemF(1:ndf*nel) = ElemF(ilist2);
        
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
        
    case 12 % energy
        
        ElemE = 0;
        ElemE1 = 0;
        ElemE2 = 0;
        
        
        % Load Guass Integration Points

        lint = nel;
        der = 1;
        bf = 1;
        ib = 0;
        
        Dmat = Emod;
        D = Dmat;
        Nmat = zeros(1,nst);
        Bmat = zeros(1,nst);
        BBmat = zeros(1,nst);
        cmat = Dmat;
        
        [tau,intb,Mp,Kp,vol] = Tau1_1d(xl,Dmat,rho,Nbeta,tstep,nel,ndm,nen);
        bave = intb/vol;
        
        sw = int1d(lint);
        [shl,shld,shls,be] = shl1d(sw(1,:),lint,nel-1);

        %Integration Loop
        for ll = 1:lint

                % Evaluate 1-D basis functions at integration points
                Wgt = sw(2,ll);
                [shg, shgs, Jdet, be] = shg1d(xl,ndm,nel,shld(ll,:),shls(ll,:),nen,bf,der,be);

                c1 = Wgt*Jdet*Aelem;

                Nmat(1:ndf:ndf*(nel-1)+1) = shl(ll,:);

                Bmat(1:ndf:ndf*(nel-1)+1) = shg';

                BBmat(1:ndf:ndf*(nel-1)+1) = shgs';
                
                %Compute F, sigma, D
                E = Bmat*reshape(ul,nst,1);
                velo = Nmat*reshape(vl,nst,1);
                
                dE = BBmat*reshape(ul,nst,1);
                divsig = D*dE;
                
                % load fine scale
                beta_n = hr(nhb-1+(ll-1)*3+1);
                betadot_n = hr(nhb-1+(ll-1)*3+2);
                
                ElemE = ElemE + c1/2*(E'*D*E + velo'*rho*velo);
                ElemE1 = ElemE1 + c1*(betadot_n'*bave*rho*velo - beta_n'*bave*divsig);
                ElemE2 = ElemE2 + c1/2*(betadot_n'*Mp/vol*betadot_n + beta_n'*Kp/vol*beta_n);
%                               + 1/2*(betadot_n'*intb*rho*velo - beta_n'*intb*divsig + betadot_n'*Mp*betadot_n + beta_n'*Kp*beta_n);
%                 ElemE = ElemE + c1/2*(E'*D*E + velo'*rho*velo + betadot_n'*intb/vol*rho*velo - beta_n'*intb/vol*divsig) ...
%                               + c1/2*(betadot_n'*intb*rho*velo - beta_n'*intb*divsig + betadot_n'*Mp*betadot_n + beta_n'*Kp*beta_n);
% %                               + 1/2*(betadot_n'*intb*rho*velo - beta_n'*intb*divsig + betadot_n'*Mp*betadot_n + beta_n'*Kp*beta_n);

        end %je
ElemE = ElemE + ElemE1 + ElemE2;
end