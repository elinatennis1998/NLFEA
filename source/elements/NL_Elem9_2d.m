
% Tim Truster
% 6/19/2013
% UIUC
% First attempt at multiscale dynamics; this one is not debugged and does
% not work. Not modified after 6/19 except to add these remarks.


% Set Material Properties

PatchE = mateprop(2);
Patchv = mateprop(3);
rho0 = mateprop(4);

Bcol1 = [1; 3];
Bcol2 = [2; 3];
col1 = [1; 2];
col2 = [2; 1];

switch isw %Task Switch
%%    
    case 1
        
        if ndf > 2
            
            for i = 3:ndf
                lie(i,1) = 0;
            end
            
        end
        
        nh1 = 4*2*3;
        
%%
    case 3 %Compute Stiffness and Residual
        
        ElemK = zeros(nst);
        ElemF = zeros(nst,1);
        ElemFn = zeros(nst,1);

        %Set integration number
        lint = IntPoint(nel);
        ib = 0;
        bf = 1;
        der = 1;

        lam = Patchv*PatchE/((1+Patchv)*(1-2*Patchv));
        mu = PatchE/(2*(1+Patchv));
        thick = 1;
        fbx = 0;
        fby = 0;
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
        
        [tau,intb,Mp,Kp] = Tau9_2d(xl,Dmat,rho0,Nbeta,tstep,nel);

        % Compute kinematic fields at intermediate time levels
        ul_n_af = (1 - Nalphaf)*ul + Nalphaf*ul_n;
        al_n_am = (1 - Nalpham)*al + Nalpham*al_n;
        ulres = reshape(ul_n_af,ndf*nel,1);
        alres = reshape(al_n_am,ndf*nel,1);
        
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
            
            Fvec = [fbx; fby];
            Tmat = tau*intb;
          
            % load fine scale
            beta_n1 = hr(nha+(ll-1)*6+1:nha+(ll-1)*6+2);
            betadot_n1 = hr(nha+(ll-1)*6+3:nha+(ll-1)*6+4);
            betaddot_n1 = hr(nha+(ll-1)*6+5:nha+(ll-1)*6+6);
            
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
            
            diagc = [cmat(1:3,1:3) zeros(3); zeros(3) cmat(1:3,1:3)];
            term8b = I4*diagc*P1; %page 8b
            BBT = term8b;
            ElemTR = Tmat*(-rho*accel + BBT*BBmat*ulres + Fvec(1:2));
            
            % update fine scales
            beta_n = ElemTR + tau*(1/(2*Nbeta*tstep^2)*Mp*(betaddot_n1*tstep^2 + 2*betadot_n1*tstep + 2*beta_n1) - Mp*betaddot_n1);
            betadot_n = (1-Ngamma/Nbeta)*betaddot_n1 + Ngamma*tstep*(1-1/(2*Nbeta))*betadot_n1 + Ngamma/(Nbeta*tstep)*(beta_n - beta_n1);
            betaddot_n = (1-1/(2*Nbeta))*betaddot_n1 - 1/(Nbeta*tstep)*betadot_n1 + 1/(Nbeta*tstep^2)*(beta_n - beta_n1);
            
            % store fine scale
            hr(nhb+(ll-1)*6+1:nhb+(ll-1)*6+2) = beta_n;
            hr(nhb+(ll-1)*6+3:nhb+(ll-1)*6+4) = betadot_n;
            hr(nhb+(ll-1)*6+5:nhb+(ll-1)*6+6) = betaddot_n;
            
            ElemF = ElemF + c1*(Nmat'*Fvec - rho*Nmat'*accel - Bmat'*sigma ...
                                - BBmat'*BBT'*be(3)*beta_n + Nmat'*rho*be(3)*betaddot_n);
            
            ElemK = ElemK + c1*(Nmat'*rho/(Nbeta*tstep^2)*Nmat + Bmat'*Dmat*Bmat ...
                                - (-rho/(Nbeta*tstep^2)*Nmat + BBT*BBmat)'*be(3)*Tmat*(-rho/(Nbeta*tstep^2)*Nmat + BBT*BBmat));

        end %je
        
    case 5 %Compute Mass
        
        ElemM = zeros(nst);
        
        thick = 1;
        
        % Load Guass Integration Points
        lint = IntPoint(nel);
        der = 0;
        bf = 1;
        ib = 0;
        Nmat = zeros(2,nst);
        
%         [tau,intb] = Tau9_2d(xl,0,rho0,1,1,nel);
%         if ~fsinit 
%             intb = 0;
%         end
        
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
                 
            end

            %Compute rho, accel
            rho = rho0; %rho = f(x,y);
            
%             ElemM = ElemM + c1*rho*(Nmat'*Nmat - intb*rho*be(3)*(Nmat'*tau*Nmat));
            ElemM = ElemM + c1*rho*(Nmat'*Nmat);

        end %je

    case 6 %Compute Residual
        
        ElemF = zeros(nst,1);
        ElemFn = zeros(nst,1);

        %Set integration number
        lint = IntPoint(nel);
        ib = 0;
        bf = 1;
        der = 1;

        lam = Patchv*PatchE/((1+Patchv)*(1-2*Patchv));
        mu = PatchE/(2*(1+Patchv));
        thick = 1;
        fbx = 0;
        fby = 0;
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
        
        [tau,intb,Mp,Kp] = Tau9_2d(xl,Dmat,rho0,Nbeta,tstep,nel);

        % Compute kinematic fields at intermediate time levels
        ul_n_af = (1 - Nalphaf)*ul + Nalphaf*ul_n;
        al_n_am = (1 - Nalpham)*al + Nalpham*al_n;
        ulres = reshape(ul_n_af,ndf*nel,1);
        alres = reshape(al_n_am,ndf*nel,1);
        
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
            
            Fvec = [fbx; fby];
            Tmat = tau*intb;
          
            % load fine scale
            beta_n1 = hr(nha+(ll-1)*6+1:nha+(ll-1)*6+2);
            betadot_n1 = hr(nha+(ll-1)*6+3:nha+(ll-1)*6+4);
            betaddot_n1 = hr(nha+(ll-1)*6+5:nha+(ll-1)*6+6);
            
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
            
            diagc = [cmat(1:3,1:3) zeros(3); zeros(3) cmat(1:3,1:3)];
            term8b = I4*diagc*P1; %page 8b
            BBT = term8b;
            ElemTR = Tmat*(-rho*accel + BBT*BBmat*ulres + Fvec(1:2));
            
            % update fine scales
            beta_n = ElemTR + tau*(1/(2*Nbeta*tstep^2)*Mp*(betaddot_n1*tstep^2 + 2*betadot_n1*tstep + 2*beta_n1) - Mp*betaddot_n1);
            betadot_n = (1-Ngamma/Nbeta)*betaddot_n1 + Ngamma*tstep*(1-1/(2*Nbeta))*betadot_n1 + Ngamma/(Nbeta*tstep)*(beta_n - beta_n1);
            betaddot_n = (1-1/(2*Nbeta))*betaddot_n1 - 1/(Nbeta*tstep)*betadot_n1 + 1/(Nbeta*tstep^2)*(beta_n - beta_n1);
            
            % store fine scale
            hr(nhb+(ll-1)*6+1:nhb+(ll-1)*6+2) = beta_n;
            hr(nhb+(ll-1)*6+3:nhb+(ll-1)*6+4) = betadot_n;
            hr(nhb+(ll-1)*6+5:nhb+(ll-1)*6+6) = betaddot_n;
            
            ElemF = ElemF + c1*(Nmat'*Fvec - rho*Nmat'*accel - Bmat'*sigma ...
                                - BBmat'*BBT'*be(3)*beta_n + Nmat'*rho*be(3)*betaddot_n);

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
        Bmat = zeros(3,nst);
        
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
                
              Bmat(Bcol1,(ie-1)*2+1) = shg(ie,col1);
              Bmat(Bcol2,(ie-1)*2+2) = shg(ie,col2);
                 
            end
            
            ElemK = ElemK + c1*(Bmat'*Dmat*Bmat);

        end %je
        
    case 40 % Initialize FS acceleration

        %Set integration number
        lint = IntPoint(nel);
        ib = 0;
        bf = 1;
        der = 1;

        lam = Patchv*PatchE/((1+Patchv)*(1-2*Patchv));
        mu = PatchE/(2*(1+Patchv));
        thick = 1;
        fbx = 0;
        fby = 0;
        Dmat = mu*diag([2 2 1]) + lam*[1; 1; 0]*[1 1 0];
        cmat = Dmat;
        Nmat = zeros(2,nst);
        Bmat = zeros(3,nst);
        BBmat = zeros(6,nst);
        I4 = [1 0 0 0 0 1
              0 0 1 0 1 0];
        P1 = [1 0 0 0 0 0
              0 0 0 0 0 1
              0 0 1 1 0 0
              0 0 1 0 0 0
              0 0 0 0 1 0
              0 1 0 0 0 1];
        Z2 = zeros(2);
        
%         if step == 0
        [tau,intb] = Tau9_2d(xl,0,rho0,1,1,nel);
%         else
%         [tau,intb] = Tau9_2d(xl,Dmat,rho0,Nbeta,tstep,nel);
%         end

        % Compute kinematic fields at intermediate time levels
        ul_n_af = ul_n;
        al_n_am = al_n;
        ulres = reshape(ul_n_af,ndf*nel,1);
        alres = reshape(al_n_am,ndf*nel,1);
        
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
            
            Fvec = [fbx; fby];
            Tmat = be(3)*tau*intb;
            
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
            
            diagc = [cmat(1:3,1:3) zeros(3); zeros(3) cmat(1:3,1:3)];
            term8b = I4*diagc*P1; %page 8b
            BBT = term8b;
            D22 = - BBT'*Tmat*BBT;
            ElemTR = Tmat*(-rho*accel + BBT*BBmat*ulres + Fvec(1:2));
            
            hr(nha+(ll-1)*6+5:nha+ll*6) = ElemTR;

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

%         %Set integration number
%         lint = IntPoint3(nel);
%         ib = 0;
%         bf = 0;
%         der = 0;
% 
%         lam = Patchv*PatchE/((1+Patchv)*(1-2*Patchv));
%         mu = PatchE/(2*(1+Patchv));
%         thick = 1;
%         fbx = 0;
%         fby = 0;
%         fbz = 0;
%         Dmat = mu*diag([2 2 2 1 1 1]) + lam*[1; 1; 1; 0; 0; 0]*[1 1 1 0 0 0];
%         Nmat = zeros(3,3*nel);
%         Bmat = zeros(6,3*nel);
% 
%         for l = 1:lint                    
% 
%             % Evaluate 1-D basis functions at integration points
%             if nel == 4 || nel == 10
%               [w,ss] =  int3d_t(l,lint,ib);
%               [shl,shld,shls,be] = shltt(ss,nel,nel,der,bf);
%               [shg,shgs,Jdet,be] = shgtt(xl,nel,shld,shls,nel,bf,der,be);
%             else
%               [w,ss] =  intpntb(l,lint,ib);
%               [shl,shld,shls,be] = shlb(ss,nel,nel,der,bf);
%               [shg,shgs,Jdet,be] = shgb(xl,nel,shld,shls,nel,bf,der,be);
%             end
%             
%             if Jdet < 0
%                 elem
%                 Jdet %#ok<NOPTS>
%             end
%             c1 = Jdet*w*thick;
%             
%             % Form B matrix
%             for ie = 1:nel
%               
%               Bmat(Bcol1,(ie-1)*3+1) = shg(ie,col1);
%               Bmat(Bcol2,(ie-1)*3+2) = shg(ie,col2);
%               Bmat(Bcol3,3*ie      ) = shg(ie,col3);
%                  
%             end
%             
%             epsil = Bmat*reshape(ul,ndf*nel,1);
%             stres = Dmat*epsil;
%             volum = c1;
%             
%             ElemSS(1:6) = ElemSS(1:6) + c1*[1; 1; 1; 1/2; 1/2; 1/2].*epsil;
%             ElemSS(7:12) = ElemSS(7:12) + c1*stres;
%             ElemSS(13) = ElemSS(13) + volum;
% 
%         end %je

    case 52 % Surface strain homogenization
        
        ElemSS = zeros(13,1);
%         t1 = zeros(3,1);
%         t2 = t1;
%         t3 = t1;
%         
%         ElemEL = matepropL(1);
%         ElemvL = matepropL(2);
%         ElemER = matepropR(1);
%         ElemvR = matepropR(2);
%         lamdaR = ElemvR*ElemER/((1+ElemvR)*(1-2*ElemvR));
%         muR = ElemER/(2*(1+ElemvR));
%         lamdaL = ElemvL*ElemEL/((1+ElemvL)*(1-2*ElemvL));
%         muL = ElemEL/(2*(1+ElemvL));
%         DmatR = muR*diag([2 2 2 1 1 1]) + lamdaR*[1; 1; 1; 0; 0; 0]*[1 1 1 0 0 0];
%         DmatL = muL*diag([2 2 2 1 1 1]) + lamdaL*[1; 1; 1; 0; 0; 0]*[1 1 1 0 0 0];
%         
%         NmatL = zeros(3,nstL);
%         BmatL = zeros(6,nstL);
%         bnAdN1 = zeros(6,nstL);
%         N1 = zeros(3,nstL);
%         NmatR = zeros(3,nstR);
%         BmatR = zeros(6,nstR);
%         bnAdN2 = zeros(6,nstR);
%         N2 = zeros(3,nstR);
%         
%         lint = 3;
%         ll = 0; % Counter for history variables
% 
%         nil = surfacesi(1);
% 
%         for intt = 1:nil %integrate on left domain
%                 
%             trinum = surfacesi(intt+1);
%             xit = zeros(ndm,3);
%             for j = 1:3
%                 node = ixt(trinum,j);
%                 for i = 1:ndm
%                     xit(i,j) = Coordinates(node,i);
%                 end
%             end
%         
%             for l = 1:lint
% 
%                 ll = ll + 1;
% 
%                 %Integration point, weight, jacobian
%                 [Wgt,litr,lits] =  intpntt(l,lint,0);
%                 [shl,shld,shls,be] = shlt(litr,lits,3,3,0,0);
%                 [shg, shgs, Jdet, be, xs] = shgt(xit,3,shld,shls,3,0,0,be);
%             
%                 %Physical location of int pt
%                 xint = xit(1,:)*shl;
%                 yint = xit(2,:)*shl;
%                 zint = xit(3,:)*shl;
% 
%                 xi = POU_Coord3(xint,yint,zint,xlL,1,nelL);
%                 rL = xi(1);
%                 sL = xi(2);
%                 tL = xi(3);
% 
%                 % Evaluate  basis functions at integration points
%                 if nelL == 4 || nelL == 10
%                   [shlL,shldL,shls,be] = shltt([rL sL tL],nelL,nel2L,0,0);
%                   PxyL = shgtt(xlL(:,1:nelL),nelL,shldL,shls,nen,0,0,be);
%                 else
%                   [shlL,shldL,shls,be] = shlb([rL sL tL],nelL,nel2L,0,0);
%                   PxyL = shgb(xlL(:,1:nelL),nelL,shldL,shls,nen,0,0,be);
%                 end
%                 QxyL = PxyL;
% 
%                 xi = POU_Coord3(xint,yint,zint,xlR,1,nelR);
%                 rR = xi(1);
%                 sR = xi(2);
%                 tR = xi(3);
% 
%                 % Evaluate  basis functions at integration points
%                 if nelR == 4 || nelR == 10
%                   [shlR,shldR,shls,be] = shltt([rR sR tR],nelR,nel2R,0,0);
%                   PxyR = shgtt(xlR(:,1:nelR),nelR,shldR,shls,nen,0,0,be);
%                 else
%                   [shlR,shldR,shls,be] = shlb([rR sR tR],nelR,nel2R,0,0);
%                   PxyR = shgb(xlR(:,1:nelR),nelR,shldR,shls,nen,0,0,be);
%                 end
%                 QxyR = PxyR;
% 
%                 %Evaluate tangent and normal vectors
%                 t1 = xs(:,1);
%                 [tm1, tu1] = VecNormalize(t1);
%                 t2 = xs(:,2);
%                 [tm2, tu2] = VecNormalize(t2);
%                 t3 = VecCrossProd(t1,t2);
%                 [tm3, tu3] = VecNormalize(t3);
%                 nLx = tu3(1);
%                 nLy = tu3(2);
%                 nLz = tu3(3);
%                 nRx = -tu3(1);
%                 nRy = -tu3(2);
%                 nRz = -tu3(3);
%                 tLx = tu1(1);
%                 tLy = tu1(2);
%                 tLz = tu1(3);
%                 nvect = [nLx 0 0 nLy 0 nLz
%                          0 nLy 0 nLx nLz 0
%                          0 0 nLz 0 nLy nLx]; %- ?
% 
%                 c1 = Wgt*tm3;
%                 
%                 for i = 1:nelL
%                     NmatL(1,(i-1)*3+1) = shlL(i);
%                     NmatL(2,(i-1)*3+2) = shlL(i);
%                     NmatL(3,3*i      ) = shlL(i);
%                     BmatL(Bcol1,(i-1)*3+1) = QxyL(i,col1);
%                     BmatL(Bcol2,(i-1)*3+2) = QxyL(i,col2);
%                     BmatL(Bcol3,3*i      ) = QxyL(i,col3);
%                 end
% 
%                 for i = 1:nelR
%                     NmatR(1,(i-1)*3+1) = shlR(i);
%                     NmatR(2,(i-1)*3+2) = shlR(i);
%                     NmatR(3,3*i      ) = shlR(i);
%                     BmatR(Bcol1,(i-1)*3+1) = QxyR(i,col1);
%                     BmatR(Bcol2,(i-1)*3+2) = QxyR(i,col2);
%                     BmatR(Bcol3,3*i      ) = QxyR(i,col3);
%                 end
%                 
% %                 dam = 1;
%                 nvec = [nLx; nLy; nLz];
%                 jumpu = NmatL*reshape(ulL,ndf*nelL,1) - NmatR*reshape(ulR,ndf*nelR,1);
%                 epsili = -1/2*(nvec*jumpu' + jumpu*nvec');
%                 epsil = [epsili(1,1); epsili(2,2); epsili(3,3); epsili(1,2); epsili(2,3); epsili(3,1)];
%                 
%                 ElemSS(1:6) = ElemSS(1:6) + c1*epsil;
%             
%             end %lint
%         
%         end %intt
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
%         t1 = zeros(3,1);
%         t2 = t1;
%         t3 = t1;
%         
%         dtol = 1e-11;
%         
% %         beta = 1;
% %         beta2 = beta^2;
% %         sigmax = 100;
% %         dc = 0.2;
% %         beta = 0.707;
% %         beta2 = beta^2;
% %         sigmax = 0.01e-3;
% %         dc = 20;
% %         
% %         Hc = sigmax/dc;
% %         rp = 100*Hc;
%         ElemEL = matepropL(1);
%         ElemvL = matepropL(2);
%         ElemER = matepropR(1);
%         ElemvR = matepropR(2);
%         lamdaR = ElemvR*ElemER/((1+ElemvR)*(1-2*ElemvR));
%         muR = ElemER/(2*(1+ElemvR));
%         lamdaL = ElemvL*ElemEL/((1+ElemvL)*(1-2*ElemvL));
%         muL = ElemEL/(2*(1+ElemvL));
%         DmatR = muR*diag([2 2 2 1 1 1]) + lamdaR*[1; 1; 1; 0; 0; 0]*[1 1 1 0 0 0];
%         DmatL = muL*diag([2 2 2 1 1 1]) + lamdaL*[1; 1; 1; 0; 0; 0]*[1 1 1 0 0 0];
%         
%         NmatL = zeros(3,nstL);
%         BmatL = zeros(6,nstL);
%         bnAdN1 = zeros(6,nstL);
%         N1 = zeros(3,nstL);
%         NmatR = zeros(3,nstR);
%         BmatR = zeros(6,nstR);
%         bnAdN2 = zeros(6,nstR);
%         N2 = zeros(3,nstR);
%         
%         lint = 4;3;
%         ll = 0; % Counter for history variables
%         
%             for l = 1:lint
% 
%                 ll = ll + 1;
% 
%                 % Evaluate  basis functions at integration points
%                 if nelL == 4 || nelL == 10
%                   [shlL,shldL,shls,be] = shltt([rL sL tL],nelL,nel2L,0,0);
%                   PxyL = shgtt(xlL(:,1:nelL),nelL,shldL,shls,nen,0,0,be);
%                 else
%                   [Wgt,ss] = intpntb(l,lint,5);
%                   [shlL,shldL,shls,be] = shlb(ss,nelL,nel2L,0,0);
%                   [PxyL,shgs,Jdet,bubble,xs] = shgb(xlL(:,1:nelL),nelL,shldL,shls,nen,0,0,be);
%                 end
%                 QxyL = PxyL;
% 
%                 %Physical location of int pt
%                 xint = xlL(1,:)*shlL;
%                 yint = xlL(2,:)*shlL;
%                 zint = xlL(3,:)*shlL;
% 
%                 xi = POU_Coord3(xint,yint,zint,xlR,1,nelR);
%                 rR = xi(1);
%                 sR = xi(2);
%                 tR = xi(3);
% 
%                 % Evaluate  basis functions at integration points
%                 if nelR == 4 || nelR == 10
%                   [shlR,shldR,shls,be] = shltt([rR sR tR],nelR,nel2R,0,0);
%                   PxyR = shgtt(xlR(:,1:nelR),nelR,shldR,shls,nen,0,0,be);
%                 else
%                   [shlR,shldR,shls,be] = shlb([rR sR tR],nelR,nel2R,0,0);
%                   PxyR = shgb(xlR(:,1:nelR),nelR,shldR,shls,nen,0,0,be);
%                 end
%                 QxyR = PxyR;
% 
%                 %Evaluate tangent and normal vectors
%                 t1 = xs(:,1);
%                 [tm1, tu1] = VecNormalize(t1);
%                 t2 = xs(:,2);
%                 [tm2, tu2] = VecNormalize(t2);
%                 t3 = VecCrossProd(t1,t2);
%                 [tm3, tu3] = VecNormalize(t3);
%                 nLx = -tu3(1);
%                 nLy = -tu3(2);
%                 nLz = -tu3(3);
%                 nRx = tu3(1);
%                 nRy = tu3(2);
%                 nRz = tu3(3);
%                 tLx = tu1(1);
%                 tLy = tu1(2);
%                 tLz = tu1(3);
%                 nvect = [nLx 0 0 nLy 0 nLz
%                          0 nLy 0 nLx nLz 0
%                          0 0 nLz 0 nLy nLx]; %- ?
%                 nvec = [nLx; nLy; nLz];
% 
%                 c1 = Wgt*tm3;
%                 
%             for i = 1:nelL
% %                 NmatL(:,3*i-2:3*i) = shlL(i)*eye(3);
% %                 BmatL(:,3*i-2:3*i) = [QxyL(i,1) 0 0 
% %                                       0 QxyL(i,2) 0 
% %                                       0 0 QxyL(i,3) 
% %                                       QxyL(i,2) QxyL(i,1) 0 
% %                                       0 QxyL(i,3) QxyL(i,2) 
% %                                       QxyL(i,3) 0 QxyL(i,1) ];
%                 NmatL(1,(i-1)*3+1) = shlL(i);
%                 NmatL(2,(i-1)*3+2) = shlL(i);
%                 NmatL(3,3*i      ) = shlL(i);
%                 BmatL(Bcol1,(i-1)*3+1) = QxyL(i,col1);
%                 BmatL(Bcol2,(i-1)*3+2) = QxyL(i,col2);
%                 BmatL(Bcol3,3*i      ) = QxyL(i,col3);
%             end
%             
%             for i = 1:nelR
% %                 NmatR(:,3*i-2:3*i) = shlR(i)*eye(3);
% %                 BmatR(:,3*i-2:3*i) = [QxyR(i,1) 0 0 
% %                                       0 QxyR(i,2) 0 
% %                                       0 0 QxyR(i,3) 
% %                                       QxyR(i,2) QxyR(i,1) 0 
% %                                       0 QxyR(i,3) QxyR(i,2) 
% %                                       QxyR(i,3) 0 QxyR(i,1)];
%                 NmatR(1,(i-1)*3+1) = shlR(i);
%                 NmatR(2,(i-1)*3+2) = shlR(i);
%                 NmatR(3,3*i      ) = shlR(i);
%                 BmatR(Bcol1,(i-1)*3+1) = QxyR(i,col1);
%                 BmatR(Bcol2,(i-1)*3+2) = QxyR(i,col2);
%                 BmatR(Bcol3,3*i      ) = QxyR(i,col3);
%             end
%             
%             % Load history
%             damhr = nh1-1+(ll-1)*7;
%             dmaxhr = nh1-1+(ll-1)*7+4;
%             dvec = [hr(damhr+1); hr(damhr+2); hr(damhr+3)];
%             dmax = hr(dmaxhr);
%                 
%             bnAdN1 = nvect*DmatL*BmatL;
%             bnAdN2 = nvect*DmatR*BmatR;
%             rhspulL = reshape(ulL,ndf*nelL,1);
%             rhspulR = reshape(ulR,ndf*nelR,1);
%             
%             tvtr = 1/2*(bnAdN1*rhspulL + bnAdN2*rhspulR);
%             sn = nvec'*tvtr;
%             jumpu = NmatR*rhspulR - NmatL*rhspulL;
%             un = jumpu'*nvec;
%             tn = sn+rp*un;
%             
%             if tn >= 0 % tension
%                 
%                 tvec = tvtr + rp*jumpu;
%                 
%             else % compression
%                 
%                 tvec = tvtr + rp*jumpu - tn*nvec;
%                 
%             end
%             
%             normtvec = sqrt(tvec'*tvec);
%             if dmax >= dc
%                 psik = 0;
%             else
%                 psik = sigmax - Hc*dmax;
%             end
%             
%             if xint > 0
%                 if yint > 0
%                     theta = atan(yint/xint);
%                 else
%                     theta = 2*pi + atan(yint/xint);
%                 end
%             else
%                 if yint > 0
%                     theta = pi + atan(yint/xint);
%                 else
%                     theta = pi + atan(yint/xint);
%                 end
%             end
%             mvec = [-nvec(2); nvec(1); 0];
%             
% %             ElemI(:,ll) = [theta
% %                            un
% %                            jumpu'*mvec
% %                            sn
% %                            tvtr'*mvec
% %                            tn
% %                            (tvtr + rp*jumpu)'*mvec
% %                            normtvec
% %                            dvec'*nvec
% %                            dvec'*mvec];
%             ElemI(:,ll) = [theta
%                            un
%                            jumpu'*mvec
%                            sn
%                            tvtr'*mvec
%                            (tvtr + rp*(jumpu-dvec))'*nvec
%                            (tvtr + rp*(jumpu-dvec))'*mvec
%                            normtvec
%                            dvec'*nvec
%                            dvec'*mvec];
%             
%             end %lint
%         
% %         end %intt
        
    case 12 % energy
        
        ElemE = 0;
        Nmat = zeros(2,nst);
        Bmat = zeros(3,nst);

        lam = Patchv*PatchE/((1+Patchv)*(1-2*Patchv));
        mu = PatchE/(2*(1+Patchv));
        
        thick = 1;
        
        % Load Guass Integration Points

        lint = IntPoint(nel);
        der = 0;
        bf = 1;
        ib = 0;
                D = [lam+2*mu lam 0 
                     lam lam+2*mu 0
                     0 0  mu];

        %Integration Loop
        for ll = 1:lint

                %Evaluate first derivatives of basis functions at int. point
                if nel == 3 || nel == 6
                  [Wgt,litr,lits] =  intpntt(ll,lint,ib);
                  [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
                  [Qxy, shgs, Jdet,be] = shgt(xl,nel,shld,shls,nel,bf,der,be);
                else
                  [Wgt,litr,lits] =  intpntq(ll,lint,ib);
                  [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
                  [Qxy, shgs, Jdet,be] = shgq(xl,nel,shld,shls,nel,bf,der,be);
                end
                
                for i = 1:nel
                
                  Nmat(1,(ie-1)*2+1) = shl(ie);
                  Nmat(2,(ie-1)*2+2) = shl(ie);
                    Bmat(:,2*i-1:2*i) = [Qxy(i,1) 0
                                         0 Qxy(i,2)
                                         Qxy(i,2) Qxy(i,1)];
                end
                
                %Compute F, sigma, D
                E = Bmat*reshape(ul,nst,1);
                velo = Nmat*reshape(vl,nst,1);
                
                % load fine scale
                beta_n = hr(nhb+(ll-1)*6+1:nhb+(ll-1)*6+2);
                betadot_n = hr(nhb+(ll-1)*6+3:nhb+(ll-1)*6+4);
                E = E + [be(1) 0
                        0 be(2)
                        be(2) be(1)]*beta_n;
                velo = velo + be(3)*betadot_n;
                
                c1 = Wgt*Jdet*thick;
                ElemE = ElemE + c1/2*(E'*D*E + velo'*rho*velo);

        end %je

end