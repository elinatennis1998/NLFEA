
% Tim Truster
% 6/24/2013
% UIUC
% Nonlinear multiscale dynamics
% First attempt; use with transient=3 for energy constraint method
% Energy is conserved with that method; example is IMSDYN_NL.m


% Set Material Properties

Patcha = mateprop(2);
Patchb = mateprop(3);
Patchc = mateprop(4);
rho0 = mateprop(5);

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
        
        nh1 = 9*2*3;
        
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

        thick = 1;
        fbx = 0;
        fby = 0;
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
        dmat = [8 0 0 0
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
                0 0 0 0]*Patchc;

        % Compute kinematic fields at intermediate time levels
        ul_n_af = (1 - Nalphaf)*ul + Nalphaf*ul_n;
        al_n_am = (1 - Nalpham)*al + Nalpham*al_n;
        ulres = reshape(ul_n_af,ndf*nel,1);
        alres = reshape(al_n_am,ndf*nel,1);
        
        [tau,intb,Mp,Kp,vol] = Tau10_2d(xl,ulres,Patcha,Patchb,Patchc,rho0,Nbeta,tstep,nel);
        
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
            E = (Bmat*ulres);
            E(3) = E(3)/2; %for extra factor of 2 on shear strain
            Ekk = E(1) + E(2);
            E2 = [E(1) E(3); E(3) E(2)];
            E2 = E2*E2;
            E2 = [E2(1,1); E2(2,2); E2(1,2)];
            sigma = Patcha*Ekk*[1; 1; 0] + Patchb*E + Patchc*E2;
            
            cmat = [Patcha+Patchb+2*Patchc*E(1) Patcha Patchc*E(3) 
                    Patcha Patcha+Patchb+2*Patchc*E(2) Patchc*E(3)
                    Patchc*E(3) Patchc*E(3)  Patchb/2+Patchc/2*Ekk];
            
            Fvec = [fbx; fby];
            Tmat = tau*intb;
            bave = intb/vol;
          
%             % load fine scale
%             beta_n1 = hr(nha-1+(ll-1)*6+1:nha-1+(ll-1)*6+2);
%             betadot_n1 = hr(nha-1+(ll-1)*6+3:nha-1+(ll-1)*6+4);
%             betaddot_n1 = hr(nha-1+(ll-1)*6+5:nha-1+(ll-1)*6+6);
%             
%             diagc = [cmat(1:3,1:3) zeros(3); zeros(3) cmat(1:3,1:3)];
%             term8b = I4*diagc*P1; %page 8b
%             BBT = term8b;
%             ElemTR = Tmat*(-rho*accel + BBT*BBmat*ulres + Fvec(1:2));
%             
%             % update fine scales
%             beta_n = ElemTR + tau*(1/(2*Nbeta*tstep^2)*Mp*(betaddot_n1*tstep^2 + 2*betadot_n1*tstep + 2*beta_n1) - Mp*betaddot_n1);
%             betadot_n = tstep*(1-Ngamma/(2*Nbeta))*betaddot_n1 + (1-Ngamma/Nbeta)*betadot_n1 + Ngamma/(Nbeta*tstep)*(beta_n - beta_n1);
%             betaddot_n = (1-1/(2*Nbeta))*betaddot_n1 - 1/(Nbeta*tstep)*betadot_n1 + 1/(Nbeta*tstep^2)*(beta_n - beta_n1);
%             
%             % store fine scale
%             hr(nhb-1+(ll-1)*6+1:nhb-1+(ll-1)*6+2) = beta_n;
%             hr(nhb-1+(ll-1)*6+3:nhb-1+(ll-1)*6+4) = betadot_n;
%             hr(nhb-1+(ll-1)*6+5:nhb-1+(ll-1)*6+6) = betaddot_n;
            
            ElemF = ElemF + c1*(Nmat'*Fvec - rho*Nmat'*accel - Bmat'*sigma);% ...
%                                 + BBmat'*BBT'*bave*beta_n - Nmat'*rho*bave*betaddot_n);
            
%             ElemK = ElemK + c1*(Nmat'*rho/(Nbeta*tstep^2)*Nmat + Bmat'*cmat*Bmat);% ...
            ElemK = ElemK + c1*(Bmat'*cmat*Bmat);% ...
%                                 - (-rho/(Nbeta*tstep^2)*Nmat + BBT*BBmat)'*bave*Tmat*(-rho/(Nbeta*tstep^2)*Nmat + BBT*BBmat));

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
        
        ulres = reshape(ul,ndf*nel,1);
        [tau,intb,Mp,Kp,vol] = Tau10_2d(xl,ulres,0,0,0,rho0,1,1,nel);
        Tmat = 0*intb^2/vol*tau;
        
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
            
            ElemM = ElemM + c1*rho*(Nmat'*Nmat - rho*(Nmat'*Tmat*Nmat));
%             ElemM = ElemM + c1*rho*(Nmat'*Nmat);

        end %je
        ElemM;

    case 6 %Compute Residual
        
        ElemF = zeros(nst,1);
        ElemFn = zeros(nst,1);

        %Set integration number
        lint = IntPoint(nel);
        ib = 0;
        bf = 1;
        der = 1;

        thick = 1;
        fbx = 0;
        fby = 0;
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
        dmat = [8 0 0 0
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
                0 0 0 0]*Patchc;

        % Compute kinematic fields at intermediate time levels
        ul_n_af = (1 - Nalphaf)*ul + Nalphaf*ul_n;
        al_n_am = (1 - Nalpham)*al + Nalpham*al_n;
        ulres = reshape(ul_n_af,ndf*nel,1);
        alres = reshape(al_n_am,ndf*nel,1);
        
        [tau,intb,Mp,Kp,vol] = Tau10_2d(xl,ulres,Patcha,Patchb,Patchc,rho0,Nbeta,tstep,nel);
        
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
            E = (Bmat*ulres);
            E(3) = E(3)/2; %for extra factor of 2 on shear strain
            Ekk = E(1) + E(2);
            E2 = [E(1) E(3); E(3) E(2)];
            E2 = E2*E2;
            E2 = [E2(1,1); E2(2,2); E2(1,2)];
            sigma = Patcha*Ekk*[1; 1; 0] + Patchb*E + Patchc*E2;
            
            cmat = [Patcha+Patchb+2*Patchc*E(1) Patcha Patchc*E(3) 
                    Patcha Patcha+Patchb+2*Patchc*E(2) Patchc*E(3)
                    Patchc*E(3) Patchc*E(3)  Patchb/2+Patchc/2*Ekk];
            
            Fvec = [fbx; fby];
            Tmat = tau*intb;
            bave = intb/vol;
          
%             % load fine scale
%             beta_n1 = hr(nha-1+(ll-1)*6+1:nha-1+(ll-1)*6+2);
%             betadot_n1 = hr(nha-1+(ll-1)*6+3:nha-1+(ll-1)*6+4);
%             betaddot_n1 = hr(nha-1+(ll-1)*6+5:nha-1+(ll-1)*6+6);
%             
%             diagc = [cmat(1:3,1:3) zeros(3); zeros(3) cmat(1:3,1:3)];
%             term8b = I4*diagc*P1; %page 8b
%             BBT = term8b;
%             ElemTR = Tmat*(-rho*accel + BBT*BBmat*ulres + Fvec(1:2));
%             
%             % update fine scales
%             beta_n = ElemTR + tau*(1/(2*Nbeta*tstep^2)*Mp*(betaddot_n1*tstep^2 + 2*betadot_n1*tstep + 2*beta_n1) - Mp*betaddot_n1);
%             betadot_n = tstep*(1-Ngamma/(2*Nbeta))*betaddot_n1 + (1-Ngamma/Nbeta)*betadot_n1 + Ngamma/(Nbeta*tstep)*(beta_n - beta_n1);
%             betaddot_n = (1-1/(2*Nbeta))*betaddot_n1 - 1/(Nbeta*tstep)*betadot_n1 + 1/(Nbeta*tstep^2)*(beta_n - beta_n1);
%             
%             % store fine scale
%             hr(nhb-1+(ll-1)*6+1:nhb-1+(ll-1)*6+2) = beta_n;
%             hr(nhb-1+(ll-1)*6+3:nhb-1+(ll-1)*6+4) = betadot_n;
%             hr(nhb-1+(ll-1)*6+5:nhb-1+(ll-1)*6+6) = betaddot_n;
            
            ElemF = ElemF + c1*(Nmat'*Fvec - rho*Nmat'*accel - Bmat'*sigma);% ...
%                                 + BBmat'*BBT'*bave*beta_n - Nmat'*rho*bave*betaddot_n);

        end %je

    case 21 %Compute Stiffness
        
        ElemK = zeros(nst);

        %Set integration number
        lint = IntPoint(nel);
        ib = 0;
        bf = 1;
        der = 1;

        thick = 1;
        fbx = 0;
        fby = 0;
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
        dmat = [8 0 0 0
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
                0 0 0 0]*Patchc;

        % Compute kinematic fields at intermediate time levels
        ul_n_af = (1 - Nalphaf)*ul + Nalphaf*ul_n;
        al_n_am = (1 - Nalpham)*al + Nalpham*al_n;
        ulres = reshape(ul_n_af,ndf*nel,1);
        alres = reshape(al_n_am,ndf*nel,1);
        
        [tau,intb,Mp,Kp,vol] = Tau10_2d(xl,ulres,Patcha,Patchb,Patchc,rho0,Nbeta,tstep,nel);
        
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
            E = (Bmat*ulres);
            E(3) = E(3)/2; %for extra factor of 2 on shear strain
            Ekk = E(1) + E(2);
            E2 = [E(1) E(3); E(3) E(2)];
            E2 = E2*E2;
            E2 = [E2(1,1); E2(2,2); E2(1,2)];
            sigma = Patcha*Ekk*[1; 1; 0] + Patchb*E + Patchc*E2;
            
            cmat = [Patcha+Patchb+2*Patchc*E(1) Patcha Patchc*E(3) 
                    Patcha Patcha+Patchb+2*Patchc*E(2) Patchc*E(3)
                    Patchc*E(3) Patchc*E(3)  Patchb/2+Patchc/2*Ekk];
            
            Fvec = [fbx; fby];
            Tmat = tau*intb;
            bave = intb/vol;
          
%             % load fine scale
%             beta_n1 = hr(nha-1+(ll-1)*6+1:nha-1+(ll-1)*6+2);
%             betadot_n1 = hr(nha-1+(ll-1)*6+3:nha-1+(ll-1)*6+4);
%             betaddot_n1 = hr(nha-1+(ll-1)*6+5:nha-1+(ll-1)*6+6);
%             
%             diagc = [cmat(1:3,1:3) zeros(3); zeros(3) cmat(1:3,1:3)];
%             term8b = I4*diagc*P1; %page 8b
%             BBT = term8b;
%             ElemTR = Tmat*(-rho*accel + BBT*BBmat*ulres + Fvec(1:2));
%             
%             % update fine scales
%             beta_n = ElemTR + tau*(1/(2*Nbeta*tstep^2)*Mp*(betaddot_n1*tstep^2 + 2*betadot_n1*tstep + 2*beta_n1) - Mp*betaddot_n1);
%             betadot_n = tstep*(1-Ngamma/(2*Nbeta))*betaddot_n1 + (1-Ngamma/Nbeta)*betadot_n1 + Ngamma/(Nbeta*tstep)*(beta_n - beta_n1);
%             betaddot_n = (1-1/(2*Nbeta))*betaddot_n1 - 1/(Nbeta*tstep)*betadot_n1 + 1/(Nbeta*tstep^2)*(beta_n - beta_n1);
%             
%             % store fine scale
%             hr(nhb-1+(ll-1)*6+1:nhb-1+(ll-1)*6+2) = beta_n;
%             hr(nhb-1+(ll-1)*6+3:nhb-1+(ll-1)*6+4) = betadot_n;
%             hr(nhb-1+(ll-1)*6+5:nhb-1+(ll-1)*6+6) = betaddot_n;
            
%             ElemK = ElemK + c1*(Nmat'*rho/(Nbeta*tstep^2)*Nmat + Bmat'*cmat*Bmat);% ...
            ElemK = ElemK + c1*(Bmat'*cmat*Bmat);% ...
%                                 - (-rho/(Nbeta*tstep^2)*Nmat + BBT*BBmat)'*bave*Tmat*(-rho/(Nbeta*tstep^2)*Nmat + BBT*BBmat));

        end %je
        
    case 40 % Initialize FS acceleration

%         %Set integration number
%         lint = IntPoint(nel);
%         ib = 0;
%         bf = 1;
%         der = 1;
% 
%         thick = 1;
%         fbx = 0;
%         fby = 0;
%         Dmat = mu*diag([2 2 1]) + lam*[1; 1; 0]*[1 1 0];
%         Nmat = zeros(2,nst);
%         Bmat = zeros(3,nst);
%         BBmat = zeros(6,nst);
%         cmat = Dmat;
%         I4 = [1 0 0 0 0 1
%               0 0 1 0 1 0];
%         P1 = [1 0 0 0 0 0
%               0 0 0 0 0 1
%               0 0 1 1 0 0
%               0 0 1 0 0 0
%               0 0 0 0 1 0
%               0 1 0 0 0 1];
% 
%         % Compute kinematic fields at intermediate time levels
%         ul_n_af = ul_n;
%         al_n_am = al_n;
%         ulres = reshape(ul_n_af,ndf*nel,1);
%         alres = reshape(al_n_am,ndf*nel,1);
%         
%         [tau,intb,Mp,Kp,vol] = Tau10_2d(xl,ulres,0,0,0,rho0,1,1,nel);
%         Tmat = intb*tau;
%         
%         for ll = 1:lint                    
% 
%             % Evaluate 1-D basis functions at integration points
%             if nel == 3 || nel == 6
%               [w,litr,lits] =  intpntt(ll,lint,ib);
%               [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
%               [shg,shgs,Jdet,be] = shgt(xl,nel,shld,shls,nel,bf,der,be);
%             else
%               [w,litr,lits] =  intpntq(ll,lint,ib);
%               [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
%               [shg,shgs,Jdet,be] = shgq(xl,nel,shld,shls,nel,bf,der,be);
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
%               Nmat(1,(ie-1)*2+1) = shl(ie);
%               Nmat(2,(ie-1)*2+2) = shl(ie);
%                 
%               Bmat(Bcol1,(ie-1)*2+1) = shg(ie,col1);
%               Bmat(Bcol2,(ie-1)*2+2) = shg(ie,col2);
%               
%               BBmat(:,(ie-1)*2+1:2*ie) = [shgs(ie,1) 0               
%                      shgs(ie,2) 0             
%                      shgs(ie,3) 0               
%                      0         shgs(ie,1)       
%                      0         shgs(ie,2)       
%                      0         shgs(ie,3)];
%                  
%             end
% 
%             %Compute rho, accel
%             rho = rho0; %rho = f(x,y);
%             accel = al_n_am*shl;
%             E = (Bmat*ulres);
%             E(3) = E(3)/2; %for extra factor of 2 on shear strain
%             Ekk = E(1) + E(2);
%             E2 = [E(1) E(3); E(3) E(2)];
%             E2 = E2*E2;
%             E2 = [E2(1,1); E2(2,2); E2(1,2)];
%             sigma = Patcha*Ekk*[1; 1; 0] + Patchb*E + Patchc*E2;
%             
%             cmat = [Patcha+Patchb+2*Patchc*E(1) Patcha Patchc*E(3) 
%                     Patcha Patcha+Patchb+2*Patchc*E(2) Patchc*E(3)
%                     Patchc*E(3) Patchc*E(3)  Patchb/2+Patchc/2*Ekk];
%             
%             Fvec = [fbx; fby];
%             Tmat = tau*intb;
%             bave = intb/vol;
% 
%             %Compute rho, accel
%             rho = rho0; %rho = f(x,y);
%             accel = al_n_am*shl;
%             
%             diagc = [cmat(1:3,1:3) zeros(3); zeros(3) cmat(1:3,1:3)];
%             term8b = I4*diagc*P1; %page 8b
%             BBT = term8b;
%             ElemTR = Tmat*(-rho*accel + BBT*BBmat*ulres + Fvec(1:2));
%             
%             hr(nha-1+(ll-1)*6+5:nha-1+ll*6) = ElemTR;
% 
%         end %je

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
        
    case 60
        
        numhr = 4;
        ElemI = zeros(10,numhr);
        
    case 12 % energy
        
        ElemE = 0;
        ElemE1 = 0;
        ElemE2 = 0;
        
        thick = 1;
        
        % Load Guass Integration Points

        lint = IntPoint(nel);
        der = 1;
        bf = 1;
        ib = 0;
        
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
        
        [tau,intb,Mp,Kp,vol] = Tau10_2d(xl,ulres,Patcha,Patchb,Patchc,rho0,Nbeta,tstep,nel);
        bave = intb/vol;

        %Integration Loop
        for ll = 1:lint

                %Evaluate first derivatives of basis functions at int. point
                if nel == 3 || nel == 6
                  [Wgt,litr,lits] =  intpntt(ll,lint,ib);
                  [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
                  [shg, shgs, Jdet,be] = shgt(xl,nel,shld,shls,nel,bf,der,be);
                else
                  [Wgt,litr,lits] =  intpntq(ll,lint,ib);
                  [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
                  [shg, shgs, Jdet,be] = shgq(xl,nel,shld,shls,nel,bf,der,be);
                end
                
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
                
                %Compute F, sigma, D
                E = Bmat*reshape(ul,nst,1);
                velo = Nmat*reshape(vl,nst,1);
                E(3) = E(3)/2; %for extra factor of 2 on shear strain
                Ekk = E(1) + E(2);
                E2 = [E(1) E(3); E(3) E(2)];
                
                U = 1/2*Patcha*Ekk^2;
                for i = 1:2
                    for j = 1:2
                        U = U + 1/2*Patchb*E2(i,j)*E2(j,i);
                        for k = 1:2
                            U = U + 1/3*Patchc*E2(i,j)*E2(j,k)*E2(k,i);
                        end
                    end
                end
            
                cmat = [Patcha+Patchb+2*Patchc*E(1) Patcha Patchc*E(3) 
                        Patcha Patcha+Patchb+2*Patchc*E(2) Patchc*E(3)
                        Patchc*E(3) Patchc*E(3)  Patchb/2+Patchc/2*Ekk];
                
                diagc = [cmat(1:3,1:3) zeros(3); zeros(3) cmat(1:3,1:3)];
                term8b = I4*diagc*P1; %page 8b
                BBT = term8b;
                divsig = BBT*BBmat*reshape(ul,nst,1);
                
                % load fine scale
                beta_n = hr(nhb-1+(ll-1)*6+1:nhb-1+(ll-1)*6+2);
                betadot_n = hr(nhb-1+(ll-1)*6+3:nhb-1+(ll-1)*6+4);
                
                c1 = Wgt*Jdet*thick;
                ElemE = ElemE + c1*U;% + c1/2*(velo'*rho*velo);
                ElemE1 = ElemE1 + c1*(betadot_n'*bave*rho*velo - beta_n'*bave*divsig);
                ElemE2 = ElemE2 + c1/2*(betadot_n'*Mp/vol*betadot_n + beta_n'*Kp/vol*beta_n);
%                               + 1/2*(betadot_n'*intb*rho*velo - beta_n'*intb*divsig + betadot_n'*Mp*betadot_n + beta_n'*Kp*beta_n);
%                 ElemE = ElemE + c1/2*(E'*D*E + velo'*rho*velo + betadot_n'*intb/vol*rho*velo - beta_n'*intb/vol*divsig) ...
%                               + c1/2*(betadot_n'*intb*rho*velo - beta_n'*intb*divsig + betadot_n'*Mp*betadot_n + beta_n'*Kp*beta_n);
% %                               + 1/2*(betadot_n'*intb*rho*velo - beta_n'*intb*divsig + betadot_n'*Mp*betadot_n + beta_n'*Kp*beta_n);

        end %je
ElemE = ElemE + ElemE1 + ElemE2;

end