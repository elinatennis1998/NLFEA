
% Tim Truster
% 2/13/2014
%
% Multiscale dynamics using either RFB or polynomial bubble
% a-form of dynamics

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

bubfun = 2;1; % 1=polynomial, 2=RFB

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
        
        if bubfun == 1
        [tau,intb,Mp,Kp,vol] = Tau1_1dad(xl,Dmat,rho,coeffm,coeffk,nel,ndm,nen);
        else
        [tau,intb,Mp,Kp,vol] = Tau1_1drfb(xl,Dmat,rho,coeffm,coeffk,Nbeta,tstep);
        end

        % Compute kinematic fields at intermediate time levels
        ul_n_af = (1 - Nalphaf)*ul + Nalphaf*ul_n;
        al_n_am = (1 - Nalpham)*al + Nalpham*al_n;
        ulres = reshape(ul_n_af,ndf*nel,1);
        alres = reshape(al_n_am,ndf*nel,1);
        
        sw = int1d(lint);
        [shl,shld,shls,be] = shl1d(sw(1,:),lint,nel-1);
        
        if(exist('iprob','var')&&iprob==1&&step==stepmax) && elem == 1
            nind = 0;
        end
        
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
            if adform == 1 %d-form
            beta_n = ElemTR + tau*(1/(2*Nbeta*tstep^2)*Mp*(betaddot_n1*tstep^2 + 2*betadot_n1*tstep + 2*beta_n1) - Mp*betaddot_n1);
            betadot_n = tstep*(1-Ngamma/(2*Nbeta))*betaddot_n1 + (1-Ngamma/Nbeta)*betadot_n1 + Ngamma/(Nbeta*tstep)*(beta_n - beta_n1);
            betaddot_n = (1-1/(2*Nbeta))*betaddot_n1 - 1/(Nbeta*tstep)*betadot_n1 + 1/(Nbeta*tstep^2)*(beta_n - beta_n1);
            else %a-form
            betaddot_n = ElemTR-tau*Kp*(beta_n1+tstep*betadot_n1+tstep^2*((1/2-Nbeta)*betaddot_n1));
            beta_n = beta_n1+tstep*betadot_n1+tstep^2*((1/2-Nbeta)*betaddot_n1+Nbeta*betaddot_n);
            betadot_n = betadot_n1+tstep*((1-Ngamma)*betaddot_n1+Ngamma*betaddot_n);
            end
            
            % store fine scale
            hr(nhb-1+(ll-1)*3+1) = beta_n;
            hr(nhb-1+(ll-1)*3+2) = betadot_n;
            hr(nhb-1+(ll-1)*3+3) = betaddot_n;
            
            if(exist('iprob','var')&&iprob==1&&step==stepmax)
                xint = Nmat*xl';
                nind = nind + 1;
                Sig1(nind,:)=[xint/L sigma/NodeLoadnp(3)];
                
            end
            
            if adform == 1 %d-form
            ElemF = ElemF + c1*(Nmat'*Fvec - rho*Nmat'*accel - Bmat'*sigma ...
                                - Nmat'*rho*bave*betaddot_n);
            
            ElemK = ElemK + c1*(Nmat'*rho/(Nbeta*tstep^2)*Nmat + Bmat'*Dmat*Bmat ...
                                - (-rho/(Nbeta*tstep^2)*Nmat)'*bave*Tmat*(-rho/(Nbeta*tstep^2)*Nmat));
            else
            ElemF = ElemF + c1*(Nmat'*Fvec - rho*Nmat'*accel - Bmat'*sigma ...
                                - Nmat'*rho*bave*betaddot_n);
%                                 - Bmat'*Dmat*bave*beta_n);
            
            ElemK = ElemK + c1*(Nmat'*rho*Nmat + Bmat'*Dmat*(Nbeta*tstep^2)*Bmat ...
                                - (-rho*Nmat)'*bave*Tmat*(-rho*Nmat));
            end

        end %je
        ElemK;
        
    case 5 %Compute Mass
        
        ElemM = zeros(nst);
        
        thick = 1;
        
        % Load Guass Integration Points
        lint = nel;
        der = 0;
        bf = 1;
        ib = 0;
        Nmat = zeros(1,nst);
        
        if bubfun == 1
        [tau,intb,Mp,Kp,vol] = Tau1_1dad(xl,0,rho,1,0,nel,ndm,nen);
        else
        [tau,intb,Mp,Kp,vol] = Tau1_1drfb(xl,Dmat,rho,1,0,Nbeta,tstep);
        end
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
        
        if bubfun == 1
        [tau,intb,Mp,Kp,vol] = Tau1_1dad(xl,Dmat,rho,coeffm,coeffk,nel,ndm,nen);
        else
        [tau,intb,Mp,Kp,vol] = Tau1_1drfb(xl,Dmat,rho,coeffm,coeffk,Nbeta,tstep);
        end

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
            if adform == 1 %d-form
            beta_n = ElemTR + tau*(1/(2*Nbeta*tstep^2)*Mp*(betaddot_n1*tstep^2 + 2*betadot_n1*tstep + 2*beta_n1) - Mp*betaddot_n1);
            betadot_n = tstep*(1-Ngamma/(2*Nbeta))*betaddot_n1 + (1-Ngamma/Nbeta)*betadot_n1 + Ngamma/(Nbeta*tstep)*(beta_n - beta_n1);
            betaddot_n = (1-1/(2*Nbeta))*betaddot_n1 - 1/(Nbeta*tstep)*betadot_n1 + 1/(Nbeta*tstep^2)*(beta_n - beta_n1);
            else %a-form
            betaddot_n = ElemTR-tau*Kp*(beta_n1+tstep*betadot_n1+tstep^2*((1/2-Nbeta)*betaddot_n1));
            beta_n = beta_n1+tstep*betadot_n1+tstep^2*((1/2-Nbeta)*betaddot_n1+Nbeta*betaddot_n);
            betadot_n = betadot_n1+tstep*((1-Ngamma)*betaddot_n1+Ngamma*betaddot_n);
            end
            
            % store fine scale
            hr(nhb-1+(ll-1)*3+1) = beta_n;
            hr(nhb-1+(ll-1)*3+2) = betadot_n;
            hr(nhb-1+(ll-1)*3+3) = betaddot_n;
            
            if adform == 1 %d-form
            ElemF = ElemF + c1*(Nmat'*Fvec - rho*Nmat'*accel - Bmat'*sigma ...
                                - Nmat'*rho*bave*betaddot_n);
            else
            ElemF = ElemF + c1*(Nmat'*Fvec - rho*Nmat'*accel - Bmat'*sigma ...
                                - Nmat'*rho*bave*betaddot_n);
%                                 - Bmat'*Dmat*bave*beta_n);
            end

        end %je

    case 21 %Compute Stiffness
        
        ElemK = zeros(ndf*nel);

        %Set integration number
        lint = nel;
        
        ib = 0;
        bf = 1;
        der = 1;
        
        fb = 0;
        Dmat = Emod;
        Bmat = zeros(1,ndf*nel);
        
        sw = int1d(lint);
        [shl,shld,shls,be] = shl1d(sw(1,:),lint,nel-1);
        
        for ll = 1:lint                    

            % Evaluate 1-D basis functions at integration points
            Wgt = sw(2,ll);
            [shg, shgs, Jdet, be] = shg1d(xl,ndm,nel,shld(ll,:),shls(ll,:),nen,bf,der,be);
            
            c1 = Wgt*Jdet*Aelem;
                
            Bmat(1:ndf:ndf*(nel-1)+1) = shg';
            
            ElemK = ElemK + c1*(Bmat'*Dmat*Bmat);

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
        
        if bubfun == 1
        [tau,intb,Mp,Kp,vol] = Tau1_1dad(xl,0,rho,1,1,nel,ndm,nen);
        else
        [tau,intb,Mp,Kp,vol] = Tau1_1drfb(xl,Dmat,rho,1,0,Nbeta,tstep);
        end
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
        
    case 51 % Volume stress/strain homogenization
        
        ElemSS = zeros(13,1);

    case 52 % Surface strain homogenization
        
        ElemSS = zeros(13,1);
%%        
    case 25 %Stress Projection2

        ElemS = zeros(nel,numstr+1);
        ElemS2 = zeros(nel,numstr);

        Dmat = Emod;
        Bmat = zeros(1,ndf*nel);
        
        % Load Guass Integration Points

        if nel == 2
            lint = 2;
            nint = 1;
        elseif nel == 3
            lint = 3;
            nint = 2;
        else
            lint = 3;
            nint = 2;
        end
        
        der = 0;
        bf = 0;
        ib = 0;
        
        ulres = reshape(ul,ndf*nel,1);

        %Stress Loop
        
        sw = int1d(nint);
        [shl,shld,shls,be] = shl1d(sw(1,:),lint,nel-1);
        
        
        for ll = 1:nint                    

            % Evaluate 1-D basis functions at integration points
            Wgt = sw(2,ll);
            [shg, shgs, Jdet, be] = shg1d(xl,ndm,nel,shld(ll,:),shls(ll,:),nen,bf,der,be);
            
            % Form B matrix
            Bmat(1:ndf:ndf*(nel-1)+1) = shg';
            
            epsil = Bmat*ulres;
            sigma = Dmat*epsil;
            
            ElemS2(ll,1) = sigma;

        end %je
        
        % interpolate stress at nodes
        if nel == 2
            plist = [-1 1];
        elseif nel == 3
            sqr3 = sqrt(3);
            plist = [-sqr3 0 sqr3];
        else
            sqr3 = sqrt(3);
            plist = [-sqr3 0 sqr3];
        end
        
        for ll = 1:lint
            
            r = plist(1,ll);
            shpS = shl1d(r,1,nint-1); % polynomial order p = nint - 1
            
%             for stres = 1:numstr
%                 
%                 sigmas = ElemS2(1:nint,stres)'*shpS;
%                 ElemS(ll,stres) = sigmas;
%                 
%             end
            ElemS(ll,1:numstr) = (shpS*ElemS2(1:nint,:))';
            
        end
        
%         %Integration Loop
%         Vol = xl(2) - xl(1);

        for i = 1:nel
        ElemS(i,numstr+1) = 1; % use this for simple average Vol; % use this for weighted average 
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
        
        if bubfun == 1
        [tau,intb,Mp,Kp,vol] = Tau1_1dad(xl,Dmat,rho,coeffm,coeffk,nel,ndm,nen);
        else
        [tau,intb,Mp,Kp,vol] = Tau1_1drfb(xl,Dmat,rho,coeffm,coeffk,Nbeta,tstep);
        end
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