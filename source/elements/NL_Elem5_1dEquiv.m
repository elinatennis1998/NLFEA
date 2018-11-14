% Tim Truster
% 05/03/2014
%
% Implementation of the VMS-stabilized wave equation with piece-wise
% constant in time FS acceleration. Method is coded such that it is
% equivalent to the mass-averaging schemes that lead to a higher-order
% mass. The stabilization parameter tau is provided in the input file as a
% material property. The final form of the equations is:
%
% M' * a_n+1 + beta*tstep^2 * K' * a_n+1 = F + F_tractions + tauF - K' * u_predictor
%
% Verified 5/3/14 using DynBarImpactVMS and DynBarImpactVMS2; the first
% input file is for linear elements and was benchmarked against NL_Elem5_1d
% by setting ModelAx_1 to zero at the beginning, because the isw=5 mass
% matrix is different between these two subroutines.
%
% Note: tau was derived considering element length = 2h, so that CFL needs
% to be defined as CFL = cwave*p*dt/h where p is the polynomial order.
% This makes a difference in the computations because M' = M +
% tau/(1-tau)*beta*dt^2*K, so the relation between tau and dt is IMPORTANT.

% Set Material Properties

cwave = mateprop(1);
tau = mateprop(2);

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
        
        iste = 2*(nen-1);
        
%%
    case 3 %Compute Stiffness and Residual
        
        ElemK = zeros(ndf*nel);
%         ElemF = zeros(ndf*nel,1);
        ElemM = ElemK;
        ElemLLT = ElemK;
        ElemT = ElemK;

        %Set integration number
        lint = nel;
        
        ib = 0;
        bf = 1;
        der = 1;
        
        fb = 0;
        Dmat = cwave^2;
        Tmat = cwave^4;
        Nmat = zeros(1,nst);
        Bmat = zeros(1,nst);
        BBmat = zeros(1,nst);
        
%         if t_on
%             % tau from Harari
%         tau = 1+6*Nbeta*CFL^2*(1-cosh(1/sqrt(Nbeta)/CFL))/(2+cosh(1/sqrt(Nbeta)/CFL));
%         else
%         tau = 0;
%         end

        % Compute kinematic fields at intermediate time levels
        ul_n_af = (1 - Nalphaf)*ul + Nalphaf*ul_n;
        al_n_am = (1 - Nalpham)*al + Nalpham*al_n;
        ulres = reshape(ul_n_af,ndf*nel,1);
        alres = reshape(al_n_am,ndf*nel,1);
        
%         ul_pred = ul_n + vl_n*tstep + tstep^2/2*(1-2*Nbeta)*al_n;
%         ulpres = reshape(ul_pred,ndf*nel,1);
        
        sw = int1d(lint);
        [shl,shld,shls,be] = shl1d(sw(1,:),lint,nel-1);
        
        if(exist('iprob','var')&&iprob==1&&step==stepmax) && elem == 1
            nind = 0;
        end
            
        for ll = 1:lint                    

            % Evaluate 1-D basis functions at integration points
            Wgt = sw(2,ll);
            [shg, shgs, Jdet, be] = shg1d(xl,ndm,nel,shld(ll,:),shls(ll,:),nen,bf,der,be);
            
            c1 = Wgt*Jdet;
            
            Nmat(1:ndf:ndf*(nel-1)+1) = shl(ll,:);
                
            Bmat(1:ndf:ndf*(nel-1)+1) = shg';
            
            BBmat(1:ndf:ndf*(nel-1)+1) = shgs';

            %Compute rho, accel
            accel = Nmat*alres;
            sigma = Dmat*(Bmat*ulres);
            
            Fvec = fb;
            
%             ElemF = ElemF + c1*(Nmat'*Fvec - (1-tau)*Bmat'*sigma);
            
            ElemK = ElemK + c1*Bmat'*Dmat*Bmat;
            ElemM = ElemM + c1*(Nmat'*Nmat);
            ElemLLT = ElemLLT + c1*(BBmat'*Dmat*Nmat + Nmat'*Dmat*BBmat);
            ElemT = ElemT + c1*BBmat'*Tmat*BBmat;
%             
%             ElemK = ElemK + c1*(Nmat'*Nmat + Bmat'*(Nbeta*tstep^2)*Dmat*Bmat ...
%                                - (Nmat - (Nbeta*tstep^2)*Dmat*BBmat)'*tau*(Nmat - (Nbeta*tstep^2)*Dmat*BBmat));
            
            if(exist('iprob','var')&&iprob==1&&step==stepmax)
                nind = nind + 1;
                xint = Nmat*xl';
                Sig1(nind,:)=[xint/L sigma/NodeLoadnp(3)];
                
            end

        end %je
        
        % Final formulas
        if tau > - 2
        ElemM = ElemM + tau/(1-tau)*Nbeta*tstep^2*(ElemLLT + ElemK); % combined mass
        ElemK = ElemK - tau/(1-tau)*Nbeta*tstep^2*ElemT; % combined stiffness
        ElemF = -ElemM*alres - ElemK*ulres;
        ElemK = ElemM + Nbeta*tstep^2*ElemK;
        else % do mass averaging from Fried et al.
        h = xl(1,nel) - xl(1,1);
        ElemM = ElemM + ElemT*h^3/16*h/15/2*(tau + 3);
        ElemF = -ElemM*alres - ElemK*ulres;
        ElemK = ElemM + Nbeta*tstep^2*ElemK;
        end
        
    case 5 %Compute Mass
        
        ElemK = zeros(ndf*nel);
        ElemF = zeros(ndf*nel,1);
        ElemM = ElemK;
        ElemLLT = ElemK;
        ElemT = ElemK;

        %Set integration number
        lint = nel;
        
        ib = 0;
        bf = 1;
        der = 1;
        
        fb = 0;
        Dmat = cwave^2;
        Tmat = cwave^4;
        Nmat = zeros(1,nst);
        Bmat = zeros(1,nst);
        BBmat = zeros(1,nst);
        
        sw = int1d(lint);
        [shl,shld,shls,be] = shl1d(sw(1,:),lint,nel-1);
            
        for ll = 1:lint                    

            % Evaluate 1-D basis functions at integration points
            Wgt = sw(2,ll);
            [shg, shgs, Jdet, be] = shg1d(xl,ndm,nel,shld(ll,:),shls(ll,:),nen,bf,der,be);
            
            c1 = Wgt*Jdet;
            
            Nmat(1:ndf:ndf*(nel-1)+1) = shl(ll,:);
                
            Bmat(1:ndf:ndf*(nel-1)+1) = shg';
            
            BBmat(1:ndf:ndf*(nel-1)+1) = shgs';
            
            ElemK = ElemK + c1*Bmat'*Dmat*Bmat;
            ElemM = ElemM + c1*(Nmat'*Nmat);
            ElemLLT = ElemLLT + c1*(BBmat'*Dmat*Nmat + Nmat'*Dmat*BBmat);
            ElemT = ElemT + c1*BBmat'*Tmat*BBmat;

        end %je
        
        % Final formulas
        if tau > - 2
        ElemM = ElemM + tau/(1-tau)*Nbeta*tstep^2*(ElemLLT + ElemK); % combined mass
        else % do mass averaging from Fried et al.
        h = xl(1,nel) - xl(1,1);
        ElemM = ElemM + ElemT*h^3/16*h/15/2*(tau + 3);
        end

    case 6 %Compute Residual
        
        ElemK = zeros(ndf*nel);
        ElemF = zeros(ndf*nel,1);
        ElemM = ElemK;
        ElemLLT = ElemK;
        ElemT = ElemK;

        %Set integration number
        lint = nel;
        
        ib = 0;
        bf = 1;
        der = 1;
        
        fb = 0;
        Dmat = cwave^2;
        Tmat = cwave^4;
        Nmat = zeros(1,nst);
        Bmat = zeros(1,nst);
        BBmat = zeros(1,nst);

        % Compute kinematic fields at intermediate time levels
        ul_n_af = (1 - Nalphaf)*ul + Nalphaf*ul_n;
        al_n_am = (1 - Nalpham)*al + Nalpham*al_n;
        ulres = reshape(ul_n_af,ndf*nel,1);
        alres = reshape(al_n_am,ndf*nel,1);
        
%         ul_pred = ul_n + vl_n*tstep + tstep^2/2*(1-2*Nbeta)*al_n;
%         ulpres = reshape(ul_pred,ndf*nel,1);
        
        sw = int1d(lint);
        [shl,shld,shls,be] = shl1d(sw(1,:),lint,nel-1);
        
            if(step==0.6*L/cwave) && elem == 1
                nind = 0;
            end
            
        for ll = 1:lint                    

            % Evaluate 1-D basis functions at integration points
            Wgt = sw(2,ll);
            [shg, shgs, Jdet, be] = shg1d(xl,ndm,nel,shld(ll,:),shls(ll,:),nen,bf,der,be);
            
            c1 = Wgt*Jdet;
            
            Nmat(1:ndf:ndf*(nel-1)+1) = shl(ll,:);
                
            Bmat(1:ndf:ndf*(nel-1)+1) = shg';
            
            BBmat(1:ndf:ndf*(nel-1)+1) = shgs';

            %Compute rho, accel
            accel = Nmat*alres;
            sigma = Dmat*(Bmat*ulres);
            
            Fvec = fb;
            
            ElemK = ElemK + c1*Bmat'*Dmat*Bmat;
            ElemM = ElemM + c1*(Nmat'*Nmat);
            ElemLLT = ElemLLT + c1*(BBmat'*Dmat*Nmat + Nmat'*Dmat*BBmat);
            ElemT = ElemT + c1*BBmat'*Tmat*BBmat;

        end %je
        
        % Final formulas
        if tau > - 2
        ElemM = ElemM + tau/(1-tau)*Nbeta*tstep^2*(ElemLLT + ElemK); % combined mass
        ElemK = ElemK - tau/(1-tau)*Nbeta*tstep^2*ElemT; % combined stiffness
        ElemF = -ElemM*alres - ElemK*ulres;
        else % do mass averaging from Fried et al.
        h = xl(1,nel) - xl(1,1);
        ElemM = ElemM + ElemT*h^3/16*h/15/2*(tau + 3);
        ElemF = -ElemM*alres - ElemK*ulres;
        end

    case 21 %Compute Stiffness
        
        ElemK = zeros(ndf*nel);
        ElemM = ElemK;
        ElemLLT = ElemK;
        ElemT = ElemK;

        %Set integration number
        lint = nel;
        
        ib = 0;
        bf = 1;
        der = 1;
        
        fb = 0;
        Dmat = cwave^2;
        Tmat = cwave^4;
        Nmat = zeros(1,nst);
        Bmat = zeros(1,nst);
        BBmat = zeros(1,nst);
        
        sw = int1d(lint);
        [shl,shld,shls,be] = shl1d(sw(1,:),lint,nel-1);
            
        for ll = 1:lint                    

            % Evaluate 1-D basis functions at integration points
            Wgt = sw(2,ll);
            [shg, shgs, Jdet, be] = shg1d(xl,ndm,nel,shld(ll,:),shls(ll,:),nen,bf,der,be);
            
            c1 = Wgt*Jdet;
            
            Nmat(1:ndf:ndf*(nel-1)+1) = shl(ll,:);
                
            Bmat(1:ndf:ndf*(nel-1)+1) = shg';
            
            BBmat(1:ndf:ndf*(nel-1)+1) = shgs';
            
            ElemK = ElemK + c1*Bmat'*Dmat*Bmat;
            ElemM = ElemM + c1*(Nmat'*Nmat);
            ElemLLT = ElemLLT + c1*(BBmat'*Dmat*Nmat + Nmat'*Dmat*BBmat);
            ElemT = ElemT + c1*BBmat'*Tmat*BBmat;

        end %je
        
        % Final formulas
        if tau > - 2
        ElemK = ElemK - tau/(1-tau)*Nbeta*tstep^2*ElemT; % combined stiffness
        else % do mass averaging from Fried et al.
        h = xl(1,nel) - xl(1,1);
        end
        
    case 40 % Initialize FS acceleration
        
        
    case 12 % energy
        
        ElemK = zeros(ndf*nel);
        ElemM = ElemK;
        ElemLLT = ElemK;
        ElemT = ElemK;

        %Set integration number
        lint = nel;
        
        ib = 0;
        bf = 1;
        der = 1;
        
        fb = 0;
        Dmat = cwave^2;
        Tmat = cwave^4;
        Nmat = zeros(1,nst);
        Bmat = zeros(1,nst);
        BBmat = zeros(1,nst);

        % Compute kinematic fields at intermediate time levels
        ul_n_af = (1 - Nalphaf)*ul + Nalphaf*ul_n;
        al_n_am = (1 - Nalpham)*al + Nalpham*al_n;
        ulres = reshape(ul_n_af,ndf*nel,1);
        alres = reshape(al_n_am,ndf*nel,1);
        vlres = reshape(vl,nst,1);
        
        sw = int1d(lint);
        [shl,shld,shls,be] = shl1d(sw(1,:),lint,nel-1);
        
            if(step==0.6*L/cwave) && elem == 1
                nind = 0;
            end
            
        for ll = 1:lint                    

            % Evaluate 1-D basis functions at integration points
            Wgt = sw(2,ll);
            [shg, shgs, Jdet, be] = shg1d(xl,ndm,nel,shld(ll,:),shls(ll,:),nen,bf,der,be);
            
            c1 = Wgt*Jdet;
            
            Nmat(1:ndf:ndf*(nel-1)+1) = shl(ll,:);
                
            Bmat(1:ndf:ndf*(nel-1)+1) = shg';
            
            BBmat(1:ndf:ndf*(nel-1)+1) = shgs';
            
            ElemK = ElemK + c1*Bmat'*Dmat*Bmat;
            ElemM = ElemM + c1*(Nmat'*Nmat);
            ElemLLT = ElemLLT + c1*(BBmat'*Dmat*Nmat + Nmat'*Dmat*BBmat);
            ElemT = ElemT + c1*BBmat'*Tmat*BBmat;

        end %je
        
        % Final formulas
        if tau > - 2
        ElemM = ElemM + tau/(1-tau)*Nbeta*tstep^2*(ElemLLT + ElemK); % combined mass
        ElemK = ElemK - tau/(1-tau)*Nbeta*tstep^2*ElemT; % combined stiffness
        else % do mass averaging from Fried et al.
        h = xl(1,nel) - xl(1,1);
        ElemM = ElemM + ElemT*h^3/16*h/15/2*(tau + 3);
        end
        
        ElemE = 1/2*(ulres'*ElemK*ulres + vlres'*ElemM*vlres);
        
%%
    case 26 % Element stresses
        
        ElemS = zeros(nestr,1);

        %Set integration point number; use the rule for integrating K
        lint = nel-1;
        
        ib = 0;
        bf = 1;
        der = 1;
        
        fb = 0;
        Dmat = cwave^2;
        Bmat = zeros(1,nst);

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
            
            c1 = Wgt*Jdet;
                
            xint = xl*shl(ll,:)';
            Bmat(1:ndf:ndf*(nel-1)+1) = shg';
            sigma = Dmat*(Bmat*ulres);
            
            ElemS(2*ll-1:2*ll) = [xint/L sigma/NodeLoadnp(3)]';%xint; sigma];

        end %je
        
end