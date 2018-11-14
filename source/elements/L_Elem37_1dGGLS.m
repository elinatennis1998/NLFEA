% Tim Truster
% 08/27/2015
%
% Combined GLS and GGLS formulations for wave equation.
% tau1 is for GLS (M), tau2 is for GGLS (K)

% Set Material Properties

cwave = mateprop(1);
tau1 = mateprop(2);
tau2 = mateprop(3);

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
        ElemL1 = ElemK;
        ElemT1 = ElemK;
        ElemL2 = ElemK;
        ElemT2 = ElemK;

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
        BBBmat = zeros(1,nst);

        % Compute kinematic fields at intermediate time levels
        ul_n_af = (1 - Nalphaf)*ul + Nalphaf*ul_n;
        al_n_am = (1 - Nalpham)*al + Nalpham*al_n;
        ulres = reshape(ul_n_af,ndf*nel,1);
        alres = reshape(al_n_am,ndf*nel,1);
        
        sw = int1d(lint);
        [shl,shld,shls,be,shlt] = shl1d(sw(1,:),lint,nel-1);
        
        if(exist('iprob','var')&&iprob==1&&step==stepmax) && elem == 1
            nind = 0;
        end
            
        for ll = 1:lint                    

            % Evaluate 1-D basis functions at integration points
            Wgt = sw(2,ll);
            [shg, shgs, Jdet, be,shgt] = shg1d(xl,ndm,nel,shld(ll,:),shls(ll,:),nen,bf,der,be,shlt(ll,:));
            
            c1 = Wgt*Jdet;
            
            Nmat(1:ndf:ndf*(nel-1)+1) = shl(ll,:);
                
            Bmat(1:ndf:ndf*(nel-1)+1) = shg';
            
            BBmat(1:ndf:ndf*(nel-1)+1) = shgs';
            
            BBBmat(1:ndf:ndf*(nel-1)+1) = shgt';

            %Compute rho, accel
            accel = Nmat*alres;
            sigma = Dmat*(Bmat*ulres);
            
            Fvec = fb;
            
            ElemK = ElemK + c1*Bmat'*Dmat*Bmat;
            ElemM = ElemM + c1*(Nmat'*Nmat);
            ElemL1 = ElemL1 + c1*BBmat'*Dmat*Nmat;
            ElemT1 = ElemT1 + c1*BBmat'*Tmat*BBmat;
            ElemL2 = ElemL2 + c1*BBBmat'*Dmat*Bmat;
            ElemT2 = ElemT2 + c1*BBBmat'*Tmat*BBBmat;
            
            if(exist('iprob','var')&&iprob==1&&step==stepmax)
                nind = nind + 1;
                xint = Nmat*xl';
                Sig1(nind,:)=[xint/L sigma/NodeLoadnp(3)];
                
            end

        end %je
        
        % Final formulas
        if tau2 > - 2
        ElemM = ElemM - tau1*(ElemM - Nbeta*tstep^2*ElemL1) + tau2*(ElemK/cwave^2 - Nbeta*tstep^2*ElemL2); % combined mass
        ElemK = ElemK + tau1*(ElemL1' - Nbeta*tstep^2*ElemT1) - tau2*(ElemL2' - Nbeta*tstep^2*ElemT2); % combined stiffness
        ElemF = -ElemM*alres - ElemK*ulres;
        ElemK = ElemM + Nbeta*tstep^2*ElemK;
        else % do mass averaging from Fried et al.
        h = xl(1,nel) - xl(1,1);
        ElemM = ElemM + ElemT1*h^3/16*h/15/2*(tau + 3);
        ElemF = -ElemM*alres - ElemK*ulres;
        ElemK = ElemM + Nbeta*tstep^2*ElemK;
        end
        
    case 5 %Compute Mass
        
        ElemK = zeros(ndf*nel);
%         ElemF = zeros(ndf*nel,1);
        ElemM = ElemK;
        ElemL1 = ElemK;
        ElemT1 = ElemK;
        ElemL2 = ElemK;
        ElemT2 = ElemK;

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
        BBBmat = zeros(1,nst);
        
        sw = int1d(lint);
        [shl,shld,shls,be,shlt] = shl1d(sw(1,:),lint,nel-1);
        
        if(exist('iprob','var')&&iprob==1&&step==stepmax) && elem == 1
            nind = 0;
        end
            
        for ll = 1:lint                    

            % Evaluate 1-D basis functions at integration points
            Wgt = sw(2,ll);
            [shg, shgs, Jdet, be,shgt] = shg1d(xl,ndm,nel,shld(ll,:),shls(ll,:),nen,bf,der,be,shlt(ll,:));
            
            c1 = Wgt*Jdet;
            
            Nmat(1:ndf:ndf*(nel-1)+1) = shl(ll,:);
                
            Bmat(1:ndf:ndf*(nel-1)+1) = shg';
            
            BBmat(1:ndf:ndf*(nel-1)+1) = shgs';
            
            BBBmat(1:ndf:ndf*(nel-1)+1) = shgt';
            
            Fvec = fb;
            
            ElemK = ElemK + c1*Bmat'*Dmat*Bmat;
            ElemM = ElemM + c1*(Nmat'*Nmat);
            ElemL1 = ElemL1 + c1*BBmat'*Dmat*Nmat;
            ElemL2 = ElemL2 + c1*BBBmat'*Dmat*Bmat;
            
            if(exist('iprob','var')&&iprob==1&&step==stepmax)
                nind = nind + 1;
                xint = Nmat*xl';
                Sig1(nind,:)=[xint/L sigma/NodeLoadnp(3)];
                
            end

        end %je
        
        % Final formulas
        if tau2 > - 2
        ElemM = ElemM - tau1*(ElemM - Nbeta*tstep^2*ElemL1) + tau2*(ElemK/cwave^2 - Nbeta*tstep^2*ElemL2); % combined mass
        else % do mass averaging from Fried et al.
        h = xl(1,nel) - xl(1,1);
        ElemM = ElemM + ElemT1*h^3/16*h/15/2*(tau + 3);
        end

    case 6 %Compute Residual
        
        ElemK = zeros(ndf*nel);
%         ElemF = zeros(ndf*nel,1);
        ElemM = ElemK;
        ElemL1 = ElemK;
        ElemT1 = ElemK;
        ElemL2 = ElemK;
        ElemT2 = ElemK;

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
        BBBmat = zeros(1,nst);

        % Compute kinematic fields at intermediate time levels
        ul_n_af = (1 - Nalphaf)*ul + Nalphaf*ul_n;
        al_n_am = (1 - Nalpham)*al + Nalpham*al_n;
        ulres = reshape(ul_n_af,ndf*nel,1);
        alres = reshape(al_n_am,ndf*nel,1);
        
        sw = int1d(lint);
        [shl,shld,shls,be,shlt] = shl1d(sw(1,:),lint,nel-1);
            
        for ll = 1:lint                    

            % Evaluate 1-D basis functions at integration points
            Wgt = sw(2,ll);
            [shg, shgs, Jdet, be,shgt] = shg1d(xl,ndm,nel,shld(ll,:),shls(ll,:),nen,bf,der,be,shlt(ll,:));
            
            c1 = Wgt*Jdet;
            
            Nmat(1:ndf:ndf*(nel-1)+1) = shl(ll,:);
                
            Bmat(1:ndf:ndf*(nel-1)+1) = shg';
            
            BBmat(1:ndf:ndf*(nel-1)+1) = shgs';
            
            BBBmat(1:ndf:ndf*(nel-1)+1) = shgt';

            %Compute rho, accel
            accel = Nmat*alres;
            sigma = Dmat*(Bmat*ulres);
            
            Fvec = fb;
            
            ElemK = ElemK + c1*Bmat'*Dmat*Bmat;
            ElemM = ElemM + c1*(Nmat'*Nmat);
            ElemL1 = ElemL1 + c1*BBmat'*Dmat*Nmat;
            ElemT1 = ElemT1 + c1*BBmat'*Tmat*BBmat;
            ElemL2 = ElemL2 + c1*BBBmat'*Dmat*Bmat;
            ElemT2 = ElemT2 + c1*BBBmat'*Tmat*BBBmat;

        end %je
        
        % Final formulas
        if tau2 > - 2
        ElemM = ElemM - tau1*(ElemM - Nbeta*tstep^2*ElemL1) + tau2*(ElemK/cwave^2 - Nbeta*tstep^2*ElemL2); % combined mass
        ElemK = ElemK + tau1*(ElemL1' - Nbeta*tstep^2*ElemT1) - tau2*(ElemL2' - Nbeta*tstep^2*ElemT2); % combined stiffness
        ElemF = -ElemM*alres - ElemK*ulres;
        else % do mass averaging from Fried et al.
        h = xl(1,nel) - xl(1,1);
        ElemM = ElemM + ElemT1*h^3/16*h/15/2*(tau + 3);
        ElemF = -ElemM*alres - ElemK*ulres;
        end

    case 21 %Compute Stiffness
        
        ElemK = zeros(ndf*nel);
        ElemL1 = ElemK;
        ElemT1 = ElemK;
        ElemL2 = ElemK;
        ElemT2 = ElemK;

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
        BBBmat = zeros(1,nst);
        
        sw = int1d(lint);
        [shl,shld,shls,be,shlt] = shl1d(sw(1,:),lint,nel-1);
            
        for ll = 1:lint                    

            % Evaluate 1-D basis functions at integration points
            Wgt = sw(2,ll);
            [shg, shgs, Jdet, be,shgt] = shg1d(xl,ndm,nel,shld(ll,:),shls(ll,:),nen,bf,der,be,shlt(ll,:));
            
            c1 = Wgt*Jdet;
            
            Nmat(1:ndf:ndf*(nel-1)+1) = shl(ll,:);
                
            Bmat(1:ndf:ndf*(nel-1)+1) = shg';
            
            BBmat(1:ndf:ndf*(nel-1)+1) = shgs';
            
            BBBmat(1:ndf:ndf*(nel-1)+1) = shgt';
            
            Fvec = fb;
            
            ElemK = ElemK + c1*Bmat'*Dmat*Bmat;
            ElemL1 = ElemL1 + c1*BBmat'*Dmat*Nmat;
            ElemT1 = ElemT1 + c1*BBmat'*Tmat*BBmat;
            ElemL2 = ElemL2 + c1*BBBmat'*Dmat*Bmat;
            ElemT2 = ElemT2 + c1*BBBmat'*Tmat*BBBmat;

        end %je
        
        % Final formulas
        if tau2 > - 2
        ElemK = ElemK + tau1*(ElemL1' - Nbeta*tstep^2*ElemT1) - tau2*(ElemL2' - Nbeta*tstep^2*ElemT2); % combined stiffness
        else % do mass averaging from Fried et al.
        h = xl(1,nel) - xl(1,1);
        end
        
    case 40 % Initialize FS acceleration
        
        
    case 12 % energy
        
        ElemK = zeros(ndf*nel);
        ElemM = ElemK;
        ElemL1 = ElemK;
        ElemT1 = ElemK;
        ElemL2 = ElemK;
        ElemT2 = ElemK;

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
        BBBmat = zeros(1,nst);

        % Compute kinematic fields at intermediate time levels
        ul_n_af = (1 - Nalphaf)*ul + Nalphaf*ul_n;
        al_n_am = (1 - Nalpham)*al + Nalpham*al_n;
        ulres = reshape(ul_n_af,ndf*nel,1);
        alres = reshape(al_n_am,ndf*nel,1);
        vlres = reshape(vl,nst,1);
        
        sw = int1d(lint);
        [shl,shld,shls,be,shlt] = shl1d(sw(1,:),lint,nel-1);
        
            if(step==0.6*L/cwave) && elem == 1
                nind = 0;
            end
            
        for ll = 1:lint                    

            % Evaluate 1-D basis functions at integration points
            Wgt = sw(2,ll);
            [shg, shgs, Jdet, be,shgt] = shg1d(xl,ndm,nel,shld(ll,:),shls(ll,:),nen,bf,der,be,shlt(ll,:));
            
            c1 = Wgt*Jdet;
            
            Nmat(1:ndf:ndf*(nel-1)+1) = shl(ll,:);
                
            Bmat(1:ndf:ndf*(nel-1)+1) = shg';
            
            BBmat(1:ndf:ndf*(nel-1)+1) = shgs';
            
            BBBmat(1:ndf:ndf*(nel-1)+1) = shgt';
            
            ElemK = ElemK + c1*Bmat'*Dmat*Bmat;
            ElemM = ElemM + c1*(Nmat'*Nmat);
            ElemL1 = ElemL1 + c1*BBmat'*Dmat*Nmat;
            ElemT1 = ElemT1 + c1*BBmat'*Tmat*BBmat;
            ElemL2 = ElemL2 + c1*BBBmat'*Dmat*Bmat;
            ElemT2 = ElemT2 + c1*BBBmat'*Tmat*BBBmat;

        end %je
        
        % Final formulas
        if tau2 > - 2
        ElemM = ElemM - tau1*(ElemM - Nbeta*tstep^2*ElemL1) + tau2*(ElemK/cwave^2 - Nbeta*tstep^2*ElemL2); % combined mass
        ElemK = ElemK + tau1*(ElemL1' - Nbeta*tstep^2*ElemT1) - tau2*(ElemL2' - Nbeta*tstep^2*ElemT2); % combined stiffness
        else % do mass averaging from Fried et al.
        h = xl(1,nel) - xl(1,1);
        ElemM = ElemM + ElemT1*h^3/16*h/15/2*(tau + 3);
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
            [shg, shgs, Jdet, be,shgt] = shg1d(xl,ndm,nel,shld(ll,:),shls(ll,:),nen,bf,der,be,shlt(ll,:));
            
            c1 = Wgt*Jdet;
                
            xint = xl*shl(ll,:)';
            Bmat(1:ndf:ndf*(nel-1)+1) = shg';
            sigma = Dmat*(Bmat*ulres);
            
            ElemS(2*ll-1:2*ll) = [xint/L sigma/NodeLoadnp(3)]';%xint; sigma];

        end %je
        
end