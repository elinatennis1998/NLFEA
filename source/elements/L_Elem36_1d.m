% Tim Truster
% 08/28/2015
%
% Combined GLS and GGLS formulations for Helmholtz equation.
% tau1 is for GLS (M), tau2 is for GGLS (K)

% Set Material Properties

cwave = mateprop(1);
kwave = mateprop(2);
tau1 = mateprop(3);
tau2 = mateprop(4);

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
        
        Nmat = zeros(1,nst);
        Bmat = zeros(1,nst);
        BBmat = zeros(1,nst);
        BBBmat = zeros(1,nst);
        
        ulres = reshape(ul,ndf*nel,1);
        
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
            
            ElemK = ElemK + c1*(Bmat'*Bmat);
            ElemM = ElemM + c1*(Nmat'*Nmat);
            ElemL1 = ElemL1 + c1*BBmat'*Nmat;
            ElemT1 = ElemT1 + c1*(BBmat'*BBmat);
            ElemL2 = ElemL2 + c1*BBBmat'*Bmat;
            ElemT2 = ElemT2 + c1*(BBBmat'*BBBmat);

        end %je
        
        % Final formulas
        ElemM = ElemM - tau1*(ElemM - ElemL1) + tau2*(ElemK - ElemL2); % combined mass
        ElemK = ElemK + tau1*(ElemL1' - ElemT1) - tau2*(ElemL2' - ElemT2); % combined stiffness
        ElemK = ElemK - kwave^2*ElemM;
        ElemF =  - ElemK*ulres;

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
        
        Nmat = zeros(1,nst);
        Bmat = zeros(1,nst);
        BBmat = zeros(1,nst);
        BBBmat = zeros(1,nst);
        
        ulres = reshape(ul,ndf*nel,1);
        
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
            
            ElemK = ElemK + c1*(Bmat'*Bmat);
            ElemM = ElemM + c1*(Nmat'*Nmat);
            ElemL1 = ElemL1 + c1*BBmat'*Nmat;
            ElemT1 = ElemT1 + c1*(BBmat'*BBmat);
            ElemL2 = ElemL2 + c1*BBBmat'*Bmat;
            ElemT2 = ElemT2 + c1*(BBBmat'*BBBmat);

        end %je
        
        % Final formulas
        ElemM = ElemM - tau1*(ElemM - ElemL1) + tau2*(ElemK - ElemL2); % combined mass
        ElemK = ElemK + tau1*(ElemL1' - ElemT1) - tau2*(ElemL2' - ElemT2); % combined stiffness
        ElemK = ElemK - kwave^2*ElemM;
        ElemF =  - ElemK*ulres;

    case 21 %Compute Stiffness
        
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
        
        Nmat = zeros(1,nst);
        Bmat = zeros(1,nst);
        BBmat = zeros(1,nst);
        BBBmat = zeros(1,nst);
        
        ulres = reshape(ul,ndf*nel,1);
        
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
            
            ElemK = ElemK + c1*(Bmat'*Bmat);
            ElemM = ElemM + c1*(Nmat'*Nmat);
            ElemL1 = ElemL1 + c1*BBmat'*Nmat;
            ElemT1 = ElemT1 + c1*(BBmat'*BBmat);
            ElemL2 = ElemL2 + c1*BBBmat'*Bmat;
            ElemT2 = ElemT2 + c1*(BBBmat'*BBBmat);

        end %je
        
        % Final formulas
        ElemM = ElemM - tau1*(ElemM - ElemL1) + tau2*(ElemK - ElemL2); % combined mass
        ElemK = ElemK + tau1*(ElemL1' - ElemT1) - tau2*(ElemL2' - ElemT2); % combined stiffness
        ElemK = ElemK - kwave^2*ElemM;
        
    case 40 % Initialize FS acceleration
        
        
        
%%
    case 26 % Element stresses
        
        ElemS = zeros(nestr,1);

        %Set integration point number; use the rule for integrating K
        lint = nel-1;
        
        ib = 0;
        bf = 1;
        der = 1;
        
        Dmat = cwave^2;
        Bmat = zeros(1,nst);
        
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