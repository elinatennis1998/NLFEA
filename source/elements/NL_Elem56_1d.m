% Tim Truster
% 1/07/2015
% 1d advection-diffusion-transient element

% Set material
mat = @simple1x1;
mprops = [mateprop(2)];

% Set Material Properties

% diffusion matrix
A11 = mateprop(1);
% A22 = mateprop(2);
% A12 = mateprop(3);
% Av = [A11 A22 2*A12];
% A = [A11 A12; A12 A22];
% % advection velocity
% av = mateprop(4:5)';
av = mateprop(2);
% transient coefficient
rho = mateprop(3);

switch isw %Task Switch
%%    
    case 1
        
        % prescribe away extra dofs greater than 1
        if ndf > 1
            
            for i = 2:ndf
                lie(i,1) = 0;
            end
            
        end
        
%%
    case 3 %Compute Stiffness and Residual
        
        ElemK = zeros(ndf*nel);
        ElemF = zeros(ndf*nel,1);

        %Set integration number
        lint = nel;
        
        ib = 0;
        bf = 1;
        der = 1;
        
        Dmat = A11;
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
            
            %                   transient/mass part               diffusion         advection
            [v_d, dv_d] = mat(Nmat*ul', mprops);
            ElemK = ElemK + c1*(coeffm*(Nmat'*rho*Nmat) + coeffk*(Bmat'*Dmat*Bmat + Nmat'*(v_d)'*Bmat + (Nmat'*(dv_d)'*Nmat)*(Bmat*ul')));
            ElemF = ElemF - c1*(Nmat'*rho*Nmat*vl'      + (Bmat'*Dmat*Bmat + Nmat'*(v_d)'*Bmat)*ul');

        end %je

        ElemF;
        
    case 6 %Compute Residual
        
        ElemF = zeros(ndf*nel,1);

        %Set integration number
        lint = nel;
        
        ib = 0;
        bf = 1;
        der = 1;
        
        Dmat = A11;
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
            
            %                   transient/mass part               diffusion      advection
            [v_d, dv_d] = mat(Nmat*ul', mprops);
            ElemF = ElemF - c1*(Nmat'*rho*Nmat*vl'      + (Bmat'*Dmat*Bmat + Nmat'*(v_d)'*Bmat)*ul');

        end %je
        
        ElemF;

    case 21 %Compute Stiffness
        
        ElemK = zeros(ndf*nel);

        %Set integration number
        lint = nel;
        
        ib = 0;
        bf = 1;
        der = 1;
        
        Dmat = A11;
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
            
            %                   transient/mass part               diffusion      advection
            [v_d, dv_d] = mat(Nmat*ul', mprops);
            ElemK = ElemK + c1*((Bmat'*Dmat*Bmat + Nmat'*(v_d)'*Bmat + (Nmat'*(dv_d)'*Nmat)*(Bmat*ul')));

        end %je
        
    case 15 % Body force

        
        
    case 5 % Mass matrix
        
        ElemM = zeros(ndf*nel);

        %Set integration number
        lint = nel;
        
        ib = 0;
        bf = 1;
        der = 1;
        
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
            
            %                   transient/mass part               diffusion      advection
            ElemM = ElemM + c1*((Nmat'*rho*Nmat));

        end %je
        
end