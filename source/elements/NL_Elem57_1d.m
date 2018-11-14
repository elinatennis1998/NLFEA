% Tim Truster
% 1/08/2015
% 1d advection-diffusion-transient element, with 3 component fields
% template for Mark to add his stuff
% Stiffness matrix looks like:
% [w1,x w2,x w3,x]*[A11 A12 A13  [u1,x]
%                   A12 A22 A23  [u2,x]
%                   A13 A23 A33]*[u3,x]
% where [u1,x]  
%       [u2,x] = B*ul
%       [u3,x]
% with ul organized like [node1 node2] and node1 = [u1_1 u1_2 u1_3]

% Set Material Properties
matmodel = @simple3x1;
mprops = [mateprop(7) mateprop(9)];

% diffusion matrix
A11 = mateprop(1);
A22 = mateprop(2);
A33 = mateprop(3);
A12 = mateprop(4);
A23 = mateprop(5);
A13 = mateprop(6);
% Av = [A11 A22 2*A12];
% A = [A11 A12; A12 A22];
% % advection velocity
% av = mateprop(4:5)';
av = mateprop(7);
% transient coefficient
rho = mateprop(8);

switch isw %Task Switch
%%    
    case 1
        
        % prescribe away extra dofs greater than 1
        if ndf > 3
            
            for i = 4:ndf
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
        
        Dmat = [A11 A12 A13; A12 A22 A23; A13 A23 A33];
        Dmat = Dmat(1:ndf,1:ndf);
        Nmat = zeros(ndf,nst);
        Bmat = zeros(ndf,nst);
        BBmat = zeros(ndf,nst);
        
        ulres = reshape(ul,nst,1);
        vlres = reshape(vl,nst,1);
        
        sw = int1d(lint);
        [shl,shld,shls,be] = shl1d(sw(1,:),lint,nel-1);
            
        for ll = 1:lint                    

            % Evaluate 1-D basis functions at integration points
            Wgt = sw(2,ll);
            [shg, shgs, Jdet, be] = shg1d(xl,ndm,nel,shld(ll,:),shls(ll,:),nen,bf,der,be);
            
            c1 = Wgt*Jdet;
            
            % Form B matrix
            for ie = 1:nel
                
              Nmat(1:ndf,(ie-1)*ndf+1:(ie-1)*ndf+ndf) = shl(ll,ie)*eye(ndf);
              Bmat(1:ndf,(ie-1)*ndf+1:(ie-1)*ndf+ndf) = shg(ie)*eye(ndf);
              BBmat(1:ndf,(ie-1)*ndf+1:(ie-1)*ndf+ndf) = shgs(ie)*eye(ndf);
                 
            end
            
            
            x = dot(Nmat(1,1:ndf:end),xl);
            [v_d, dv_d] = matmodel(x, Nmat*ulres , mprops);
            
            %                   transient/mass part               diffusion         advection
            ElemK = ElemK + c1*(coeffm*(Nmat'*rho*Nmat) + coeffk*(Bmat'*Dmat*Bmat + Nmat'*v_d*Bmat + (Nmat'*Bmat)*(ulres*(dv_d'*Nmat))));
            ElemF = ElemF - c1*(Nmat'*rho*Nmat*vlres      + (Bmat'*Dmat*Bmat + Nmat'*v_d*Bmat)*ulres);

        end %je

        ElemF;
        
    case 6 %Compute Residual
        
        ElemF = zeros(ndf*nel,1);

        %Set integration number
        lint = nel;
        
        ib = 0;
        bf = 1;
        der = 1;
        
        Dmat = [A11 A12 A13; A12 A22 A23; A13 A23 A33];
        Dmat = Dmat(1:ndf,1:ndf);
        Nmat = zeros(ndf,nst);
        Bmat = zeros(ndf,nst);
        BBmat = zeros(ndf,nst);
        
        ulres = reshape(ul,nst,1);
        vlres = reshape(vl,nst,1);
        
        sw = int1d(lint);
        [shl,shld,shls,be] = shl1d(sw(1,:),lint,nel-1);
            
        for ll = 1:lint                    

            % Evaluate 1-D basis functions at integration points
            Wgt = sw(2,ll);
            [shg, shgs, Jdet, be] = shg1d(xl,ndm,nel,shld(ll,:),shls(ll,:),nen,bf,der,be);
            
            c1 = Wgt*Jdet;
            
            % Form B matrix
            for ie = 1:nel
                
              Nmat(1:ndf,(ie-1)*ndf+1:(ie-1)*ndf+ndf) = shl(ll,ie)*eye(ndf);
              Bmat(1:ndf,(ie-1)*ndf+1:(ie-1)*ndf+ndf) = shg(ie)*eye(ndf);
              BBmat(1:ndf,(ie-1)*ndf+1:(ie-1)*ndf+ndf) = shgs(ie)*eye(ndf);
                 
            end
            
            x = dot(Nmat(1,1:ndf:end),xl);
            [v_d, dv_d] = matmodel(x, Nmat*ulres , mprops);         
            %                   transient/mass part               diffusion         advection
            ElemF = ElemF - c1*(Nmat'*rho*Nmat*vlres      + (Bmat'*Dmat*Bmat + Nmat'*v_d*Bmat)*ulres);

        end %je
        
        ElemF;

    case 21 %Compute Stiffness
        
        ElemK = zeros(ndf*nel);

        %Set integration number
        lint = nel;
        
        ib = 0;
        bf = 1;
        der = 1;
        
        Dmat = [A11 A12 A13; A12 A22 A23; A13 A23 A33];
        Dmat = Dmat(1:ndf,1:ndf);
        Nmat = zeros(ndf,nst);
        Bmat = zeros(ndf,nst);
        BBmat = zeros(ndf,nst);
        
        ulres = reshape(ul,nst,1);
        vlres = reshape(vl,nst,1);
        
        sw = int1d(lint);
        [shl,shld,shls,be] = shl1d(sw(1,:),lint,nel-1);
            
        for ll = 1:lint                    

            % Evaluate 1-D basis functions at integration points
            Wgt = sw(2,ll);
            [shg, shgs, Jdet, be] = shg1d(xl,ndm,nel,shld(ll,:),shls(ll,:),nen,bf,der,be);
            
            c1 = Wgt*Jdet;
            
            % Form B matrix
            for ie = 1:nel
                
              Nmat(1:ndf,(ie-1)*ndf+1:(ie-1)*ndf+ndf) = shl(ll,ie)*eye(ndf);
              Bmat(1:ndf,(ie-1)*ndf+1:(ie-1)*ndf+ndf) = shg(ie)*eye(ndf);
              BBmat(1:ndf,(ie-1)*ndf+1:(ie-1)*ndf+ndf) = shgs(ie)*eye(ndf);
                 
            end

            x = dot(Nmat(1,1:ndf:end),xl);
            [v_d, dv_d] = matmodel(x, Nmat*ulres , mprops);          
            ElemK = ElemK + c1*((Bmat'*Dmat*Bmat + Nmat'*av*Bmat + (Nmat'*Bmat)*(ulres*(dv_d'*Nmat))));

        end %je
        
    case 15 % Body force

        
        
    case 5 % Mass matrix
        
        ElemM = zeros(ndf*nel);

        %Set integration number
        lint = nel;
        
        ib = 0;
        bf = 1;
        der = 1;
        
        Nmat = zeros(ndf,nst);
        Bmat = zeros(ndf,nst);
        BBmat = zeros(ndf,nst);
        
        sw = int1d(lint);
        [shl,shld,shls,be] = shl1d(sw(1,:),lint,nel-1);
            
        for ll = 1:lint                    

            % Evaluate 1-D basis functions at integration points
            Wgt = sw(2,ll);
            [shg, shgs, Jdet, be] = shg1d(xl,ndm,nel,shld(ll,:),shls(ll,:),nen,bf,der,be);
            
            c1 = Wgt*Jdet;
            
            % Form B matrix
            for ie = 1:nel
                
              Nmat(1:ndf,(ie-1)*ndf+1:(ie-1)*ndf+ndf) = shl(ll,ie)*eye(ndf);
                 
            end
            
            %                   transient/mass part               diffusion      advection
            ElemM = ElemM + c1*((Nmat'*rho*Nmat));

        end %je
        
end