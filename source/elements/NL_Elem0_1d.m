% Tim Truster
% 12/10/2013
% 1-D dynamic bar element
% Tested using dynamic_bar1d.m

Emod = mateprop(1);
Aelem = mateprop(2);
rho = mateprop(3);
Et = mateprop(4);
eps_y = mateprop(5);
sig_y = Emod*eps_y;

switch isw %Task Switch
    
    case 1
        
        if ndf > 1
            
            for i = 2:ndf
                lie(i,1) = 0;
            end
            
        end
%%
    case 3

        ElemK = zeros(ndf*nel);
        ElemF = zeros(ndf*nel,1);
        Nmat = zeros(1,ndf*nel);
        Bmat = zeros(1,ndf*nel);
        
        ulres = reshape(ul,ndf*nen,1);
        alres = reshape(al,ndf*nen,1);

        % Load Gauss Points for quadrature
        lint = nel;

        der = 0;
        bf = 0;
        ib = 0;
        sw = int1d(lint);
        [shl,shld,shls,be] = shl1d(sw(1,:),lint,nel-1);
        
        if(exist('iprob','var')&&iprob==1&&step==stepmax) && elem == 1
            nind = 0;
        end
                
        for je = 1:lint

            Wgt = sw(2,je);
            [shg, shgs, Jdet, be] = shg1d(xl,ndm,nel,shld(je,:),shls(je,:),nen,bf,der,be);
                    
            c1 = Wgt*Jdet*Aelem;
                
            Nmat(1:ndf:ndf*(nel-1)+1) = shl(je,:);
                
            Bmat(1:ndf:ndf*(nel-1)+1) = shg';
            
            accel = Nmat*alres;
            E = Bmat*ulres;
            if abs(E) < eps_y % linear
            Dmat = Emod;
            sigma = Dmat*E;
            else % second branch
            Dmat = Et;
            sigma = Dmat*(E-sign(E)*eps_y) + sign(E)*sig_y; % very important to have sign function for compression
            end
            
            if(exist('iprob','var')&&iprob==1&&step==stepmax)
                xint = Nmat*xl';
                nind = nind + 1;
                Sig1(nind,:)=[xint/L sigma/NodeLoadnp(3)];
                
            end
            
            ElemF = ElemF + c1*(- rho*Nmat'*accel - Bmat'*sigma);
            if (transient <= 4  && transient >= 3) || transient == 13 % mass is included in NR_Loop3
            ElemK = ElemK + c1*(Bmat'*Dmat*Bmat);
            else
            ElemK = ElemK + c1*(Nmat'*rho*coeffm*Nmat + Bmat'*Dmat*coeffk*Bmat);
            end
            
        end %je
ElemK;
%%
    case 5 %Compute Mass
        
        ElemM = zeros(nst);
        Nmat = zeros(1,ndf*nel);
        
        % Load Gauss Points for quadrature
        lint = nel;

        der = 0;
        bf = 0;
        ib = 0;
        sw = int1d(lint);
        [shl,shld,shls,be] = shl1d(sw(1,:),lint,nel-1);
        
        for je = 1:lint

            Wgt = sw(2,je);
            [shg, shgs, Jdet, be] = shg1d(xl,ndm,nel,shld(je,:),shls(je,:),nen,bf,der,be);
                    
            c1 = Wgt*Jdet*Aelem;
                
            Nmat(1:ndf:ndf*(nel-1)+1) = shl(je,:);
            
            ElemM = ElemM + c1*(Nmat'*rho*Nmat);
            
        end %je
        
%%
    case 6

        ElemF = zeros(ndf*nel,1);
        Nmat = zeros(1,ndf*nel);
        Bmat = zeros(1,ndf*nel);
        
        ulres = reshape(ul,ndf*nen,1);
        alres = reshape(al,ndf*nen,1);

        Dmat = Emod;
        
        % Load Gauss Points for quadrature
        lint = nel;

        der = 0;
        bf = 0;
        ib = 0;
        sw = int1d(lint);
        [shl,shld,shls,be] = shl1d(sw(1,:),lint,nel-1);
        
        for je = 1:lint

            Wgt = sw(2,je);
            [shg, shgs, Jdet, be] = shg1d(xl,ndm,nel,shld(je,:),shls(je,:),nen,bf,der,be);
                    
            c1 = Wgt*Jdet*Aelem;
                
            Nmat(1:ndf:ndf*(nel-1)+1) = shl(je,:);
                
            Bmat(1:ndf:ndf*(nel-1)+1) = shg';
            
            accel = Nmat*alres;
            E = Bmat*ulres;
            if abs(E) < eps_y % linear
            Dmat = Emod;
            sigma = Dmat*E;
            else % second branch
            Dmat = Et;
            sigma = Dmat*(E-sign(E)*eps_y) + sign(E)*sig_y; % very important to have sign function for compression
            end
            
            ElemF = ElemF + c1*(- rho*Nmat'*accel - Bmat'*sigma);
            
        end %je
%%
    case 15 % body force for linear element
        
%         ElemF = zeros(nst,1); In FormFE
        
        Nmat = zeros(1,nst);
        
        % Load Gauss Points for quadrature
        lint = nel - 1;
%         if nel == 3
%             lint = lintt3;%13;
%         elseif nel == 4
%             lint = lintq4;
%         elseif nel == 6
%             lint = lintt6;%13;
%         elseif nel == 9
%             lint = lintq9;
%         end

        der = 0;
        bf = 0;
        
        if exist('iprob','var') && iprob == 1
            fb = 10;
        else
            fb = bodyf;
        end
        
        sw = int1d(lint);
        [shl,shld,shls,be] = shl1d(sw(1,:),lint,nel-1);
        
        for je = 1:lint

            Wgt = sw(2,je);
            [shg, shgs, Jdet, be] = shg1d(xl,ndm,nel,shld(je,:),shls(je,:),nen,bf,der,be);
                    
            c1 = Wgt*Jdet*Aelem;
                
            Nmat(1:ndf:ndf*(nel-1)+1) = shl(je,:);
            
            ElemF = ElemF + c1*(Nmat'*fb);
                
        end %je
        
    case 12 % energy
        
        ElemE = 0;
        Nmat = zeros(1,nst);
        Bmat = zeros(1,nst);
        
        D = Emod;
        ulres = reshape(ul,ndf*nen,1);
        vlres = reshape(vl,ndf*nen,1);
        
        % Load Guass Integration Points
        lint = nel;

        der = 0;
        bf = 0;
        ib = 0;
        sw = int1d(lint);
        [shl,shld,shls,be] = shl1d(sw(1,:),lint,nel-1);

        %Integration Loop
        for je = 1:lint

            Wgt = sw(2,je);
            [shg, shgs, Jdet, be] = shg1d(xl,ndm,nel,shld(je,:),shls(je,:),nen,bf,der,be);
                    
            c1 = Wgt*Jdet*Aelem;
                
            Nmat(1:ndf:ndf*(nel-1)+1) = shl(je,:);
                
            Bmat(1:ndf:ndf*(nel-1)+1) = shg';
            
            %Compute F, sigma, D
            E = Bmat*ulres;
            velo = Nmat*vlres;
            
            % Strain energy
            if abs(E) < eps_y % linear
            U = 1/2*E'*Emod*E;
            else % second branch
            E2 = (E-sign(E)*eps_y); % very important to have sign function for compression
            U = 1/2*E2'*Et*E2 + sign(E)*sig_y*E - 1/2*eps_y'*Emod*eps_y; % computed so that U(eps_y) = 1/2*eps_y'*Emod*eps_y
            end
                
            if transient <= 5  && transient >= 4 % mass is included in NR_Loop3
            ElemE = ElemE + c1*U;
            else
            ElemE = ElemE + c1*(U + 1/2*velo'*rho*velo);
            end

        end %je
%%
    case 21

        ElemK = zeros(ndf*nel);
        Bmat = zeros(1,ndf*nel);
        
        ulres = reshape(ul,ndf*nen,1);

        % Load Gauss Points for quadrature
        lint = nel;

        der = 0;
        bf = 0;
        ib = 0;
        sw = int1d(lint);
        [shl,shld,shls,be] = shl1d(sw(1,:),lint,nel-1);
                
        for je = 1:lint

            Wgt = sw(2,je);
            [shg, shgs, Jdet, be] = shg1d(xl,ndm,nel,shld(je,:),shls(je,:),nen,bf,der,be);
                    
            c1 = Wgt*Jdet*Aelem;
                
            Nmat(1:ndf:ndf*(nel-1)+1) = shl(je,:);
                
            Bmat(1:ndf:ndf*(nel-1)+1) = shg';
            
            E = Bmat*ulres;
            if abs(E) < eps_y % linear
            Dmat = Emod;
            sigma = Dmat*E;
            else % second branch
            Dmat = Et;
            sigma = Dmat*(E-sign(E)*eps_y) + sign(E)*sig_y; % very important to have sign function for compression
            end
            
            ElemK = ElemK + c1*(Bmat'*Dmat*Bmat);
            
        end %je
        
end %Task Switch
