% Tim Truster
% 09/30/2015
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
        
        nh3 = 15;
%%
    case 3

        ElemK = zeros(ndf*nel);
        ElemF = zeros(ndf*nel,1);
        b1 = ElemF;
        b3 = ElemF;
        Nmat = zeros(1,ndf*nel);
        Bmat = zeros(1,ndf*nel);
        ElemM = zeros(nst);
        
        litk = zeros(ndf*nel,1);
        litkk = zeros(ndf*nel);
        litf = zeros(ndf*nel,1);
        kbar = zeros(ndf*nel);
        A13 = kbar;
        A31 = kbar;
        A32 = ElemF;
        A33 = kbar;
        A11 = kbar;
        
        ulres = reshape(ul,ndf*nen,1);
        alres = reshape(al,ndf*nen,1);

        % uncondense
        if initia == 0
        uldd = zeros(ndf, nen);
        uldd(ELDOFTa) = del_ModelDx(EGDOFTa)';
        tau11 = hr(nhc+0:nhc+1);
        A33A31 = reshape(hr(nhc+4:nhc+7),2,2);
        A33b3 = hr(nhc+7:nhc+8);
        A33A32 = hr(nhc+9:nhc+10);
        del_tau11 = -A33A31*uldd' + A33b3 - del_lener*A33A32;
%         y1 = zeros(ndf, nen);
%         y1(ELDOFTa) = y(EGDOFTa)';
%         z1 = zeros(ndf, nen);
%         z1(ELDOFTa) = z(EGDOFTa)';
%         y2 = -A33A31*y1'+A33b3;
%         z2 = -A33A31*z1'+A33A32;
%         dd2 = y2 - del_lener*z2;
        tau11 = tau11 + del_tau11;
%         if elem == 10
%             del_tau11
%             tau11
%         end
        else
            tau11 = ElemF;
            del_tau11 = ElemF;
            hr(nhc+0:nhc+1) = tau11;
        end
        
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
            
            % NOTE: these formulas assume linear finite elements (constant
            % derivatives)
b1 = b1 -(((1-lener)*c1*(Nmat'*rho*Nmat) + (Nmat*tau11)*c1*(Bmat'*Dmat*Bmat))*(4/tstep^2*(ul' - ul_n') - 4/tstep*vl_n' - al_n') - c1*( - Bmat'*sigma));
A11 = A11  + (1-lener)*(4/tstep^2*c1*(Nmat'*rho*Nmat)) + (Nmat*tau11)*4/tstep^2*c1*(Bmat'*Bmat) + c1*(Bmat'*Dmat*Bmat);
%             ElemM = ElemM + c1*(Nmat'*rho*Nmat);
%             ElemK = ElemK + c1*(Bmat'*Dmat*Bmat);
            litkk(:,je) = c1*(4/tstep^2*(Bmat'*Bmat)*ulres - Bmat'*Bmat*al_n' - 4/tstep^2*Bmat'*Bmat*ul_n' + 4/tstep*Bmat'*Bmat*vl_n');
            litf = c1*(- Bmat'*sigma - Nmat'*Nmat*al_n' - 2/tstep*Nmat'*Nmat*vl_n');
            kbar = kbar + c1*(Bmat'*Bmat);
            A13(:,je) = c1*(Bmat'*Bmat)*(4/tstep^2*(ul' - ul_n') - 4/tstep*vl_n' - al_n');
b3 = b1 -(lener*litf - (Nmat*tau11)*c1*(4/tstep^2*(Bmat'*Bmat)*ulres - Bmat'*Bmat*al_n' - 4/tstep^2*Bmat'*Bmat*ul_n' + 4/tstep*Bmat'*Bmat*vl_n'));
A31 =  A31 + lener*c1*(Bmat'*Dmat*Bmat) - c1*(4/tstep^2*(Bmat'*Bmat));
A32 = A32 + litf;
A33 = A33 -c1*(4/tstep^2*(Bmat'*Bmat)*ulres - Bmat'*Bmat*al_n' - 4/tstep^2*Bmat'*Bmat*ul_n' + 4/tstep*Bmat'*Bmat*vl_n')*Nmat;
            
        end %je
if initia == 0
% apply static condensation

% j = 1;
% Untitled3
% A31(j) = -A31j;
% A11(:,j) = -A11j;
% j = 2;
% Untitled3
% A31(j) = -A31j;
% A11(:,j) = -A11j;
% Untitled4
% A33 = -A33j;
% A13 = -A13j;
% Untitled7
% A32 = -A32j;

ElemK = A11 - (A13/A33)*A31;
ElemF = b1 - (A13/A33)*b3;
% store
hr(nhc+0:nhc+1) = tau11;
hr(nhc+2:nhc+3) = del_tau11;
hr(nhc+4:nhc+7) = A33\A31;
hr(nhc+7:nhc+8) = A33\b3;
hr(nhc+9:nhc+10) = A33\A32;
else
% b1 = -(((1-lener)*ElemM + tau11*kbar)*(4/tstep^2*(ul' - ul_n') - 4/tstep*vl_n' - al_n') - ElemF);
% A11 = (1-lener)*(4/tstep^2*ElemM) + tau11*4/tstep^2*kbar + ElemK;
ElemK = A11;
ElemF = b1;
end

% % apply static condensation
% b1 = -(((1-lener)*ElemM + tau11*kbar)*(4/tstep^2*(ul' - ul_n') - 4/tstep*vl_n' - al_n') - ElemF);
% A11 = (1-lener)*(4/tstep^2*ElemM) + tau11*4/tstep^2*kbar + ElemK;
% % A13 = kbar*(4/tstep^2*(ul' - ul_n') - 4/tstep*vl_n' - al_n');
% % b3 = -(lener/(litk'*litk)*litk'*litf - tau11);
% % A31 = lener/(litk'*litk)*litf'*4/tstep^2*kbar - lener/(litk'*litk)^2*8/tstep^2*(kbar*litk*litk'*litf)' - lener/(litk'*litk)*litk'*ElemK;
% % A32 = 1/(litk'*litk)*litk'*litf;
% % A33 = -1;
% 
% ElemK = A11;% - A13/A33*A31;
% ElemF = b1;% - A13/A33*b3;
% % % store
% % hr(nhc+0) = tau11;
% % hr(nhc+1) = del_tau11;
% % hr(nhc+2:nhc+3) = 1/A33*A31;
% % hr(nhc+4) = 1/A33*b3;
% % hr(nhc+5) = 1/A33*A32;

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
            
            ElemF = ElemF + c1*(- Bmat'*sigma);
            
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
        
    case 80 % compute A12 for energy constraint method

        ElemF = zeros(ndf*nel,1);
        Nmat = zeros(1,ndf*nel);
        Bmat = zeros(1,ndf*nel);
        ElemM = zeros(nst);
        
        litk = zeros(ndf*nel,1);
        litf = zeros(ndf*nel,1);
        kbar = zeros(ndf*nel);
        
        ulres = reshape(ul,ndf*nen,1);
        alres = reshape(al,ndf*nen,1);
        
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
            
            ElemM = ElemM + c1*(Nmat'*rho*Nmat);
            litk = litk + c1*(4/tstep^2*(Bmat'*Bmat)*ulres - Bmat'*Bmat*al_n' - 4/tstep^2*Bmat'*Bmat*ul_n' + 4/tstep*Bmat'*Bmat*vl_n');
            litf = litf + c1*(- Bmat'*sigma - Nmat'*Nmat*al_n' - 2/tstep*Nmat'*Nmat*vl_n');
            kbar = kbar + c1*(Bmat'*Bmat);
            
        end %je
A12 = -ElemM*(4/tstep^2*(ul' - ul_n') - 4/tstep*vl_n' - al_n');
A13 = kbar*(4/tstep^2*(ul' - ul_n') - 4/tstep*vl_n' - al_n');
A32 = 1/(litk'*litk)*litk'*litf;
A33 = -1;

Untitled4
A33 = -A33j;
A13 = -A13j;
Untitled7
A32 = -A32j;
A12 = -A12j;
ElemF = A12 - A13/A33*A32;

        
end %Task Switch
