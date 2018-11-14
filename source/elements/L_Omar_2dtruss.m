% Tim Truster
% 07/08/2013
% 1-D bar element
% Tested using patchbar.m and patchframe.m
% Verified for p=1,2,3 with equal spacing; I don't think unequal spacing or
% p>3 gives the correct results; second derivatives and bubble function not
% yet verified.

PatchE = mateprop(1);
PatchA = mateprop(2);

switch isw %Task Switch
    
    case 1
        
        if ndf > 3
            
            for i = 2:ndf
                lie(i,1) = 0;
            end
            
        end
        
        istv = 1;
        iste = nen;
%%
    case 3

        ElemK = zeros(ndf*nel);
        ElemF = zeros(ndf*nel,1);
        ElemL =norm(diff(xl'));
        
        EA_L=PatchE*PatchA/ElemL; % Get EA/L 
        
        c2=(diff(xl(1,:))/ElemL)^2; %cos^2
        cs=diff(xl(1,:))*diff(xl(2,:))/(ElemL^2); %cos*sin
        s2=1-c2; %sin^2
        % Convert from local to global coordinates
        Local2Globalcoor=[c2,  cs, -c2, -cs
                          cs,  s2, -cs, -s2
                         -c2, -cs,  c2,  cs
                         -cs, -s2,  cs,  s2];
                     
        ElemK=EA_L.*Local2Globalcoor;
        
%         Nmat = zeros(1,ndf*nel);
%         Bmat = zeros(1,ndf*nel);
%         
%         Dmat = PatchE;
%         ulres = reshape(ul,ndf*nen,1);
% 
%         % Load Gauss Points for quadrature
%         lint = nel - 1;
% %         if nel == 2
% %             lint = lintt3;%13;
% %         elseif nel == 3
% %             lint = lintq4;
% %         elseif nel == 4
% %             lint = lintt6;%13;
% %         else
% %             lint = lintq9;
% %         end
% 
%         der = 0;
%         bf = 0;
%         ib = 0;
%         sw = int1d(lint);
%         [shl,shld,shls,be] = shl1d(sw(1,:),lint,nel-1);
%         
%         for je = 1:lint
% 
%             Wgt = sw(2,je);
%             [shg, shgs, Jdet, be] = shg1d(xl,ndm,nel,shld(je,:),shls(je,:),nen,bf,der,be);
%                     
%             c1 = Wgt*Jdet*PatchA;
%                 
%             Nmat(1:ndf:ndf*(nel-1)+1) = shl(je,:);
%                 
%             Bmat(1:ndf:ndf*(nel-1)+1) = shg';
%             
%             sigma = Dmat*Bmat*ulres;
%             
%             ElemF = ElemF - c1*Bmat'*sigma;
%             ElemK = ElemK + c1*Bmat'*Dmat*Bmat;
%             
%         end %je
ElemK
        
%%
    case 6

        ElemF = zeros(ndf*nel,1);
        Nmat = zeros(1,ndf*nel);
        Bmat = zeros(1,ndf*nel);
        
        Dmat = PatchE;
        ulres = reshape(ul,ndf*nen,1);

        % Load Gauss Points for quadrature
        lint = nel - 1;
%         if nel == 2
%             lint = lintt3;%13;
%         elseif nel == 3
%             lint = lintq4;
%         elseif nel == 4
%             lint = lintt6;%13;
%         else
%             lint = lintq9;
%         end

        der = 0;
        bf = 0;
        ib = 0;
        sw = int1d(lint);
        [shl,shld,shls,be] = shl1d(sw(1,:),lint,nel-1);
        
        for je = 1:lint

            Wgt = sw(2,je);
            [shg, shgs, Jdet, be] = shg1d(xl,ndm,nel,shld(je,:),shls(je,:),nen,bf,der,be);
                    
            c1 = Wgt*Jdet*PatchA;
                
            Nmat(1:ndf:ndf*(nel-1)+1) = shl(je,:);
                
            Bmat(1:ndf:ndf*(nel-1)+1) = shg';
            
            sigma = Dmat*Bmat*ulres;
            
            ElemF = ElemF - c1*Bmat'*sigma;
            
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
                    
            c1 = Wgt*Jdet*PatchA;
                
            Nmat(1:ndf:ndf*(nel-1)+1) = shl(je,:);
            
            ElemF = ElemF + c1*(Nmat'*fb);
                
        end %je

%%
    case 11

        ElemE = zeros(numEn,1);
        
        % Load Gauss Points for quadrature
        if nel == 3
            lint = lintt3;%13;
        elseif nel == 4
            lint = lintq4;
        elseif nel == 6
            lint = lintt6;%13;
        elseif nel == 9
            lint = lintq9;
        end

        fb = 0;
        ib = 0;
        bf = 0;
        der = 0;

        el2el = zeros(2,1);
        eprixel = zeros(2,1);
        epriyel = zeros(2,1);
        ue = zeros(2,1);
        duex = ue;
        duey = ue;

        for je = 1:lint

            if nel == 3 || nel == 6
                [Wgt,litr,lits] = intpntt(je,lint,ib);
                [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
                [shg, shgs, Jdet, be, xs] = shgt(xl,nel,shld,shls,nen,bf,der,be);
            elseif nel == 4 || nel == 9
                [Wgt,litr,lits] = intpntq(je,lint,ib);
                [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
                [shg, shgs, Jdet, be, xs] = shgq(xl,nel,shld,shls,nen,bf,der,be);
            end
            
            c1 = Wgt*Jdet*thick;

            xint = xl(1,1:nel)*shl;
            yint = xl(2,1:nel)*shl;
                
            dux = ul(:,1:nel)*shg(:,1);
            duy = ul(:,1:nel)*shg(:,2);
            u = ul(:,1:nel)*shl;

            %Evaluate residual

            %Compute value of exact fields at int. point
            if iprob == 2
                [ue,duex,duey] = uexactbb(xint,yint,PatchE,Patchv,Dbeam);
            elseif iprob == 3
                [ue,duex,duey] = uexacttb(xint,yint,PatchE,Patchv,alph,kT,rhoT,cT,q,qtau,h);
            elseif iprob == 5
                [ue,duex,duey] = uexact_selfw(xint,yint,PatchE,Patchv,rho);
            end

            %Add standard int. point error to element standard error
            for in = 1:2
                un   = c1 * ( (u(in)-ue(in))^2 );
                upnx   = c1 * ( (dux(in)-duex(in))^2 );
                upny   = c1 * ( (duy(in)-duey(in))^2 );
                el2el(in)   = el2el(in)   + un;
                eprixel(in) = eprixel(in) + upnx;
                epriyel(in) = epriyel(in) + upny;
            end

        end %je

        for in= 1:2
            ElemE(in) = el2el(in);
            ElemE(in+3) = eprixel(in);
            ElemE(in+6) = epriyel(in);
        end
        
        H1u = eprixel(1)+epriyel(1)+eprixel(2)+epriyel(2);
        Ieffvals(elem,:) = [1 0 H1u 0 0];

%%        
    case 25 %Stress Projection2

        ElemS = zeros(nel,npstr+1);
        ElemS2 = zeros(nel,npstr);

        Dmat = PatchE;
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
            
%             for stres = 1:npstr
%                 
%                 sigmas = ElemS2(1:nint,stres)'*shpS;
%                 ElemS(ll,stres) = sigmas;
%                 
%             end
            ElemS(ll,1:npstr) = (shpS*ElemS2(1:nint,:))';
            
        end
        
%         %Integration Loop
%         Vol = xl(2) - xl(1);

        for i = 1:nel
        ElemS(i,npstr+1) = 1; % use this for simple average Vol; % use this for weighted average 
        end
        
end %Task Switch
