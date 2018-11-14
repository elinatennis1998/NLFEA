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
        
        if ndf > 2
            
            for i = 3:ndf
                lie(i,1) = 0;
            end
            
        end
        
        istv = 1;
        iste = nen;
%%
    case 3

        ElemK = zeros(ndf*nel);
        ElemF = zeros(ndf*nel,1);
        Nmat = zeros(1,ndf*nel);
        Bmat = zeros(1,ndf*nel);
        
        % Transformation
        d_ex = (xl(1,2)-xl(1,1));
        d_ey = (xl(2,2)-xl(2,1));
        Lelem = sqrt((d_ex)^2 + (d_ey)^2);
        lx = d_ex/Lelem;
        ly = d_ey/Lelem;
        if ndf == 3
            Tmat = [lx ly 0   0   0 0
                   -ly lx 0   0   0 0
                    0  0  1   0   0 0
                    0  0  0  lx  ly 0
                    0  0  0 -ly  lx 0
                    0  0  0   0   0 1];
        else
            Tmat = [lx ly   0   0
                   -ly lx   0   0
                    0  0   lx  ly
                    0  0  -ly  lx];
        end
        
        Dmat = PatchE;
        ulres = Tmat*reshape(ul,ndf*nen,1);

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
            
            sigma = Dmat*Bmat*ulres(1:ndf*nel);
            
            ElemF = ElemF - c1*Bmat'*sigma;
            ElemK = ElemK + c1*Bmat'*Dmat*Bmat;
            
        end %je
ElemF = Tmat'*ElemF;
ElemK = Tmat'*ElemK*Tmat;
        
%%
    case 6

        ElemF = zeros(ndf*nel,1);
        Nmat = zeros(1,ndf*nel);
        Bmat = zeros(1,ndf*nel);
        
        % Transformation
        d_ex = (xl(1,2)-xl(1,1));
        d_ey = (xl(2,2)-xl(2,1));
        Lelem = sqrt((d_ex)^2 + (d_ey)^2);
        lx = d_ex/Lelem;
        ly = d_ey/Lelem;
        if ndf == 3
            Tmat = [lx ly 0   0   0 0
                   -ly lx 0   0   0 0
                    0  0  1   0   0 0
                    0  0  0  lx  ly 0
                    0  0  0 -ly  lx 0
                    0  0  0   0   0 1];
        else
            Tmat = [lx ly   0   0
                   -ly lx   0   0
                    0  0   lx  ly
                    0  0  -ly  lx];
        end
        
        Dmat = PatchE;
        ulres = Tmat*reshape(ul,ndf*nen,1);

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
ElemF = Tmat'*ElemF;
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
        
        ulres = reshape(ul,ndf*nen,1);

        %Stress Loop
        
        sw = int1d(nint);
        [shl,shld,shls,be] = shl1d(sw(1,:),lint,nel-1);
        
        
        for ll = 1:nint                    

            % Evaluate 1-D basis functions at integration points
            Wgt = sw(2,ll);
            [shg, shgs, Jdet, be] = shg1d(xl,ndm,nel,shld(ll,:),shls(ll,:),nen,bf,der,be);
            
            % Form B matrix
            Bmat(1:ndf:ndf*(nel-1)+1) = shg';
            
            epsil = Bmat*ulres(1:ndf*nel);
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
%%
    case 26
        
        ElemS = zeros(nestr,1);

        ElemF = zeros(ndf*nel,1);
        Nmat = zeros(1,ndf*nel);
        Bmat = zeros(1,ndf*nel);
        
        % Transformation
        d_ex = (xl(1,2)-xl(1,1));
        d_ey = (xl(2,2)-xl(2,1));
        Lelem = sqrt((d_ex)^2 + (d_ey)^2);
        lx = d_ex/Lelem;
        ly = d_ey/Lelem;
        if ndf == 3
            Tmat = [lx ly 0   0   0 0
                   -ly lx 0   0   0 0
                    0  0  1   0   0 0
                    0  0  0  lx  ly 0
                    0  0  0 -ly  lx 0
                    0  0  0   0   0 1];
        else
            Tmat = [lx ly   0   0
                   -ly lx   0   0
                    0  0   lx  ly
                    0  0  -ly  lx];
        end
        
        Dmat = PatchE;
        ulres = Tmat*reshape(ul,ndf*nen,1);

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
        
        ElemS(1:2*ndf) = ElemF;
        
end %Task Switch
