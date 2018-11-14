% Tim Truster
% 10/13/2013
% Stabilized DG for Poisson problem
% Revised implementation according to VMS derivations, using MVT to pull
% terms outside integrals
% Adapted from IVMS folder

switch isw %Task Switch
    
    case 1
        
        if ndf > 1
            
            for i = 2:ndf
                lie(i,1) = 0;
            end
            
        end

%%
    case 3
 
        A11 = mateprop(1);
        A22 = mateprop(2);
        A12 = mateprop(3);
        Av = [A11 A22 2*A12];
        A = [A11 A12; A12 A22];
        thick = 1;

        ElemK = zeros(nst);
        ElemF = zeros(nst,1);
        Nvec = zeros(1,nst);
        Bmat = zeros(2,nst);

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

        ideriv = 0;
        fb = 0;
        der = 0;
        bf = 0;
        ib = 0;
        
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

            xint = xl(1,:)*shl;
            yint = xl(2,:)*shl;

            if iprob == 1

                L = 1; % CHANGE THIS IF OTHERWISE
                pL = pi/L;
                ddu = [-pL^2*sin(pL*xint)*sin(pL*yint)
                       -pL^2*sin(pL*xint)*sin(pL*yint)
                        pL^2*cos(pL*xint)*cos(pL*yint)];
                fb = -Av*ddu;

            elseif iprob == 2

                L = 1; % CHANGE THIS IF OTHERWISE
                yg = 0.25;
                pL = pi/L;
                ddu = [-pL^2*sin(pL*xint)*sin(pL*yint)
                       -pL^2*sin(pL*xint)*sin(pL*yint)
                        pL^2*cos(pL*xint)*cos(pL*yint)];
                if yint > yg
                    ddu(2) = ddu(2) + (pi^3*cos((pi*yg)/L)*sin((pi*xint)/L)*(MateT(1,1) - MateT(2,1))/MateT(2,1)*(yg - yint))/L^3;
                end
                fb = -Av*ddu;

            end
            
            for i = 1:nel
                Nvec(i) = shl(i);
                Bmat(1,i) = shg(i,1);
                Bmat(2,i) = shg(i,2);
            end
            
            ElemF = ElemF + c1*(fb*Nvec' - Bmat'*A*Bmat*ul');
            ElemK = ElemK + c1*Bmat'*A*Bmat;
            
        end %je
% ElemK

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

        ideriv = 0;
        fb = 0;
        ib = 0;
        bf = 1;
        der = 1;

        el2el = zeros(1,1);
        eprixel = zeros(1,1);
        epriyel = zeros(1,1);
        ue = zeros(1,1);
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

            xint = xl(1,:)*shl;
            yint = xl(2,:)*shl;
                
            if iprob == 1

                L = 1; % CHANGE THIS IF OTHERWISE
                pL = pi/L;
                ddu = [-pL^2*sin(pL*xint)*sin(pL*yint)
                       -pL^2*sin(pL*xint)*sin(pL*yint)
                        pL^2*cos(pL*xint)*cos(pL*yint)];
                fb = -Av*ddu;

            elseif iprob == 2

                L = 1; % CHANGE THIS IF OTHERWISE
                yg = 0.25;
                pL = pi/L;
                ddu = [-pL^2*sin(pL*xint)*sin(pL*yint)
                       -pL^2*sin(pL*xint)*sin(pL*yint)
                        pL^2*cos(pL*xint)*cos(pL*yint)];
                if yint > yg
                    ddu(2) = ddu(2) + (pi^3*cos((pi*yg)/L)*sin((pi*xint)/L)*(MateT(1,1) - MateT(2,1))/MateT(2,1)*(yg - yint))/L^3;
                end
                fb = -Av*ddu;

            end
                
            dux = ul*shg(:,1);
            duy = ul*shg(:,2);
            u = ul*shl;

            %Evaluate residual

            %Compute value of exact fields at int. point
            if iprob == 1
                [ue,duex,duey] = uexact_pois(xint,yint,L);
            elseif iprob == 2
                [ue,duex,duey] = uexact_bimat(xint,yint,L,yg,MateT(1,1),MateT(2,1),1);
            end

            %Add standard int. point error to element standard error
            for in = 1:1
                un   = c1 * ( (u(in)-ue(in))^2 );
                upnx   = c1 * ( (dux(in)-duex(in))^2 );
                upny   = c1 * ( (duy(in)-duey(in))^2 );
                el2el(in)   = el2el(in)   + un;
                eprixel(in) = eprixel(in) + upnx;
                epriyel(in) = epriyel(in) + upny;
            end

        end %je

        for in= 1:1
            ElemE(in) = el2el(in);
            ElemE(in+1) = eprixel(in);
            ElemE(in+2) = epriyel(in);
        end
        
        H1u = eprixel(1)+epriyel(1);
        Ieffvals(elem,:) = [1 0 H1u];
%%
    case -1 % boundary tractions
        
        ElemF = zeros(nst,1);
        
        % Reorder nodes for corresponding face of integration
        SurfOrientEtoS

        % Determine bounds of integration
        
        if nel == 4 || nel == 9
            
        dr = 2;
        ro = -1;
        
        % Upper Limit
        if nodeA == ElemFlag(2)
            eR2 = 1;
        elseif nodeA == ElemFlag(nel2)
            eR2 = ep;
        elseif nodeA == ElemFlag(5) %pressure midway through edge
            eR2 = ep;
        else %
            eR2 = 0;
        end
        % Lower Limit
        if nodeB == ElemFlag(1)
            eR1 = -1;
        elseif nodeB == ElemFlag(nel2)
            eR1 = ep;
        elseif nodeA == ElemFlag(5) %pressure midway through edge
            eR1 = ep;
        else %nodeA == ElemFlag(5)
            eR1 = 0;
        end
        
        elseif nel == 3 || nel == 6
            
        dr = 1;
        ro = 0;
            
        % Upper Limit
        if nodeA == ElemFlag(2)
            eR2 = 1;
        elseif nodeA == ElemFlag(nel2)
            eR2 = ep;
        elseif nodeA == ElemFlag(4) %pressure midway through edge
            eR2 = ep;
        else %nodeA == ElemFlagR(5)
            eR2 = 1/2;
        end
        % Lower Limit
        if nodeB == ElemFlag(1)
            eR1 = 0;
        elseif nodeB == ElemFlag(nel2)
            eR1 = ep;
        elseif nodeB == ElemFlag(4) %pressure midway through edge
            eR1 = ep;
        else %nodeA == ElemFlag(5)
            eR1 = 1/2;
        end
        
        end
        
        % Set jacobian for integration space
        drdr = (eR2 - eR1)/dr;
        
        lint = 4;
%         % Load Gauss Points for quadrature
%         if enrich == 1
%             [rlist, rWgts, rnum] = GaussPoints(pr+2);
%             slist = -1;
%             sWgts = 1;
%             snum = 1;
%         else
%             [rlist, rWgts, rnum] = GaussPoints(pr+1);
%             slist = -1;
%             sWgts = 1;
%             snum = 1;
%         end

%         lamda = Patchv*PatchE/((1+Patchv)*(1-2*Patchv));
%         mu = PatchE/(2*(1+Patchv));
%         thick = 1;
        ideriv = 0;
        der = 0;
        bf = 0;

        Nmat = zeros(2,ndf*nel);
        
        for je = 1:lint

            if nel == 3 || nel == 6
                [Wgt,r,s] = intpntt(je,lint,1);
            else %if nel == 4
                [Wgt,r,s] = intpntq(je,lint,1);
            end
            r = drdr*(r-ro)+eR1;

            if nel == 3 || nel == 6
                [shl,shld,shls,be] = shlt(r,s,nel,nel,der,bf);
                [shg, shgs, Jdet, be, xs] = shgt(xl,nel,shld,shls,nen,bf,der,be);
            elseif nel == 4 || nel == 9
                [shl,shld,shls,be] = shlq(r,s,nel,nel,der,bf);
                [shg, shgs, Jdet, be, xs] = shgq(xl,nel,shld,shls,nen,bf,der,be);
            end
                    
            %Evaluate tangent and normal vectors
            t1 = [xs(:,1); 0];
            [tm1, tu1] = VecNormalize(t1);
            t2 = [0; 0; 1];
            tm2 = 1;
            tu2 = t2';
            t3 = VecCrossProd(t1,t2);
            [tm3, tu3] = VecNormalize(t3);
            
            Traction = traction(1);
            
            c1 = Wgt*tm3*drdr;
            
            Nmat = shl';
            
            ElemF = ElemF + c1*Nmat'*Traction';

        end %ie
        
        % Reorder nodes back to the orientation of the element
        SurfOrientStoE
        ElemF(1:ndf*nel) = ElemF(ilist2);

%%
    case 25 %Stress Projection2

        ElemS = zeros(nel,npstr+1);
        ElemS2 = zeros(nel,npstr);

        A11 = mateprop(1);
        A22 = mateprop(2);
        A12 = mateprop(3);
        Av = [A11 A22 2*A12];
        A = [A11 A12; A12 A22];
        thick = 1;
%         Nmat = zeros(2,2*nel);
%         Bmat = zeros(3,2*nel);
        
        % Load Guass Integration Points

        if nel == 3
            lint = 1;
            nint = 1;
        elseif nel == 4
%             lint = 4;
            lint = 4;
            nint = 1;
        elseif nel == 6
            lint = 7;
            nint = 3;
        else
            lint = 9;
            nint = 4;
        end
        
        der = 0;
        bf = 0;
        ib = 0;
        
        ulres = reshape(ul,ndf*nel,1);

        %Stress Loop
        for ll = 1:nint

            %Evaluate first derivatives of basis functions at int. point
            if nel == 3 || nel == 6
              [w,litr,lits] =  intpntt(ll,nint,ib);
              [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
              [shg,shgs,Jdet,be] = shgt(xl,nel,shld,shls,nel,bf,der,be);
            else
              [w,litr,lits] =  intpntq(ll,nint,ib);
              [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
              [shg,shgs,Jdet,be] = shgq(xl,nel,shld,shls,nel,bf,der,be);
            end
            
            Bmat = shg';
            
            epsil = Bmat*ulres;
            sigma = A*epsil;
            
            ElemS2(ll,1:2) = sigma';

        end %je
        
        % interpolate stress at nodes
        if nel == 3
            plist = [0 1 0
                     0 0 1];
        elseif nel == 4
            plist = [-1 1 1 -1
                     -1 -1 1 1];
        elseif nel == 6
            plist = [-1/3 5/3 -1/3 2/3 2/3 -1/3
                     -1/3 -1/3 5/3 -1/3 2/3 2/3];
        else
            sqr3 = sqrt(3);
            plist = [-sqr3 sqr3 sqr3 -sqr3 0 sqr3 0 -sqr3 0
                     -sqr3 -sqr3 sqr3 sqr3 -sqr3 0 sqr3 0 0];
        end
        
        for ll = 1:nelS
            
            r = plist(1,ll);
            s = plist(2,ll);
            shpS = sshp2d(r,s,nint);
            
%             for stres = 1:npstr
%                 
%                 sigmas = ElemS2(1:nint,stres)'*shpS;
%                 ElemS(ll,stres) = sigmas;
%                 
%             end
            ElemS(ll,1:npstr) = (ElemS2(1:nint,:)'*shpS)';
            
        end
        
%         %Integration Loop
%         Vol = 0;
%         for ll = 1:lint
% 
%             %Evaluate first derivatives of basis functions at int. point
%             if nel == 3 || nel == 6
%               [Wgt,litr,lits] =  intpntt(ll,lint,ib);
%               [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
%               [Qxy, shgs, Jdet,be] = shgt(xl+ul(1:2,:),nel,shld,shls,nel,bf,der,be);
%             else
%               [Wgt,litr,lits] =  intpntq(ll,lint,ib);
%               [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
%               [Qxy, shgs, Jdet,be] = shgq(xl+ul(1:2,:),nel,shld,shls,nel,bf,der,be);
%             end
% 
%             w = Wgt*Jdet*thick;
%             
%             Vol = Vol + w;
% 
%         end %je

        for i = 1:nel
        ElemS(i,npstr+1) = 1; % use this for simple average Vol; % use this for weighted average 
        end

end %Task Switch
