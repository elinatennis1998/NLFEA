% Tim Truster
% 09/06/2012
% Mixed, stabilized velocity-pressure element for Darcy flow,
% CG or DG treatment. Part of Inter_FEA_Program.
% T3, Q4, T6, and Q9 supported
% Based off of original in Darcy Flow folder
% Updated 9/19/2012 with average bubble definitions

% Sign convention uses pressure positive in compression, as adopted in the
% IJNMF paper.

kappaE = mateprop(1);
muE = mateprop(2);
gcE = mateprop(3);
rhoE = mateprop(4);
thick = 1;
w = 1.05;6;
alp = 0.5;

nelV = nel;
alpha = 0;1;

Bcol1 = [4; 6];
Bcol2 = [5; 7];
col1 = [1; 2];
col2 = [1; 2];

switch isw %Task Switch
%%
    case 3
 
        ElemK = zeros(nst);
        ElemF = zeros(nst,1);
        Bmat = zeros(9,nst);

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

        ib = 0;
        der = 0;
        bf = 0;
        gx = 0;
        gy = 0;
        psi = 0;
        phi = 0;
        h = xl(1,2) - xl(1,1);
        
        if iprob == -1 %Real Darcy Problem
            if numel > 13200
                % For refined mesh, find out which coarse element the fine
                % element is in and grab those material properties
                if mod(elem,60*hDS) == 0
                    ie = 60*hDS;
                else
                    ie = mod(elem,60*hDS);
                end
                je = ceil(elem/(60*hDS));
                Ie = ceil(ie/hDS);
                Je = ceil(je/hDS);
                ElemDS = (Je-1)*60 + Ie;
                kappaE = PermTable(ElemDS,1); % Smooth Data
    %             kappaE = PermTable(ElemDS+58*220*60,1); % Channelized Data
                muE = 0.3;
            else
                kappaE = PermTable(elem,1); % Smooth Data
    %             kappaE = PermTable(elem+58*220*60,1); % Channelized Data
                muE = 0.3;
            end
        end
        muEkappaE = muE/kappaE;
        kappaEmuE = kappaE/muE;
        rhoEgcE = rhoE/gcE;
        gvec = -1/2*rhoEgcE*[gx; gy; 0];
            
        Dmat = [1/2*muEkappaE*eye(2) zeros(2,1) zeros(2,4) -1/2*eye(2,2)
                zeros(1,2) 0 [1 0 0 1] zeros(1,2)
                zeros(4,2) -[1; 0; 0; 1] alpha/2*muEkappaE*h^2*([1; 0; 0; 1]*[1 0 0 1]) zeros(4,2)
                1/2*eye(2) zeros(2,1) zeros(2,4) 1/2*kappaEmuE*eye(2)];
        Dvec = [1/2*rhoEgcE*eye(2) zeros(2,1)
                zeros(1,2) 1
                zeros(4,2) alpha/2*muEkappaE*h^2*[1; 0; 0; 1]
                -1/2*kappaEmuE*eye(2) zeros(2,1)];
        
        ulres = reshape(ul,ndf*nen,1);
        
        for je = 1:lint

            if nel == 3 || nel == 6
                [Wgt,litr,lits] = intpntt(je,lint,ib);
                [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
                [shg, shgs, Jdet, be, xs] = shgt(xl,nel,shld,shls,nen,bf,der,be);
                [shlp,shld,shls,be] = shlt(litr,lits,nelP,nel,der,bf);
                shp = shgt(xl,nelP,shld,shls,nen,bf,der,be);
            elseif nel == 4 || nel == 9
                [Wgt,litr,lits] = intpntq(je,lint,ib);
                [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
                [shg, shgs, Jdet, be, xs] = shgq(xl,nel,shld,shls,nen,bf,der,be);
                [shlp,shld,shls,be] = shlq(litr,lits,nelP,nel,der,bf);
                shp = shgq(xl,nelP,shld,shls,nen,bf,der,be);
            end
            
            c1 = Wgt*Jdet*thick;

            xint = xl(1,1:nelV)*shl;
            yint = xl(2,1:nelV)*shl;

            if iprob == 1

                phi = 2*kappaEmuE*(2*pi)^2*sin(2*pi*xint)*sin(2*pi*yint);
                gvec(3) = phi;
                
            elseif iprob == 2

                phi = 2*kappaEmuE*(2*pi)^2*cos(2*pi*xint)*cos(2*pi*yint);
                gvec(3) = phi;
                
            elseif iprob == 4
                
                % NOTE: never verified relation between mu and K in the
                % formulas
                if alp == 0
                gvec = [ -((w*yint*cos(w*xint)*(muE - 1)))
                         -(- ((muE - 1)*(2*yint - 48*sin(w*xint) + 1))/(48) - (((muE - 1)*(16*alp + 32*alp*yint - 48*alp*sin(w*xint)))/48 - (alp*(muE - 1)*(2*yint - 48*sin(w*xint) + 1))/48)/((alp + (muE)^(1/2))))
                         (-(5/8)*1/2 - w^2*yint*sin(w*xint) - 17/48)];
                else
                gvec = [ ((w*yint*cos(w*xint)*(muE - 1)))
                         (- ((muE - 1)*(2*yint - 48*sin(w*xint) + 1))/(48) - (((muE - 1)*(16*alp + 32*alp*yint - 48*alp*sin(w*xint)))/48 - (alp*(muE - 1)*(2*yint - 48*sin(w*xint) + 1))/48)/((alp + (muE)^(1/2))))
                         ((5*((muE)^(1/2)/alp - 1))/(8*((2*(muE)^(1/2))/alp + 2)) - w^2*yint*sin(w*xint) - 17/48)];
                end
                
            elseif iprob == 5
                
                w = 1.05;
                G = sqrt(1)/alp;
                gvec = [  (exp(yint/G)*sin(w + xint/G)*(muE - 1))/1
                         -(exp(yint/G)*cos(w + xint/G)*(muE - 1))/1
                         0];
                
            elseif iprob == 6 % Final Stokes-Darcy problem used in paper
                
                G = sqrt(1)/alp;
                gvec = [ 0
                         0
                         w^2*cos(w*xint)];

            end

            % Form B matrix
            for ie = 1:nel
                
              Bmat(1,(ie-1)*3+1) = shl(ie);
              Bmat(2,(ie-1)*3+2) = shl(ie);
                
              Bmat(Bcol1,(ie-1)*3+1) = shg(ie,col1);
              Bmat(Bcol2,(ie-1)*3+2) = shg(ie,col2);
                 
            end
            for ie = 1:nelP
                
              Bmat(3,(ie-1)*3+3) = shlp(ie);
                
              Bmat(8,(ie-1)*3+3) = shp(ie,1);
              Bmat(9,(ie-1)*3+3) = shp(ie,2);
                 
            end
            
            vpfields = Bmat*ulres;%(1:ndf*nel);
            
            ElemF = ElemF + c1*(Bmat'*Dvec*gvec - Bmat'*Dmat*vpfields);
            ElemK = ElemK + c1*(Bmat'*Dmat*Bmat);
            
        end
ElemK;
%%
    case 15 %Body force
 
        if ismember(iprob,[1 2 4 5 6])
            
%         ElemF = zeros(nst,1);
        Bmat = zeros(9,nst);

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

        ib = 0;
        der = 0;
        bf = 0;
        gx = 0;
        gy = 0;
        psi = 0;
        phi = 0;
        h = xl(1,2) - xl(1,1);
        
        muEkappaE = muE/kappaE;
        kappaEmuE = kappaE/muE;
        rhoEgcE = rhoE/gcE;
        gvec = -1/2*rhoEgcE*[gx; gy; 0];
            
        Dmat = [1/2*muEkappaE*eye(2) zeros(2,1) zeros(2,4) -1/2*eye(2,2)
                zeros(1,2) 0 [1 0 0 1] zeros(1,2)
                zeros(4,2) -[1; 0; 0; 1] alpha/2*muEkappaE*h^2*([1; 0; 0; 1]*[1 0 0 1]) zeros(4,2)
                1/2*eye(2) zeros(2,1) zeros(2,4) 1/2*kappaEmuE*eye(2)];
        Dvec = [1/2*rhoEgcE*eye(2) zeros(2,1)
                zeros(1,2) 1
                zeros(4,2) alpha/2*muEkappaE*h^2*[1; 0; 0; 1]
                -1/2*kappaEmuE*eye(2) zeros(2,1)];
        
        for je = 1:lint

            if nel == 3 || nel == 6
                [Wgt,litr,lits] = intpntt(je,lint,ib);
                [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
                [shg, shgs, Jdet, be, xs] = shgt(xl,nel,shld,shls,nen,bf,der,be);
                [shlp,shld,shls,be] = shlt(litr,lits,nelP,nel,der,bf);
                shp = shgt(xl,nelP,shld,shls,nen,bf,der,be);
            elseif nel == 4 || nel == 9
                [Wgt,litr,lits] = intpntq(je,lint,ib);
                [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
                [shg, shgs, Jdet, be, xs] = shgq(xl,nel,shld,shls,nen,bf,der,be);
                [shlp,shld,shls,be] = shlq(litr,lits,nelP,nel,der,bf);
                shp = shgq(xl,nelP,shld,shls,nen,bf,der,be);
            end
            
            c1 = Wgt*Jdet*thick;

            xint = xl(1,1:nelV)*shl;
            yint = xl(2,1:nelV)*shl;

            if iprob == 1

                phi = 2*kappaEmuE*(2*pi)^2*sin(2*pi*xint)*sin(2*pi*yint);
                gvec(3) = phi;
                
            elseif iprob == 2

                phi = 2*kappaEmuE*(2*pi)^2*cos(2*pi*xint)*cos(2*pi*yint);
                gvec(3) = phi;
                
            elseif iprob == 4
                
                % NOTE: never verified relation between mu and K in the
                % formulas
                if alp == 0
                gvec = [ -((w*yint*cos(w*xint)*(muE - 1)))
                         -(- ((muE - 1)*(2*yint - 48*sin(w*xint) + 1))/(48) - (((muE - 1)*(16*alp + 32*alp*yint - 48*alp*sin(w*xint)))/48 - (alp*(muE - 1)*(2*yint - 48*sin(w*xint) + 1))/48)/((alp + (muE)^(1/2))))
                         (-(5/8)*1/2 - w^2*yint*sin(w*xint) - 17/48)];
                else
                gvec = [ ((w*yint*cos(w*xint)*(muE - 1)))
                         (- ((muE - 1)*(2*yint - 48*sin(w*xint) + 1))/(48) - (((muE - 1)*(16*alp + 32*alp*yint - 48*alp*sin(w*xint)))/48 - (alp*(muE - 1)*(2*yint - 48*sin(w*xint) + 1))/48)/((alp + (muE)^(1/2))))
                         ((5*((muE)^(1/2)/alp - 1))/(8*((2*(muE)^(1/2))/alp + 2)) - w^2*yint*sin(w*xint) - 17/48)];
                end
                
            elseif iprob == 5
                
                G = sqrt(1)/alp;
                gvec = [  (exp(yint/G)*sin(w + xint/G)*(muE - 1))/1
                         -(exp(yint/G)*cos(w + xint/G)*(muE - 1))/1
                         0];
                
            elseif iprob == 6 % Final Stokes-Darcy problem used in paper
                
                G = sqrt(1)/alp;
                gvec = [ 0
                         0
                         w^2*cos(w*xint)];

            end

            % Form B matrix
            for ie = 1:nel
                
              Bmat(1,(ie-1)*3+1) = shl(ie);
              Bmat(2,(ie-1)*3+2) = shl(ie);
                
              Bmat(Bcol1,(ie-1)*3+1) = shg(ie,col1);
              Bmat(Bcol2,(ie-1)*3+2) = shg(ie,col2);
                 
            end
            for ie = 1:nelP
                
              Bmat(3,(ie-1)*3+3) = shlp(ie);
                
              Bmat(8,(ie-1)*3+3) = shp(ie,1);
              Bmat(9,(ie-1)*3+3) = shp(ie,2);
                 
            end
            
            ElemF = ElemF + c1*(Bmat'*Dvec*gvec);
            
        end
ElemF;
        end
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
        
        muEkappaE = muE/kappaE;
        kappaEmuE = kappaE/muE;
        rhoEgcE = rhoE/gcE;

        el2el = zeros(3,1);
        eprixel = zeros(3,1);
        epriyel = zeros(3,1);
        ue = zeros(3,1);
        duex = ue;
        duey = ue;

        for je = 1:lint

            if nel == 3 || nel == 6
                [Wgt,litr,lits] = intpntt(je,lint,ib);
                [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
                [shg, shgs, Jdet, be, xs] = shgt(xl,nel,shld,shls,nen,bf,der,be);
                [shlp,shld,shls,be] = shlt(litr,lits,nelP,nel,der,bf);
                shp = shgt(xl,nelP,shld,shls,nen,bf,der,be);
            elseif nel == 4 || nel == 9
                [Wgt,litr,lits] = intpntq(je,lint,ib);
                [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
                [shg, shgs, Jdet, be, xs] = shgq(xl,nel,shld,shls,nen,bf,der,be);
                [shlp,shld,shls,be] = shlq(litr,lits,nelP,nel,der,bf);
                shp = shgq(xl,nelP,shld,shls,nen,bf,der,be);
            end
            c1 = Wgt*Jdet;

            xint = xl(1,1:nelV)*shl;
            yint = xl(2,1:nelV)*shl;
            dux = zeros(3,1);
            duy = zeros(3,1);
            u = zeros(3,1);
            dux(1:2) = ul(1:2,1:nelV)*shg(:,1);
            duy(1:2) = ul(1:2,1:nelV)*shg(:,2);
            u(1:2) = ul(1:2,1:nelV)*shl;
            dux(3) = ul(3,1:nelP)*shp(:,1);
            duy(3) = ul(3,1:nelP)*shp(:,2);
            u(3) = ul(3,1:nelP)*shlp;

            %Compute value of exact fields at int. point
            if iprob == 1
                [ue,duex,duey] = uexact_cs(xint,yint,kappaEmuE);
            elseif iprob == 2
                [ue,duex,duey] = uexact_csb(xint,yint,kappaEmuE);
            elseif iprob == 4
                [ue,duex,duey] = uexact_vasil(xint,yint,muE,1,0.5,w);
            elseif iprob == 5
                [ue,duex,duey] = uexact_vasil1(xint,yint,muE,1,0.5,w);
            elseif iprob == 6
                [ue,duex,duey] = uexact_vasiltim(xint,yint,muE,1,0.5,w);
            elseif iprob == 7
                [ue,duex,duey] = uexact_DPT(xint,yint,muEkappaE);
            end

            %Add standard int. point error to element standard error
            for in = 1:3
                un   = c1 * ( (u(in)-ue(in))^2 );
                upnx   = c1 * ( (dux(in)-duex(in))^2 );
                upny   = c1 * ( (duy(in)-duey(in))^2 );
                el2el(in)   = el2el(in)   + un;
                eprixel(in) = eprixel(in) + upnx;
                epriyel(in) = epriyel(in) + upny;
            end

        end %je

        for in= 1:3
            ElemE(in) = el2el(in);
            ElemE(in+3) = eprixel(in);
            ElemE(in+6) = epriyel(in);
        end
        
%%
    case -1 % boundary terms
        
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
        else %
            eR2 = 0;
        end
        % Lower Limit
        if nodeB == ElemFlag(1)
            eR1 = -1;
        else %nodeA == ElemFlag(5)
            eR1 = 0;
        end
        
        elseif nel == 3 || nel == 6
            
        dr = 1;
        ro = 0;
            
        % Upper Limit
        if nodeA == ElemFlag(2)
            eR2 = 1;
        else %nodeA == ElemFlagR(5)
            eR2 = 1/2;
        end
        % Lower Limit
        if nodeB == ElemFlag(1)
            eR1 = 0;
        else %nodeA == ElemFlag(5)
            eR1 = 1/2;
        end
        
        end
        
        % Set jacobian for integration space
        drdr = (eR2 - eR1)/dr;
        
        lint = 4;
        ib  = 1;
        der = 0;
        bf = 0;
        
        for je = 1:lint

            if nel == 3 || nel == 6
                [Wgt,r,s] = intpntt(je,lint,1);
            else %if nel == 4
                [Wgt,r,s] = intpntq(je,lint,1);
            end
            r = drdr*(r-ro)+eR1;

            if nel == 3 || nel == 6
                [shl,shld,shls,be] = shlt(r,s,nel,nel,der,bf);
                [Qxy, cartd2, Jdet, be, xs] = shgt(xl,nel,shld,shls,nen,bf,der,be);
            elseif nel == 4 || nel == 9
                [shl,shld,shls,be] = shlq(r,s,nel,nel,der,bf);
                [Qxy, cartd2, Jdet, be, xs] = shgq(xl,nel,shld,shls,nen,bf,der,be);
            end
                    
            %Evaluate tangent and normal vectors
            t1 = [xs(:,1); 0];
            [tm1, tu1] = VecNormalize(t1);
            t2 = [0; 0; 1];
            tm2 = 1;
            tu2 = t2';
            t3 = VecCrossProd(t1,t2);
            [tm3, tu3] = VecNormalize(t3);
            
            if iprob == 5
                
                alp = 0.5;
                w = 1.05;
                xint = xl(1,1:nel)*shl;
                yint = xl(2,1:nel)*shl;
                G = sqrt(1)/alp;
                pd = G/1*cos(xint/G+w)*exp(yint/G);
                Traction(1) = -pd*tu3(1);
                Traction(2) = -pd*tu3(2); % minus pressure is important
                
            elseif iprob == 6
                
                alp = 0.5;
                w = 1.05;
                mu = 1;
                xint = xl(1,1:nel)*shl;
                yint = xl(2,1:nel)*shl;
                G = sqrt(1)/alp;
                pd = mu*G/1*cos(xint/G+w)*exp(yint/G) + mu/1*cos(w*xint);
                Traction(1) = -pd*tu3(1);
                Traction(2) = -pd*tu3(2); % minus pressure is important
                
            elseif iprob == 7
                
                xint = xl(1,1:nel)*shl;
                yint = xl(2,1:nel)*shl;
                pd = 2*xint*(yint-1);
                Traction(1) = -pd*tu3(1);
                Traction(2) = -pd*tu3(2); % minus pressure is important
                
            else
                Traction = traction;
            end
            
            c1 = Wgt*tm3*drdr*thick;
            
            for o=1:nel
                
                don = shl(o);
                
                F = don*Traction';

                ElemF(ndf*(o-1)+1) = ElemF(ndf*(o-1)+1) + F(1)*c1;

                ElemF(ndf*(o-1)+2) = ElemF(ndf*(o-1)+2) + F(2)*c1;

            end %o

        end %ie
        ElemF;
        
        % Reorder nodes back to the orientation of the element
        SurfOrientStoE
        ElemF(1:ndf*nel) = ElemF(ilist2);

end %Task Switch
