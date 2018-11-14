
% Stokes element for IJNMF paper
% 05/2013
%
% Sign convention uses pressure positive in compression, as adopted in the
% paper.

PatchE = mateprop(1);
Patchv = mateprop(2);
thick = mateprop(3);
w = 1.05;6;

Bcol1 = [1; 3];
Bcol2 = [2; 3];
col1 = [1; 2];
col2 = [2; 1];

switch isw %Task Switch
    
    case 1
        
        if ndf > 3
            
            for i = 4:ndf
                lie(i,1) = 0;
            end
            
        end
%%
    case 3 %stiffness and force
        
        ElemK = zeros(nst);
        ElemF = zeros(nst,1);
        Nmat = zeros(3,nst);
        Bmat = zeros(4,nst);
        BBmat = zeros(2,nst);
        
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

        lam = Patchv*PatchE/((1+Patchv)*(1-2*Patchv));
        mu = PatchE/(2*(1+Patchv));
        Dmat = mu*diag([2 2 1]);
        Dmat = [Dmat  [1; 1; 0]
                [1 1 0] -1/lam];
        
        fbx = 0;
        fby = 0;

        der = 1;
        bf = 1;
        
        [t11,t12,t21,t22] = Tau3_2d(xl,mu,nel,nen,lint);
        
        ulres = reshape(ul,ndf*nen,1);
        
        for je = 1:lint

            if nel == 3 || nel == 6
                [Wgt,r,s] = intpntt(je,lint,0);
            else %if nel == 4 || nel == 9
                [Wgt,r,s] = intpntq(je,lint,0);
            end
            
            if nel == 3 || nel == 6
                [shl,shld,shls,be] = shlt(r,s,nel,nel,der,bf);
                [shg, shgs, Jdet, be, xs] = shgt(xl,nel,shld,shls,nen,bf,der,be);
            elseif nel == 4 || nel == 9
                [shl,shld,shls,be] = shlq(r,s,nel,nel,der,bf);
                [shg, shgs, Jdet, be, xs] = shgq(xl,nel,shld,shls,nen,bf,der,be);
            end
            c1 = Wgt*Jdet*thick;

            if iprob == 4
                alp = 0;%.5;
                xint = xl(1,1:nel)*shl;
                yint = xl(2,1:nel)*shl;
                fbx = - 2*mu*(2*xint - 4) - (w*yint*cos(w*xint))/(2);
                fby = 2*mu*(w^2*sin(w*xint) - 2*xint + 4) - pi*sin(pi*yint);
                
            elseif iprob == 5
                
                alp = 0.5;
                w = 1.05;
                xint = xl(1,1:nel)*shl;
                yint = xl(2,1:nel)*shl;
                G = sqrt(1)/alp;
                fbx = (exp(1/(2*G))*sin(w + xint/G)*(2*mu - G^2))/(G^2);
                fby = 1;
                
            elseif iprob == 6 % Final Stokes-Darcy problem used in paper
                
                alp = 0.5;
                w = 1.05;
                xint = xl(1,1:nel)*shl;
                yint = xl(2,1:nel)*shl;
                G = sqrt(1)/alp;                 
                fbx =  - (exp(1/(2*G))*sin(w + xint/G)*(G/1 - (2*mu)/G))/G - (2*mu*w*yint*sin(w*xint))/1;
                fby = (2*mu*cos(w*xint))/1 + 1;
                
            end
            fb = [fbx; fby; 0];
            
            Tmat = be(3)*[t11 t12; t21 t22];
            
            % Form B matrix
            for i = 1:nel
                
              Nmat(1,(i-1)*ndf+1) = shl(i);
              Nmat(2,(i-1)*ndf+2) = shl(i);
              Nmat(3,(i-1)*ndf+3) = shl(i);
                
              Bmat(Bcol1,(i-1)*ndf+1) = shg(i,col1);
              Bmat(Bcol2,(i-1)*ndf+2) = shg(i,col2);
              Bmat(4    ,(i-1)*ndf+3) = -shl(i);
                                   
              BBmat(:,(i-1)*ndf+1:(i-1)*ndf+3) = [mu*(2*shgs(i,1)+shgs(i,2)) mu*shgs(i,3) -shg(i,1)
                                          mu*shgs(i,3) mu*(shgs(i,1)+2*shgs(i,2)) -shg(i,2)];
                                    
            end
            Bmatw = diag([1 1 1 -1])*Bmat;
            BBmatw = BBmat;
            BBmatw(:,3:3:3*nel) = -BBmatw(:,3:3:3*nel);
            
            upfields = Bmat*ulres;%(1:ndf*nel);
            
            ElemF = ElemF + c1*(Nmat'*fb + BBmatw'*Tmat*fb(1:2) - Bmatw'*Dmat*upfields + BBmatw'*Tmat*BBmat*ulres);
            ElemK = ElemK + c1*(Bmatw'*Dmat*Bmat - BBmatw'*Tmat*BBmat);
                
        end %je
        
        ElemK;
        
%%
    case 15 % body force for linear element
        
%         ElemF = zeros(nst,1); In FormFE
        
        if ismember(iprob,[4 5 6])
            
        Nmat = zeros(3,nst);
        BBmat = zeros(2,nst);
        
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

        lam = Patchv*PatchE/((1+Patchv)*(1-2*Patchv));
        mu = PatchE/(2*(1+Patchv));
        Dmat = mu*diag([2 2 1]);
        Dmat = [Dmat  [1; 1; 0]
                [1 1 0] -1/lam];

        der = 1;
        bf = 1;
        
        fbx = 0;
        fby = 0;

        [t11,t12,t21,t22] = Tau3_2d(xl,mu,nel,nen,lint);
        
        for je = 1:lint

            if nel == 3 || nel == 6
                [Wgt,r,s] = intpntt(je,lint,0);
            else %if nel == 4 || nel == 9
                [Wgt,r,s] = intpntq(je,lint,0);
            end
            
            if nel == 3 || nel == 6
                [shl,shld,shls,be] = shlt(r,s,nel,nel,der,bf);
                [shg, shgs, Jdet, be, xs] = shgt(xl,nel,shld,shls,nen,bf,der,be);
            elseif nel == 4 || nel == 9
                [shl,shld,shls,be] = shlq(r,s,nel,nel,der,bf);
                [shg, shgs, Jdet, be, xs] = shgq(xl,nel,shld,shls,nen,bf,der,be);
            end
            c1 = Wgt*Jdet*thick;

            if iprob == 4
                alp = 0;%.5;
                xint = xl(1,1:nel)*shl;
                yint = xl(2,1:nel)*shl;
                fbx = - 2*mu*(2*xint - 4) - (w*yint*cos(w*xint))/(2);
                fby = 2*mu*(w^2*sin(w*xint) - 2*xint + 4) - pi*sin(pi*yint);
            elseif iprob == 5
                
                alp = 0.5;
                w = 1.05;
                xint = xl(1,1:nel)*shl;
                yint = xl(2,1:nel)*shl;
                G = sqrt(1)/alp;
                fbx = (exp(1/(2*G))*sin(w + xint/G)*(2*mu - G^2))/(G^2);
                fby = 1;
                
            elseif iprob == 6 % Final Stokes-Darcy problem used in paper
                
                alp = 0.5;
                w = 1.05;
                xint = xl(1,1:nel)*shl;
                yint = xl(2,1:nel)*shl;
                G = sqrt(1)/alp;                 
                fbx =  - (exp(1/(2*G))*sin(w + xint/G)*(G/1 - (2*mu)/G))/G - (2*mu*w*yint*sin(w*xint))/1;
                fby = (2*mu*cos(w*xint))/1 + 1;
                
            end
            fb = [fbx; fby; 0];
            
            Tmat = be(3)*[t11 t12; t21 t22];
            
            % Form B matrix
            for i = 1:nel
                
              Nmat(1,(i-1)*ndf+1) = shl(i);
              Nmat(2,(i-1)*ndf+2) = shl(i);
              Nmat(3,(i-1)*ndf+3) = shl(i);
                                   
              BBmat(:,(i-1)*ndf+1:(i-1)*ndf+3) = [mu*(2*shgs(i,1)+shgs(i,2)) mu*shgs(i,3) shg(i,1)
                                          mu*shgs(i,3) mu*(shgs(i,1)+2*shgs(i,2)) shg(i,2)];
                                    
            end
            
            ElemF = ElemF + c1*(Nmat'*fb + BBmat'*Tmat*fb(1:2));
                
        end %je
        
        end
        
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
                mu = 1;
                xint = xl(1,1:nel)*shl;
                yint = xl(2,1:nel)*shl;
                G = sqrt(1)/alp;
                Traction(1) = (2*mu*exp(1/G)*sin(w + xint/G))/G;
                Traction(2) =  - exp(1/(2*G))*cos(w + xint/G)*(G/1 - (2*mu)/G) - (2*mu*exp(1/G)*cos(w + xint/G))/G - 1/2;
                
            elseif iprob == 6
                
                alp = 0.5;
                w = 1.05;
                mu = 1;
                xint = xl(1,1:nel)*shl;
                yint = xl(2,1:nel)*shl;
                G = sqrt(1)/alp;
                cosw = cos(w + xint/G);
                sinw = sin(w + xint/G);
                eyG = exp(yint/G);
                sigma = [ exp(1/(2*G))*cosw*((2*mu)/G - (G*mu)/1) - yint + (2*mu*eyG*cosw)/G - (2*mu*yint*cos(w*xint))/1 + 1/2, (2*mu*eyG*sinw)/G
                          (2*mu*eyG*sinw)/G, exp(1/(2*G))*cosw*((2*mu)/G - (G*mu)/1) - yint - (2*mu*eyG*cosw)/G - (2*mu*yint*cos(w*xint))/1 + 1/2];
                Traction(1) = sigma(1,:)*tu3(1:2)';
                Traction(2) = sigma(2,:)*tu3(1:2)';
                
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
%%       
    case 11
        
        if elem > 1 && darcystokeserror == 0 % Compute Darcy error separate from Stokes
            
            ErrD = Energy(1:9);
            Fine = Energy(10:15);
            L2uerrD = sqrt(ErrD(1) + ErrD(2));
            L2uerrD = log(L2uerrD)/log(10);
            L2perrD = sqrt(ErrD(3));
            L2perrD = log(L2perrD)/log(10);
            H1uerrD = sqrt(ErrD(4) + ErrD(5) + ErrD(7) + ErrD(8));
            H1uerrD = log(H1uerrD)/log(10);
            L2divuerrD = sqrt(ErrD(4) + ErrD(8));
            L2divuerrD = log(L2divuerrD)/log(10);
            H1perrD = sqrt(ErrD(6) + ErrD(9));
            H1perrD = log(H1perrD)/log(10);
            L2Fine = log(sqrt(Fine(1)+Fine(2)))/log(10);
            H1Fine = log(sqrt(Fine(3)+Fine(4)+Fine(5)+Fine(6)))/log(10);

            fprintf('%i  %1.7e  %1.7e  %1.7e  %1.7e  %i\n',numel,L2uerrD,L2perrD,H1uerrD,H1perrD,iprob)
            fprintf('%i  %1.7e  %1.7e \n',numel,L2Fine,H1Fine)
            
            Energy = 0*Energy;
            darcystokeserror = 1;
    
        end
        
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

        lam = Patchv*PatchE/((1+Patchv)*(1-2*Patchv));
        mu = PatchE/(2*(1+Patchv));
        ideriv = 1;
        fbx = 0;
        fby = 0;
        ib = 0;
        bf = 1;
        der = 1;

        el2el = zeros(3,1);
        eprixel = zeros(3,1);
        epriyel = zeros(3,1);
        el2fine = zeros(7,1);
        ue = zeros(3,1);
        duex = ue;
        duey = ue;

        [t11,t12,t21,t22] = Tau3_2d(xl,mu,nel,nen,lint);
        
        if nel == 4
            r = 0;
            s = 0;
        elseif nel == 3
            r = 1/3;
            s = 1/3;
        end
        
        if nel == 3 || nel == 6
            [shl,shld,shls,be] = shlt(r,s,nel,nel,der,bf);
            [shg, shgs, Jdet, be, xs] = shgt(xl,nel,shld,shls,nen,bf,der,be);
        elseif nel == 4 || nel == 9
            [shl,shld,shls,be] = shlq(r,s,nel,nel,der,bf);
            [shg, shgs, Jdet, be, xs] = shgq(xl,nel,shld,shls,nen,bf,der,be);
        end
        c1 = Wgt*Jdet*thick;
        b = be(3);

        if nelP == 3 || nelP == 6
            [shlp,shld,shls,bub] = shlt(r,s,nelP,nel,0,0);
            [shp, shp2, Jdet, be, xs] = shgt(xl,nelP,shld,shls,nen,0,0,bub);
        elseif nelP == 4 || nelP == 9
            [shlp,shld,shls,bub] = shlq(r,s,nelP,nel,0,0);
            [shp, shp2, Jdet, bub, xs] = shgq(xl,nelP,shld,shls,nen,0,0,bub);
        end

        xint = xl(1,1:nel)*shl;
        yint = xl(2,1:nel)*shl;
        dux = ul(1:2,1:nel)*shg(:,1);
        duy = ul(1:2,1:nel)*shg(:,2);
        u = ul(1:2,1:nel)*shl;
        dux(3) = ul(3,1:nelP)*shp(:,1);
        duy(3) = ul(3,1:nelP)*shp(:,2);
        u(3) = ul(3,1:nelP)*shlp;
        px = dux(3);
        py = duy(3);
        ux_xx = ul(1,1:nel)*shgs(:,1);
        ux_yy = ul(1,1:nel)*shgs(:,2);
        ux_xy = ul(1,1:nel)*shgs(:,3);
        uy_xx = ul(2,1:nel)*shgs(:,1);
        uy_yy = ul(2,1:nel)*shgs(:,2);
        uy_xy = ul(2,1:nel)*shgs(:,3);
        
        if iprob == 4
            alp = 0;%.5;
            xint = xl(1,1:nel)*shl;
            yint = xl(2,1:nel)*shl;
            fbx = - 2*mu*(2*xint - 4) - (w*yint*cos(w*xint))/(2);
            fby = 2*mu*(w^2*sin(w*xint) - 2*xint + 4) - pi*sin(pi*yint);
        elseif iprob == 5

            alp = 0.5;
            w = 1.05;
            xint = xl(1,1:nel)*shl;
            yint = xl(2,1:nel)*shl;
            G = sqrt(1)/alp;
            fbx = (exp(1/(2*G))*sin(w + xint/G)*(2*mu - G^2))/(G^2);
            fby = 1;

        elseif iprob == 6 % Final Stokes-Darcy problem used in paper

            alp = 0.5;
            w = 1.05;
            xint = xl(1,1:nel)*shl;
            yint = xl(2,1:nel)*shl;
            G = sqrt(1)/alp;                 
            fbx =  - (exp(1/(2*G))*sin(w + xint/G)*(G/1 - (2*mu)/G))/G - (2*mu*w*yint*sin(w*xint))/1;
            fby = (2*mu*cos(w*xint))/1 + 1;

        end

        %Evaluate residual of equilibrium equation
        rx = fbx + -px + mu*(2*ux_xx + ux_yy + uy_xy);
        ry = fby + -py + mu*(uy_xx + ux_xy + 2*uy_yy);

        %Evaluate explicit fine scale
        bubblevals(elem,1) = (t11*rx+t12*ry)*b;
        bubblevals(elem,2) = (t21*rx+t22*ry)*b;

        for je = 1:lint

            if nel == 3 || nel == 6
                [Wgt,r,s] = intpntt(je,lint,0);
            elseif nel == 4 || nel == 9
                [Wgt,r,s] = intpntq(je,lint,0);
            end

            if nel == 3 || nel == 6
                [shl,shld,shls,be] = shlt(r,s,nel,nel,der,bf);
                [shg, shgs, Jdet, be, xs] = shgt(xl,nel,shld,shls,nen,bf,der,be);
            elseif nel == 4 || nel == 9
                [shl,shld,shls,be] = shlq(r,s,nel,nel,der,bf);
                [shg, shgs, Jdet, be, xs] = shgq(xl,nel,shld,shls,nen,bf,der,be);
            end
            c1 = Wgt*Jdet*thick;
            b = be(3);

            if nelP == 3 || nelP == 6
                [shlp,shld,shls,bub] = shlt(r,s,nelP,nel,0,0);
                [shp, shp2, Jdet, bub, xs] = shgt(xl,nelP,shld,shls,nen,0,0,bub);
            elseif nelP == 4 || nelP == 9
                [shlp,shld,shls,bub] = shlq(r,s,nelP,nel,0,0);
                [shp, shp2, Jdet, bub, xs] = shgq(xl,nelP,shld,shls,nen,0,0,bub);
            end

            xint = xl(1,1:nel)*shl;
            yint = xl(2,1:nel)*shl;
            dux = ul(1:2,1:nel)*shg(:,1);
            duy = ul(1:2,1:nel)*shg(:,2);
            u = ul(1:2,1:nel)*shl;
            dux(3) = ul(3,1:nelP)*shp(:,1);
            duy(3) = ul(3,1:nelP)*shp(:,2);
            u(3) = ul(3,1:nelP)*shlp;
            px = dux(3);
            py = duy(3);
            ux_xx = ul(1,1:nel)*shgs(:,1);
            ux_yy = ul(1,1:nel)*shgs(:,2);
            ux_xy = ul(1,1:nel)*shgs(:,3);
            uy_xx = ul(2,1:nel)*shgs(:,1);
            uy_yy = ul(2,1:nel)*shgs(:,2);
            uy_xy = ul(2,1:nel)*shgs(:,3);
            
            if iprob == 4
                alp = 0;%.5;
                xint = xl(1,1:nel)*shl;
                yint = xl(2,1:nel)*shl;
                fbx = - 2*mu*(2*xint - 4) - (w*yint*cos(w*xint))/(2);
                fby = 2*mu*(w^2*sin(w*xint) - 2*xint + 4) - pi*sin(pi*yint);
            elseif iprob == 5
                
                alp = 0.5;
                w = 1.05;
                xint = xl(1,1:nel)*shl;
                yint = xl(2,1:nel)*shl;
                G = sqrt(1)/alp;
                fbx = (exp(1/(2*G))*sin(w + xint/G)*(2*mu - G^2))/(G^2);
                fby = 1;
                
            elseif iprob == 6 % Final Stokes-Darcy problem used in paper
                
                alp = 0.5;
                w = 1.05;
                xint = xl(1,1:nel)*shl;
                yint = xl(2,1:nel)*shl;
                G = sqrt(1)/alp;                 
                fbx =  - (exp(1/(2*G))*sin(w + xint/G)*(G/1 - (2*mu)/G))/G - (2*mu*w*yint*sin(w*xint))/1;
                fby = (2*mu*cos(w*xint))/1 + 1;
                
            end

            %Evaluate residual of equilibrium equation
            rx = fbx + -px + mu*(2*ux_xx + ux_yy + uy_xy);
            ry = fby + -py + mu*(uy_xx + ux_xy + 2*uy_yy);
                
            %Evaluate explicit fine scale
            ufinex = (t11*rx+t12*ry)*b;
            ufiney = (t21*rx+t22*ry)*b;
            ufinex_x = (t11*rx+t12*ry)*be(1);
            ufiney_x = (t21*rx+t22*ry)*be(1);
            ufinex_y = (t11*rx+t12*ry)*be(2);
            ufiney_y = (t21*rx+t22*ry)*be(2);

            %Compute value of exact fields at int. point
            if iprob == 4
                [ue,duex,duey] = uexact_vasil(xint,yint,mu,1,0.5,w);
            elseif iprob == 5
                [ue,duex,duey] = uexact_vasil1(xint,yint,mu,1,0.5,w);
            elseif iprob == 6
                [ue,duex,duey] = uexact_vasiltim(xint,yint,mu,1,0.5,w);
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

            %Add explicit int. point error to element explicit error
            el2fine(1) = el2fine(1) + c1*ufinex^2;
            el2fine(2) = el2fine(2) + c1*ufiney^2;
            el2fine(3) = el2fine(3) + c1*ufinex_x^2;
            el2fine(4) = el2fine(4) + c1*ufiney_x^2;
            el2fine(5) = el2fine(5) + c1*ufinex_y^2;
            el2fine(6) = el2fine(6) + c1*ufiney_y^2;

        end %je

        for in= 1:3
            ElemE(in) = el2el(in);
            ElemE(in+3) = eprixel(in);
            ElemE(in+6) = epriyel(in);
            ElemE(in+9) = el2fine(in);
            ElemE(in+12) = el2fine(in+3);
        end
        
        H1up = el2fine(3)+el2fine(4)+el2fine(5)+el2fine(6);
        H1u = eprixel(1)+eprixel(2)+epriyel(1)+epriyel(2);
        Ieffvals(elem,:) = [sqrt(H1up/H1u) H1up H1u];
        
end %Task Switch
