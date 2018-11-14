
% Tim Truster
% 09/05/2012
% Stabilized DG for elasticity problem
% Revised implementation according to VMS derivations, using MVT to pull
% terms outside integrals
% Part of Inter_FEA_Program.
% Temperature added

PatchE = mateprop(1);
Patchv = mateprop(2);
thick = mateprop(3);
if ~exist('PSPS','var')
PSPS = 'n'; %plane stresS/straiN flag
end
if ~exist('iprob','var')
iprob = 0;
end

Bcol1 = [1; 3];
Bcol2 = [2; 3];
col1 = [1; 2];
col2 = [2; 1];

% Temperature parameters
alph = 1;
if iprob == 3
    kT = 1;
    rhoT = 1;
    cT = 1;
    q = 1;
    qtau = 1/2;
    h = 1/2;
end

switch isw %Task Switch
    
    case 1
        
        if ndf > 2
            
            for i = 3:ndf
                lie(i,1) = 0;
            end
            
        end
        
        istv = 7; % number of stresses per node for post-processing
        iste = 7; % number of stresses per element
%%
    case 3

        ElemK = zeros(ndf*nel);
        ElemF = zeros(ndf*nel,1);
        Nmat = zeros(2,ndf*nel);
        Bmat = zeros(3,ndf*nel);

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
        
        cvec = zeros(3,1); %temperature vector, see Hughes book
        if PSPS == 'n'
        Dmat = mu*diag([2 2 1]) + lam*[1; 1; 0]*[1 1 0];
        else
        Dmat = PatchE/(1-Patchv^2)*[1      Patchv  0
                                  Patchv  1      0
                                  0      0      (1-Patchv)/2];
        end

        der = 0;
        bf = 0;
        ib = 0;
        
        ulres = reshape(ul,ndf*nen,1);
        
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
            
            % Form B matrix
            for ie = 1:nel
                
              Nmat(1,(ie-1)*ndf+1) = shl(ie);
              Nmat(2,(ie-1)*ndf+2) = shl(ie);
                
              Bmat(Bcol1,(ie-1)*ndf+1) = shg(ie,col1);
              Bmat(Bcol2,(ie-1)*ndf+2) = shg(ie,col2);
                 
            end
            
            xint = xl(:,1:nel)*shl;
            
            if iprob == 3
                T1 = h*q/12/kT*(6*xint(2)/h - 1) + 1/(2*h*rhoT*cT)*qtau;
                T2 = q/(4*kT*h)*xint(2)^2;
                T = T1 + T2;
                cvec = alph*T*[1; 1; 0];
            elseif iprob == 4
                cvec = alph*1*[1; 1; 0];
            else
                cvec = zeros(3,1);
            end
            
            sigma = Dmat*(Bmat*ulres(1:nel*ndf) - cvec);
            
            ElemF = ElemF - c1*Bmat'*sigma;
            ElemK = ElemK + c1*Bmat'*Dmat*Bmat;
            
        end %je
ElemK
        
%%
    case 15 % body force for linear element
        
%         ElemF = zeros(nst,1); In FormFE
        
        if iprob == 5
            
        Nmat = zeros(2,nst);
        
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
        
%         thick = 1;

        der = 0;
        bf = 0;
        
        if iprob == 5
            fbx = 0;
            grav = 9.81;
            fby = -rho*grav;
            fb = [fbx; fby];
        else
            fb = bodyf;
        end
        
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
%             % Evaluate 1-D basis functions at integration points
%             [shp,shp2,be] = shpl_2d(r,s,nel,ideriv,ep,enrich);
%             %Evaluate first derivatives of basis functions at int. point
%             [Qxy, cartd2, Jdet] = shpg_2d(shp,shp2,xl,nel2,ideriv,be);
            c1 = Wgt*Jdet*thick;
            
            % Form B matrix
            for ie = 1:nel
                
              Nmat(1,(ie-1)*ndf+1) = shl(ie);
              Nmat(2,(ie-1)*ndf+2) = shl(ie);
                                    
            end
            
            ElemF = ElemF + c1*(Nmat'*fb);
                
        end %je
        
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
            if exist('iprob','var') && iprob == 1
%                 [ue,duex,duey] = uexact_thickcylE(xint,yint,0,ElemYM,Elemv,5,15,.32*lamda,0);
                [ue,duex,duey] = uexact_thickcylEP(xint,yint,0,PatchE,Patchv,PatchE*100,Ri,Ro,Ro,pPPC*lamda,0);
                duey(1) = 0;
                dux(2) = 1/2*(dux(2) + duy(1)); % strain component
                duy(1) = 0;
            elseif exist('iprob','var') && iprob == 2
                [ue,duex,duey] = uexactbb(xint,yint,PatchE,Patchv,Dbeam);
            elseif exist('iprob','var') && iprob == 3
                [ue,duex,duey] = uexacttb(xint,yint,PatchE,Patchv,alph,kT,rhoT,cT,q,qtau,h);
            elseif exist('iprob','var') && iprob == 5
                [ue,duex,duey] = uexact_selfw(xint,yint,PatchE,Patchv,rho);
            elseif exist('iprob','var') && iprob == 6
                [ue,duex,duey] = uexact_beam(xint,yint,PatchE,Patchv,Load);
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
    case -1 % boundary tractions
        
        ElemF = zeros(nst,1);
        
        % Reorder nodes for corresponding face of integration
        SurfOrientEtoS

        % Determine bounds of integration
        
        if nel == 4 || nel == 9
            
        dr = 2;
        ro = -1;
        
        % Upper Limit
%         if nodeA == ElemFlag(2)
            eR2 = 1;
%         else %
%             eR2 = 0;
%         end
        % Lower Limit
%         if nodeB == ElemFlag(1)
            eR1 = -1;
%         else %nodeA == ElemFlag(5)
%             eR1 = 0;
%         end
        
        elseif nel == 3 || nel == 6
            
        dr = 1;
        ro = 0;
            
        % Upper Limit
%         if nodeA == ElemFlag(2)
            eR2 = 1;
%         else %nodeA == ElemFlagR(5)
%             eR2 = 1/2;
%         end
        % Lower Limit
%         if nodeB == ElemFlag(1)
            eR1 = 0;
%         else %nodeA == ElemFlag(5)
%             eR1 = 1/2;
%         end
        
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
            Traction = zeros(1,2);
            
            if iprob == 3
                
                z = xl(2,1:nel)*shl;
                T1 = h*q/12/kT*(6*z/h - 1) + 1/(2*h*rhoT*cT)*qtau;
                T2 = q/(4*kT*h)*z^2;
                T = T1 + T2;
                NT = alph*PatchE/(rhoT*cT)*qtau;
                MT = alph*PatchE*q/(3*kT)*h^3;
                A = 2*h;
                I = 2*h^3/3;
                s_xx = (-alph*PatchE*T + NT/A + MT/I*z);
                s_xy = 0;
                s_yy = 0;
                Traction(1) = s_xx*tu3(1) + s_xy*tu3(2);
                Traction(2) = s_xy*tu3(1) + s_yy*tu3(2);
            elseif iprob == 4
                Traction = zeros(2,1);
            elseif iprob == 2
                x = 0;
                y = 0;

                for j = 1:nel
                   x = x + xl(1,j)*shl(j);
                   y = y + xl(2,j)*shl(j);
                end

                s_xx = bMoment*y/bInertia;
                s_xy = 0;
                s_yy = 0;
                Traction(1) = s_xx*tu3(1) + s_xy*tu3(2);
                Traction(2) = s_xy*tu3(1) + s_yy*tu3(2);
            elseif iprob == 5
                
                x = xl(1,1:nel)*shl;
                y = xl(2,1:nel)*shl;
                grav = 9.81d0;
                Cwid = D;
                Len = L/2;
                x2 = x-Len;
                s_xx = rho*grav*y - 1.d0/10.d0*rho*grav*y*(15.d0*Len^2-15.d0*x2^2 +4.d0*Cwid^2 +10.d0*y^2)/Cwid^2;
                s_xy = 3.d0/2.d0*rho*grav*x2*(Cwid^2-y^2)/Cwid^2;
                s_yy = rho*grav*y - 1.d0/2.d0*rho*grav*y*(3.d0*Cwid^2-y^2)/Cwid^2;
                Traction(1) = s_xx*tu3(1) + s_xy*tu3(2);
                Traction(2) = s_xy*tu3(1) + s_yy*tu3(2);
            elseif iprob == 6
                x = 0;
                y = 0;

                for j = 1:nel
                   x = x + xl(1,j)*shl(j);
                   y = y + xl(2,j)*shl(j);
                end
                Cwid = D;
                Len = L;
%                 Load = 2560;
                s_xx = 3.d0*(x-Len)*y*Load/(2*Cwid^3);
                s_xy = 3.d0/4.d0*(Cwid^2-y^2)*Load/Cwid^3;
                s_yy = 0;
                Traction(1) = s_xx*tu3(1) + s_xy*tu3(2);
                Traction(2) = s_xy*tu3(1) + s_yy*tu3(2);
            elseif edge < 0 % negative edge denotes pressure
                pressur = -traction(1); % pressure is defined positive in compression
                Traction = pressur*tu3(1:2);
            else
                Traction = traction(1:2);
            end
            
            c1 = Wgt*tm3*drdr*thick;
            
            for ie = 1:nel
                
              Nmat(1,(ie-1)*ndf+1) = shl(ie);
              Nmat(2,(ie-1)*ndf+2) = shl(ie);
                 
            end
            
            ElemF = ElemF + c1*Nmat'*Traction';

        end %ie
        
        % Reorder nodes back to the orientation of the element
        SurfOrientStoE
        ElemF(1:ndf*nel) = ElemF(ilist2);


%%    
    case 8 %interface stiffness
        
        ElemKRR = zeros(nstR,nstR);
        ElemKRL = zeros(nstR,nstL);
        ElemKLR = zeros(nstL,nstR);
        ElemKLL = zeros(nstL,nstL);
        ElemFR = zeros(nstR,1);
        ElemFL = zeros(nstL,1);
        t1 = zeros(3,1);
        t2 = t1;
        t3 = t1;
        thick = 1;
        
        bf1 = 0;
        bf2 = 0;
        
        ElemEL = matepropL(1);
        ElemvL = matepropL(2);
        ElemER = matepropR(1);
        ElemvR = matepropR(2);
        lamdaR = ElemvR*ElemER/((1+ElemvR)*(1-2*ElemvR));
        muR = ElemER/(2*(1+ElemvR));
        lamdaL = ElemvL*ElemEL/((1+ElemvL)*(1-2*ElemvL));
        muL = ElemEL/(2*(1+ElemvL));
        if PSPS == 'n'
        DmatR = muR*diag([2 2 1]) + lamdaR*[1; 1; 0]*[1 1 0];
        DmatL = muL*diag([2 2 1]) + lamdaL*[1; 1; 0]*[1 1 0];
        else
        DmatL = ElemEL/(1-ElemvL^2)*[1      ElemvL  0
                                  ElemvL  1      0
                                  0      0      (1-ElemvL)/2];
        DmatR = ElemER/(1-ElemvR^2)*[1      ElemvR  0
                                  ElemvR  1      0
                                  0      0      (1-ElemvR)/2];
        end
        
        NmatL = zeros(2,nstL);
        BmatL = zeros(3,nstL);
        bnAdN1 = zeros(3,nstL);
        N1 = zeros(2,nstL);
        NmatR = zeros(2,nstR);
        BmatR = zeros(3,nstR);
        bnAdN2 = zeros(3,nstR);
        N2 = zeros(2,nstR);
        
        ulresL = reshape(ulL,ndf*nelL,1);
        ulresR = reshape(ulR,ndf*nelR,1);
        
        % Determin bounds of integration segment
        InterConn2D2 % InterConn2DT % 
        
%         h = 2/(hR + hL);

        % Set jacobian for integration space
        drdr = (eL2 - eL1)/drL;
        
        m = (eR2-eR1)/(eL1-eL2);
        
        % Perform integration of various matrices
        
%         [tauL,intbL] = TauE1_2dT(xlintL,DmatL,lintt6,nelL);
%         [tauR,intbR] = TauE1_2dT(xlintR,DmatR,lintt6,nelR);
% For separate bubble types on T and Q
        [tauL,intbL] = TauE1_2d(xlintL,DmatL,nelL,lintt6,lintq9);
        [tauR,intbR] = TauE1_2d(xlintR,DmatR,nelR,lintt6,lintq9);
%         % lumping
%         tauL(1,1) = sum(tauL(:,1));
%         tauL(2,2) = sum(tauL(:,2));
%         tauL(1,2) = 0;
%         tauL(2,1) = 0;
%         tauR(1,1) = sum(tauR(:,1));
%         tauR(2,2) = sum(tauR(:,2));
%         tauR(1,2) = 0;
%         tauR(2,1) = 0;
        tau = tauL + tauR;
        
        lint = 3;
        ib = 0;
        der = 0;
        bf = 0;
        ebL = 0;
        ebR = 0;
        intedge = 0;
        
        % THIS LOOP COULD BE EVALUATED IN CLOSED-FORM
        for ie = 1:lint

% % For triangle bubbles on T and Q
%             if nelL == 3 || nelL == 6
%                 [Wgt,r,s] = intpntt(ie,lint,1);
%                 rT = r;
%                 sT = 0;
%             elseif nelL == 4 || nelL == 9
%                 [Wgt,r,s] = intpntq(ie,lint,1);
%                 rT = (r + 1)/2;
%                 sT = 0;
%             end
%                     
%             ebeL = edgebubble(rT,sT);
%             ebeR = ebeL;
            
% For separate bubble types on T and Q
            if nelL == 3 || nelL == 6
                [Wgt,r,s] = intpntt(ie,lint,1);
                ebeL = edgebubble(r,s,nelL);
            elseif nelL == 4 || nelL == 9
                [Wgt,r,s] = intpntq(ie,lint,1);
                ebeL = edgebubbleQ(r,s,nelL);
            end
                    
            if nelL == 3 || nelL == 6
                [Wgt,rR,sR] = intpntt(ie,lint,1);
                ebeR = edgebubble(rR,sR,nelR);
            elseif nelL == 4 || nelL == 9
                [Wgt,rR,sR] = intpntq(ie,lint,1);
                ebeR = edgebubbleQ(rR,sR,nelR);
            end
            
            r = drdr*(r-roL)+eL1;
            
            if nelL == 3 || nelL == 6
                [shlL,shld,shls,be] = shlt(r,s,nelL,nelL,der,bf);
                [QxyL, shgs, Jdet, be, xsL] = shgt(xlL(:,1:nelL),nelL,shld,shls,nen,bf,der,be);
            elseif nelL == 4 || nelL == 9
                [shlL,shld,shls,be] = shlq(r,s,nelL,nelL,der,bf);
                [QxyL, shgs, Jdet, be, xsL] = shgq(xlL(:,1:nelL),nelL,shld,shls,nen,bf,der,be);
            end
                    
            %Evaluate tangent and normal vectors
            t1 = [xsL(:,1); 0];
            [tm1, tu1] = VecNormalize(t1);
            t2 = [0; 0; 1];
            tm2 = 1;
            tu2 = t2';
            t3 = VecCrossProd(t1,t2);
            [tm3, tu3] = VecNormalize(t3);
            nLx = tu3(1);
            nLy = tu3(2);
            nvect = [nLx nLy];
                    
            c1 = Wgt*tm3*drdr*thick;
            
            ebL = ebL + c1*ebeL;
            ebR = ebR + c1*ebeR;
            intedge = intedge + c1;
            
        end
        
        if nitvms == 1
        % VMS
%         volL = getvol(xlL,nelL);
%         volR = getvol(xlR,nelR);
%         volbL = getvol(xlintL,nelL);
%         volbR = getvol(xlintR,nelR);
%         tauL = tauL*(volL/volbL)^-1;
%         tauR = tauR*(volR/volbR)^-1;
        edgeK = tauL*ebL^2 + tauR*ebR^2;
        gamL = ebL^2*(edgeK\tauL);
        gamR = ebR^2*(edgeK\tauR);
        ep = pencoeff*intedge*inv(edgeK); %#ok<MINV>
        elseif nitvms == 2
        % Nitsche
        volL = getvol(xlL,nelL);
        volR = getvol(xlR,nelR);
        h = 2/(intedge/volR + intedge/volL);
        gamL = 1/2*eye(2);
        gamR = 1/2*eye(2);
        ep = pencoeff*eye(2)*max(muL,muR)/h;
        else
        % RFB
        tauL = RFBEfunc(xlintL,muL,lamdaL,8);
        tauR = RFBEfunc(xlintR,muR,lamdaR,8);
%         if inter == 65
%             tauL = 1.0e-008 *[
%    0.315216684586325   0.060542382391646
%    0.060542382391646   0.389745645413441];
%             tauR = 1.0e-009 *[
%    0.498227303287427   0.043693780768305
%    0.043693780768297   0.972566080948220];
%         end
%         tauL = RFBtaus(2*D/m1,1/((n1/m1)/(L/D)));
%         tauR = RFBtaus(2*D/m2,1/((n2/m2)/(L/D)));
        edgeK = tauL + tauR;
        gamL = (edgeK\tauL);
        gamR = (edgeK\tauR);
        ep = pencoeff*intedge*inv(edgeK); %#ok<MINV>
        end
        
        for ie = 1:lint
            
            if nelL == 3 || nelL == 6
                [Wgt,r,s] = intpntt(ie,lint,1);
            elseif nelL == 4 || nelL == 9
                [Wgt,r,s] = intpntq(ie,lint,1);
            end
            
            r = drdr*(r-roL)+eL1;
            
            if nelL == 3 || nelL == 6
                [shlL,shld,shls,be] = shlt(r,s,nelL,nelL,der,bf);
                [QxyL, shgs, Jdet, be, xsL] = shgt(xlL(:,1:nelL),nelL,shld,shls,nen,bf,der,be);
            elseif nelL == 4 || nelL == 9
                [shlL,shld,shls,be] = shlq(r,s,nelL,nelL,der,bf);
                [QxyL, shgs, Jdet, be, xsL] = shgq(xlL(:,1:nelL),nelL,shld,shls,nen,bf,der,be);
            end
            
            rR = m*(r-eL2) + eR1;
            
            if nelR == 3 || nelR == 6
                s = 0;
            else %if nelR == 4
                s = -1;
            end
            
            if nelR == 3 || nelR == 6
                [shlR,shld,shls,be] = shlt(rR,s,nelR,nelR,der,bf);
                [QxyR, shgs, Jdet, be, xsR] = shgt(xlR(:,1:nelR),nelR,shld,shls,nen,bf,der,be);
            elseif nelR == 4 || nelR == 9
                [shlR,shld,shls,be] = shlq(rR,s,nelR,nelR,der,bf);
                [QxyR, shgs, Jdet, be, xsR] = shgq(xlR(:,1:nelR),nelR,shld,shls,nen,bf,der,be);
            end
            
            %Evaluate tangent and normal vectors
            t1 = [xsL(:,1); 0];
            [tm1, tu1] = VecNormalize(t1);
            t2 = [0; 0; 1];
            tm2 = 1;
            tu2 = t2';
            t3 = VecCrossProd(t1,t2);
            [tm3, tu3] = VecNormalize(t3);
            nLx = tu3(1);
            nLy = tu3(2);
            nvect = [nLx 0 nLy
                     0 nLy nLx];
            nvec = [nLx; nLy];
                    
            c1 = Wgt*tm3*drdr*thick;
            
            for i = 1:nelL
                NmatL(1,(i-1)*2+1) = shlL(i);
                NmatL(2,(i-1)*2+2) = shlL(i);
                BmatL(Bcol1,(i-1)*2+1) = QxyL(i,col1);
                BmatL(Bcol2,(i-1)*2+2) = QxyL(i,col2);
            end
            
            for i = 1:nelR
                NmatR(1,(i-1)*2+1) = shlR(i);
                NmatR(2,(i-1)*2+2) = shlR(i);
                BmatR(Bcol1,(i-1)*2+1) = QxyR(i,col1);
                BmatR(Bcol2,(i-1)*2+2) = QxyR(i,col2);
            end
            
            bnAdN1 = gamL*nvect*DmatL*BmatL;
            bnAdN2 = gamR*nvect*DmatR*BmatR;
        
            xint = xlL*shlL;
            
            if iprob == 3
                T1 = h*q/12/kT*(6*xint(2)/h - 1) + 1/(2*h*rhoT*cT)*qtau;
                T2 = q/(4*kT*h)*xint(2)^2;
                T = T1 + T2;
                cvec = alph*T*[1; 1; 0];
            elseif iprob == 4
                cvec = alph*1*[1; 1; 0];
            end
            
            tvtr = (bnAdN1*ulresL + bnAdN2*ulresR) - (gamL*nvect*DmatL + gamR*nvect*DmatR)*cvec; %average stress
            jumpu = NmatR*ulresR - NmatL*ulresL;        %displacement jump

            % Not needed for linear problems
            ElemFL = ElemFL - c1*( - NmatL'*(tvtr + ep*jumpu) + bnAdN1'*jumpu);
            ElemFR = ElemFR - c1*( + NmatR'*(tvtr + ep*jumpu) + bnAdN2'*jumpu);

            ElemKLL = ElemKLL - c1*NmatL'*bnAdN1 - c1*bnAdN1'*NmatL;
            ElemKLR = ElemKLR - c1*NmatL'*bnAdN2 + c1*bnAdN1'*NmatR;
            ElemKRL = ElemKRL + c1*NmatR'*bnAdN1 - c1*bnAdN2'*NmatL;
            ElemKRR = ElemKRR + c1*NmatR'*bnAdN2 + c1*bnAdN2'*NmatR;

            ElemKLL = ElemKLL + c1*(NmatL'*ep*NmatL);
            ElemKLR = ElemKLR - c1*(NmatL'*ep*NmatR);
            ElemKRL = ElemKRL - c1*(NmatR'*ep*NmatL);
            ElemKRR = ElemKRR + c1*(NmatR'*ep*NmatR);
            
        end
% ElemKLL

%%
    case 22 %Stress Projection
        
        ElemM = zeros(nst);
        ElemF = zeros(nst,1);
        
        % Load Gauss Points for quadrature
        if enrich == 1
            if nel == 3
                lint = lintt3;%13;
            elseif nel == 4
                lint = lintq4;
            elseif nel == 6
                lint = lintt6;%13;
            elseif nel == 9
                lint = lintq9;
            end
%             [rlist, rWgts, rnum] = GaussPoints(pr+2);
%             [slist, sWgts, snum] = GaussPoints(pr+1);
        else
            if nel == 3
                lint = lintt3;%13;
            elseif nel == 4
                lint = lintq4;
            elseif nel == 6
                lint = lintt6;%13;
            elseif nel == 9
                lint = lintq9;
            end
%             [rlist, rWgts, rnum] = GaussPoints(pr+1);
%             [slist, sWgts, snum] = GaussPoints(pr+1);
        end

%         if PSPS == 's'
%             PatchE = PatchE*(1 + 2*Patchv)/(1 + Patchv)^2;
%             Patchv = Patchv/(1 + Patchv);
%         end
        lam = Patchv*PatchE/((1+Patchv)*(1-2*Patchv));
        mu = PatchE/(2*(1+Patchv));
        vis = mu;
%         thick = 1;
        ideriv = 0;
        fbx = 0;
        fby = 0;
        hx = 0;
        hy = 0;
        II = eye(2,2);
        inds = [1 2 1
                1 2 2];

        for je = 1:lint

            if nel == 3 || nel == 6
                [Wgt,litr,lits] = intpntt(je,lint,ib);
                [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
                [shg, shgs, Jdet, be, xs] = shgt(xl,nel,shld,shls,nen,bf,der,be);
                shp = shlt(litr,lits,nelP,nelP,0,0);
                shpS = shlt(litr,lits,nelS,nelS,0,0);
            elseif nel == 4 || nel == 9
                [Wgt,litr,lits] = intpntq(je,lint,ib);
                [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
                [shg, shgs, Jdet, be, xs] = shgq(xl,nel,shld,shls,nen,bf,der,be);
                shp = shlq(litr,lits,nelP,nelP,0,0);
                shpS = shlq(litr,lits,nelS,nelS,0,0);
            end
            
            xint = xl(1,:)*shl;
            yint = xl(2,:)*shl;
            dux = ul(1:2,:)*shg(:,1);
            duy = ul(1:2,:)*shg(:,2);
            u = ul(1:2,:)*shl;
            strains = [2*dux(1) 2*duy(2) dux(2)+duy(1)];
            px = ul(3,1:nelP)*shp(:,1)*pfact;
            py = ul(3,1:nelP)*shp(:,2)*pfact;
            p = ul(3,1:nelP)*shl*pfact;
            
            sigmas = p*II(inds(1,stres),inds(2,stres)) + mu*strains(stres);

            for jn = 1:nelS

                djn=shpS(3,jn);

                ElemF(jn) = ElemF(jn) + c1*djn*sigmas;

                for in=1:nelS

                    din=shpS(3,in);

                    ElemM(in,jn)=ElemM(in,jn) + c1*din*djn;

                end %in

            end %jn

        end %je
%%
    case 23
        
                ElemKRR = zeros(nstR,nstR);
        ElemKRL = zeros(nstR,nstL);
        ElemKLR = zeros(nstL,nstR);
        ElemKLL = zeros(nstL,nstL);
        t1 = zeros(3,1);
        t2 = t1;
        t3 = t1;
        ElemEL = matepropL(1);
        ElemvL = matepropL(2);
        ElemER = matepropR(1);
        ElemvR = matepropR(2);
        lamdaR = ElemvR*ElemER/((1+ElemvR)*(1-2*ElemvR));
        muR = ElemER/(2*(1+ElemvR));
        lamdaL = ElemvL*ElemEL/((1+ElemvL)*(1-2*ElemvL));
        muL = ElemEL/(2*(1+ElemvL));
        
        if neqc > 0 %Lagrange Multipliers
            ElemKcL = zeros(nstC,nstL);
            ElemKcR = zeros(nstC,nstR);
            ElemKLc = zeros(nstL,nstC);
            ElemKRc = zeros(nstR,nstC);
        end

        % Determine bounds of integration, right
        
        if nelR == 4 || nelR == 9
            
            drR = 2;
            roR = -1;

            % Upper Limit
            if nodeAR == ElemFlagR(2)
                eR2 = 1;
            elseif nodeAR == ElemFlagR(nel2R)
                eR2 = epR;
            elseif nodeAR == -1 %no enrichment but DG instead
                xy = xlL(:,1);
                POUxi = POU_Coord(xy(1),xy(2),xlR,1,nelR);
                eR2 = POUxi(1);
            elseif nelR == 9 && nodeAR == ElemFlagR(5)
                eR2 = 0;
            end
            % Lower Limit
            if nodeBR == ElemFlagR(1)
                eR1 = -1;
            elseif nodeBR == ElemFlagR(nel2R)
                eR1 = epR;
            elseif nodeBR == -1 %no enrichment but DG instead
                xy = xlL(:,2);
                POUxi = POU_Coord(xy(1),xy(2),xlR,1,nelR);
                eR1 = POUxi(1);
            elseif nelR == 9 && nodeBR == ElemFlagR(5)
                eR1 = 0;
            end
        
        elseif nelR == 3 || nelR == 6
            
            drR = 1;
            roR = 0;

            % Upper Limit
            if nodeAR == ElemFlagR(2)
                eR2 = 1;
            elseif nodeAR == ElemFlagR(nel2R)
                eR2 = epR;
            elseif nodeAR == -1 %no enrichment but DG instead
                xy = xlL(:,1);
                POUxi = POU_Coord(xy(1),xy(2),xlR,1,nelR);
                eR2 = POUxi(1);
            elseif nelR == 6 && nodeAR == ElemFlagR(4)
                eR2 = 1/2;
            end
            % Lower Limit
            if nodeBR == ElemFlagR(1)
                eR1 = 0;
            elseif nodeBR == ElemFlagR(nel2R)
                eR1 = epR;
            elseif nodeBR == -1 %no enrichment but DG instead
                xy = xlL(:,2);
                POUxi = POU_Coord(xy(1),xy(2),xlR,1,nelR);
                eR1 = POUxi(1);
            elseif nelR == 6 && nodeBR == ElemFlagR(4)
                eR1 = 1/2;
            end
        
        end
        
        % Determine bounds of integration, left
        
        if nelL == 4 || nelL == 9
            
            drL = 2;
            roL = -1;

            % Upper Limit
            if nodeAL == ElemFlagL(1)
                eL1 = -1;
            elseif nodeAL == ElemFlagL(nel2L)
                eL1 = epL;
            elseif nodeAL == -1 %no enrichment but DG instead
                xy = xlR(:,2);
                POUxi = POU_Coord(xy(1),xy(2),xlL,1,nelL);
                eL1 = POUxi(1);
            elseif nelL == 9 && nodeAL == ElemFlagL(5)
                eL1 = 0;
            end
            % Lower Limit
            if nodeBL == ElemFlagL(2)
                eL2 = 1;
            elseif nodeBL == ElemFlagL(nel2L)
                eL2 = epL;
            elseif nodeBL == -1 %no enrichment but DG instead
                xy = xlR(:,1);
                POUxi = POU_Coord(xy(1),xy(2),xlL,1,nelL);
                eL2 = POUxi(1);
            elseif nelL == 9 && nodeB == ElemFlagL(5)
                eL2 = 0;
            end
        
        elseif nelL == 3 || nelL == 6
            
            drL = 1;
            roL = 0;

            % Upper Limit
            if nodeAL == ElemFlagL(1)
                eL1 = 0;
            elseif nodeAL == ElemFlagL(nel2L)
                eL1 = epL;
            elseif nodeAL == -1 %no enrichment but DG instead
                xy = xlR(:,2);
                POUxi = POU_Coord(xy(1),xy(2),xlL,1,nelL);
                eL1 = POUxi(1);
            elseif nelL == 6 && nodeAL == ElemFlagL(4)
                eL1 = 1/2;
            end
            % Lower Limit
            if nodeBL == ElemFlagL(2)
                eL2 = 1;
            elseif nodeBL == ElemFlagL(nel2L)
                eL2 = epL;
            elseif nodeBL == -1 %no enrichment but DG instead
                xy = xlR(:,1);
                POUxi = POU_Coord(xy(1),xy(2),xlL,1,nelL);
                eL2 = POUxi(1);
            elseif nelL == 6 && nodeBL == ElemFlagL(4)
                eL2 = 1/2;
            end
        
        end
        
        % Set jacobian for integration space
        drdr = (eL2 - eL1)/drL;
        
        m = (eR2-eR1)/(eL1-eL2);
        
        lint = 4;
        ideriv = 0;

        s = 0;
        
        for ie = 1:lint

            if nelL == 3 || nelL == 6
                [Wgt,r,s] = intpntt(ie,lint,1);
            elseif nelL == 4 || nel == 9
                [Wgt,r,s] = intpntq(ie,lint,1);
            end
            
            r = drdr*(r-roL)+eL1;
            
            if nelL == 3 || nelL == 6
                [shlL,shld,shls,be] = shlt(r,s,nelL,nelL,der,bf);
                [QxyL, shgs, Jdet, be, xsL] = shgt(xlL,nelL,shld,shls,nen,bf,der,be);
                shpL = shlt(r,s,nelLP,nelLP,der,bf);
                shpLS = shlt(r,s,nelLS,nelLS,der,bf);
            elseif nelL == 4 || nelL == 9
                [shlL,shld,shls,be] = shlq(r,s,nelL,nelL,der,bf);
                [QxyL, shgs, Jdet, be, xsL] = shgq(xlL,nelL,shld,shls,nen,bf,der,be);
                shpL = shlq(r,s,nelLP,nelLP,der,bf);
                shpLS = shlq(r,s,nelLS,nelLS,der,bf);
            end
            
            rR = m*(r-eL2) + eR1;
            
            if nelR == 3 || nelR == 6
                s = 0;
            else %if nelR == 4
                s = -1;
            end
            
            if nelR == 3 || nelR == 6
                [shlR,shld,shls,be] = shlt(rR,s,nelR,nelR,der,bf);
                [QxyR, shgs, Jdet, be, xsR] = shgt(xlR,nelR,shld,shls,nen,bf,der,be);
                shpR = shlt(rR,s,nelRP,nelRP,der,bf);
                shpRS = shlt(rR,s,nelRS,nelRS,der,bf);
            elseif nelR == 4 || nelR == 9
                [shlR,shld,shls,be] = shlq(rR,s,nelR,nelR,der,bf);
                [QxyR, shgs, Jdet, be, xsR] = shgq(xlR,nelR,shld,shls,nen,bf,der,be);
                shpR = shlq(rR,s,nelRP,nelRP,der,bf);
                shpRS = shlq(rR,s,nelRS,nelRS,der,bf);
            end
                    
            %Evaluate tangent and normal vectors
            t1 = [xsL(:,1); 0];
            [tm1, tu1] = VecNormalize(t1);
            t2 = [0; 0; 1];
            tm2 = 1;
            tu2 = t2';
            t3 = VecCrossProd(t1,t2);
            [tm3, tu3] = VecNormalize(t3);
            nLx = tu3(1);
            nLy = tu3(2);
            nRx = -tu3(1);
            nRy = -tu3(2);
            tLx = tu1(1);
            tLy = tu1(2);
            
            dux = ulL(1:2,:)*QxyL(1,:);
            duy = ulL(1:2,:)*QxyL(2,:);
            strainsL = [2*dux(1) 2*duy(2) dux(2)+duy(1)];
            pL = ulL(3,1:nelLP)*shpL*pfact;
            sigmasL = pL*II(inds(1,stres),inds(2,stres)) + mu*strainsL(stres);
            dux = ulR(1:2,:)*QxyR(1,:);
            duy = ulR(1:2,:)*QxyR(2,:);
            pR = ulR(3,1:nelRP)*shpR*pfact;
            strainsR = [2*dux(1) 2*duy(2) dux(2)+duy(1)];
            sigmasR = pR*II(inds(1,stres),inds(2,stres)) + mu*strainsR(stres);
            
            sigmas = sigmasR-sigmasL;
                    
            c1 = Wgt*tm3*drdr*thick;
                
            for a = 1:nelLS
            
                NLa = shpLS(a)*c1;

                for c = 1:nelLS

                    NLc = shpLS(c);

                    ElemKLL(c,a) = ElemKLL(c,a) - ...
                                             10*(NLc*NLa);
                                         
                end
                                         
                for d = 1:nelRS

                    NRd = shpRS(d);

                    ElemKRL(d,a) = ElemKRL(d,a) + ...
                                             10*(NRd*NLa);
                                         
                end

            end
            
            for b = 1:nelRS
            
                NRb = shpRS(b)*c1;

                for c = 1:nelLS

                    NLc = shpLS(c);

                    ElemKLR(c,b) = ElemKLR(c,b) + ...
                                             10*(NLc*NRb);
                                         
                end
                                         
                for d = 1:nelRS

                    NRd = shpRS(d);

                    ElemKRR(d,b) = ElemKRR(d,b) - ...
                                             10*(NRd*NRb);
                                         
                end

            end
            
        end %lint
        
        ElemKLL = 1/10*ElemKLL;
        ElemKRL = 1/10*ElemKRL;
        ElemKLR = 1/10*ElemKLR;
        ElemKRR = 1/10*ElemKRR;

%%        
    case 25 %Stress Projection2

        ElemS = zeros(nel,npstr+1);
        ElemS2 = zeros(nel,npstr);
        stresID = [1 2 3 0 1 3];

        lam = Patchv*PatchE/((1+Patchv)*(1-2*Patchv));
        mu = PatchE/(2*(1+Patchv));
        thick = 1;
        fbx = 0;
        fby = 0;
        fbz = 0;
        if PSPS == 'n'
        Dmat = mu*diag([2 2 1]) + lam*[1; 1; 0]*[1 1 0];
        else
        Dmat = PatchE/(1-Patchv^2)*[1      Patchv  0
                                  Patchv  1      0
                                  0      0      (1-Patchv)/2];
        end
        Nmat = zeros(2,ndf*nel);
        Bmat = zeros(3,ndf*nel);
        I1 = [1; 1; 0];
        
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
        
        ulres = reshape(ul,nst,1);

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
            
            % Form B matrix
            for ie = 1:nel
              
              Bmat(Bcol1,(ie-1)*ndf+1) = shg(ie,col1);
              Bmat(Bcol2,(ie-1)*ndf+2) = shg(ie,col2);
                 
            end
            
            epsil = Bmat*ulres(1:ndf*nel);
            sigma = Dmat*epsil;
            
            for stres = 1:npstr
            
            if stres <= 3 % stress components
                sigmas = sigma(stresID(stres));
            elseif stres >= 5
                if stres <= 6 % principal stresses
                    if PSPS == 'n'
                        sigz = lam*(epsil(1)+epsil(2));
                    else
                        sigz = 0;
                    end
                    sigma2 = [sigma(1) sigma(3) 0; sigma(3) sigma(2) 0; 0 0 sigz];
                    psig = eig(sigma2);
                    sigmas = psig(stresID(stres));
                else % hydrostatic stress
                    if PSPS == 'n'
                        sigz = lam*(epsil(1)+epsil(2));
                    else
                        sigz = 0;
                    end
                    sigmas = 1/3*(sigma'*I1 + sigz);
                end
            else % von Mises stress
                if PSPS == 'n'
                    sigz = lam*(epsil(1)+epsil(2));
                else
                    sigz = 0;
                end
                trs = sigma'*I1 + sigz;
                dsig = sigma - 1/3*trs*I1;
                sigmas = sqrt(3/2*(dsig'*dsig));
            end
            
            ElemS2(ll,stres) = sigmas;
            
            end

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
        
    case 26
        
        ElemS = zeros(npstr,1);
        stresID = [1 2 3 0 1 3];

        lam = Patchv*PatchE/((1+Patchv)*(1-2*Patchv));
        mu = PatchE/(2*(1+Patchv));
        thick = 1;
        fbx = 0;
        fby = 0;
        fbz = 0;
        if PSPS == 'n'
        Dmat = mu*diag([2 2 1]) + lam*[1; 1; 0]*[1 1 0];
        else
        Dmat = PatchE/(1-Patchv^2)*[1      Patchv  0
                                  Patchv  1      0
                                  0      0      (1-Patchv)/2];
        end
        Nmat = zeros(2,ndf*nel);
        Bmat = zeros(3,ndf*nel);
        I1 = [1; 1; 0];
        
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
        
        ulres = reshape(ul,nst,1);

        %Stress Loop
        ll = 1;
%         for ll = 1:nint

            %Evaluate first derivatives of basis functions at int. point
            if nel == 3 || nel == 6
              [w,litr,lits] =  intpntt(ll,ll,ib);
              [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
              [shg,shgs,Jdet,be] = shgt(xl,nel,shld,shls,nel,bf,der,be);
            else
              [w,litr,lits] =  intpntq(ll,ll,ib);
              [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
              [shg,shgs,Jdet,be] = shgq(xl,nel,shld,shls,nel,bf,der,be);
            end
            
            % Form B matrix
            for ie = 1:nel
              
              Bmat(Bcol1,(ie-1)*ndf+1) = shg(ie,col1);
              Bmat(Bcol2,(ie-1)*ndf+2) = shg(ie,col2);
                 
            end
            
            epsil = Bmat*ulres(1:ndf*nel);
            sigma = Dmat*epsil;
            
            for stres = 1:npstr
            
            if stres <= 3 % stress components
                sigmas = sigma(stresID(stres));
            elseif stres >= 5
                if stres <= 6 % principal stresses
                    if PSPS == 'n'
                        sigz = lam*(epsil(1)+epsil(2));
                    else
                        sigz = 0;
                    end
                    sigma2 = [sigma(1) sigma(3) 0; sigma(3) sigma(2) 0; 0 0 sigz];
                    psig = eig(sigma2);
                    sigmas = psig(stresID(stres));
                else % hydrostatic stress
                    if PSPS == 'n'
                        sigz = lam*(epsil(1)+epsil(2));
                    else
                        sigz = 0;
                    end
                    sigmas = 1/3*(sigma'*I1 + sigz);
                end
            else % von Mises stress
                if PSPS == 'n'
                    sigz = lam*(epsil(1)+epsil(2));
                else
                    sigz = 0;
                end
                trs = sigma'*I1 + sigz;
                dsig = sigma - 1/3*trs*I1;
                sigmas = sqrt(3/2*(dsig'*dsig));
            end
            
            ElemS(stres) = sigmas;
            
            end

%         end %je
        
    case 60
        
        numhr = 3;
        ElemI = zeros(14,numhr);
        t1 = zeros(3,1);
        t2 = t1;
        t3 = t1;
        thick = 1;
        
        bf1 = 0;
        bf2 = 0;
        
        ElemEL = matepropL(1);
        ElemvL = matepropL(2);
        ElemER = matepropR(1);
        ElemvR = matepropR(2);
        lamdaR = ElemvR*ElemER/((1+ElemvR)*(1-2*ElemvR));
        muR = ElemER/(2*(1+ElemvR));
        lamdaL = ElemvL*ElemEL/((1+ElemvL)*(1-2*ElemvL));
        muL = ElemEL/(2*(1+ElemvL));
        if PSPS == 'n'
        DmatR = muR*diag([2 2 1]) + lamdaR*[1; 1; 0]*[1 1 0];
        DmatL = muL*diag([2 2 1]) + lamdaL*[1; 1; 0]*[1 1 0];
        else
        DmatL = ElemEL/(1-ElemvL^2)*[1      ElemvL  0
                                  ElemvL  1      0
                                  0      0      (1-ElemvL)/2];
        DmatR = ElemER/(1-ElemvR^2)*[1      ElemvR  0
                                  ElemvR  1      0
                                  0      0      (1-ElemvR)/2];
        end
        
        NmatL = zeros(2,nstL);
        BmatL = zeros(3,nstL);
        bnAdN1 = zeros(3,nstL);
        N1 = zeros(2,nstL);
        NmatR = zeros(2,nstR);
        BmatR = zeros(3,nstR);
        bnAdN2 = zeros(3,nstR);
        N2 = zeros(2,nstR);
        
        ulresL = reshape(ulL,ndf*nelL,1);
        ulresR = reshape(ulR,ndf*nelR,1);
        
        % Determin bounds of integration segment
        InterConn2D2 % InterConn2DT % 
        
%         h = 2/(hR + hL);

        % Set jacobian for integration space
        drdr = (eL2 - eL1)/drL;
        
        m = (eR2-eR1)/(eL1-eL2);
        
        % Perform integration of various matrices
        
%         [tauL,intbL] = TauE1_2dT(xlintL,DmatL,lintt6,nelL);
%         [tauR,intbR] = TauE1_2dT(xlintR,DmatR,lintt6,nelR);
% For separate bubble types on T and Q
        [tauL,intbL] = TauE1_2d(xlintL,DmatL,nelL,lintt6,lintq9);
        [tauR,intbR] = TauE1_2d(xlintR,DmatR,nelR,lintt6,lintq9);
%         % lumping
%         tauL(1,1) = sum(tauL(:,1));
%         tauL(2,2) = sum(tauL(:,2));
%         tauL(1,2) = 0;
%         tauL(2,1) = 0;
%         tauR(1,1) = sum(tauR(:,1));
%         tauR(2,2) = sum(tauR(:,2));
%         tauR(1,2) = 0;
%         tauR(2,1) = 0;
        tau = tauL + tauR;
        
        lint = 3;
        ib = 0;
        der = 0;
        bf = 0;
        ebL = 0;
        ebR = 0;
        intedge = 0;
        
        % THIS LOOP COULD BE EVALUATED IN CLOSED-FORM
        for ie = 1:lint

% % For triangle bubbles on T and Q
%             if nelL == 3 || nelL == 6
%                 [Wgt,r,s] = intpntt(ie,lint,1);
%                 rT = r;
%                 sT = 0;
%             elseif nelL == 4 || nelL == 9
%                 [Wgt,r,s] = intpntq(ie,lint,1);
%                 rT = (r + 1)/2;
%                 sT = 0;
%             end
%                     
%             ebeL = edgebubble(rT,sT);
%             ebeR = ebeL;
            
% For separate bubble types on T and Q
            if nelL == 3 || nelL == 6
                [Wgt,r,s] = intpntt(ie,lint,1);
                ebeL = edgebubble(r,s,nelL);
            elseif nelL == 4 || nelL == 9
                [Wgt,r,s] = intpntq(ie,lint,1);
                ebeL = edgebubbleQ(r,s,nelL);
            end
                    
            if nelL == 3 || nelL == 6
                [Wgt,rR,sR] = intpntt(ie,lint,1);
                ebeR = edgebubble(rR,sR,nelR);
            elseif nelL == 4 || nelL == 9
                [Wgt,rR,sR] = intpntq(ie,lint,1);
                ebeR = edgebubbleQ(rR,sR,nelR);
            end
            
            r = drdr*(r-roL)+eL1;
            
            if nelL == 3 || nelL == 6
                [shlL,shld,shls,be] = shlt(r,s,nelL,nelL,der,bf);
                [QxyL, shgs, Jdet, be, xsL] = shgt(xlL(:,1:nelL),nelL,shld,shls,nen,bf,der,be);
            elseif nelL == 4 || nelL == 9
                [shlL,shld,shls,be] = shlq(r,s,nelL,nelL,der,bf);
                [QxyL, shgs, Jdet, be, xsL] = shgq(xlL(:,1:nelL),nelL,shld,shls,nen,bf,der,be);
            end
                    
            %Evaluate tangent and normal vectors
            t1 = [xsL(:,1); 0];
            [tm1, tu1] = VecNormalize(t1);
            t2 = [0; 0; 1];
            tm2 = 1;
            tu2 = t2';
            t3 = VecCrossProd(t1,t2);
            [tm3, tu3] = VecNormalize(t3);
            nLx = tu3(1);
            nLy = tu3(2);
            nvect = [nLx nLy];
                    
            c1 = Wgt*tm3*drdr*thick;
            
            ebL = ebL + c1*ebeL;
            ebR = ebR + c1*ebeR;
            intedge = intedge + c1;
            
        end
        
        if nitvms == 1
        % VMS
%         volL = getvol(xlL,nelL);
%         volR = getvol(xlR,nelR);
%         volbL = getvol(xlintL,nelL);
%         volbR = getvol(xlintR,nelR);
%         tauL = tauL*(volL/volbL)^-1;
%         tauR = tauR*(volR/volbR)^-1;
        edgeK = tauL*ebL^2 + tauR*ebR^2;
        gamL = ebL^2*(edgeK\tauL);
        gamR = ebR^2*(edgeK\tauR);
        ep = pencoeff*intedge*inv(edgeK); %#ok<MINV>
        elseif nitvms == 2
        % Nitsche
        volL = getvol(xlL,nelL);
        volR = getvol(xlR,nelR);
        h = 2/(intedge/volR + intedge/volL);
        gamL = 1/2*eye(2);
        gamR = 1/2*eye(2);
        ep = pencoeff*eye(2)*max(muL,muR)/h;
        else
        % RFB
        tauL = RFBEfunc(xlintL,muL,lamdaL,8);
        tauR = RFBEfunc(xlintR,muR,lamdaR,8);
%         if inter == 65
%             tauL = 1.0e-008 *[
%    0.315216684586325   0.060542382391646
%    0.060542382391646   0.389745645413441];
%             tauR = 1.0e-009 *[
%    0.498227303287427   0.043693780768305
%    0.043693780768297   0.972566080948220];
%         end
%         tauL = RFBtaus(2*D/m1,1/((n1/m1)/(L/D)));
%         tauR = RFBtaus(2*D/m2,1/((n2/m2)/(L/D)));
        edgeK = tauL + tauR;
        gamL = (edgeK\tauL);
        gamR = (edgeK\tauR);
        ep = pencoeff*intedge*inv(edgeK); %#ok<MINV>
        end
        
        for ie = 1:lint
            
            if nelL == 3 || nelL == 6
                [Wgt,r,s] = intpntt(ie,lint,1);
            elseif nelL == 4 || nelL == 9
                [Wgt,r,s] = intpntq(ie,lint,1);
            end
            
            r = drdr*(r-roL)+eL1;
            
            if nelL == 3 || nelL == 6
                [shlL,shld,shls,be] = shlt(r,s,nelL,nelL,der,bf);
                [QxyL, shgs, Jdet, be, xsL] = shgt(xlL(:,1:nelL),nelL,shld,shls,nen,bf,der,be);
            elseif nelL == 4 || nelL == 9
                [shlL,shld,shls,be] = shlq(r,s,nelL,nelL,der,bf);
                [QxyL, shgs, Jdet, be, xsL] = shgq(xlL(:,1:nelL),nelL,shld,shls,nen,bf,der,be);
            end
            
            rR = m*(r-eL2) + eR1;
            xint = xlL(1,1:nelL)*shlL;
            yint = xlL(2,1:nelL)*shlL;
            
            if nelR == 3 || nelR == 6
                s = 0;
            else %if nelR == 4
                s = -1;
            end
            
            if nelR == 3 || nelR == 6
                [shlR,shld,shls,be] = shlt(rR,s,nelR,nelR,der,bf);
                [QxyR, shgs, Jdet, be, xsR] = shgt(xlR(:,1:nelR),nelR,shld,shls,nen,bf,der,be);
            elseif nelR == 4 || nelR == 9
                [shlR,shld,shls,be] = shlq(rR,s,nelR,nelR,der,bf);
                [QxyR, shgs, Jdet, be, xsR] = shgq(xlR(:,1:nelR),nelR,shld,shls,nen,bf,der,be);
            end
            
            %Evaluate tangent and normal vectors
            t1 = [xsL(:,1); 0];
            [tm1, tu1] = VecNormalize(t1);
            t2 = [0; 0; 1];
            tm2 = 1;
            tu2 = t2';
            t3 = VecCrossProd(t1,t2);
            [tm3, tu3] = VecNormalize(t3);
            nLx = tu3(1);
            nLy = tu3(2);
            nvect = [nLx 0 nLy
                     0 nLy nLx];
            nvec = [nLx; nLy];
                    
            c1 = Wgt*tm3*drdr*thick;
            
            for i = 1:nelL
                NmatL(1,(i-1)*2+1) = shlL(i);
                NmatL(2,(i-1)*2+2) = shlL(i);
                BmatL(Bcol1,(i-1)*2+1) = QxyL(i,col1);
                BmatL(Bcol2,(i-1)*2+2) = QxyL(i,col2);
            end
            
            for i = 1:nelR
                NmatR(1,(i-1)*2+1) = shlR(i);
                NmatR(2,(i-1)*2+2) = shlR(i);
                BmatR(Bcol1,(i-1)*2+1) = QxyR(i,col1);
                BmatR(Bcol2,(i-1)*2+2) = QxyR(i,col2);
            end
            
            bnAdN1 = gamL*nvect*DmatL*BmatL;
            bnAdN2 = gamR*nvect*DmatR*BmatR;
        
            xint = xlL*shlL;
            
            if iprob == 3
                T1 = h*q/12/kT*(6*xint(2)/h - 1) + 1/(2*h*rhoT*cT)*qtau;
                T2 = q/(4*kT*h)*xint(2)^2;
                T = T1 + T2;
                cvec = alph*T*[1; 1; 0];
            elseif iprob == 4
                cvec = alph*1*[1; 1; 0];
            end
            
            tvtr = (bnAdN1*ulresL + bnAdN2*ulresR) - (gamL*nvect*DmatL + gamR*nvect*DmatR)*cvec; %average stress
            jumpu = NmatR*ulresR - NmatL*ulresL;        %displacement jump

        %Compute value of exact fields at int. point
        if iprob == 1
        elseif iprob == 2
            [ue,duex,duey] = uexactbb(xint(1),xint(2),ElemEL,ElemvL,D);
            strvec = [duex(1); duey(2); duex(2)+duey(1)];
        elseif iprob == 3
            [ue,duex,duey] = uexacttb(xint(1),xint(2),ElemEL,ElemvL,alph,kT,rhoT,cT,q,qtau,h);
            T1 = h*q/12/kT*(6*xint(2)/h - 1) + 1/(2*h*rhoT*cT)*qtau;
            T2 = q/(4*kT*h)*xint(2)^2;
            T = T1 + T2;
            cvec = alph*T*[1; 1; 0];
            strvec = [duex(1); duey(2); duex(2)+duey(1)] - cvec;
            sxx2 = alph*ElemEL/(1-ElemvL)*(q/4/h/kT)*(h^2/3-xint(2)^2);
        end
        t_exact = nvect*DmatL*strvec;

        ElemI(:,ie) = [xint
                       jumpu
                       tvtr
                       (tvtr + ep*jumpu)
                       t_exact
                       (tvtr + ep*jumpu)-t_exact
                       (tvtr)-t_exact];
        end %lint
        
    case 16 % J-Integral
        
        ElemJ = zeros(ndm,1);
        Bmat = zeros(3,ndf*nel);
        enmom = zeros(2,2);

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
        
        cvec = zeros(3,1); %temperature vector, see Hughes book
        if PSPS == 'n'
        Dmat = mu*diag([2 2 1]) + lam*[1; 1; 0]*[1 1 0];
        else
        Dmat = PatchE/(1-Patchv^2)*[1      Patchv  0
                                  Patchv  1      0
                                  0      0      (1-Patchv)/2];
        end

        der = 0;
        bf = 0;
        ib = 0;
        
        ulres = reshape(ul,ndf*nen,1);
        
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
            
            % Form B matrix
            for ie = 1:nel
                
              Bmat(Bcol1,(ie-1)*ndf+1) = shg(ie,col1);
              Bmat(Bcol2,(ie-1)*ndf+2) = shg(ie,col2);
                 
            end
            
            xint = xl(:,1:nel)*shl;
            
            epsil = Bmat*ulres(1:nel*ndf);
            sigma = Dmat*epsil;
            dqdx = ElemQ*shg;
            gradu = ul*shg;
            
            % from FEAP
            Weng = 1/2*sigma'*epsil;

        enmom(1,1) = sigma(1)*gradu(1,1) ...
                   + sigma(3)*gradu(2,1) - Weng;
        enmom(1,2) = sigma(1)*gradu(1,2) ...
                   + sigma(3)*gradu(2,2);
        enmom(2,1) = sigma(3)*gradu(1,1) ...
                   + sigma(2)*gradu(2,1);
        enmom(2,2) = sigma(3)*gradu(1,2) ...
                   + sigma(2)*gradu(2,2) - Weng;
%         enmom(3,3) = sig(3,l)*gradu(3,3) - Weng;
            
% c       Accumulate integral

%         for k = 1:nel
%           for i = 1:ndm
%             thi = shp(i,k,l)*dvol(l);
%             thv = thi*th(k);
%             for j = 1:ndm
              ElemJ = ElemJ + enmom'*dqdx'*c1;
%             end %do ! j
%           end %do ! i
%         end %do ! k
    % Verified the above in-line to be the same as the loop or adding
    % separate numbers together.
            
        end %je
        
        ElemJ = -ElemJ;
    case 26
        
        ElemS = zeros(npstr,1);
        stresID = [1 2 3 0 1 3];

        lam = Patchv*PatchE/((1+Patchv)*(1-2*Patchv));
        mu = PatchE/(2*(1+Patchv));
        thick = 1;
        fbx = 0;
        fby = 0;
        fbz = 0;
        if PSPS == 'n'
        Dmat = mu*diag([2 2 1]) + lam*[1; 1; 0]*[1 1 0];
        else
        Dmat = PatchE/(1-Patchv^2)*[1      Patchv  0
                                  Patchv  1      0
                                  0      0      (1-Patchv)/2];
        end
        Nmat = zeros(2,ndf*nel);
        Bmat = zeros(3,ndf*nel);
        I1 = [1; 1; 0];
        
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
        
        ulres = reshape(ul,nst,1);

        %Stress Loop
        ll = 1;
%         for ll = 1:nint

            %Evaluate first derivatives of basis functions at int. point
            if nel == 3 || nel == 6
              [w,litr,lits] =  intpntt(ll,ll,ib);
              [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
              [shg,shgs,Jdet,be] = shgt(xl,nel,shld,shls,nel,bf,der,be);
            else
              [w,litr,lits] =  intpntq(ll,ll,ib);
              [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
              [shg,shgs,Jdet,be] = shgq(xl,nel,shld,shls,nel,bf,der,be);
            end
            
            % Form B matrix
            for ie = 1:nel
              
              Bmat(Bcol1,(ie-1)*ndf+1) = shg(ie,col1);
              Bmat(Bcol2,(ie-1)*ndf+2) = shg(ie,col2);
                 
            end
            
            epsil = Bmat*ulres(1:ndf*nel);
            sigma = Dmat*epsil;
            
            for stres = 1:npstr
            
            if stres <= 3 % stress components
                sigmas = sigma(stresID(stres));
            elseif stres >= 5
                if stres <= 6 % principal stresses
                    if PSPS == 'n'
                        sigz = lam*(epsil(1)+epsil(2));
                    else
                        sigz = 0;
                    end
                    sigma2 = [sigma(1) sigma(3) 0; sigma(3) sigma(2) 0; 0 0 sigz];
                    psig = eig(sigma2);
                    sigmas = psig(stresID(stres));
                else % hydrostatic stress
                    if PSPS == 'n'
                        sigz = lam*(epsil(1)+epsil(2));
                    else
                        sigz = 0;
                    end
                    sigmas = 1/3*(sigma'*I1 + sigz);
                end
            else % von Mises stress
                if PSPS == 'n'
                    sigz = lam*(epsil(1)+epsil(2));
                else
                    sigz = 0;
                end
                trs = sigma'*I1 + sigz;
                dsig = sigma - 1/3*trs*I1;
                sigmas = sqrt(3/2*(dsig'*dsig));
            end
            
            ElemS(stres) = sigmas;
            
            end

%         end %je
        
end %Task Switch
