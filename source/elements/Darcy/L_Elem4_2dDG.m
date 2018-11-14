% Tim Truster
% 09/06/2012
% Mixed, stabilized velocity-pressure element for Darcy flow,
% CG or DG treatment. Part of Inter_FEA_Program.
% T3, Q4, T6, and Q9 supported
% Based off of original in Darcy Flow folder
% Updated 9/19/2012 with average bubble definitions

% Signs flipped to fit with elasticity for Darcy-Stokes


if isw ~= 1
CGtoDGarrays

inter = elem - (numel - numSI);
nodeAR = SurfacesI(inter,1);
nodeBR = SurfacesI(inter,2);
nodeAL = SurfacesI(inter,3);
nodeBL = SurfacesI(inter,4);

nelLP = nelL;
nelRP = nelR;
end

nitvms = 1;
if nitvms == 1 %VMS
pencoeffv = 1;1/40;
pencoeffp = 1/pencoeffv;
elseif nitvms == 2 %Nitsche
pencoeff = 1/4;1/80;1;2;
else %RFB
pencoeff = 1;4;2;
end

nelV = nel;
alpha = 0;1;

Bcol1 = [4; 6];
Bcol2 = [5; 7];
col1 = [1; 2];
col2 = [1; 2];

switch isw %Task Switch
    
    case 1
        
        if ndf > 3
            
            for i = 4:ndf
                lie(i,1) = 0;
            end
            
        end
%%    
    case 3 %interface stiffness
        
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
        
        kappaL = matepropL(1);
        muL = matepropL(2);
        gcL = matepropL(3);
        rhoL = matepropL(4);
        kappaR = matepropR(1);
        muR = matepropR(2);
        gcR = matepropR(3);
        rhoR = matepropR(4);
        mukappaL = muL/kappaL;
        kappamuL = kappaL/muL;
        rhogcL = rhoL/gcL;
        mukappaR = muR/kappaR;
        kappamuR = kappaR/muR;
        rhogcR = rhoR/gcR;
        
        Dmatv = diag([1 1 0]);
        Dmatp = diag([0 0 1]);
        
        NmatL = zeros(3,nstL);
        bnAdN1 = zeros(3,nstL);
        NmatR = zeros(3,nstR);
        bnAdN2 = zeros(3,nstR);
        
        ulresL = reshape(ulL,ndf*nelL,1);
        ulresR = reshape(ulR,ndf*nelR,1);
        
        % Determin bounds of integration segment
        InterConn2D2 % InterConn2DT % 
        
%         h = 2/(hR + hL);

        % Set jacobian for integration space
        drdr = (eL2 - eL1)/drL;
        
        m = (eR2-eR1)/(eL1-eL2);
        
        % Perform integration of various matrices
        
        [tauL,intbL] = TauE4_2d(xlintL,mukappaL,nelL,lintt6,lintq9);
        [tauR,intbR] = TauE4_2d(xlintR,mukappaR,nelR,lintt6,lintq9);
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
        
        if nitvms == 1 %VMS
        bLave = ebL/intedge;
        bRave = ebR/intedge;
        taubL = tauL*ebL*bLave;
        taubR = tauR*ebR*bRave;
        edgeK = taubL + taubR;
        gamL = taubL/edgeK;
        gamR = taubR/edgeK;
        ep = pencoeffv*inv(edgeK);
%         ep2 = -pencoeffp*taubL*taubR/edgeK; % 5/4/13
        ep2 = pencoeffp*taubL*taubR/edgeK;
        elseif nitvms == 2 %Nitsche
        volL = getvol2D(xlL,[lintt3,lintq4,lintt6,lintq9],nelL);
        volR = getvol2D(xlR,[lintt3,lintq4,lintt6,lintq9],nelR);
        h = 2/(intedge/volR + intedge/volL);
        gamL = 1/2;
        gamR = 1/2;
%         ep = 10*max(mukappaL,mukappaR)/h;
        ep = pencoeff*max(mukappaL,mukappaR)*h;
        ep2 = -1/ep;
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
                shlpL = shlt(r,s,nelLP,nelL,der,bf);
            elseif nelL == 4 || nelL == 9
                [shlL,shld,shls,be] = shlq(r,s,nelL,nelL,der,bf);
                [QxyL, shgs, Jdet, be, xsL] = shgq(xlL(:,1:nelL),nelL,shld,shls,nen,bf,der,be);
                shlpL = shlq(r,s,nelLP,nelL,der,bf);
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
                shlpR = shlt(rR,s,nelRP,nelR,der,bf);
            elseif nelR == 4 || nelR == 9
                [shlR,shld,shls,be] = shlq(rR,s,nelR,nelR,der,bf);
                [QxyR, shgs, Jdet, be, xsR] = shgq(xlR(:,1:nelR),nelR,shld,shls,nen,bf,der,be);
                shlpR = shlq(rR,s,nelRP,nelR,der,bf);
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
            nvect = [0 0 nLx
                     0 0 nLy
                     0 0 0];
            nvec = [nLx; nLy];
                    
            c1 = Wgt*tm3*drdr*thick;
            
            for i = 1:nelL
                NmatL(1,(i-1)*ndf+1) = shlL(i);
                NmatL(2,(i-1)*ndf+2) = shlL(i);
            end
            for i = 1:nelLP
                NmatL(3,(i-1)*ndf+3) = shlpL(i);
            end
            
            for i = 1:nelR
                NmatR(1,(i-1)*ndf+1) = shlR(i);
                NmatR(2,(i-1)*ndf+2) = shlR(i);
            end
            for i = 1:nelRP
                NmatR(3,(i-1)*ndf+3) = shlpR(i);
            end
            
            bnAdN1 = gamL*nvect*NmatL;
            bnAdN2 = gamR*nvect*NmatR;
            
            tvtr = (bnAdN1*ulresL + bnAdN2*ulresR); %average pressure
            jump = NmatR*ulresR - NmatL*ulresL;        %v-p jump
%             jumpv = jump(1:2);
%             jumpp = jump(3);

            % Not needed for linear problems
            ElemFL = ElemFL - c1*( - NmatL'*(-tvtr + ep*Dmatv*jump - ep2*Dmatp*jump) + bnAdN1'*Dmatv*jump);
            ElemFR = ElemFR - c1*( + NmatR'*(-tvtr + ep*Dmatv*jump - ep2*Dmatp*jump) + bnAdN2'*Dmatv*jump);

%             ElemKLL = ElemKLL + c1*NmatL'*Dmatv*bnAdN1 - c1*bnAdN1'*Dmatv*NmatL;
%             ElemKLR = ElemKLR + c1*NmatL'*Dmatv*bnAdN2 + c1*bnAdN1'*Dmatv*NmatR;
%             ElemKRL = ElemKRL - c1*NmatR'*Dmatv*bnAdN1 - c1*bnAdN2'*Dmatv*NmatL;
%             ElemKRR = ElemKRR - c1*NmatR'*Dmatv*bnAdN2 + c1*bnAdN2'*Dmatv*NmatR;
            ElemKLL = ElemKLL + c1*NmatL'*Dmatv*bnAdN1 + c1*bnAdN1'*Dmatv*NmatL; % 5/4/13
            ElemKLR = ElemKLR + c1*NmatL'*Dmatv*bnAdN2 - c1*bnAdN1'*Dmatv*NmatR;
            ElemKRL = ElemKRL - c1*NmatR'*Dmatv*bnAdN1 + c1*bnAdN2'*Dmatv*NmatL;
            ElemKRR = ElemKRR - c1*NmatR'*Dmatv*bnAdN2 - c1*bnAdN2'*Dmatv*NmatR;

%             ElemKLL = ElemKLL + c1*(NmatL'*ep*Dmatv*NmatL) - c1*(NmatL'*1/ep*Dmatp*NmatL);
%             ElemKLR = ElemKLR - c1*(NmatL'*ep*Dmatv*NmatR) + c1*(NmatL'*1/ep*Dmatp*NmatR);
%             ElemKRL = ElemKRL - c1*(NmatR'*ep*Dmatv*NmatL) + c1*(NmatR'*1/ep*Dmatp*NmatL);
%             ElemKRR = ElemKRR + c1*(NmatR'*ep*Dmatv*NmatR) - c1*(NmatR'*1/ep*Dmatp*NmatR);
            ElemKLL = ElemKLL + c1*(NmatL'*ep*Dmatv*NmatL) - c1*(NmatL'*ep2*Dmatp*NmatL);
            ElemKLR = ElemKLR - c1*(NmatL'*ep*Dmatv*NmatR) + c1*(NmatL'*ep2*Dmatp*NmatR);
            ElemKRL = ElemKRL - c1*(NmatR'*ep*Dmatv*NmatL) + c1*(NmatR'*ep2*Dmatp*NmatL);
            ElemKRR = ElemKRR + c1*(NmatR'*ep*Dmatv*NmatR) - c1*(NmatR'*ep2*Dmatp*NmatR);
            
        end
        
        ElemKLL;
            ElemK = [ElemKLL ElemKLR; ElemKRL ElemKRR];
            ElemF = [ElemFL; ElemFR];
%%
    case 11
        
        ElemE = zeros(numEn,1);
        
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
                [Wgt,r,s] = intpntt(je,lint,0);
            elseif nel == 4 || nel == 9
                [Wgt,r,s] = intpntq(je,lint,0);
            end

            % Evaluate 1-D basis functions at integration points
            [shp,shp2,be] = shpl_2d(r,s,nel,ideriv,ep,enrich);
            %Evaluate first derivatives of basis functions at int. point
            [shg, shgs, Jdet, xs, be] = shpg_2d(shp,shp2,xl,nel2,ideriv,be);
            c1 = Wgt*Jdet*thick;
            shp = shpl_2d(r,s,nelP,ideriv,0,0);
            shpS = shpl_2d(r,s,nelS,ideriv,0,0);
            
            xint = shg(3,:)*xl(1,:)';
            yint = shg(3,:)*xl(2,:)';
            dux = ul(1:2,:)*shg(1,:)';
            duy = ul(1:2,:)*shg(2,:)';
            u = ul(1:2,:)*shg(3,:)';
            strains = [2*dux(1) 2*duy(2) dux(2)+duy(1)];
            px = ul(3,1:nelP)*shp(1,:)'*pfact;
            py = ul(3,1:nelP)*shp(2,:)'*pfact;
            p = ul(3,1:nelP)*shp(3,:)'*pfact;
            
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
            
            % Evaluate 1-D basis functions at integration points
            [shpL,shp2,be] = shpl_2d(r,s,nelL,ideriv,epL,enrichL);
            %Evaluate first derivatives of basis functions at int. point
            [QxyL, cartd2, Jdet, xsL] = shpg_2d(shpL,shp2,xlL,nel2L,ideriv,be); %#ok<ASGLU>
            shpL = shpl_2d(r,s,nelLP,ideriv,0,0);
            shpLS = shpl_2d(r,s,nelLS,ideriv,0,0);
            
            rR = m*(r-eL2) + eR1;
            
            if nelR == 3 || nelR == 6
                s = 0;
            else %if nelR == 4
                s = -1;
            end
            
            % Evaluate 1-D basis functions at integration points
            [shpR,shp2,be] = shpl_2d(rR,s,nelR,ideriv,epR,enrichR);
            %Evaluate first derivatives of basis functions at int. point
            [QxyR, cartd2, Jdet, xsR] = shpg_2d(shpR,shp2,xlR,nel2R,ideriv,be);
            % Evaluate 1-D basis functions at integration points
            shpR = shpl_2d(rR,s,nelRP,ideriv,0,0);
            shpRS = shpl_2d(rR,s,nelRS,ideriv,0,0);
                    
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
            
            dux = ulL(1:2,:)*QxyL(1,:)';
            duy = ulL(1:2,:)*QxyL(2,:)';
            strainsL = [2*dux(1) 2*duy(2) dux(2)+duy(1)];
            pL = ulL(3,1:nelLP)*shpL(3,:)'*pfact;
            sigmasL = pL*II(inds(1,stres),inds(2,stres)) + mu*strainsL(stres);
            dux = ulR(1:2,:)*QxyR(1,:)';
            duy = ulR(1:2,:)*QxyR(2,:)';
            pR = ulR(3,1:nelRP)*shpR(3,:)'*pfact;
            strainsR = [2*dux(1) 2*duy(2) dux(2)+duy(1)];
            sigmasR = pR*II(inds(1,stres),inds(2,stres)) + mu*strainsR(stres);
            
            sigmas = sigmasR-sigmasL;
                    
            c1 = Wgt*tm3*drdr*thick;
                
            for a = 1:nelLS
            
                NLa = shpLS(3,a)*c1;

                for c = 1:nelLS

                    NLc = shpLS(3,c);

                    ElemKLL(c,a) = ElemKLL(c,a) - ...
                                             10*(NLc*NLa);
                                         
                end
                                         
                for d = 1:nelRS

                    NRd = shpRS(3,d);

                    ElemKRL(d,a) = ElemKRL(d,a) + ...
                                             10*(NRd*NLa);
                                         
                end

            end
            
            for b = 1:nelRS
            
                NRb = shpRS(3,b)*c1;

                for c = 1:nelLS

                    NLc = shpLS(3,c);

                    ElemKLR(c,b) = ElemKLR(c,b) + ...
                                             10*(NLc*NRb);
                                         
                end
                                         
                for d = 1:nelRS

                    NRd = shpRS(3,d);

                    ElemKRR(d,b) = ElemKRR(d,b) - ...
                                             10*(NRd*NRb);
                                         
                end

            end
            
        end %lint
        
        ElemKLL = 1/10*ElemKLL;
        ElemKRL = 1/10*ElemKRL;
        ElemKLR = 1/10*ElemKLR;
        ElemKRR = 1/10*ElemKRR;
        
    case 60
        
        numhr = 3;
        ElemI = zeros(13,numhr);
        t1 = zeros(3,1);
        t2 = t1;
        t3 = t1;
        thick = 1;
        
        bf1 = 0;
        bf2 = 0;
        
        kappaL = matepropL(1);
        muL = matepropL(2);
        gcL = matepropL(3);
        rhoL = matepropL(4);
        kappaR = matepropR(1);
        muR = matepropR(2);
        gcR = matepropR(3);
        rhoR = matepropR(4);
        mukappaL = muL/kappaL;
        kappamuL = kappaL/muL;
        rhogcL = rhoL/gcL;
        mukappaR = muR/kappaR;
        kappamuR = kappaR/muR;
        rhogcR = rhoR/gcR;
        
        Dmatv = diag([1 1 0]);
        Dmatp = diag([0 0 1]);
        
        NmatL = zeros(3,nstL);
        bnAdN1 = zeros(3,nstL);
        NmatR = zeros(3,nstR);
        bnAdN2 = zeros(3,nstR);
        
        ulresL = reshape(ulL,ndf*nelL,1);
        ulresR = reshape(ulR,ndf*nelR,1);
        
        % Determin bounds of integration segment
        InterConn2D2 % InterConn2DT % 
        
%         h = 2/(hR + hL);

        % Set jacobian for integration space
        drdr = (eL2 - eL1)/drL;
        
        m = (eR2-eR1)/(eL1-eL2);
        
        % Perform integration of various matrices
        
        [tauL,intbL] = TauE4_2d(xlintL,mukappaL,nelL,lintt6,lintq9);
        [tauR,intbR] = TauE4_2d(xlintR,mukappaR,nelR,lintt6,lintq9);
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
                ebeL = edgebubble(r,s,nel);
            elseif nelL == 4 || nelL == 9
                [Wgt,r,s] = intpntq(ie,lint,1);
                ebeL = edgebubbleQ(r,s,nel);
            end
                    
            if nelL == 3 || nelL == 6
                [Wgt,rR,sR] = intpntt(ie,lint,1);
                ebeR = edgebubble(rR,sR,nel);
            elseif nelL == 4 || nelL == 9
                [Wgt,rR,sR] = intpntq(ie,lint,1);
                ebeR = edgebubbleQ(rR,sR,nel);
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
        
        
        if nitvms == 1 %VMS
        bLave = ebL/intedge;
        bRave = ebR/intedge;
        taubL = tauL*ebL*bLave;
        taubR = tauR*ebR*bRave;
        edgeK = taubL + taubR;
        gamL = taubL/edgeK;
        gamR = taubR/edgeK;
        ep = pencoeff*inv(edgeK);
        ep2 = -pencoeff*taubL*taubR/edgeK;
        elseif nitvms == 2 %Nitsche
        volL = getvol2D(xlL,[lintt3,lintq4,lintt6,lintq9],nelL);
        volR = getvol2D(xlR,[lintt3,lintq4,lintt6,lintq9],nelR);
        h = 2/(intedge/volR + intedge/volL);
        gamL = 1/2;
        gamR = 1/2;
%         ep = 10*max(mukappaL,mukappaR)/h;
        ep = pencoeff*max(mukappaL,mukappaR)*h;
        ep2 = -1/ep;
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
                shlpL = shlt(r,s,nelLP,nelL,der,bf);
            elseif nelL == 4 || nelL == 9
                [shlL,shld,shls,be] = shlq(r,s,nelL,nelL,der,bf);
                [QxyL, shgs, Jdet, be, xsL] = shgq(xlL(:,1:nelL),nelL,shld,shls,nen,bf,der,be);
                shlpL = shlq(r,s,nelLP,nelL,der,bf);
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
                shlpR = shlt(rR,s,nelRP,nelR,der,bf);
            elseif nelR == 4 || nelR == 9
                [shlR,shld,shls,be] = shlq(rR,s,nelR,nelR,der,bf);
                [QxyR, shgs, Jdet, be, xsR] = shgq(xlR(:,1:nelR),nelR,shld,shls,nen,bf,der,be);
                shlpR = shlq(rR,s,nelRP,nelR,der,bf);
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
            nvect = [0 0 nLx
                     0 0 nLy
                     0 0 0];
            nvec = [nLx; nLy];
                    
            c1 = Wgt*tm3*drdr*thick;
            
            for i = 1:nelL
                NmatL(1,(i-1)*ndf+1) = shlL(i);
                NmatL(2,(i-1)*ndf+2) = shlL(i);
            end
            for i = 1:nelLP
                NmatL(3,(i-1)*ndf+3) = shlpL(i);
            end
            
            for i = 1:nelR
                NmatR(1,(i-1)*ndf+1) = shlR(i);
                NmatR(2,(i-1)*ndf+2) = shlR(i);
            end
            for i = 1:nelRP
                NmatR(3,(i-1)*ndf+3) = shlpR(i);
            end
            
            bnAdN1 = gamL*nvect*NmatL;
            bnAdN2 = gamR*nvect*NmatR;
            
            tvtr = (bnAdN1*ulresL + bnAdN2*ulresR); %average pressure
            jump = NmatR*ulresR - NmatL*ulresL;        %v-p jump
%             jumpv = jump(1:2);
%             jumpp = jump(3);

        ElemI(:,ie) = [xint
                       yint
                       0
                       jump
                       tvtr
                       (tvtr + [ep; ep; -ep2].*jump)
                       ep];
        end %lint
%%        
end %Task Switch
