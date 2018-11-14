
% 05/2013
% DG Darcy element

% Sign convention uses pressure positive in compression, as adopted in the
% IJNMF paper.


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
pencoeffp = -1/pencoeffv;
elseif nitvms == 2 %Nitsche
pencoeff = 1/4;1/80;1;2;
else %RFB
pencoeff = 1;4;2;
end

nelV = nel;
alpha = 0;1;
thick = 1;

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
            Dmatv(1:2,1:2) = nvec*nvec'; %5/22/2013, only normal velocity is penalized
                    
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
            bnAdNL = nvect*NmatL;
            bnAdNR = nvect*NmatR;
            
            tvtr = (bnAdN1*ulresL + bnAdN2*ulresR); %average pressure
            jumpflux = (bnAdNL*ulresL - bnAdNR*ulresR); %flux jump
            jump = NmatR*ulresR - NmatL*ulresL;        %v-p jump
%             jumpv = jump(1:2);
%             jumpp = jump(3);

            % Not needed for linear problems
%             ElemFL = ElemFL - c1*( - NmatL'*(-tvtr + ep*Dmatv*jump - ep2*Dmatp*jump) + bnAdN1'*Dmatv*jump);
%             ElemFR = ElemFR - c1*( + NmatR'*(-tvtr + ep*Dmatv*jump - ep2*Dmatp*jump) + bnAdN2'*Dmatv*jump);
            ElemFL = ElemFL - c1*( - NmatL'*(-tvtr + ep*Dmatv*jump) + bnAdN1'*(Dmatv*jump) - bnAdNL'*(ep2*jumpflux));
            ElemFR = ElemFR - c1*( + NmatR'*(-tvtr + ep*Dmatv*jump) + bnAdN2'*(Dmatv*jump) + bnAdNR'*(ep2*jumpflux));

            ElemKLL = ElemKLL + c1*NmatL'*Dmatv*bnAdN1 - c1*bnAdN1'*Dmatv*NmatL;
            ElemKLR = ElemKLR + c1*NmatL'*Dmatv*bnAdN2 + c1*bnAdN1'*Dmatv*NmatR;
            ElemKRL = ElemKRL - c1*NmatR'*Dmatv*bnAdN1 - c1*bnAdN2'*Dmatv*NmatL;
            ElemKRR = ElemKRR - c1*NmatR'*Dmatv*bnAdN2 + c1*bnAdN2'*Dmatv*NmatR;

%             ElemKLL = ElemKLL + c1*(NmatL'*ep*Dmatv*NmatL) - c1*(NmatL'*1/ep*Dmatp*NmatL);
%             ElemKLR = ElemKLR - c1*(NmatL'*ep*Dmatv*NmatR) + c1*(NmatL'*1/ep*Dmatp*NmatR);
%             ElemKRL = ElemKRL - c1*(NmatR'*ep*Dmatv*NmatL) + c1*(NmatR'*1/ep*Dmatp*NmatL);
%             ElemKRR = ElemKRR + c1*(NmatR'*ep*Dmatv*NmatR) - c1*(NmatR'*1/ep*Dmatp*NmatR);
%             ElemKLL = ElemKLL + c1*(NmatL'*ep*Dmatv*NmatL) - c1*(NmatL'*ep2*Dmatp*NmatL); % changed 5/22/2013, using flux instead of jump in p
%             ElemKLR = ElemKLR - c1*(NmatL'*ep*Dmatv*NmatR) + c1*(NmatL'*ep2*Dmatp*NmatR); % this makes a difference of a factor of 4 when gam = 1/2
%             ElemKRL = ElemKRL - c1*(NmatR'*ep*Dmatv*NmatL) + c1*(NmatR'*ep2*Dmatp*NmatL);
%             ElemKRR = ElemKRR + c1*(NmatR'*ep*Dmatv*NmatR) - c1*(NmatR'*ep2*Dmatp*NmatR);
            ElemKLL = ElemKLL + c1*(NmatL'*ep*Dmatv*NmatL) - c1*(bnAdNL'*ep2*bnAdNL);
            ElemKLR = ElemKLR - c1*(NmatL'*ep*Dmatv*NmatR) + c1*(bnAdNL'*ep2*bnAdNR);
            ElemKRL = ElemKRL - c1*(NmatR'*ep*Dmatv*NmatL) + c1*(bnAdNR'*ep2*bnAdNL);
            ElemKRR = ElemKRR + c1*(NmatR'*ep*Dmatv*NmatR) - c1*(bnAdNR'*ep2*bnAdNR);
            
        end
        
        ElemKLL;
            ElemK = [ElemKLL ElemKLR; ElemKRL ElemKRR];
            ElemF = [ElemFL; ElemFR];
%%
    case 11
        
        ElemE = zeros(numEn,1);
        
end %Task Switch
