% Tim Truster
% 07/16/2013
%
% Routine to plot fine scales along interface for SS13N49U1.m with 2 Q9
% elements through the thickness.
% Run the input file and NL_FEA_Program first to generate the solution

%Background of Q9
plotModelContP(Coordinates(:,[2 1]), 0*Node_U_V(:,2), ix([4 8],:), 2, ndm, nen, 2, 1, 1, '')

for elem = numel-numSI+1:numel-numSI+8
    
    %Extract patch nodal coordinates
    
    ElemFlag = ix(elem,1:nen);
    actnode = find(ElemFlag>0);
    xl = zeros(ndm,nen);
    xl(1:ndm,actnode) = Coordinates(ElemFlag(actnode),1:ndm)';

    [EGDOFTa, EGDOFTi, ELDOFTa, ELDOFTi] = plocal(NDOFT, ElemFlag, squeeze(iedof(:,:,ma)), actnode, nen, ndf, neq);
    
    %Extract patch solution values
    ul = zeros(ndf, nen);
    ul_n = zeros(ndf, nen);
    uld = zeros(ndf, nen);
    ul(ELDOFTa) = ModelDx(EGDOFTa)';
    ul(ELDOFTi) = gBC(EGDOFTi)';
    ul_n(ELDOFTa) = ModelDxn_1(EGDOFTa)';
    ul_n(ELDOFTi) = gBC_n(EGDOFTi)';
    uld(ELDOFTa) = s_del_ModelDx(EGDOFTa)';
    uld(ELDOFTi) = (gBC(EGDOFTi) - gBC_n(EGDOFTi))';
    
    CGtoDGarrays

    inter = elem - (numel - numSI);
    nodeAR = SurfacesI(inter,1);
    nodeBR = SurfacesI(inter,2);
    nodeAL = SurfacesI(inter,3);
    nodeBL = SurfacesI(inter,4);

    nelLP = nelL;
    nelRP = nelR;
    
        ElemER = matepropR(1);
        ElemvR = matepropR(2);
        if MatTypeTable(2,maR) == 1
        lamdaR = ElemvR*ElemER/((1+ElemvR)*(1-2*ElemvR));
        muR = ElemER/(2*(1+ElemvR));
        DmatR = muR*diag([2 2 1]) + lamdaR*[1; 1; 0]*[1 1 0];
        elseif MatTypeTable(2,maR) == 3
        lamdaR = ElemvR*ElemER/((1+ElemvR)*(1-2*ElemvR));
        muR = ElemER/(2*(1+ElemvR));
        DmatR = muR*diag([2 2 1]);
        elseif MatTypeTable(2,maR) == 7 % Darcy
        muR = 0;
        DmatR = muR*diag([2 2 1]);
        end
        
        ElemEL = matepropL(1);
        ElemvL = matepropL(2);
        if MatTypeTable(2,maL) == 1
        lamdaL = ElemvL*ElemEL/((1+ElemvL)*(1-2*ElemvL));
        muL = ElemEL/(2*(1+ElemvL));
        DmatL = muL*diag([2 2 1]) + lamdaL*[1; 1; 0]*[1 1 0];
        elseif MatTypeTable(2,maL) == 3
        lamdaL = ElemvL*ElemEL/((1+ElemvL)*(1-2*ElemvL));
        muL = ElemEL/(2*(1+ElemvL));
        DmatL = muL*diag([2 2 1]);
        elseif MatTypeTable(2,maL) == 7 % Darcy
        muL = 0;
        DmatL = muL*diag([2 2 1]);
        end
        
        ulresL = reshape(ulL,ndf*nelL,1);
        ulresR = reshape(ulR,ndf*nelR,1);
        
        NmatL = zeros(3,nstL);
        BmatL = zeros(4,nstL);
        bnAdN1 = zeros(4,nstL);
        N1 = zeros(3,nstL);
        NmatR = zeros(3,nstR);
        BmatR = zeros(4,nstR);
        bnAdN2 = zeros(4,nstR);
        N2 = zeros(3,nstR);
        
        % Determin bounds of integration segment
        InterConn2D2 % InterConn2DT % 
        
        h = 2/(hR + hL);
        
        % Set jacobian for integration space
        drdr = (eL2 - eL1)/drL;
        
        m = (eR2-eR1)/(eL1-eL2);

        s = -1;
        
% For separate bubble types on T and Q
        if MatTypeTable(2,maL) ~= 7 % Darcy
        [tauL,intbL] = TauE1_2d(xlintL,DmatL,nelL,lintt6,lintq9);
        else
            tauL = 0;
        end
        if MatTypeTable(2,maR) ~= 7 % Darcy
        [tauR,intbR] = TauE1_2d(xlintR,DmatR,nelR,lintt6,lintq9);
        else
            tauR = 0;
        end
        
        lint = 3;
        ib = 0;
        der = 0;
        bf = 0;
        ebL = 0;
        ebR = 0;
        intedge = 0;
        
        % THIS LOOP COULD BE EVALUATED IN CLOSED-FORM
        for ie = 1:lint
            
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
        Kinv = [ep zeros(2,1); zeros(1,2) 0];
        elseif nitvms == 2
        % Nitsche
        volL = getvol(xlL,nelL);
        volR = getvol(xlR,nelR);
%         h = 2/(intedge/volR + intedge/volL);
        gamL = 1/2*eye(2);
        gamR = 1/2*eye(2);
        if MatTypeTable(2,maL) == 7 || MatTypeTable(2,maR) == 7 % Darcy
        ep = pencoeff*([nLx nLy]'*[nLx nLy])*max(muL,muR)/h;
        else
        ep = pencoeff*eye(2)*max(muL,muR)/h;
        end
        ePP = -pencoeff*2/max(muL,muR)*h;
        Kinv = [ep zeros(2,1); zeros(1,2) ePP];
        end

        resL = zeros(2,1);
        resR = resL;
        
        for ie = 1:lint
            
% For separate bubble types on T and Q
            if nelL == 3 || nelL == 6
                [Wgt,r,s] = intpntt(ie,lint,1);
                ebeL = edgebubble(r,s,nelL);
            elseif nelL == 4 || nelL == 9
                [Wgt,r,s] = intpntq(ie,lint,1);
                ebeL = edgebubbleQ(r,s,nelL);
            end
            % Okay for now because the bubble is symmetric; otherwise need
            % to map it over
            if nelL == 3 || nelL == 6
                [Wgt,rR,sR] = intpntt(ie,lint,1);
                ebeR = edgebubble(rR,sR,nelR);
            elseif nelL == 4 || nelL == 9
                [Wgt,rR,sR] = intpntq(ie,lint,1);
                ebeR = edgebubbleQ(rR,sR,nelR);
            end
            
            if nelL == 3 || nelL == 6
                [Wgt,r,s] = intpntt(ie,lint,1);
            elseif nelL == 4 || nelL == 9
                [Wgt,r,s] = intpntq(ie,lint,1);
            end
            
            r = drdr*(r-roL)+eL1;

            if nelLP == 3 || nelLP == 6
                [shpL,shld,shls,be] = shlt(r,s,nelLP,nelLP,0,0);
                [PxyL, shgs, Jdet] = shgt(xlL(:,1:nelLP),nelLP,shld,shls,nen,0,0,be);
            elseif nelLP == 4 || nelLP == 9
                [shpL,shld,shls,be] = shlq(r,s,nelLP,nelLP,0,0);
                [PxyL, shgs, Jdet] = shgq(xlL(:,1:nelLP),nelLP,shld,shls,nen,0,0,be);
            end
            
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
            
            if nelRP == 3 || nelRP == 6
                [shpR,shld,shls,be] = shlt(rR,s,nelRP,nelRP,0,0);
                [PxyR, shgs, Jdet, be, xsR] = shgt(xlR(:,1:nelRP),nelRP,shld,shls,nen,0,0,be);
            elseif nelRP == 4 || nelRP == 9
                [shpR,shld,shls,be] = shlq(rR,s,nelRP,nelRP,0,0);
                [PxyR, shgs, Jdet, be, xsR] = shgq(xlR(:,1:nelRP),nelRP,shld,shls,nen,0,0,be);
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
                      0 nLy nLx]; %- ?
            
            c1 = Wgt*tm3*drdr*thick;
                
            for i = 1:nelL
%                 NmatL(:,3*i-2:3*i) = [shlL(i)*eye(2) zeros(2,1)];
                NmatL(:,(i-1)*ndf+1:(i-1)*ndf+3) = shlL(i)*eye(3);
                BmatL(:,(i-1)*ndf+1:(i-1)*ndf+3) = [QxyL(i,1) 0 0
                                      0 QxyL(i,2) 0
                                      QxyL(i,2) QxyL(i,1) 0
                                      0 0 shlL(i)];
            end
            
            for i = 1:nelR
%                 NmatR(:,3*i-2:3*i) = [shlR(i)*eye(2) zeros(2,1)];
                NmatR(:,(i-1)*ndf+1:(i-1)*ndf+3) = shlR(i)*eye(3);
                BmatR(:,(i-1)*ndf+1:(i-1)*ndf+3) = [QxyR(i,1) 0 0
                                      0 QxyR(i,2) 0
                                      QxyR(i,2) QxyR(i,1) 0
                                      0 0 shpR(i)];
            end
               
            
            N1 = NmatL;
            N2 = NmatR;
            
            bnAdN1 = gamL*nvect*[DmatL [1; 1; 0]]*BmatL;
            bnAdN2 = gamR*nvect*[DmatR [1; 1; 0]]*BmatR;
            jumpu = NmatL*ulresL - NmatR*ulresR;
            tvtr = bnAdN1*ulresL + bnAdN2*ulresR;
            lammult = tvtr - Kinv(1:2,1:2)*jumpu(1:2);
            
            resL = resL + c1*ebeL*(lammult - nvect*[DmatL [1; 1; 0]]*BmatL*ulresL);
            resR = resR + c1*ebeR*(-lammult + nvect*[DmatR [1; 1; 0]]*BmatR*ulresR);

        end %ie
        
        betaL = tauL*resL;
        betaR = tauR*resR;
        
        plotFSelembub(xlintL([2 1],:)', betaL(1), ndm, nelL, 2, 1, 1, '')
        plotFSelembub(xlintR([2 1],:)', betaR(1), ndm, nelR, 2, 1, 1, '')
    
end

%Mesh lines
plotModelDefoP(Coordinates(:,[2 1]), ix([4 8 (17:8:17+8*8)],:), 10, nen, 2, 1, 1, '', 'f', 'n')

axis on
axis equal
colorbar('EastOutside','FontSize',14,'FontName','Times New Roman','LineWidth',1)
set(gca,'FontName','Times New Roman','FontSize',14)
xlabel('y','FontName','Times New Roman','FontSize',14)
ylabel('x','FontName','Times New Roman','FontSize',14)