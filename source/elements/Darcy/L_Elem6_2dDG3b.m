
% 05/2013
% DG Stokes-Darcy element

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
alp = mateprop(5);

% Mixture of Darcy and mixed elasticity (Stokes) DG element
% Signs flipped on pressure

% Set Material Properties

nitvms = 1;
if nitvms == 1 %VMS
pencoeff = 1;
elseif nitvms == 2 %Nitsche
pencoeff = 10;1;2;
else %RFB
pencoeff = 1;4;2;
end

switch isw %Task Switch
%%
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
        
        ElemER = matepropR(1);
        ElemvR = matepropR(2);
        if MatTypeTable(2,maR) == 1 % pure-displacement elasticity
        lamdaR = ElemvR*ElemER/((1+ElemvR)*(1-2*ElemvR));
        muR = ElemER/(2*(1+ElemvR));
        DmatR = muR*diag([2 2 1]) + lamdaR*[1; 1; 0]*[1 1 0];
        elseif MatTypeTable(2,maR) == 3 % Stokes
        lamdaR = ElemvR*ElemER/((1+ElemvR)*(1-2*ElemvR));
        muR = ElemER/(2*(1+ElemvR));
        DmatR = muR*diag([2 2 1]);
        elseif MatTypeTable(2,maR) == 7 % Darcy
        kappaR = ElemER;
        muR = ElemvR;
        mukappaR = muR/kappaR;
        DmatR = 0*diag([2 2 1]);
        end
        
        ElemEL = matepropL(1);
        ElemvL = matepropL(2);
        if MatTypeTable(2,maL) == 1 % pure-displacement elasticity
        lamdaL = ElemvL*ElemEL/((1+ElemvL)*(1-2*ElemvL));
        muL = ElemEL/(2*(1+ElemvL));
        DmatL = muL*diag([2 2 1]) + lamdaL*[1; 1; 0]*[1 1 0];
        elseif MatTypeTable(2,maL) == 3 % Stokes
        lamdaL = ElemvL*ElemEL/((1+ElemvL)*(1-2*ElemvL));
        muL = ElemEL/(2*(1+ElemvL));
        DmatL = muL*diag([2 2 1]);
        elseif MatTypeTable(2,maL) == 7 % Darcy
        kappaL = ElemEL;
        muL = ElemvL;
        mukappaL = muL/kappaL;
        DmatL = 0*diag([2 2 1]);
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
%         h = 0.125;
        
        % Set jacobian for integration space
        drdr = (eL2 - eL1)/drL;
        
        m = (eR2-eR1)/(eL1-eL2);

        s = -1;
        
        % For separate bubble types on T and Q
        if MatTypeTable(2,maL) == 3 % Stokes
        [tauL,intbL] = TauE1_2d(xlintL,DmatL,nelL,lintt6,lintq9);
        elseif MatTypeTable(2,maL) == 7 % Darcy
        [tauL,intbL] = TauE4_2d(xlintL,mukappaL,nelL,lintt6,lintq9);
        else
            tauL = 0;
        end
        if MatTypeTable(2,maR) == 3 % Stokes
        [tauR,intbR] = TauE1_2d(xlintR,DmatR,nelR,lintt6,lintq9);
        elseif MatTypeTable(2,maR) == 7 % Darcy
        [tauR,intbR] = TauE4_2d(xlintR,mukappaR,nelR,lintt6,lintq9);
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
        
        nmat = nvect'*nvect;
        if MatTypeTable(2,maL) == 3 % Stokes
        tauL = nvect*tauL*nvect';
        end
        if MatTypeTable(2,maR) ~= 3 % Stokes
        tauR = nvect*tauR*nvect';
        end        
        
        if nitvms == 1
        % VMS
% %         volL = getvol(xlL,nelL);
% %         volR = getvol(xlR,nelR);
% %         volbL = getvol(xlintL,nelL);
% %         volbR = getvol(xlintR,nelR);
% %         tauL = tauL*(volL/volbL)^-1;
% %         tauR = tauR*(volR/volbR)^-1;
%         edgeK = tauL*ebL^2 + tauR*ebR^2;
%         gamL = ebL^2*(edgeK\tauL);
%         gamR = ebR^2*(edgeK\tauR);
%         ep = pencoeff*intedge*inv(edgeK); %#ok<MINV>
%         Kinv = [ep*nmat zeros(2,1); zeros(1,2) 0];
%         ep2 = -1/ep;
%         Kinvp = ep2*nmat;
        bLave = ebL/intedge;
        bRave = ebR/intedge;
        taubL = tauL*ebL*bLave;
        taubR = tauR*ebR*bRave;
        edgeK = taubL + taubR;
        gamL = taubL/edgeK;
        gamR = taubR/edgeK;
        ep = pencoeff/edgeK;
        ep2 = pencoeff*taubL*taubR/edgeK;
        Kinv = [ep*nmat zeros(2,1); zeros(1,2) 0];
        Kinvp = ep2*nmat;
        elseif nitvms == 2
        % Nitsche
        volL = getvol(xlL,nelL);
        volR = getvol(xlR,nelR);
%         h = 2/(intedge/volR + intedge/volL);
        gamL = 1/2*eye(2);
        gamR = 1/2*eye(2);
        if MatTypeTable(2,maL) == 7 || MatTypeTable(2,maR) == 7 % Darcy
        ep = pencoeff*max(muL,muR)/h;
        Kinv = [ep*nmat zeros(2,1); zeros(1,2) 0];
        ep2 = -1/ep;
        Kinvp = ep2*nmat;
        else
        ep = pencoeff*eye(2)*max(muL,muR)/h;
        ePP = -pencoeff*2/max(muL,muR)*h;
        Kinv = [ep zeros(2,1); zeros(1,2) ePP];
        end
        end

        for ie = 1:lint

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
                                      0 0 -shlL(i)];
            end
            
            for i = 1:nelR
%                 NmatR(:,3*i-2:3*i) = [shlR(i)*eye(2) zeros(2,1)];
                NmatR(:,(i-1)*ndf+1:(i-1)*ndf+3) = shlR(i)*eye(3);
                BmatR(:,(i-1)*ndf+1:(i-1)*ndf+3) = [QxyR(i,1) 0 0
                                      0 QxyR(i,2) 0
                                      QxyR(i,2) QxyR(i,1) 0
                                      0 0 -shlR(i)];
            end
               
%%       
            if numSI > 0
            
            N1 = NmatL;
            N2 = NmatR;
            
            bnAdN1 = nmat*gamL*nvect*[DmatL [1; 1; 0]]*BmatL;
            bnAdN2 = nmat*gamR*nvect*[DmatR [1; 1; 0]]*BmatR;
            bnAdNL = nmat*nvect*[DmatL [1; 1; 0]]*BmatL;
            bnAdNR = nmat*nvect*[DmatR [1; 1; 0]]*BmatR;
            bnAdN1w = nmat*gamL*nvect*[DmatL [1; 1; 0]]*diag([1 1 1 -1])*BmatL;
            bnAdN2w = nmat*gamR*nvect*[DmatR [1; 1; 0]]*diag([1 1 1 -1])*BmatR;
            bnAdNLw = nmat*nvect*[DmatL [1; 1; 0]]*diag([1 1 1 -1])*BmatL;
            bnAdNRw = nmat*nvect*[DmatR [1; 1; 0]]*diag([1 1 1 -1])*BmatR;
            
            tvtr = (bnAdN1*ulresL + bnAdN2*ulresR); %average pressure
            jumpflux = (bnAdNL*ulresL - bnAdNR*ulresR); %flux jump
            jump = NmatR*ulresR - NmatL*ulresL;        %v-p jump
%             jumpv = jump(1:2);
%             jumpp = jump(3);

            % Not needed for linear problems
            ElemFL = ElemFL - c1*( - NmatL(1:2,:)'*(tvtr + Kinv(1:2,1:2)*jump(1:2)) + bnAdN1w'*(jump(1:2)) - bnAdNLw'*(Kinvp*jumpflux));
            ElemFR = ElemFR - c1*( + NmatR(1:2,:)'*(tvtr + Kinv(1:2,1:2)*jump(1:2)) + bnAdN2w'*(jump(1:2)) + bnAdNRw'*(Kinvp*jumpflux));
        
            % Signs set to agree with L_Elem3_2d, original implementation
            ElemKLL = ElemKLL - c1*NmatL(1:2,:)'*bnAdN1 - c1*bnAdN1w'*NmatL(1:2,:);
            ElemKLR = ElemKLR - c1*NmatL(1:2,:)'*bnAdN2 + c1*bnAdN1w'*NmatR(1:2,:);
            ElemKRL = ElemKRL + c1*NmatR(1:2,:)'*bnAdN1 - c1*bnAdN2w'*NmatL(1:2,:);
            ElemKRR = ElemKRR + c1*NmatR(1:2,:)'*bnAdN2 + c1*bnAdN2w'*NmatR(1:2,:);
        
            ElemKLL = ElemKLL + c1*N1'*Kinv*N1 - c1*bnAdNLw'*Kinvp*bnAdNL;
            ElemKLR = ElemKLR - c1*N1'*Kinv*N2 + c1*bnAdNLw'*Kinvp*bnAdNR;
            ElemKRL = ElemKRL - c1*N2'*Kinv*N1 + c1*bnAdNRw'*Kinvp*bnAdNL;
            ElemKRR = ElemKRR + c1*N2'*Kinv*N2 - c1*bnAdNRw'*Kinvp*bnAdNR;
            
            if MatTypeTable(2,maL) == 3 % Beavers-Joseph-Saffman law
            ElemKLL = ElemKLL + c1*2*N1'*alp*muL/sqrt(kappaR)*[[(eye(2)-nmat) zeros(2,1)]; zeros(1,3)]*N1;
            ElemFL = ElemFL - c1*2*N1'*alp*muL/sqrt(kappaR)*[[(eye(2)-nmat) zeros(2,1)]; zeros(1,3)]*N1*ulresL;
            end
            
            if MatTypeTable(2,maR) == 3 % Beavers-Joseph-Saffman law
            ElemKRR = ElemKRR + c1*2*N2'*alp*muR/sqrt(kappaL)*[[(eye(2)-nmat) zeros(2,1)]; zeros(1,3)]*N2;
            ElemFR = ElemFR - c1*2*N2'*alp*muR/sqrt(kappaL)*[[(eye(2)-nmat) zeros(2,1)]; zeros(1,3)]*N2*ulresR;
            end
        
            end %numSI

        end %ie
        
        ElemKLL;
            ElemK = [ElemKLL ElemKLR; ElemKRL ElemKRR];
            ElemF = [ElemFL; ElemFR];
%%
    case 11
        
        ElemE = zeros(numEn,1);

end %Task Switch
