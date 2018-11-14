
if isw ~= 1
CGtoDGarrays
nelLP = nelL;
nelRP = nelR;

inter = elem - (numel - numSI);
nodeAR = SurfacesI(inter,1);
nodeBR = SurfacesI(inter,2);
nodeAL = SurfacesI(inter,3);
nodeBL = SurfacesI(inter,4);
end


% Tim Truster
% 08/10/2014
% DG debonding element for 2D problems
% Built from L_Elem1_2dDG and NL_Elem3_3d_bubble
% Uses scalar value for stability terms

nitvms = 1;2;
if nitvms == 1 %VMS
pencoeff = 1;
elseif nitvms == 2 %Nitsche
pencoeff = 1;4;2;
else %RFB
pencoeff = 1;4;2;
end

Bcol1 = [1; 3];
Bcol2 = [2; 3];
col1 = [1; 2];
col2 = [2; 1];

diagt = 1; % make taus diagonal
scalt = 1; % make taus scalar
minmaxt = 1;

switch isw %Task Switch
    
    case 1
        
        if ndf > 2
            
            for i = 3:ndf
                lie(i,1) = 0;
            end
            
        end
        
        nh1 = 7*4*1; % 8 variables, 4 int points, 1
        nh3 = 1*4*1; % dstate variable
%%
    case {3,6} %interface stiffness
        
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
        
        dtol = 1e-11;
        itchat = 7;100;6;12;8; %flag for iteration number to suppress chatter
        
        bf1 = 0;
        bf2 = 0;
        
        sigmax = mateprop(5);
        dc = mateprop(6);
        beta = mateprop(7);
        
        Hc = sigmax/dc;
        
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
%         DmatR = muR*diag([2 2 1]) + 0*[1; 1; 0]*[1 1 0];
%         DmatL = muL*diag([2 2 1]) + 0*[1; 1; 0]*[1 1 0];
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
        
        if inter == 1
            numdbond = 0;
        end
        
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

            % Modifications to tau
            if exist('diagt','var') && diagt == 1 % make tauL and tauR diagonal
                if exist('scalt','var') && scalt == 1 % use a scalar of tau times the identity
                    if exist('minmaxt','var') && minmaxt == -1 % minimum entry
                        tau = diag(inv(tauL));
                        tau = min(tau);%*ones(ndm,1);
                        tauL = inv(diag(tau));
                        tau = diag(inv(tauR));
                        tau = min(tau);%*ones(ndm,1);
                        tauR = inv(diag(tau));
                    else %if minmaxt == 1 % maximum entry
                        tau = diag(inv(tauL));
                        tau = max(tau);%*ones(ndm,1);
                        tauL = inv(diag(tau));
                        tau = diag(inv(tauR));
                        tau = max(tau);%*ones(ndm,1);
                        tauR = inv(diag(tau));
                    end
                else % diagonal only
                    tauL = inv(diag(diag(inv(tauL))));
                    tauR = inv(diag(diag(inv(tauR))));
                end
            end

            if exist('equat','var') && equat == 1 % make tauL = tauR
                if exist('scalt','var') && scalt == 1 % use a scalar of tau times the identity
                    if exist('minmaxt','var') && minmaxt == -1 % minimum entry
                        tau = diag(inv(tauL));
                        tau = min(tau);%*ones(ndm,1);
                        tauL = inv(diag(tau));
                        tau = diag(inv(tauR));
                        tau = min(tau);%*ones(ndm,1);
                        tauR = inv(diag(tau));
                    else %if minmaxt == 1 % maximum entry
                        tau = diag(inv(tauL));
                        tau = max(tau);%*ones(ndm,1);
                        tauL = inv(diag(tau));
                        tau = diag(inv(tauR));
                        tau = max(tau);%*ones(ndm,1);
                        tauR = inv(diag(tau));
                    end
                end
    %             else % equal only
                    tau = inv(tauL + tauR);
                    tauL = 1/2*(tauL*tau*tauR + tauR*tau*tauL);
                    tauR = tauL;
    %             end
            end
        
        maxnel = max(nelL,nelR);
        if maxnel == 3 || maxnel == 4
            lint = 2;
        else
            lint = 3;
        end
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
                    
            if nelR == 3 || nelR == 6
                [WgtR,rR,sR] = intpntt(ie,lint,1);
                ebeR = edgebubble(rR,sR,nelR);
            elseif nelR == 4 || nelR == 9
                [WgtR,rR,sR] = intpntq(ie,lint,1);
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
        ep = pencoeff*(eye(2)*max(muL,muR)/h + (nvect'*nvect)*max(lamdaL,lamdaR)/h); %Larson CMAME version
        %1000;pencoeff*eye(2)*max(muL,muR)/h;
        else
        % RFB
        tauL = RFBEfunc(xlintL,muL,lamdaL,8);
        tauR = RFBEfunc(xlintR,muR,lamdaR,8);
        edgeK = tauL + tauR;
        gamL = (edgeK\tauL);
        gamR = (edgeK\tauR);
        ep = pencoeff*intedge*inv(edgeK); %#ok<MINV>
        end
        rp = ep;
        
        for ll = 1:lint
            
            if nelL == 3 || nelL == 6
                [Wgt,r,s] = intpntt(ll,lint,1);
            elseif nelL == 4 || nelL == 9
                [Wgt,r,s] = intpntq(ll,lint,1);
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
                NmatL(1,(i-1)*ndf+1) = shlL(i);
                NmatL(2,(i-1)*ndf+2) = shlL(i);
                BmatL(Bcol1,(i-1)*ndf+1) = QxyL(i,col1);
                BmatL(Bcol2,(i-1)*ndf+2) = QxyL(i,col2);
            end
            
            for i = 1:nelR
                NmatR(1,(i-1)*ndf+1) = shlR(i);
                NmatR(2,(i-1)*ndf+2) = shlR(i);
                BmatR(Bcol1,(i-1)*ndf+1) = QxyR(i,col1);
                BmatR(Bcol2,(i-1)*ndf+2) = QxyR(i,col2);
            end
            
            bnAdN1 = gamL*nvect*DmatL*BmatL;
            bnAdN2 = gamR*nvect*DmatR*BmatR;
        
            xint = xlL*shlL;
            
            % Load history
            damhr = nh1-1+(ll-1)*7;
            dmaxhr = nh1-1+(ll-1)*7+4;
            chat1hr = nh2-1+(ll-1)*7+6;
            chat2hr = nh2-1+(ll-1)*7+7;
            dvec = [hr(damhr+1); hr(damhr+2)];
            dmax = hr(dmaxhr);
            bonddbond1 = hr(chat1hr);
            chatter = hr(chat2hr);
            dstate = hr(nh3-1+(ll-1)*1+1);
            
            tvtr = (bnAdN1*ulresL + bnAdN2*ulresR); %average stress
            jumpu = NmatR*ulresR - NmatL*ulresL;    %displacement jump          
            
            
            % Evaluate debonding interface model
            elast_init = transient == 0 && initia;
            [dvec,dddt,dmax,dstate,bonddbond2,dbond,chatter] = RadialCZM(tvtr,jumpu,dmax, ...
                dvec,nvec,sigmax,Hc,dc,rp,dstate,iter,elast_init,itchat,chatter,bonddbond1,dtol,ndm);

            
            % Force terms
            ElemFL = ElemFL - c1*( - NmatL'*(tvtr + rp*(jumpu - dvec)) + bnAdN1'*(jumpu - dvec));
            ElemFR = ElemFR - c1*( + NmatR'*(tvtr + rp*(jumpu - dvec)) + bnAdN2'*(jumpu - dvec));

            % Nitsche stiffness terms
            ElemKLL = ElemKLL - c1*NmatL'*bnAdN1 - c1*bnAdN1'*NmatL;
            ElemKLR = ElemKLR - c1*NmatL'*bnAdN2 + c1*bnAdN1'*NmatR;
            ElemKRL = ElemKRL + c1*NmatR'*bnAdN1 - c1*bnAdN2'*NmatL;
            ElemKRR = ElemKRR + c1*NmatR'*bnAdN2 + c1*bnAdN2'*NmatR;

            % Penalty terms
            ElemKLL = ElemKLL + c1*rp*(NmatL'*NmatL);
            ElemKLR = ElemKLR - c1*rp*(NmatL'*NmatR);
            ElemKRL = ElemKRL - c1*rp*(NmatR'*NmatL);
            ElemKRR = ElemKRR + c1*rp*(NmatR'*NmatR);

            % Damage terms
            ElemKLL = ElemKLL - c1*(bnAdN1'-rp*NmatL')*dddt*(bnAdN1-rp*NmatL);
            ElemKLR = ElemKLR - c1*(bnAdN1'-rp*NmatL')*dddt*(bnAdN2+rp*NmatR);
            ElemKRL = ElemKRL - c1*(bnAdN2'+rp*NmatR')*dddt*(bnAdN1-rp*NmatL);
            ElemKRR = ElemKRR - c1*(bnAdN2'+rp*NmatR')*dddt*(bnAdN2+rp*NmatR);
                    
            
            % Store history
            damhr = nh2-1+(ll-1)*7;
            dmaxhr = nh2-1+(ll-1)*7+4;
            dbonhr = nh2-1+(ll-1)*7+5;
            chat1hr = nh2-1+(ll-1)*7+6;
            chat2hr = nh2-1+(ll-1)*7+7;
            hr(damhr+1) = dvec(1);
            hr(damhr+2) = dvec(2);
%             hr(damhr+3) = dvec(3);
            hr(dmaxhr) = dmax;
            hr(dbonhr) = dbond;
            hr(chat1hr) = bonddbond2;
            hr(chat2hr) = chatter;
            hr(nh3-1+(ll-1)*1+1) = dstate;
            numdbond = numdbond + abs(dbond);
            
        end %lint
        
        % Print number of debonded integration points
        if inter == numSI
%             numdbond
        end
            
% ElemKLL
            ElemK = [ElemKLL ElemKLR; ElemKRL ElemKRR];
            ElemF = [ElemFL; ElemFR];
%%
    case 11
        
        ElemE = zeros(numEn,1);
        
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
        
end %Task Switch
