% MRDG PBC subroutine 
% Elina Geut 
% Created 5/23/2019 
% Last Modified 5/23/2019

if isw ~= 1
    if exist('flag_PBC','var') && flag_PBC ==1
        CGtoDGPBCarrays
    else
        CGtoDGarrays
    end
nelLP = nelL;
nelRP = nelR;

inter = elem - (numel - numSI);
nodeAR = -1;%SurfacesI(inter,1);
nodeBR = -1;%SurfacesI(inter,2);
nodeAL = -1;%SurfacesI(inter,3);
nodeBL = -1;%SurfacesI(inter,4);
end

dx =  xlR(1,2)- xlL(1,1);
dy =  xlR(2,2) - xlL(2,1);

% Tim Truster
% Stabilized DG for elasticity problem, sectors are still portions of the
% element, since the Tau functions were already general enough.
% 05/16/15

nitvms = 1;2;
if nitvms == 1 %VMS
pencoeff = 4;1;
elseif nitvms == 2 %Nitsche
pencoeff = 1;4;2;
else %RFB
pencoeff = 1;4;2;
end

Bcol1 = [1; 3];
Bcol2 = [2; 3];
col1 = [1; 2];
col2 = [2; 1];

switch isw %Task Switch
    
    case 1
        
        if ndf > 2
            
            for i = 3:ndf
                lie(i,1) = 0;
            end
            
        end
        iste = 20; % number of stresses per element
%%
    case {3,6,21} %interface stiffness
        
        nelLg = (ndm+1);
        nelRg = (ndm+1);
        nstLg = ndf*nelLg;
        nstRg = ndf*nelRg;
        
        ElemK = zeros(nst,nst);
        ElemF = zeros(nst,1);
        ElemKRR = zeros(nstR,nstR);
        ElemKRL = zeros(nstR,nstL);
        ElemKLR = zeros(nstL,nstR);
        ElemKLL = zeros(nstL,nstL);
        ElemFR = zeros(nstR,1);
        ElemFL = zeros(nstL,1);
       
        ElemKRRg = zeros(nstR,ndf*(ndm+1));
        ElemKRLg = zeros(nstR,ndf*(ndm+1));
        ElemKLRg = zeros(nstL,ndf*(ndm+1));
        ElemKLLg = zeros(nstL,ndf*(ndm+1));
        ElemKRgR = zeros(ndf*(ndm+1),nstR);
        ElemKRgL = zeros(ndf*(ndm+1),nstL);
        ElemKLgR = zeros(ndf*(ndm+1),nstR);
        ElemKLgL = zeros(ndf*(ndm+1),nstL);
        ElemKRgRg = zeros(ndf*(ndm+1),ndf*(ndm+1));
        ElemKRgLg = zeros(ndf*(ndm+1),ndf*(ndm+1));
        ElemKLgRg = zeros(ndf*(ndm+1),ndf*(ndm+1));
        ElemKLgLg = zeros(ndf*(ndm+1),ndf*(ndm+1));
        ElemFRg = zeros(ndf*(ndm+1),1);
        ElemFLg = zeros(ndf*(ndm+1),1);
        t1 = zeros(3,1);
        t2 = t1;
        t3 = t1;
        thick = 1;
        
        xBg_mid = BoundXmid(:,ma-nummatCG); % midpoint of grain boundary; formed by InterIntegV.m
        
        % Flags for using mesoscale or microscale; these are used only in
        % the stiffness and force vectors; 0 turns off, 1 turns on
        Lmeso = mateprop(5);
        Rmeso = mateprop(6);
        Lmicr = mateprop(7);
        Rmicr = mateprop(8);
        
        %Modification of 3/21/19 (EG)
         if all(ElemFlagL==ElemFlagR)
            DiricBC = 1; % Need to flip normal vector t3 for Dirichlet edges
        else
            DiricBC = 0;
            interI = ma-MaterTypeNum(2)+1;
        end
        %End of modification of 3/21/19 (EG)
        
        ElemFlagLg = ElemFlag(nelL+nelR+1:nelL+nelR+ndm+1);
        ElemFlagRg = ElemFlag(nelL+nelR+ndm+1+1:nelL+nelR+ndm+1+ndm+1);
        xlLg = xl(1:ndm,nelL+nelR+1:nelL+nelR+ndm+1);
        xlRg = xl(1:ndm,nelL+nelR+ndm+1+1:nelL+nelR+ndm+1+ndm+1);
        ulLg = ul(1:ndf,nelL+nelR+1:nelL+nelR+ndm+1);
        ulRg = ul(1:ndf,nelL+nelR+ndm+1+1:nelL+nelR+ndm+1+ndm+1);
        
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
        
        ulresLg = reshape(ulLg,ndf*nelLg,1);
        ulresRg = reshape(ulRg,ndf*nelRg,1);
        ulresM = reshape(ulM,ndf*ndm,1);
        xlresM = reshape(xlM,ndf*ndm,1);
        
        if nel == 3
            lint = lintt3;%13;
        elseif nel == 4
            lint = lintq4;
        elseif nel == 6
            lint = lintt6;%13;
        elseif nel == 9
            lint = lintq9;
        end
        
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
            elem
            xint = xl(:,1:4)*shl;
            ElemE(1:2) = ElemE(1:2) + c1*xint;
            
        end %je
        
%         zeta = ElemE;
%         zetaX = [dx
%                  dy];
%         NmatM = zeros(ndm,nstM);  % Need to correct this. currently giving wrong results. 

        
        [zetaX,zeta,NmatM,dx,dy] = compute_zeta3(xlR,xlL,ulresM,xlresM);
        elem
        %            [zetaX,zeta, NmatM] = compute_zeta(xlR,xlL,ulresM);
        %             [zetaX1,zetaX2,zeta, NmatM] = compute_zeta1(xlR,xlL,ulresM);
%         [zetaX,zeta, NmatM,dx,dy] = compute_zeta3(xlR,xlL,ulresM,xlresM);

        InterConn2D2PBC
        
        % Set jacobian for integration space
        drdr = (eL2 - eL1)/drL;
        
        m = (eR2-eR1)/(eL1-eL2);
        
        % Perform integration of various matrices
        
        % For separate bubble types on T and Q
        [tauL,intbL] = TauE1_2d(xlintL,DmatL,nelL,lintt6,lintq9);
        [tauR,intbR] = TauE1_2d(xlintR,DmatR,nelR,lintt6,lintq9);
        
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
                NmatL(1,(i-1)*ndf+1) = shlL(i);
                NmatL(2,(i-1)*ndf+2) = shlL(i);
            end
            
            for i = 1:nelR
                NmatR(1,(i-1)*ndf+1) = shlR(i);
                NmatR(2,(i-1)*ndf+2) = shlR(i);
            end
        
            xint = xlL*shlL;
            xeL = xint - xlLg(:,1);
            xvectL = [xeL(1) 0 xeL(2)/2
                     0 xeL(2) xeL(1)/2];
            xeR = xint - xlRg(:,1);
            xvectR = [xeR(1) 0 xeR(2)/2
                     0 xeR(2) xeR(1)/2];
            
            trjump = reshape(nvect*DmatR - nvect*DmatL,6,1);
            trjumpL = reshape(inv(DmatL),9,1);
            trjumpR = reshape(inv(DmatR),9,1);
            
            ElemE(1:3+9+9) = ElemE(1:3+9+9) + c1*[xint; 0; trjumpL; trjumpR];
            ElemS(1:7) = ElemS(1:7) + c1*[1; trjump];
            
          end
         
         %Part from TaylorSet (EG Modification)
             if ndm == 2
                 trjump = reshape(trjump(1:6),2,3);
             else
                 error('nmd=3')
             end
             Load = NodeLoad'; %Modified of 5/25/19
%              Load = [NodeLoad(1,3); NodeLoad(4,3); NodeLoad(2,3)];
             fluxT = factorT*trjump*Load;
             %Part from TaylorSet end
             
             %Part from Sachs set
             jumpSu = zeros(2,1);
             jumpSe = zeros(3,1);

            
             
    %End part from Sachs model
    
         for ie = 1:lint
         if Lmeso == 1 || Rmeso == 1
             Rmeso;
        end
        if all(ElemFlagL==ElemFlagR)
            DiricBC = 1; % Need to flip normal vector t3 for Dirichlet edges
            fluxjc = zeros(2,1);
            fluxjx = zeros(2,1);
            fluxjy = zeros(2,1);
            jumpic = zeros(2,1);
            jumpix = zeros(2,1);
            jumpiy = zeros(2,1);
            jumpiu = zeros(2,1);
            jumpie = zeros(3,1); 
        else
            DiricBC = 0;
            fluxjc = fluxT;
            fluxjx = [0;0];
            fluxjy = [0;0];
%             jumpic = zeros(2,1);
%             jumpix = zeros(2,1);
%             jumpiy = zeros(2,1);
            jumpiu = jumpSu;
            jumpie = jumpSe; [jumpSe(3) 0];
%             jumpic = CoordinatesI(2*(ndm+1)*(interI-1)+4,:)';
%             jumpix = CoordinatesI(2*(ndm+1)*(interI-1)+5,:)';
%             jumpiy = CoordinatesI(2*(ndm+1)*(interI-1)+6,:)';
        end
            

     %End of modified part 3/21/19 EG
         end
        
     if Rmicr + Lmicr > 0
        if nitvms == 1
            if maL == maR % Dirichlet BC assignment
                % VMS
                edgeK = tauL*ebL^2;
                gamL = eye(2);
                gamR = zeros(2,2);
                ep = pencoeff*intedge*inv(edgeK); %#ok<MINV>
            else
                % VMS
                edgeK = tauL*ebL^2 + tauR*ebR^2;
                gamL = ebL^2*(edgeK\tauL);
                gamR = ebR^2*(edgeK\tauR);
                ep = pencoeff*intedge*inv(edgeK); %#ok<MINV>
            end
        elseif nitvms == 2
        % Nitsche
        volL = getvol(xlL,nelL);
        volR = getvol(xlR,nelR);
        h = 2/(intedge/volR + intedge/volL);
        gamL = 1/2*eye(2);
        gamR = 1/2*eye(2);
        ep = pencoeff*(eye(2)*max(muL,muR)/h + (nvect'*nvect)*max(lamdaL,lamdaR)/h); %Larson CMAME version
        else
        % RFB
        tauL = RFBEfunc(xlintL,muL,lamdaL,8);
        tauR = RFBEfunc(xlintR,muR,lamdaR,8);
        edgeK = tauL + tauR;
        gamL = (edgeK\tauL);
        gamR = (edgeK\tauR);
        ep = pencoeff*intedge*inv(edgeK); %#ok<MINV>
        end
        else
        % DECIDE what to do for meso-meso case
        edgeK = tauL*ebL^2 + tauR*ebR^2;
        gamL = ebL^2*(edgeK\tauL);
        gamR = ebR^2*(edgeK\tauR);
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
            if DiricBC
                t3 = -t3;
            end
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
            % Micro quantities
            bnAdN1 = gamL*nvect*DmatL*BmatL;
            bnAdN2 = gamR*nvect*DmatR*BmatR;
            xint = xlL*shlL; % coordinates of current integration point on edge
            xe = xint - xBg_mid; % relative coordinate to the grain edge center
            xeL = xint - xlLg(:,1); % relative coordinate to the L grain; set within GrainIntegV
            xeR = xint - xlRg(:,1); % relative coordinate to the R grain; set within GrainIntegV
            xvectL = [xeL(1) 0 xeL(2)/2
                     0 xeL(2) xeL(1)/2];
            xvectR = [xeR(1) 0 xeR(2)/2
                     0 xeR(2) xeR(1)/2];
            
%             zeta = NmatM*ulresM; %Added for PBC 4/23/19
            tvtr = (bnAdN1*ulresL + bnAdN2*ulresR); %average stress
            jumpu = NmatR*ulresR - NmatL*ulresL;        %displacement jump
            jumpi = jumpiu + (xvectR*inv(DmatR) - xvectL*inv(DmatL))*jumpie;
%             jumpi = jumpic + xe(1)*jumpix + xe(2)*jumpiy;
            fluxj = fluxjc + xe(1)*fluxjx + xe(2)*fluxjy;
%             
%             tvtr = (bnAdN1*ulresL + bnAdN2*ulresR); %average stress
%             jumpu = NmatR*ulresR - NmatL*ulresL;        %displacement jump

             
            ElemKLL = ElemKLL - c1*NmatL'*bnAdN1;
            ElemKLR = ElemKLR - c1*NmatL'*bnAdN2;
            ElemKRL = ElemKRL + c1*NmatR'*bnAdN1;
            ElemKRR = ElemKRR + c1*NmatR'*bnAdN2;

            ElemKLL = ElemKLL - c1*bnAdN1'*NmatL;
            ElemKLR = ElemKLR + c1*bnAdN1'*NmatR;
            ElemKRL = ElemKRL - c1*bnAdN2'*NmatL;
            ElemKRR = ElemKRR + c1*bnAdN2'*NmatR;

            ElemKLL = ElemKLL + c1*(NmatL'*ep*NmatL);
            ElemKLR = ElemKLR - c1*(NmatL'*ep*NmatR);
            ElemKRL = ElemKRL - c1*(NmatR'*ep*NmatL);
            ElemKRR = ElemKRR + c1*(NmatR'*ep*NmatR);
            
            % add interface prescribed flux
            ElemFL = ElemFL - c1*( + NmatL'*gamR'*fluxj);
            ElemFR = ElemFR - c1*( + NmatR'*gamL'*fluxj);
            ElemFL = ElemFL + c1*( - NmatL'*ep*jumpi);
            ElemFR = ElemFR + c1*( + NmatR'*ep*jumpi);
            ElemFL = ElemFL + c1*( + bnAdN1'*jumpi);
            ElemFR = ElemFR + c1*( + bnAdN2'*jumpi);
            
            % Meso quantities
            BmatLg = [0 0 1 0 0 0
                      0 0 0 1 0 0
                      0 0 0 0 1 0];
            BmatRg = BmatLg;
            bnAdN1g = gamL*nvect*DmatL*BmatLg;
            bnAdN2g = gamR*nvect*DmatR*BmatRg;
            NmatLg = [1 0 xeL(1) 0 0.5*xeL(2) 0.5*xeL(2)
                      0 1 0 xeL(2) 0.5*xeL(1) -0.5*xeL(1)];
            NmatRg = [1 0 xeR(1) 0 0.5*xeR(2) 0.5*xeR(2)
                      0 1 0 xeR(2) 0.5*xeR(1) -0.5*xeR(1)];
            tvtrg = (bnAdN1g*ulresLg + bnAdN2g*ulresRg); %average stress
            jumpug = NmatRg*ulresRg - NmatLg*ulresLg;        %displacement jump
            
            ElemKLgLg = ElemKLgLg - c1*NmatLg'*bnAdN1g;
            ElemKLgRg = ElemKLgRg - c1*NmatLg'*bnAdN2g;
            ElemKRgLg = ElemKRgLg + c1*NmatRg'*bnAdN1g;
            ElemKRgRg = ElemKRgRg + c1*NmatRg'*bnAdN2g;

            ElemKLgLg = ElemKLgLg - c1*bnAdN1g'*NmatLg;
            ElemKLgRg = ElemKLgRg + c1*bnAdN1g'*NmatRg;
            ElemKRgLg = ElemKRgLg - c1*bnAdN2g'*NmatLg;
            ElemKRgRg = ElemKRgRg + c1*bnAdN2g'*NmatRg;

            ElemKLgLg = ElemKLgLg + c1*(NmatLg'*ep*NmatLg);
            ElemKLgRg = ElemKLgRg - c1*(NmatLg'*ep*NmatRg);
            ElemKRgLg = ElemKRgLg - c1*(NmatRg'*ep*NmatLg);
            ElemKRgRg = ElemKRgRg + c1*(NmatRg'*ep*NmatRg);
            
            % Meso-micro quantities

            ElemKLgL = ElemKLgL - c1*NmatLg'*bnAdN1;
            ElemKLgR = ElemKLgR - c1*NmatLg'*bnAdN2;
            ElemKRgL = ElemKRgL + c1*NmatRg'*bnAdN1;
            ElemKRgR = ElemKRgR + c1*NmatRg'*bnAdN2;

            ElemKLgL = ElemKLgL - c1*bnAdN1g'*NmatL;
            ElemKLgR = ElemKLgR + c1*bnAdN1g'*NmatR;
            ElemKRgL = ElemKRgL - c1*bnAdN2g'*NmatL;
            ElemKRgR = ElemKRgR + c1*bnAdN2g'*NmatR;

            ElemKLgL = ElemKLgL + c1*(NmatLg'*ep*NmatL);
            ElemKLgR = ElemKLgR - c1*(NmatLg'*ep*NmatR);
            ElemKRgL = ElemKRgL - c1*(NmatRg'*ep*NmatL);
            ElemKRgR = ElemKRgR + c1*(NmatRg'*ep*NmatR);
            
            % Micro-meso quantities

            ElemKLLg = ElemKLLg - c1*NmatL'*bnAdN1g;
            ElemKLRg = ElemKLRg - c1*NmatL'*bnAdN2g;
            ElemKRLg = ElemKRLg + c1*NmatR'*bnAdN1g;
            ElemKRRg = ElemKRRg + c1*NmatR'*bnAdN2g;

            ElemKLLg = ElemKLLg - c1*bnAdN1'*NmatLg;
            ElemKLRg = ElemKLRg + c1*bnAdN1'*NmatRg;
            ElemKRLg = ElemKRLg - c1*bnAdN2'*NmatLg;
            ElemKRRg = ElemKRRg + c1*bnAdN2'*NmatRg;

            ElemKLLg = ElemKLLg + c1*(NmatL'*ep*NmatLg);
            ElemKLRg = ElemKLRg - c1*(NmatL'*ep*NmatRg);
            ElemKRLg = ElemKRLg - c1*(NmatR'*ep*NmatLg);
            ElemKRRg = ElemKRRg + c1*(NmatR'*ep*NmatRg);
            
            %Forces
            ElemFLg = ElemFLg - c1*( + NmatLg'*gamR'*fluxj);
            ElemFRg = ElemFRg - c1*( + NmatRg'*gamL'*fluxj);
            ElemFLg = ElemFLg + c1*( - NmatLg'*ep*jumpi);
            ElemFRg = ElemFRg + c1*( + NmatRg'*ep*jumpi);
            ElemFLg = ElemFLg + c1*( + bnAdN1g'*jumpi);
            ElemFRg = ElemFRg + c1*( + bnAdN2g'*jumpi);
            
            %PBC Modification 4/16/2019
                %PBC  modification of 4/16/2019
                ElemFM = zeros(nstM,1);
                ElemKMR = zeros(nstM,nstR); %4x8
                ElemKML = zeros(nstM,nstL); %4x8
                ElemKRM = zeros(nstR,nstM); %8x4
                ElemKLM = zeros(nstL,nstM); %8x4
                ElemKMM = zeros(nstM,nstM); %4x4
                
                ElemKLMg = zeros(nstL,nstM); %8x4
                ElemKRMg = zeros(nstR,nstM); %8x4
                ElemKMLg = zeros(nstM,nstM+2);%4x6
                ElemKMRg = zeros(nstM,nstM+2);%4x6
                ElemKMMg = zeros(nstM,nstM);%4x4
                
                ElemKRgMg = zeros(nstM+2,nstM);%6x4
                ElemKLgMg = zeros(nstM+2,nstM);%6x4
                ElemKMgRg = zeros(nstM,nstM+2);%4x6
                ElemKMgLg = zeros(nstM,nstM+2);%4x6
                ElemKMgMg = zeros(nstM,nstM); %4x4
                
                ElemKLgM = zeros(nstM+2,nstM); %6x4
                ElemKRgM = zeros(nstM+2,nstM); %6x4
                ElemKMgL = zeros(nstM,nstL); %4x8
                ElemKMgR = zeros(nstM,nstR); %4x8
                ElemKMgM = zeros(nstM,nstM); %4x4
                
                NmatMg = zeros(ndm,nstM);
                
                TmatL = bnAdN1 + ep*NmatL*Lmicr;
                TmatR = bnAdN2 - ep*NmatR*Lmicr;
                TmatLg = bnAdN1g + ep*NmatLg*Lmeso;
                TmatRg = bnAdN2g - ep*NmatRg*Lmeso;
            %End of modification
                %Forces
                ElemFLg = ElemFLg - c1*( - NmatLg'*ep*jumpi - bnAdN1g'*jumpi + NmatLg'*gamR'*fluxj);
                ElemFRg = ElemFRg - c1*( + NmatRg'*ep*jumpi - bnAdN2g'*jumpi + NmatRg'*gamL'*fluxj);
%                 ElemFM = ElemFM - c1*( - NmatM'*(tvtr + ep*(jumpu-zeta)));
                
                %Modification of 5/15/19
                ElemKLM = ElemKLM - c1*TmatL'*NmatM;
                ElemKRM = ElemKRM - c1*TmatR'*NmatM;
                ElemKLgMg = ElemKLgMg - c1*TmatLg'*NmatMg;
                ElemKRgMg = ElemKRgMg - c1*TmatRg'*NmatMg;
                ElemKLMg = ElemKLMg - c1*TmatL'*NmatMg;
                ElemKRMg = ElemKRMg - c1*TmatR'*NmatMg;
                ElemKLgM = ElemKLgM - c1*TmatLg'*NmatM;
                ElemKRgM = ElemKRgM - c1*TmatRg'*NmatM;
                
                ElemKML = ElemKML - c1*NmatM'*TmatL;
                ElemKMR = ElemKMR - c1*NmatM'*TmatR;
                ElemKMgLg = ElemKMgLg - c1*NmatMg'*TmatLg;
                ElemKMgRg = ElemKMgRg - c1*NmatMg'*TmatRg;
                ElemKMLg = ElemKMLg - c1*NmatM'*TmatLg;
                ElemKMRg = ElemKMRg - c1*NmatM'*TmatRg;
                ElemKMgL = ElemKMgL - c1*NmatMg'*TmatL;
                ElemKMgR = ElemKMgR - c1*NmatMg'*TmatR;
                
                ElemKMM = ElemKMM + c1*NmatM'*ep*NmatM; 
                ElemKMMg = ElemKMMg + c1*NmatM'*ep*NmatMg;
                ElemKMgM = ElemKMgM + c1*NmatMg'*ep*NmatM;
                ElemKMgMg = ElemKMgMg + c1*NmatMg'*ep*NmatMg;
                %End of modification of 5/15/19

            %End of modification 4/16/2019
         end
         
%         ElemK = [Lmicr*Lmicr* ElemKLL Lmicr*Rmicr* ElemKLR ElemKLM Lmicr*Lmeso* ElemKLLg Lmicr*Rmeso* ElemKLRg ElemKLMg
%                 Rmicr*Lmicr* ElemKRL Rmicr*Rmicr* ElemKRR ElemKRM Rmicr*Lmeso* ElemKRLg Rmicr*Rmeso* ElemKRRg ElemKRMg
%                 ElemKML ElemKMR ElemKMM ElemKMLg ElemKMRg ElemKMMg
%                 Lmeso*Lmicr*ElemKLgL Lmeso*Rmicr*ElemKLgR ElemKLgM Lmeso*Lmeso*ElemKLgLg Lmeso*Rmeso*ElemKLgRg ElemKLgMg
%                 Rmeso*Lmicr*ElemKRgL Rmeso*Rmicr*ElemKRgR ElemKRgM Rmeso*Lmeso*ElemKRgLg Rmeso*Rmeso*ElemKRgRg ElemKRgMg
%                 ElemKMgL ElemKMgR ElemKMgM ElemKMgLg ElemKMgRg ElemKMgMg];
%         ElemK = [Lmicr*Lmicr* ElemKLL Lmicr*Rmicr* ElemKLR Lmicr*Lmeso* ElemKLLg Lmicr*Rmeso* ElemKLRg ElemKLM
%                 Rmicr*Lmicr* ElemKRL Rmicr*Rmicr* ElemKRR Rmicr*Lmeso* ElemKRLg Rmicr*Rmeso* ElemKRRg  ElemKRM
%                 Lmeso*Lmicr*ElemKLgL Lmeso*Rmicr*ElemKLgR Lmeso*Lmeso*ElemKLgLg Lmeso*Rmeso*ElemKLgRg ElemKLgMg
%                 Rmeso*Lmicr*ElemKRgL Rmeso*Rmicr*ElemKRgR Rmeso*Lmeso*ElemKRgLg Rmeso*Rmeso*ElemKRgRg ElemKRgMg
%                 ElemKML ElemKMR ElemKMgLg ElemKMgRg ElemKMMg];
%              

        ElemK;
        
        case 51 
        
        ElemS = zeros(20+6,1);
        ElemE = zeros(20,1);
        t1 = zeros(3,1);
        t2 = t1;
        t3 = t1;
        thick = 1;
        
        bf1 = 0;
        bf2 = 0;
        
        ElemFlagLg = ElemFlag(nelL+nelR+1:nelL+nelR+ndm+1);
        ElemFlagRg = ElemFlag(nelL+nelR+ndm+1+1:nelL+nelR+ndm+1+ndm+1);
        xlLg = xl(1:ndm,nelL+nelR+1:nelL+nelR+ndm+1);
        xlRg = xl(1:ndm,nelL+nelR+ndm+1+1:nelL+nelR+ndm+1+ndm+1);
        ulLg = ul(1:ndf,nelL+nelR+1:nelL+nelR+ndm+1);
        ulRg = ul(1:ndf,nelL+nelR+ndm+1+1:nelL+nelR+ndm+1+ndm+1);
        
        xBg_mid = BoundXmid(:,ma-nummatCG); % midpoint of grain boundary
        
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

        % Set jacobian for integration space
        drdr = (eL2 - eL1)/drL;
        
        m = (eR2-eR1)/(eL1-eL2);
        
        % Perform integration of various matrices
        
        % For separate bubble types on T and Q
        [tauL,intbL] = TauE1_2d(xlintL,DmatL,nelL,lintt6,lintq9);
        [tauR,intbR] = TauE1_2d(xlintR,DmatR,nelR,lintt6,lintq9);
        
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
        else
        % RFB
        tauL = RFBEfunc(xlintL,muL,lamdaL,8);
        tauR = RFBEfunc(xlintR,muR,lamdaR,8);
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
            xe = xint - xBg_mid;
            xeL = xint - xlLg(:,1);
            xeR = xint - xlRg(:,1);
            xvectL = [xeL(1) 0 xeL(2)/2
                     0 xeL(2) xeL(1)/2];
            xvectR = [xeR(1) 0 xeR(2)/2
                     0 xeR(2) xeR(1)/2];
            xvect = [xint(1) 0 xint(2)/2
                     0 xint(2) xint(1)/2];
            tracLR = nvect*sigLR./[volLR; volLR];
            epsV = epsLR./[volLR; volLR; volLR];
            jumpug = (disLR(1:2,2) + xvectR*epsV(1:3,2)) - (disLR(1:2,1) + xvectL*epsV(1:3,1));        %displacement jump
            
            tvtr = (bnAdN1*ulresL + bnAdN2*ulresR); %average stress
            jumpu = NmatR*ulresR - NmatL*ulresL;        %displacement jump
            
            straiLR = nvect*epsV;
            jumpsig = nvect*DmatR*BmatR*ulresR - nvect*DmatL*BmatL*ulresL;        %traction jump
            tracAve = gamL*tracLR(:,1) + gamR*tracLR(:,2);
            tracJump = tracLR(:,2) - tracLR(:,1);
            straiJump = straiLR(:,2) - straiLR(:,1);
            straiAve = gamL*straiLR(:,1) + gamR*straiLR(:,2);
            sL = DmatL*BmatL*ulresL; sR = DmatR*BmatR*ulresR;
            uL = NmatL*ulresL; uR = NmatR*ulresR;
%             HL1 = ((sL(1)-GrainSig(1,5)/9)*nLx+(sL(3)-GrainSig(3,5)/9)*nLy)*(uL(1) - GrainEps(1,5)/9*xint(1) - GrainEps(3,5)/9*xint(2));
%             HL2 = ((sL(3)-GrainSig(3,5)/9)*nLx+(sL(2)-GrainSig(2,5)/9)*nLy)*(uL(2) - GrainEps(3,5)/9*xint(1) - GrainEps(2,5)/9*xint(2));
%             HR1 = ((sR(1)-GrainSig(1,5)/9)*-nLx+(sR(3)-GrainSig(3,5)/9)*-nLy)*(uR(1) - GrainEps(1,5)/9*xint(1) - GrainEps(3,5)/9*xint(2));
%             HR2 = ((sR(3)-GrainSig(3,5)/9)*-nLx+(sR(2)-GrainSig(2,5)/9)*-nLy)*(uR(2) - GrainEps(3,5)/9*xint(1) - GrainEps(2,5)/9*xint(2));
            HL1 = ((sL(1)-GrainSig(1,maL)/9)*nLx+(sL(3)-GrainSig(3,maL)/9)*nLy)*(uL(1) - GrainEps(1,maL)/9*xint(1) - GrainEps(3,maL)/9*xint(2));
            HL2 = ((sL(3)-GrainSig(3,maL)/9)*nLx+(sL(2)-GrainSig(2,maL)/9)*nLy)*(uL(2) - GrainEps(3,maL)/9*xint(1) - GrainEps(2,maL)/9*xint(2));
            HR1 = ((sR(1)-GrainSig(1,maR)/9)*-nLx+(sR(3)-GrainSig(3,maR)/9)*-nLy)*(uR(1) - GrainEps(1,maR)/9*xint(1) - GrainEps(3,maR)/9*xint(2));
            HR2 = ((sR(3)-GrainSig(3,maR)/9)*-nLx+(sR(2)-GrainSig(2,maR)/9)*-nLy)*(uR(2) - GrainEps(3,maR)/9*xint(1) - GrainEps(2,maR)/9*xint(2));
            HillL = (nvect*DmatL*BmatL*ulresL - tracLR(:,1))'*(NmatL*ulresL - xvect*epsV(:,1));
            HillR = (-nvect*DmatR*BmatR*ulresR + tracLR(:,2))'*(NmatR*ulresR - xvect*epsV(:,2));
            HillL1 = (nvect*DmatL*BmatL*ulresL - tracLR(:,1));
            HillR1 = (-nvect*DmatR*BmatR*ulresR + tracLR(:,2));
            HillL2 = (NmatL*ulresL - xvect*epsV(:,1));
            HillR2 = (NmatR*ulresR - xvect*epsV(:,2));
            
            ElemS(1:18+8) = ElemS(1:18+8) + c1*[tvtr; (jumpu-jumpug); jumpsig; tracAve; tracJump; straiAve; straiJump; diag(gamL); diag(gamR);...
                xe(1)*(jumpu-jumpug); xe(2)*(jumpu-jumpug); xe(1)*(jumpsig-tracJump); xe(2)*(jumpsig-tracJump);];
%             ElemE(1:11) = ElemE(1:11) + c1*[1; HillL; HillR; HillL1; HillR1; HillL2; HillR2];
            ElemE(1:11) = ElemE(1:11) + c1*[1; HL1+HL2; HR1+HR2; HillL1; HillR1; HillL2; HillR2];
            
            
        end
        %%
        case 53 
        
        ElemS = zeros(7,1);
        ElemE = zeros(3+9+9,1);
        t1 = zeros(3,1);
        t2 = t1;
        t3 = t1;
        thick = 1;
        
        xlLg = xl(1:ndm,nelL+nelR+1:nelL+nelR+ndm+1);
        xlRg = xl(1:ndm,nelL+nelR+ndm+1+1:nelL+nelR+ndm+1+ndm+1);
        
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
        NmatR = zeros(2,nstR);
        
        % Determin bounds of integration segment
        InterConn2D2 % InterConn2DT % 

        % Set jacobian for integration space
        drdr = (eL2 - eL1)/drL;
        
        m = (eR2-eR1)/(eL1-eL2);
        
        % Perform integration of various matrices
        
        % For separate bubble types on T and Q
        [tauL,intbL] = TauE1_2d(xlintL,DmatL,nelL,lintt6,lintq9);
        [tauR,intbR] = TauE1_2d(xlintR,DmatR,nelR,lintt6,lintq9);
        
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
        else
        % RFB
        tauL = RFBEfunc(xlintL,muL,lamdaL,8);
        tauR = RFBEfunc(xlintR,muR,lamdaR,8);
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
                NmatL(1,(i-1)*ndf+1) = shlL(i);
                NmatL(2,(i-1)*ndf+2) = shlL(i);
            end
            
            for i = 1:nelR
                NmatR(1,(i-1)*ndf+1) = shlR(i);
                NmatR(2,(i-1)*ndf+2) = shlR(i);
            end
        
            xint = xlL*shlL;
            xeL = xint - xlLg(:,1);
            xvectL = [xeL(1) 0 xeL(2)/2
                     0 xeL(2) xeL(1)/2];
            xeR = xint - xlRg(:,1);
            xvectR = [xeR(1) 0 xeR(2)/2
                     0 xeR(2) xeR(1)/2];
            
            trjump = reshape(nvect*DmatR - nvect*DmatL,6,1);
            trjumpL = reshape(inv(DmatL),9,1);
            trjumpR = reshape(inv(DmatR),9,1);
            
            ElemE(1:3+9+9) = ElemE(1:3+9+9) + c1*[xint; 0; trjumpL; trjumpR];
            ElemS(1:7) = ElemS(1:7) + c1*[1; trjump];
            
        end

%%
    case 26 % Element stresses = fluxes
        
        ElemS = zeros(20,1);
        t1 = zeros(3,1);
        t2 = t1;
        t3 = t1;
        thick = 1;
        
        bf1 = 0;
        bf2 = 0;
        
        % Flags for using mesoscale or microscale
        Lmeso = mateprop(5);
        Rmeso = mateprop(6);
        Lmicr = mateprop(7);
        Rmicr = mateprop(8);
        if Lmeso == 1 || Rmeso == 1
            Rmeso;
        end
        if all(ElemFlagL==ElemFlagR)
            DiricBC = 1; % Need to flip normal vector for Dirichlet edges
            fluxjc = zeros(2,1);
            fluxjx = zeros(2,1);
            fluxjy = zeros(2,1);
            jumpic = zeros(2,1);
            jumpix = zeros(2,1);
            jumpiy = zeros(2,1);
        else
            DiricBC = 0;
            interI = ma-MaterTypeNum(2)+1;
            fluxjc = CoordinatesI(2*(ndm+1)*(interI-1)+1,:)';
            fluxjx = CoordinatesI(2*(ndm+1)*(interI-1)+2,:)';
            fluxjy = CoordinatesI(2*(ndm+1)*(interI-1)+3,:)';
            jumpic = CoordinatesI(2*(ndm+1)*(interI-1)+4,:)';
            jumpix = CoordinatesI(2*(ndm+1)*(interI-1)+5,:)';
            jumpiy = CoordinatesI(2*(ndm+1)*(interI-1)+6,:)';
        end
        % Grain stresses
        matLR = mateprop(1:2);
        sigLR = GrainSig(1:3,matLR);
        epsLR = GrainEps(1:3,matLR);
        volLR = GrainVol(1,matLR);
        
        ElemFlagLg = ElemFlag(nelL+nelR+1:nelL+nelR+ndm+1);
        ElemFlagRg = ElemFlag(nelL+nelR+ndm+1+1:nelL+nelR+ndm+1+ndm+1);
        xlLg = xl(1:ndm,nelL+nelR+1:nelL+nelR+ndm+1);
        xlRg = xl(1:ndm,nelL+nelR+ndm+1+1:nelL+nelR+ndm+1+ndm+1);
        ulLg = ul(1:ndf,nelL+nelR+1:nelL+nelR+ndm+1);
        ulRg = ul(1:ndf,nelL+nelR+ndm+1+1:nelL+nelR+ndm+1+ndm+1);
        
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

        % Set jacobian for integration space
        drdr = (eL2 - eL1)/drL;
        
        m = (eR2-eR1)/(eL1-eL2);
        
        % Perform integration of various matrices
        
        % For separate bubble types on T and Q
        [tauL,intbL] = TauE1_2d(xlintL,DmatL,nelL,lintt6,lintq9);
        [tauR,intbR] = TauE1_2d(xlintR,DmatR,nelR,lintt6,lintq9);
        
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
        else
        % RFB
        tauL = RFBEfunc(xlintL,muL,lamdaL,8);
        tauR = RFBEfunc(xlintR,muR,lamdaR,8);
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
            xvect = [xint(1) 0 xint(2)/2
                     0 xint(2) xint(1)/2];
            xint = xlL*shlL;
            xe = xint - xBg_mid;
            xeL = xint - xlLg(:,1);
            xeR = xint - xlRg(:,1);
            xvectL = [xeL(1) 0 xeL(2)/2
                     0 xeL(2) xeL(1)/2];
            xvectR = [xeR(1) 0 xeR(2)/2
                     0 xeR(2) xeR(1)/2];
            
            tvtr = (bnAdN1*ulresL + bnAdN2*ulresR); %average stress
            jumpu = NmatR*ulresR - NmatL*ulresL;        %displacement jump
            jumpi = jumpic + xe(1)*jumpix + xe(2)*jumpiy;
            fluxj = fluxjc + xe(1)*fluxjx + xe(2)*fluxjy;
            
            tracLR = nvect*sigLR./[volLR; volLR];
            epsV = epsLR./[volLR; volLR; volLR];
            straiLR = nvect*epsV;
            jumpsig = nvect*DmatR*BmatR*ulresR - nvect*DmatL*BmatL*ulresL;        %traction jump
            tracAve = gamL*tracLR(:,1) + gamR*tracLR(:,2);
            tracJump = tracLR(:,1) - tracLR(:,2);
            straiJump = straiLR(:,1) - straiLR(:,2);
            straiAve = gamL*straiLR(:,1) + gamR*straiLR(:,2);
            
            ElemS(1:18-4) = ElemS(1:18-4) + c1*[tvtr; jumpu; jumpsig; tracAve; tracJump; straiAve; straiJump];
            
        end
        ElemS(19-4:20) = [jumpic; fluxjc; nvec];
        
end %Task Switch
