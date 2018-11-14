% 07/27/2014
% Tim Truster
%
% Mixed dimensional coupling element by Nitsche method for axial loads.
% See Rabinovich et al. - 2014 - The Nitsche method applied to a class .pdf

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

pencoeff = 10;00;

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
%%
    case {3,6} %interface stiffness
        
        % For now, assume bar on right, continuum on left
        
        ElemKRR = zeros(nstR,nstR);
        ElemKRL = zeros(nstR,nstL);
        ElemKLR = zeros(nstL,nstR);
        ElemKLL = zeros(nstL,nstL);
        ElemFR = zeros(nstR,1);
        ElemFL = zeros(nstL,1);
        t1 = zeros(3,1);
        t2 = t1;
        t3 = t1;
        
        ElemEL = matepropL(1);
        ElemvL = matepropL(2);
        thick = matepropL(3);
        lamdaL = ElemvL*ElemEL/((1+ElemvL)*(1-2*ElemvL));
        muL = ElemEL/(2*(1+ElemvL));
        if PSPS == 'n'
        DmatL = muL*diag([2 2 1]) + lamdaL*[1; 1; 0]*[1 1 0];
        else
        DmatL = ElemEL/(1-ElemvL^2)*[1      ElemvL  0
                                  ElemvL  1      0
                                  0      0      (1-ElemvL)/2];
        end
        
        PatchER = matepropR(1);
        PatchAR = matepropR(2);
        DmatR = PatchER; % Area is taken care of through integration (Jdet*thick)
        
        NmatL = zeros(2,nstL);
        BmatL = zeros(3,nstL);
        bnAdN1 = zeros(3,nstL);
        N1 = zeros(1,nstL);
        NmatR = zeros(1,nstR);
        BmatR = zeros(1,nstR);
        bnAdN2 = zeros(1,nstR);
        N2 = zeros(1,nstR);
        
        ulresL = reshape(ulL,ndf*nelL,1);
        ulresR = reshape(ulR,ndf*nelR,1);
        
%         % Determin bounds of integration segment
%         InterConn2D2 % InterConn2DT % 
% 
%         % Set jacobian for integration space
        drdr = 1;%(eL2 - eL1)/drL;
%         
%         m = (eR2-eR1)/(eL1-eL2);

        gamL = 0;1/2;
        gamR = 1;1/2;
        ep = pencoeff;
        if nelL == 3 || nelL == 4
            lint = 2;
        else
            lint = 3;
        end
        der = 0;
        bf = 0;
            
            [shlR,shld,shls,be] = shl1d(-1,1,nelR-1);
            [shgR, shgsR] = shg1d(xlR,ndm,nelR,shld(1,:),shls(1,:),nen,bf,der,be);
                
            NmatR(1:ndf:ndf*(nelR-1)+1) = shlR(1,:);
                
            BmatR(1:ndf:ndf*(nelR-1)+1) = shgR';
        
        for ie = 1:lint
            
            if nelL == 3 || nelL == 6
                [Wgt,r,s] = intpntt(ie,lint,1);
            elseif nelL == 4 || nelL == 9
                [Wgt,r,s] = intpntq(ie,lint,1);
            end
            
%             r = drdr*(r-roL)+eL1;
            
            if nelL == 3 || nelL == 6
                [shlL,shld,shls,be] = shlt(r,s,nelL,nelL,der,bf);
                [QxyL, shgs, Jdet, be, xsL] = shgt(xlL(:,1:nelL),nelL,shld,shls,nen,bf,der,be);
            elseif nelL == 4 || nelL == 9
                [shlL,shld,shls,be] = shlq(r,s,nelL,nelL,der,bf);
                [QxyL, shgs, Jdet, be, xsL] = shgq(xlL(:,1:nelL),nelL,shld,shls,nen,bf,der,be);
            end
            
%             rR = m*(r-eL2) + eR1;
            
            if nelR == 3 || nelR == 6
                s = 0;
            else %if nelR == 4
                s = -1;
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
            
            bnAdN1 = nvec'*(gamL*nvect*DmatL*BmatL);
            bnAdN2 = gamR*DmatR*BmatR;
            N1 = nvec'*NmatL;
            N2 = NmatR;
        
            xint = xlL*shlL;
            
            tvtr = (bnAdN1*ulresL + bnAdN2*ulresR); %average stress
            jumpu = N2*ulresR - N1*ulresL;        %displacement jump

            % Not needed for linear problems
            ElemFL = ElemFL - c1*( - N1'*(tvtr + ep*jumpu) + bnAdN1'*jumpu);
            ElemFR = ElemFR - c1*( + N2'*(tvtr + ep*jumpu) + bnAdN2'*jumpu);

            ElemKLL = ElemKLL - c1*N1'*bnAdN1;
            ElemKLR = ElemKLR - c1*N1'*bnAdN2;
            ElemKRL = ElemKRL + c1*N2'*bnAdN1;
            ElemKRR = ElemKRR + c1*N2'*bnAdN2;

            ElemKLL = ElemKLL - c1*bnAdN1'*N1;
            ElemKLR = ElemKLR + c1*bnAdN1'*N2;
            ElemKRL = ElemKRL - c1*bnAdN2'*N1;
            ElemKRR = ElemKRR + c1*bnAdN2'*N2;

            ElemKLL = ElemKLL + c1*(N1'*ep*N1);
            ElemKLR = ElemKLR - c1*(N1'*ep*N2);
            ElemKRL = ElemKRL - c1*(N2'*ep*N1);
            ElemKRR = ElemKRR + c1*(N2'*ep*N2);
            
        end
% ElemKLL
            ElemK = [ElemKLL ElemKLR; ElemKRL ElemKRR];
            ElemF = [ElemFL; ElemFR];
%%
    case 11
        
        ElemE = zeros(numEn,1);
        
end %Task Switch
