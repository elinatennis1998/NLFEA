% 07/27/2014
% Tim Truster
%
% Mixed dimensional coupling element by Nitsche method for flexural loads.
% See nitsche1.pdf

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

pencoeff = 10;

Bcol1 = [1; 3];
Bcol2 = [2; 3];
col1 = [1; 2];
col2 = [2; 1];

switch isw %Task Switch
    
    case 1
            
        if ndf == 3 % beam and bar; hard-coded for Q9 right now

            for j = 2:10
            for i = ndf
                lie(i,j) = 0;
            end
            end
            
            for j = 11:12
            for i = 1
                lie(i,j) = 0;
            end
            end

        elseif ndf > 3 % beam only version

            for i = 2:ndf
                lie(i,1) = 0;
            end

        end
%%
    case {3,6} %interface stiffness
        
        % For now, assume bar on right, continuum on left
        
%         if elem == n1*m1/4+3 % Works only for BeamCouple8.m; node numbers need adjusted otherwise
%             VMforce = [0 0 0 0];
%         end
        
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
        PatchIR = matepropR(2);
        DmatR = PatchER; % Inertia is taken care of through integration (Jdet*thick)
        
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
        
%         % Determin bounds of integration segment
%         InterConn2D2 % InterConn2DT % 
% 
%         % Set jacobian for integration space
        drdr = 1;%(eL2 - eL1)/drL;
%         
%         m = (eR2-eR1)/(eL1-eL2);

        gamL = diag([0 1]);1/2;
        gamR = diag([1 0]);1/2;
        ep = pencoeff*diag([0 1]);
        lint = 3;
        der = 0;
        bf = 0;
            
        if ndf == 3
        bdofs = [2 5 3 6];
        else
        bdofs = [1 3 2 4];
        end
        lenR = xlR(1,2) - xlR(1,1);
        [shpw,shpt] = shp1dh(-1,lenR);
        
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
        
            xint = xlL*shlL;
            yint = xint(2);
                
            NmatR(1:2,bdofs) = [-yint*shpw(1,:) -yint*shpt(1,:)
                            shpw(4,:) shpt(4,:)];
                
            BmatR(1:3,bdofs) = [-yint*shpw(2,:) -yint*shpt(2,:)
                            zeros(1,4)
                            -shpw(3,:)/thick*bInertia -shpt(3,:)/thick*bInertia];
              % For a beam element with length 2, i.e. from BeamCouple7.m,
              % the line below is equivalent to the line above
%             BmatR(3,:) = -[3.750000000000001   3.750000000000001  -3.750000000000001   3.750000000000001]/10/thick;
            
            bnAdN1 = (gamL*nvect*DmatL*BmatL);
            bnAdN2 = (gamR*nvect*DmatR*BmatR);
            N1 = NmatL;
            N2 = NmatR;
            
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
%             VMforce = VMforce + [sum(ElemFL(2:2:18)) ElemFR(1) 0 ElemFR(2)];
%%
    case 11
        
        ElemE = zeros(numEn,1);
        
end %Task Switch
