% Tim Truster
% 06/11/2014
%
% Nonlinear DG element for small deformation von Mises plasticity
% Coupling method for mixed and displacement solid mechanics
% Implementation copied from NL_Elem54_2d

if isw ~= 1
CGtoDGarrays
nelLP = nelL;
nelRP = nelR;

inter = elem - (numel - numSI);
nodeAR = SurfacesI(inter,1);
nodeBR = SurfacesI(inter,2);
nodeAL = SurfacesI(inter,3);
nodeBL = SurfacesI(inter,4);
ielmatL = MatTypeTable(2,maL); % bulk element type for neighboring solid element
ielmatR = MatTypeTable(2,maR);
if length(matepropL)>6
    plasmodelL = matepropL(7);
    if length(matepropL)>7
        plasversionL = matepropL(8);
    else
        plasversionL = 1;2;
    end
else
    plasmodelL = 1;
    plasversionL = 1;2;
end
if length(matepropR)>6
    plasmodelR = matepropR(7);
    if length(matepropR)>7
        plasversionR = matepropR(8);
    else
        plasversionR = 1;2;
    end
else
    plasmodelR = 1;
    plasversionR = 1;2;
end
end

nitvms = 1;2;
if nitvms == 1 %VMS parameter for the stability tensor rp
pencoeff = 4;20;6;3;1;1;
elseif nitvms == 2 %Nitsche
pencoeff = 1;4;2;
else %RFB
pencoeff = 1;4;2;
end

if ~exist('modifDG','var') % allows the value to be changed in input file
modifDG = 3;1;4;
end
forceelast = 1; % force elastic value of tau


switch isw %Task Switch
%%
    case 1
        
        if ndf > 3
            
            for i = 4:ndf
                lie(i,1) = 0;
            end
            
        end
        
        % History is stored as follows:
        % 1-7 is for left element, interior (tau calculation)
        % 8-14 is for right element, interior
        % 15-35 is for left element, interface/DG integration points
        % 36-56 is for right element, interface
        nh1 = 7*2 + 7*3*2;
        % Store the values of penalty, weighting tensors
        nh3 = 12;
        
    case {3,6}%interface stiffness

        % Allocate arrays
        ElemKRR = zeros(nstR,nstR);
        ElemKRL = zeros(nstR,nstL);
        ElemKLR = zeros(nstL,nstR);
        ElemKLL = zeros(nstL,nstL);
        ElemFR = zeros(nstR,1);
        ElemFL = zeros(nstL,1);
        t1 = zeros(3,1);
        t2 = t1;
        t3 = t1;           
        
        % Set Material Properties from left and right side
        ElemYML = matepropL(2);
        ElemvL = matepropL(3);
        thickL = matepropL(1);
        KhardL = matepropL(4);
        HhardL = matepropL(5);
        sigyL = matepropL(6);
        ElemYMR = matepropR(2);
        ElemvR = matepropR(3);
        thickR = matepropR(1);
        KhardR = matepropR(4);
        HhardR = matepropR(5);
        sigyR = matepropR(6);
        
        bulkL = ElemYML/(3*(1-2*ElemvL));
        muL = ElemYML/(2*(1+ElemvL));
        bulkR = ElemYMR/(3*(1-2*ElemvR));
        muR = ElemYMR/(2*(1+ElemvR));
        One = [1; 1; 0];
        ind3to2 = [1 2 4];
        OneOne = ([1.0d0; 1.0d0; 1.0d0; 0.0d0; 0.0d0; 0.0d0]*[1.0d0; 1.0d0; 1.0d0; 0.0d0; 0.0d0; 0.0d0]');
       
        NmatL = zeros(2,nstL);
        BmatL = zeros(4,nstL);
        NmatR = zeros(2,nstR);
        BmatR = zeros(4,nstR);
        dmatL = zeros(9,3);
        dmatR = zeros(9,3);
        
        ul_nL = ul_n(1:ndf,1:nelL);
        ul_nR = ul_n(1:ndf,nelL+1:nelL+nelR);
       

        % Determin bounds of integration segment
        InterConn2D2 % InterConn2DT %
        %%%%%%%%%% NOTE: having to integrate over the interior to get tau,
        %%%%%%%%%% on a truncated sector (i.e. different integration
        %%%%%%%%%% points) will be BAD for crystal plasticity; think about
        %%%%%%%%%% using a constant value or the mid-point of the element
        %%%%%%%%%% only for computing tau. Do this in the future as a study
        
        % Set jacobian for integration space
        drdr = (eL2 - eL1)/drL;
        drdrR = (eR2 - eR1)/drR; % This causes errors for partial interface segments
        m = (eR2-eR1)/(eL1-eL2);
%         if nelL == 3 || nelL == 6  
%         lint = 3;
%         else
        lint = 3;
%         end
        ideriv = 0;

        
        %Form tau matrix: compute on first iteration and store; this gives
        %quadratic convergence when tau isn't recalculated (no
        %linearization of tau required).
        if nitvms == 1
            
        iterset = 1;
        if iter  <=iterset % == 0 %
        
        % Extract history from left element interior
        history = hr(nh1:nh1-1+7);
        % Compute stability tensor for left side
        if ielmatL == 51
        [tauL,history1,intb] = Tau2d53(xlintL,xlL,ulL,ul_nL,history,matepropL,plasmodelL,plasversionL,nelL,nen,ndf,forceelast);
        else
        [tauL,history1,intb] = Tau2d54(xlintL,xlL,ulL,ul_nL,history,matepropL,plasmodelL,plasversionL,nelL,nen,ndf,forceelast);
        end
        % Store history from left element interior
        hr(nh2:nh2-1+7) = history1;
        % Extract history from right element interior
        history = hr(nh1+7:nh1-1+7*2);
        % Compute stability tensor for right side
        if ielmatR == 51
        [tauR,history1,intb] = Tau2d53(xlintR,xlR,ulR,ul_nR,history,matepropR,plasmodelR,plasversionR,nelR,nen,ndf,forceelast);
        else
        [tauR,history1,intb] = Tau2d54(xlintR,xlR,ulR,ul_nR,history,matepropR,plasmodelR,plasversionR,nelR,nen,ndf,forceelast);
        end
        % Store history from right element interior
        hr(nh2+7:nh2-1+7*2) = history1;


        % Integrate edge bubble functions over interface segment
        ebL = 0;
        ebR = 0;
        intedge = 0;
        
        for ie = 1:lint            
            % For separate bubble types on T and Q
            if nelL == 3 || nelL == 6
                [Wgt,litr,lits] = intpntt(ie,lint,1); 
                 ebeL = edgebubble(litr,lits,nelL);
            else
                [Wgt,litr,lits] = intpntq(ie,lint,1);
                 ebeL = edgebubbleQ(litr,lits,nelL);
            end
                            
            if nelR == 3 || nelR == 6
                [WgtR,litrR,litsR] = intpntt(ie,lint,1); 
                 ebeR = edgebubble(litrR,litsR,nelR);
            else
                [WgtR,litrR,litsR] = intpntq(ie,lint,1);
                 ebeR = edgebubbleQ(litrR,litsR,nelR);
            end         
                    
            if nelL == 3 || nelL == 6
                rL = drdr*(litr-roL)+eL1;
                [shlL,shldL,shlsL,be] = shlt(rL,lits,nelL,nelL,0,0);
                [PxyL,shgsL,JdetL,bubbleL,XsL] = shgt(xlL(:,1:nelL),nelL,shldL,shlsL,nen,0,0,be); 
            else    
                rL = drdr*(litr-roL)+eL1;
                [shlL,shldL,shlsL,be] = shlq(rL,lits,nelL,nelL,0,0);
                [PxyL,shgsL,JdetL,bubbleL,XsL] = shgq(xlL(:,1:nelL),nelL,shldL,shlsL,nen,0,0,be); 
            end

            %Evaluate tangent and normal vectors
            T1L = [XsL(:,1); 0];
            [Tm1L, Tu1L] = VecNormalize(T1L);
            T2L = [0; 0; 1];
            Tm2L = 1;
            Tu2L = T2L';
            T3L = VecCrossProd(T1L,T2L);
            [Tm3L, Tu3L] = VecNormalize(T3L);            
            C1L = drdr*Wgt*Tm3L;
            
            ebL = ebL + C1L*ebeL;
            ebR = ebR + C1L*ebeR;
            intedge = intedge + C1L;

        end   
        
        % Compute penalty and weighting tensors, store for later iterations
        edgeK = tauL*ebL^2 + tauR*ebR^2;
        gamL = ebL^2*(edgeK\tauL);
        hr(nh3:nh3+3) = [gamL(1,1),gamL(1,2),gamL(2,1),gamL(2,2)];
        gamR = ebR^2*(edgeK\tauR);
        hr(nh3+4:nh3+7) = [gamR(1,1),gamR(1,2),gamR(2,1),gamR(2,2)];
        ep = pencoeff*intedge*inv(edgeK); 
        hr(nh3+8:nh3+11) = [ep(1,1),ep(1,2),ep(2,1),ep(2,2)];        
        else % load taus from previous iterations
        gamL = [hr(nh3) hr(nh3+1)
                hr(nh3+2) hr(nh3+3)];
        gamR = [hr(nh3+4) hr(nh3+5)
                hr(nh3+6) hr(nh3+7)];
        ep =   [hr(nh3+8) hr(nh3+9)
                hr(nh3+10) hr(nh3+11)];
        end %iterset
        
        else % Nitsche method, hard-coded penalty terms
            
        gamL = [0.5 0
                0 0.5];
        gamR = [0.5 0
                0 0.5];     
        ep = 20000*pencoeff*eye(2);
        
        end %nitvms
        

        rhspulL = reshape(ulL,ndf*nelL,1);
        rhspulR = reshape(ulR,ndf*nelR,1);
        rhspdulL = reshape(ulL-ul_nL,ndf*nelL,1);
        rhspdulR = reshape(ulR-ul_nR,ndf*nelR,1);
        % Main integration loop for DG residual/stiffness terms
        ll=0;           
        for l = 1:lint

            ll = ll + 1;
            if nelL == 3 || nelL == 6
                [Wgt,litr,lits] = intpntt(ll,lint,1);
            elseif nelL == 4 || nelL == 9
                [Wgt,litr,lits] = intpntq(ll,lint,1);
            end

            rL = drdr*(litr-roL)+eL1;
            if nelL == 3 || nelL == 6
            if ielmatL == 52
              shlpL = shlt(rL,lits,3,nelL,0,0);
            end
            [shlL,shldL,shlsL,be] = shlt(rL,lits,nelL,nelL,0,0);
            [QxyL,shgsL,JdetL,bubbleL,XsL] = shgt(xlL(:,1:nelL),nelL,shldL,shlsL,nen,0,0,be);
            else
            if ielmatL == 52
              shlpL = shlq(rL,lits,4,nelL,0,0);
            end
            [shlL,shldL,shlsL,be] = shlq(rL,lits,nelL,nelL,0,0);
            [QxyL,shgsL,JdetL,bubbleL,XsL] = shgq(xlL(:,1:nelL),nelL,shldL,shlsL,nen,0,0,be);     
            end
            
            rR = m*(rL-eL2) + eR1;
            if nelR == 3 || nelR == 6
                sR = 0;
            else %if nelR == 4
                sR = -1;
            end

            if nelR == 3 || nelR == 6  
            if ielmatL == 52
              shlpR = shlt(rR,sR,3,nelR,0,0);
            end          
            [shlR,shldR,shlsR,be] = shlt(rR,sR,nelR,nelR,0,0);
            [QxyR,shgsR,JdetR,bubbleR,XsR] = shgt(xlR(:,1:nelR),nelR,shldR,shlsR,nen,0,0,be);
            else
            if ielmatR == 52
              shlpR = shlq(rR,sR,4,nelR,0,0);
            end
            [shlR,shldR,shlsR,be] = shlq(rR,sR,nelR,nelR,0,0);
            [QxyR,shgsR,JdetR,bubbleR,XsR] = shgq(xlR(:,1:nelR),nelR,shldR,shlsR,nen,0,0,be);
            end
             
            % Evaluate tangent and normal vectors
            T1L = [XsL(:,1); 0];
            [Tm1L, Tu1L] = VecNormalize(T1L);
            T2L = [0; 0; 1];
            Tm2L = 1;
            Tu2L = T2L';
            T3L = VecCrossProd(T1L,T2L);
            [Tm3L, Tu3L] = VecNormalize(T3L);
            C1L = drdr*Wgt*Tm3L;
            C1R = C1L;
            
            % Copy into "current" configuration placeholders
            tu3L = Tu3L;
            tu3R = -Tu3L;
            c1L = C1L;
            c1R = C1L;

            %normal vectors
            nLx = tu3L(1);
            nLy = tu3L(2);
            nRx = tu3R(1);
            nRy = tu3R(2);   
            NLx = Tu3L(1);
            NLy = Tu3L(2);         
            nvectL1= [nLx 0   nLy   
                      0  nLy  nLx ];
            nvectR1 = [nRx 0  nRy  
                      0  nRy  nRx ];                   
            nvecL = [nLx; nLy];
            nvecR = [nRx; nRy];

            
            %%%%%%%%%%% Start Left side

            % Method/formulation type branch
            if ielmatL == 51
                
                % Shape function arrays
                for ie = 1:nelL
                  NmatL(1,(ie-1)*ndf+1) = shlL(ie);
                  NmatL(2,(ie-1)*ndf+2) = shlL(ie);

                  BmatL(Bcol1,(ie-1)*ndf+1) = QxyL(ie,col1);
                  BmatL(Bcol2,(ie-1)*ndf+2) = QxyL(ie,col2);
                end
                
            %% Material Models
            % Extract history for left element along interface segment
            if plasversionL == 1
                % Compute input for Radial Return
                eps2d = BmatL*rhspulL; %2D strain
                pn_1 = bulkL*One'*eps2d(1:3);
                eps3d = [eps2d(1); eps2d(2); 0; eps2d(3); 0; 0];
                ephr = nh1-1+14+(l-1)*7; %pointer for plastic strain at pt l
                betahr = nh1-1+14+(l-1)*7+3;
                ahr = nh1-1+14+l*7;
                ep_n = [hr(ephr+1); hr(ephr+2); 0; hr(ephr+3); 0; 0];
                ep_n(3) = -(ep_n(1) + ep_n(2));
                beta_n = [hr(betahr+1); hr(betahr+2); 0; hr(betahr+3); 0; 0];
                beta_n(3) = -(beta_n(1) + beta_n(2));
                a_n = hr(ahr);
                
                [sig_n1,C_n1,ep_n1,beta_n1,a_n1] = J2RadialReturn0(eps3d,ep_n,beta_n,a_n,muL,bulkL,KhardL,HhardL,sigyL);
                 
                if modifDG == 0
                    if a_n1 > 0
                    [dmatL1]=dmat2_ep_numer(du,ep_n,beta_n,a_n,ElemYM,Elemv,Khard,Hhard,sigy,Cdev_n1,plasversion,plasmodel);
                    else
                        dmatL1 = zeros(9,3);
                    end
                end
                
                % Convert output from 3D to 2D
                sigma2L = sig_n1(ind3to2);
                cmatL = C_n1(ind3to2,ind3to2);
                
                % Store history variables
                ephr = nh2-1+14+(l-1)*7; %pointer for plastic strain at pt l
                betahr = nh2-1+14+(l-1)*7+3;
                ahr = nh2-1+14+l*7;
                hr(ephr+1) = ep_n1(1);
                hr(ephr+2) = ep_n1(2);
                hr(ephr+3) = ep_n1(4);
                hr(betahr+1) = beta_n1(1);
                hr(betahr+2) = beta_n1(2);
                hr(betahr+3) = beta_n1(4);
                hr(ahr) = a_n1;
                
            else
                % Compute input for Radial Return
                eps2d = BmatL*rhspulL; %2D strain
                pn_1 = bulkL*One'*eps2d(1:3);
                deps2d = BmatL*rhspdulL; %3D enhanced strain
                deps_n1 = [deps2d(1); deps2d(2); deps2d(3); 0];
                ephr = nh1-1+14+(l-1)*7; %pointer for plastic strain at pt l
%                 betahr = nh1-1+(l-1)*7+3;
                ahr = nh1-1+14+l*7;
                ee_n = [hr(ephr+1); hr(ephr+2); hr(ephr+3); hr(ephr+4)];
%                 beta_n = [hr(betahr+1); hr(betahr+2); 0; hr(betahr+3); 0; 0];
%                 beta_n(3) = -(beta_n(1) + beta_n(2));
                a_n = hr(ahr);
                eps_n1 = ee_n + deps_n1;
                
                if plasmodelL == 1
                [DGAMA,LALGVA,RSTAVA,STRES] = SUVM(2,2,[ElemYML,ElemvL],[0,100;sigyL,KhardL*100+sigyL],[ee_n; a_n],eps_n1);
                C_n1 = CTVM(DGAMA,LALGVA(1),2,2,[ElemYML,ElemvL],[0,100;sigyL,KhardL*100+sigyL],RSTAVA,STRES);
                elseif plasmodelL == 2
                [DGAM,LALGVA,RSTAVA,STRES] = SUTR(2,2,[ElemYML,ElemvL],[0,100;sigy,KhardL*100+sigyL],[ee_n; a_n],eps_n1);
                C_n1 = CTTR(LALGVA(1),2,LALGVA,2,[ElemYML,ElemvL],[0,100;sigyL,KhardL*100+sigyL],RSTAVA,eps_n1,STRES);
                end
                
                if modifDG == 0
                    if RSTAVA(5) > 0
                    [dmatL1]=dmat2_ep_numer(deps2d,ee_n,0,a_n,ElemYM,Elemv,Khard,Hhard,sigy,C_n1(1:3,1:3) - bulkL*[1 1 0]'*[1 1 0],plasversion,plasmodel);
                    else
                        dmatL1 = zeros(9,3);
                    end
                end
                
                % Convert output from 3D to 2D
                sig_n1 = [STRES(1); STRES(2); STRES(4); STRES(3); 0; 0];
                sigma2L = sig_n1(ind3to2);
                cmatL = C_n1(1:3,1:3);
                
                % Store history variables
                ephr = nh2-1+14+(l-1)*7; %pointer for plastic strain at pt l
%                 betahr = nh2-1+(l-1)*7+3;
                ahr = nh2-1+14+l*7;
                hr(ephr+1) = RSTAVA(1);
                hr(ephr+2) = RSTAVA(2);
                hr(ephr+3) = RSTAVA(3);
                hr(ephr+4) = RSTAVA(4);
%                 hr(betahr+1) = beta_n1(1);
%                 hr(betahr+2) = beta_n1(2);
%                 hr(betahr+3) = beta_n1(4);
                hr(ahr) = RSTAVA(5);
            end
          
            
            % Cauchy stress tensors for both sides
            SmatL1=[sigma2L(1), sigma2L(3)
                    sigma2L(3), sigma2L(2)];


            % stiffness matrix terms for weighted stress and material moduli
            term18L=(gamL*nvectL1*[cmatL(1:3,1:3) 0*One])'; % One is a placeholder
            
            else
                
                % Shape function arrays
                for mm = 1:nelL    
                NmatL(:,3*mm-2:3*mm) = [shlL(mm,1)     0    0       
                                        0        shlL(mm,1) 0]; 

                BmatL(:,3*mm-2:3*mm) = [QxyL(mm,1) 0         0
                                        0         QxyL(mm,2) 0        
                                        QxyL(mm,2) QxyL(mm,1) 0
                                        0 0 shlpL(mm)];    
                end 
           
            % Material Models
            % Extract history for left element along interface segment
            if plasversionL == 1
                % Compute input for Radial Return
                du = BmatL*rhspulL; %2D strain
                pn_1L = du(4);
                eps2d = du(1:3); %2D enhanced strain
                eps3d = [eps2d(1); eps2d(2); 0; eps2d(3); 0; 0];
                ephr = nh1-1+14+(l-1)*7; %pointer for plastic strain at pt l
                betahr = nh1-1+14+(l-1)*7+3;
                ahr = nh1-1+14+l*7;
                ep_n = [hr(ephr+1); hr(ephr+2); 0; hr(ephr+3); 0; 0];
                ep_n(3) = -(ep_n(1) + ep_n(2));
                beta_n = [hr(betahr+1); hr(betahr+2); 0; hr(betahr+3); 0; 0];
                beta_n(3) = -(beta_n(1) + beta_n(2));
                a_n = hr(ahr);
                
                [sig_n1,C_n1,ep_n1,beta_n1,a_n1,sdev3,Cdev_n1] = J2RadialReturn0(eps3d,ep_n,beta_n,a_n,muL,bulkL,KhardL,HhardL,sigyL);
                 
                if modifDG == 0
                    if a_n1 > 0
                    [dmatL1]=dmat2_ep_numer(du,ep_n,beta_n,a_n,ElemYM,Elemv,Khard,Hhard,sigy,Cdev_n1,plasversion,plasmodel);
                    else
                        dmatL1 = zeros(9,3);
                    end
                end
                
                % Convert output from 3D to 2D
                sigma2L = sdev3(ind3to2);
                cmatL = Cdev_n1(ind3to2,ind3to2);
%                 sigma2L = sig_n1(ind3to2);
%                 cmatL = C_n1(ind3to2,ind3to2);
                
                % Store history variables
                ephr = nh2-1+14+(l-1)*7; %pointer for plastic strain at pt l
                betahr = nh2-1+14+(l-1)*7+3;
                ahr = nh2-1+14+l*7;
                hr(ephr+1) = ep_n1(1);
                hr(ephr+2) = ep_n1(2);
                hr(ephr+3) = ep_n1(4);
                hr(betahr+1) = beta_n1(1);
                hr(betahr+2) = beta_n1(2);
                hr(betahr+3) = beta_n1(4);
                hr(ahr) = a_n1;
                
            else
                % Compute input for Radial Return
                du = BmatL*rhspulL; %2D strain
                pn_1L = du(4);
                eps2d = du(1:3); %2D enhanced strain
                deps2d = BmatL*rhspdulL; %3D enhanced strain
                deps_n1 = [deps2d(1); deps2d(2); deps2d(3); 0];
                ephr = nh1-1+14+(l-1)*7; %pointer for plastic strain at pt l
%                 betahr = nh1-1+(l-1)*7+3;
                ahr = nh1-1+14+l*7;
                ee_n = [hr(ephr+1); hr(ephr+2); hr(ephr+3); hr(ephr+4)];
%                 beta_n = [hr(betahr+1); hr(betahr+2); 0; hr(betahr+3); 0; 0];
%                 beta_n(3) = -(beta_n(1) + beta_n(2));
                a_n = hr(ahr);
                eps_n1 = ee_n + deps_n1;
                
                if plasmodelL == 1
                [DGAMA,LALGVA,RSTAVA,STRES] = SUVM(2,2,[ElemYML,ElemvL],[0,100;sigyL,KhardL*100+sigyL],[ee_n; a_n],eps_n1);
                C_n1 = CTVM(DGAMA,LALGVA(1),2,2,[ElemYML,ElemvL],[0,100;sigyL,KhardL*100+sigyL],RSTAVA,STRES);
                elseif plasmodelL == 2
                [DGAM,LALGVA,RSTAVA,STRES] = SUTR(2,2,[ElemYML,ElemvL],[0,100;sigy,KhardL*100+sigyL],[ee_n; a_n],eps_n1);
                C_n1 = CTTR(LALGVA(1),2,LALGVA,2,[ElemYML,ElemvL],[0,100;sigyL,KhardL*100+sigyL],RSTAVA,eps_n1,STRES);
                end
                
                if modifDG == 0
                    if RSTAVA(5) > 0
                    [dmatL1]=dmat2_ep_numer(deps2d,ee_n,0,a_n,ElemYM,Elemv,Khard,Hhard,sigy,C_n1(1:3,1:3) - bulkL*[1 1 0]'*[1 1 0],plasversion,plasmodel);
                    else
                        dmatL1 = zeros(9,3);
                    end
                end
                
                % Convert output from 3D to 2D
                sig_n1 = [STRES(1); STRES(2); STRES(4); STRES(3); 0; 0];
                sdev3 = sig_n1 - 1/3*OneOne*sig_n1;
                sigma2L = sdev3(ind3to2);
                cmatL = C_n1(1:3,1:3) - bulkL*[1 1 0]'*[1 1 0];
                
                % Store history variables
                ephr = nh2-1+14+(l-1)*7; %pointer for plastic strain at pt l
%                 betahr = nh2-1+(l-1)*7+3;
                ahr = nh2-1+14+l*7;
                hr(ephr+1) = RSTAVA(1);
                hr(ephr+2) = RSTAVA(2);
                hr(ephr+3) = RSTAVA(3);
                hr(ephr+4) = RSTAVA(4);
%                 hr(betahr+1) = beta_n1(1);
%                 hr(betahr+2) = beta_n1(2);
%                 hr(betahr+3) = beta_n1(4);
                hr(ahr) = RSTAVA(5);
            end
          
            
            % Cauchy stress tensors for both sides
            SmatL1=[sigma2L(1), sigma2L(3)
                    sigma2L(3), sigma2L(2)]+pn_1L*eye(2);


            % stiffness matrix terms for weighted stress and material moduli
            term18L=(gamL*nvectL1*[cmatL(1:3,1:3) One])';
            
            end % method type
            
            %%%%%%%%%%% End Left side
            
            
            
            %%%%%%%%%%% Start Right side

            % Method/formulation type branch
            if ielmatR == 51
                
                % Shape function arrays
                for ie = 1:nelR
                  NmatR(1,(ie-1)*ndf+1) = shlR(ie);
                  NmatR(2,(ie-1)*ndf+2) = shlR(ie);

                  BmatR(Bcol1,(ie-1)*ndf+1) = QxyR(ie,col1);
                  BmatR(Bcol2,(ie-1)*ndf+2) = QxyR(ie,col2);
                end 
                
            % Extract history for right element along interface segment
            if plasversionR == 1
                % Compute input for Radial Return
                eps2d = BmatR*rhspulR; %2D strain
                pn_1 = bulkR*One'*eps2d(1:3);
                eps3d = [eps2d(1); eps2d(2); 0; eps2d(3); 0; 0];
                ephr = nh1-1+35+(l-1)*7; %pointer for plastic strain at pt l
                betahr = nh1-1+35+(l-1)*7+3;
                ahr = nh1-1+35+l*7;
                ep_n = [hr(ephr+1); hr(ephr+2); 0; hr(ephr+3); 0; 0];
                ep_n(3) = -(ep_n(1) + ep_n(2));
                beta_n = [hr(betahr+1); hr(betahr+2); 0; hr(betahr+3); 0; 0];
                beta_n(3) = -(beta_n(1) + beta_n(2));
                a_n = hr(ahr);
                
                [sig_n1,C_n1,ep_n1,beta_n1,a_n1] = J2RadialReturn0(eps3d,ep_n,beta_n,a_n,muR,bulkR,KhardR,HhardR,sigyR);
                
                if modifDG == 0 
                    if a_n1 > 0
                    [dmatR1]=dmat2_ep_numer(du,ep_n,beta_n,a_n,ElemYM,Elemv,Khard,Hhard,sigy,Cdev_n1,plasversion,plasmodel);
                    else
                        dmatR1 = zeros(9,3);
                    end
                end
                
                % Convert output from 3D to 2D
                sigma2R = sig_n1(ind3to2);
                cmatR = C_n1(ind3to2,ind3to2);
                
                % Store history variables
                ephr = nh2-1+35+(l-1)*7; %pointer for plastic strain at pt l
                betahr = nh2-1+35+(l-1)*7+3;
                ahr = nh2-1+35+l*7;
                hr(ephr+1) = ep_n1(1);
                hr(ephr+2) = ep_n1(2);
                hr(ephr+3) = ep_n1(4);
                hr(betahr+1) = beta_n1(1);
                hr(betahr+2) = beta_n1(2);
                hr(betahr+3) = beta_n1(4);
                hr(ahr) = a_n1;
                
            else
                % Compute input for Radial Return
                eps2d = BmatR*rhspulR; %2D strain
                pn_1 = bulkR*One'*eps2d(1:3);
                deps2d = BmatR*rhspdulR; %3D enhanced strain
                deps_n1 = [deps2d(1); deps2d(2); deps2d(3); 0];
                ephr = nh1-1+35+(l-1)*7; %pointer for plastic strain at pt l
%                 betahr = nh1-1+(l-1)*7+3;
                ahr = nh1-1+35+l*7;
                ee_n = [hr(ephr+1); hr(ephr+2); hr(ephr+3); hr(ephr+4)];
%                 beta_n = [hr(betahr+1); hr(betahr+2); 0; hr(betahr+3); 0; 0];
%                 beta_n(3) = -(beta_n(1) + beta_n(2));
                a_n = hr(ahr);
                eps_n1 = ee_n + deps_n1;
                
                if plasmodelR == 1
                [DGAMA,LALGVA,RSTAVA,STRES] = SUVM(2,2,[ElemYMR,ElemvR],[0,100;sigyR,KhardR*100+sigyR],[ee_n; a_n],eps_n1);
                C_n1 = CTVM(DGAMA,LALGVA(1),2,2,[ElemYMR,ElemvR],[0,100;sigyR,KhardR*100+sigyR],RSTAVA,STRES);
                elseif plasmodelR == 2
                [DGAM,LALGVA,RSTAVA,STRES] = SUTR(2,2,[ElemYMR,ElemvR],[0,100;sigy,KhardR*100+sigyR],[ee_n; a_n],eps_n1);
                C_n1 = CTTR(LALGVA(1),2,LALGVA,2,[ElemYMR,ElemvR],[0,100;sigyR,KhardR*100+sigyR],RSTAVA,eps_n1,STRES);
                end
                 
                if modifDG == 0 
                    if RSTAVA(5) > 0
                    [dmatR1]=dmat2_ep_numer(deps2d,ee_n,0,a_n,ElemYM,Elemv,Khard,Hhard,sigy,C_n1(1:3,1:3) - bulkR*[1 1 0]'*[1 1 0],plasversion,plasmodel);
                    else
                        dmatR1 = zeros(9,3);
                    end
                end
                
                % Convert output from 3D to 2D
                sig_n1 = [STRES(1); STRES(2); STRES(4); STRES(3); 0; 0];
                sigma2R = sig_n1(ind3to2);
                cmatR = C_n1(1:3,1:3);
                
                % Store history variables
                ephr = nh2-1+35+(l-1)*7; %pointer for plastic strain at pt l
%                 betahr = nh2-1+(l-1)*7+3;
                ahr = nh2-1+35+l*7;
                hr(ephr+1) = RSTAVA(1);
                hr(ephr+2) = RSTAVA(2);
                hr(ephr+3) = RSTAVA(3);
                hr(ephr+4) = RSTAVA(4);
%                 hr(betahr+1) = beta_n1(1);
%                 hr(betahr+2) = beta_n1(2);
%                 hr(betahr+3) = beta_n1(4);
                hr(ahr) = RSTAVA(5);
            end
            %% End Material Models
          
            
            % Cauchy stress tensors for both sides
            SmatR1=[sigma2R(1), sigma2R(3)
                    sigma2R(3), sigma2R(2)];


            % stiffness matrix terms for weighted stress and material moduli
            term18R=(gamR*nvectR1*[cmatR(1:3,1:3) 0*One])'; % One is a placeholder
                
            else
                
                % Shape function arrays
                for mm = 1:nelR    
                NmatR(:,3*mm-2:3*mm) = [shlR(mm,1)     0     0      
                                         0        shlR(mm,1) 0]; 

                BmatR(:,3*mm-2:3*mm) = [QxyR(mm,1) 0          0    
                                        0         QxyR(mm,2)  0
                                        QxyR(mm,2) QxyR(mm,1) 0
                                        0 0 shlpR(mm)];                     
                end  
                
            % Extract history for right element along interface segment
            if plasversionR == 1
                % Compute input for Radial Return
                du = BmatR*rhspulR; %2D strain
                pn_1R = du(4);
                eps2d = du(1:3); %2D enhanced strain
                eps3d = [eps2d(1); eps2d(2); 0; eps2d(3); 0; 0];
                ephr = nh1-1+35+(l-1)*7; %pointer for plastic strain at pt l
                betahr = nh1-1+35+(l-1)*7+3;
                ahr = nh1-1+35+l*7;
                ep_n = [hr(ephr+1); hr(ephr+2); 0; hr(ephr+3); 0; 0];
                ep_n(3) = -(ep_n(1) + ep_n(2));
                beta_n = [hr(betahr+1); hr(betahr+2); 0; hr(betahr+3); 0; 0];
                beta_n(3) = -(beta_n(1) + beta_n(2));
                a_n = hr(ahr);
                
                [sig_n1,C_n1,ep_n1,beta_n1,a_n1,sdev3,Cdev_n1] = J2RadialReturn0(eps3d,ep_n,beta_n,a_n,muR,bulkR,KhardR,HhardR,sigyR);
                
                if modifDG == 0 
                    if a_n1 > 0
                    [dmatR1]=dmat2_ep_numer(du,ep_n,beta_n,a_n,ElemYM,Elemv,Khard,Hhard,sigy,Cdev_n1,plasversion,plasmodel);
                    else
                        dmatR1 = zeros(9,3);
                    end
                end
                
                % Convert output from 3D to 2D
                sigma2R = sdev3(ind3to2);
                cmatR = Cdev_n1(ind3to2,ind3to2);
%                 sigma2R = sig_n1(ind3to2);
%                 cmatR = C_n1(ind3to2,ind3to2);
                
                % Store history variables
                ephr = nh2-1+35+(l-1)*7; %pointer for plastic strain at pt l
                betahr = nh2-1+35+(l-1)*7+3;
                ahr = nh2-1+35+l*7;
                hr(ephr+1) = ep_n1(1);
                hr(ephr+2) = ep_n1(2);
                hr(ephr+3) = ep_n1(4);
                hr(betahr+1) = beta_n1(1);
                hr(betahr+2) = beta_n1(2);
                hr(betahr+3) = beta_n1(4);
                hr(ahr) = a_n1;
                
            else
                % Compute input for Radial Return
                du = BmatR*rhspulR; %2D strain
                pn_1R = du(4);
                eps2d = du(1:3); %2D enhanced strain
                deps2d = BmatR*rhspdulR; %3D enhanced strain
                deps_n1 = [deps2d(1); deps2d(2); deps2d(3); 0];
                ephr = nh1-1+35+(l-1)*7; %pointer for plastic strain at pt l
%                 betahr = nh1-1+(l-1)*7+3;
                ahr = nh1-1+35+l*7;
                ee_n = [hr(ephr+1); hr(ephr+2); hr(ephr+3); hr(ephr+4)];
%                 beta_n = [hr(betahr+1); hr(betahr+2); 0; hr(betahr+3); 0; 0];
%                 beta_n(3) = -(beta_n(1) + beta_n(2));
                a_n = hr(ahr);
                eps_n1 = ee_n + deps_n1;
                
                if plasmodelR == 1
                [DGAMA,LALGVA,RSTAVA,STRES] = SUVM(2,2,[ElemYMR,ElemvR],[0,100;sigyR,KhardR*100+sigyR],[ee_n; a_n],eps_n1);
                C_n1 = CTVM(DGAMA,LALGVA(1),2,2,[ElemYMR,ElemvR],[0,100;sigyR,KhardR*100+sigyR],RSTAVA,STRES);
                elseif plasmodelR == 2
                [DGAM,LALGVA,RSTAVA,STRES] = SUTR(2,2,[ElemYMR,ElemvR],[0,100;sigy,KhardR*100+sigyR],[ee_n; a_n],eps_n1);
                C_n1 = CTTR(LALGVA(1),2,LALGVA,2,[ElemYMR,ElemvR],[0,100;sigyR,KhardR*100+sigyR],RSTAVA,eps_n1,STRES);
                end
                 
                if modifDG == 0 
                    if RSTAVA(5) > 0
                    [dmatR1]=dmat2_ep_numer(deps2d,ee_n,0,a_n,ElemYM,Elemv,Khard,Hhard,sigy,C_n1(1:3,1:3) - bulkR*[1 1 0]'*[1 1 0],plasversion,plasmodel);
                    else
                        dmatR1 = zeros(9,3);
                    end
                end
                
                % Convert output from 3D to 2D
                sig_n1 = [STRES(1); STRES(2); STRES(4); STRES(3); 0; 0];
                sdev3 = sig_n1 - 1/3*OneOne*sig_n1;
                sigma2R = sdev3(ind3to2);
                cmatR = C_n1(1:3,1:3) - bulkR*[1 1 0]'*[1 1 0];
                
                % Store history variables
                ephr = nh2-1+35+(l-1)*7; %pointer for plastic strain at pt l
%                 betahr = nh2-1+(l-1)*7+3;
                ahr = nh2-1+35+l*7;
                hr(ephr+1) = RSTAVA(1);
                hr(ephr+2) = RSTAVA(2);
                hr(ephr+3) = RSTAVA(3);
                hr(ephr+4) = RSTAVA(4);
%                 hr(betahr+1) = beta_n1(1);
%                 hr(betahr+2) = beta_n1(2);
%                 hr(betahr+3) = beta_n1(4);
                hr(ahr) = RSTAVA(5);
            end
            % End Material Models
          
            
            % Cauchy stress tensors for both sides
            SmatR1=[sigma2R(1), sigma2R(3)
                    sigma2R(3), sigma2R(2)]+pn_1R*eye(2);


            % stiffness matrix terms for weighted stress and material moduli
            term18R=(gamR*nvectR1*[cmatR(1:3,1:3) One])';
            
            end % method type
            
            %%%%%%%%%%% End Right side
            
            
            % average stress term for force vector
            term28L=NmatL'*(c1L*gamL*SmatL1*nvecL-c1R*gamR*SmatR1*nvecR);
            term28R=NmatR'*(c1L*gamL*SmatL1*nvecL-c1R*gamR*SmatR1*nvecR);
            
            % penalty terms
            jumpu = NmatL*rhspulL - NmatR*rhspulR;  %jumpu=uL-uR

            term30L=NmatL'*ep*jumpu;
            term30R=NmatR'*ep*jumpu;

            if modifDG == 0
            % Unexpected terms, involves jumpu and 6th order tensor          
            
             % Revised by TJT to use reshape function and smaller
             % multiplications
             dmatL2 = dmatL1*nvectL1'*gamL'*jumpu;
             dmatL2 = [reshape(dmatL2,3,3) zeros(3,1); zeros(1,4)];
             dmatR2 = dmatR1*nvectR1'*gamR'*jumpu;
             dmatR2 = [reshape(dmatR2,3,3) zeros(3,1); zeros(1,4)];
             
             term7L = BmatL'*dmatL2*BmatL;
             term7R = BmatR'*dmatR2*BmatR;
            end

             
            % Combine contributions into element force vector and stiffness
            % matrix
            if exist('modifDG','var') && modifDG > 0 % modify the DG terms used in the formulation

                % Penalty terms
                 
                ElemFL = ElemFL-(+C1L*term30L);
                ElemFR = ElemFR-(-C1R*term30R);

                ElemKLL = ElemKLL + C1L*(NmatL'*ep*NmatL);
                ElemKLR = ElemKLR - C1R*(NmatL'*ep*NmatR);
                ElemKRL = ElemKRL - C1L*(NmatR'*ep*NmatL);
                ElemKRR = ElemKRR + C1R*(NmatR'*ep*NmatR);  
             
            if modifDG >= 1 
                 
                % average stress terms

                ElemFL = ElemFL-(-term28L);
                ElemFR = ElemFR-(+term28R);

                ElemKLL = ElemKLL  - c1L*NmatL'*(term18L')*BmatL;
                ElemKLR = ElemKLR  + c1R*NmatL'*(term18R')*BmatR;
                ElemKRL = ElemKRL  + c1L*NmatR'*(term18L')*BmatL;
                ElemKRR = ElemKRR  - c1R*NmatR'*(term18R')*BmatR;
             
            if modifDG == 3 
                 
                % nonsymmetric terms

                ElemFL = ElemFL-(-c1L*BmatL'*(term18L)*jumpu);
                ElemFR = ElemFR-(c1R* BmatR'*(term18R)*jumpu);

                ElemKLL = ElemKLL - c1L*BmatL'*(term18L)*NmatL;
                ElemKLR = ElemKLR + c1L*BmatL'*(term18L)*NmatR;
                ElemKRL = ElemKRL + c1R*BmatR'*(term18R)*NmatL;
                ElemKRR = ElemKRR - c1R*BmatR'*(term18R)*NmatR;
             
            elseif modifDG == 4 
                 
                % nonsymmetric terms

                ElemKLL = ElemKLL - c1L*BmatL'*(term18L)*NmatL;
                ElemKLR = ElemKLR + c1L*BmatL'*(term18L)*NmatR;
                ElemKRL = ElemKRL + c1R*BmatR'*(term18R)*NmatL;
                ElemKRR = ElemKRR - c1R*BmatR'*(term18R)*NmatR;
             
            end
            end

            else % full method
             
                ElemFL = ElemFL-(-term28L+C1L*term30L);
                ElemFR = ElemFR-(+term28R-C1R*term30R);  

                ElemFL = ElemFL-(-c1L*BmatL'*(term18L)*jumpu);
                ElemFR = ElemFR-(c1R* BmatR'*(term18R)*jumpu); 

                ElemKLL = ElemKLL + C1L*(NmatL'*ep*NmatL);
                ElemKLR = ElemKLR - C1R*(NmatL'*ep*NmatR);
                ElemKRL = ElemKRL - C1L*(NmatR'*ep*NmatL);
                ElemKRR = ElemKRR + C1R*(NmatR'*ep*NmatR);   

                ElemKLL = ElemKLL  - c1L*NmatL'*(term18L')*BmatL;
                ElemKLR = ElemKLR  + c1R*NmatL'*(term18R')*BmatR;
                ElemKRL = ElemKRL  + c1L*NmatR'*(term18L')*BmatL;
                ElemKRR = ElemKRR  - c1R*NmatR'*(term18R')*BmatR;             

                ElemKLL = ElemKLL - c1L*BmatL'*(term18L)*NmatL;
                ElemKLR = ElemKLR + c1L*BmatL'*(term18L)*NmatR;
                ElemKRL = ElemKRL + c1R*BmatR'*(term18R)*NmatL;
                ElemKRR = ElemKRR - c1R*BmatR'*(term18R)*NmatR;           

                %relate to jumpu terms additional terms
                ElemKLL = ElemKLL - c1L*(term7L);  %                     
                ElemKRR = ElemKRR + c1R*(term7R);  %

             end
        end %lint
% ElemKLL
            ElemK = [ElemKLL ElemKLR; ElemKRL ElemKRR];
            ElemF = [ElemFL; ElemFR];
            
%%

    case -1 % boundary tractions
        
        ElemF = zeros(nst,1);
        
    case 61 % form data structure for interface segments
        
        % NOT FUNCTIONAL - needs to handle partial segments and Q4 elements

        % Set the number of interface quantities per node to be stored
        if ~exist('numIQ','var')
            numIQ = 6; % 3 disp-jump, 3 traction, 3 numerical-flux
        else
            numIQ = max(numIQ,6);
        end
        
        segment = segment + 1;
        IElemSeg(1,inter) = segment;
        IElemSeg(2,inter) = segment;
        
        % Always use face of left element for the nodes
        % nel == 9 ONLY
            Iix(segment,1:3) = (segnode+1:segnode+3);
            ICoordinates(segnode+1:segnode+3,1:ndm) = xlL(1:ndm,[1 5 2])';
            segnode = segnode + 3;
            
        Iix(segment,nenseg+1) = ma; % good habit to copy material ID too
        
    case 60 % output interface quantities for plotting
        
        % NOT FUNCTIONAL - needs to handle partial segments and Q4 elements
        
        % get segment number for DG element (usually segment=inter)
        segment = IElemSeg(1,inter);
    
        % get nodes on the interface segment
        ElemFlagI = Iix(segment,1:nenseg);
        xlI = ICoordinates(ElemFlagI,1:ndm)';
        

        % Determin bounds of integration segment
        InterConn2D2 % InterConn2DT %
        
        % Set jacobian for integration space
        drdr = (eL2 - eL1)/drL;
        drdrR = (eR2 - eR1)/drR; % This causes errors for partial interface segments
        m = (eR2-eR1)/(eL1-eL2);
        
        
        % interface traction, extrapolated from integration points
        % pull out weighting tensors and penalty parameter for this segment
        gamL = [hr(nh3) hr(nh3+1)
                hr(nh3+2) hr(nh3+3)];
        gamR = [hr(nh3+4) hr(nh3+5)
                hr(nh3+6) hr(nh3+7)];
        ep =   [hr(nh3+8) hr(nh3+9)
                hr(nh3+10) hr(nh3+11)];
        NmatL = zeros(2,nstL);
        BmatL = zeros(4,nstL);
        NmatR = zeros(2,nstR);
        BmatR = zeros(4,nstR);
            
            nelseg = 3;
        jumps = zeros(2,nelseg);
        tractions_int = zeros(4,nelseg); % tractions at integration points
        tractions = zeros(6,nelseg); % tractions at segment nodes
        ll=0;           
        for l = 1:nelseg

            ll = ll + 1;
            if nelL == 3 || nelL == 6
                [Wgt,litr,lits] = intpntt(ll,nelseg,1);
            elseif nelL == 4 || nelL == 9
                [Wgt,litr,lits] = intpntq(ll,nelseg,1);
            end
            rL = drdr*(litr-roL)+eL1;
            if nelL == 3 || nelL == 6
            [shlL,shldL,shlsL,be] = shlt(rL,lits,nelL,nelL,0,0);
            [QxyL,shgsL,JdetL,bubbleL,XsL] = shgt(xlL(:,1:nelL),nelL,shldL,shlsL,nen,0,0,be);
            else
            [shlL,shldL,shlsL,be] = shlq(rL,lits,nelL,nelL,0,0);
            [QxyL,shgsL,JdetL,bubbleL,XsL] = shgq(xlL(:,1:nelL),nelL,shldL,shlsL,nen,0,0,be);     
            end 
%                 if nelP == 3 || nelP == 6
%                   [shlp,shld,shls,bub] = shlt(litr,lits,nelP,nel,0,0);
% %                   Pxy = shgt(xl,nelP,shld,shls,nel,bf,der,bub);
%                 else
            if ielmatL == 52
              shlpL = shlq(rL,lits,4,nelL,0,0);
            end
%                   Pxy = shgq(xl,nelP,shld,shls,nel,bf,der,bub);
%                 end
            
            rR = m*(rL-eL2) + eR1;
            if nelR == 3 || nelR == 6
                sR = 0;
            else %if nelR == 4
                sR = -1;
            end

            if nelR == 3 || nelR == 6            
            [shlR,shldR,shlsR,be] = shlt(rR,sR,nelR,nelR,0,0);
            [QxyR,shgsR,JdetR,bubbleR,XsR] = shgt(xlR(:,1:nelR),nelR,shldR,shlsR,nen,0,0,be);
            else
            [shlR,shldR,shlsR,be] = shlq(rR,sR,nelR,nelR,0,0);
            [QxyR,shgsR,JdetR,bubbleR,XsR] = shgq(xlR(:,1:nelR),nelR,shldR,shlsR,nen,0,0,be);
            end
%                 if nelP == 3 || nelP == 6
%                   [shlp,shld,shls,bub] = shlt(litr,lits,nelP,nel,0,0);
% %                   Pxy = shgt(xl,nelP,shld,shls,nel,bf,der,bub);
%                 else
            if ielmatR == 52
              shlpR = shlq(rR,sR,4,nelR,0,0);
            end              
          
            
             
            % Evaluate tangent and normal vectors
            T1L = [XsL(:,1); 0];
            [Tm1L, Tu1L] = VecNormalize(T1L);
            T2L = [0; 0; 1];
            Tm2L = 1;
            Tu2L = T2L';
            T3L = VecCrossProd(T1L,T2L);
            [Tm3L, Tu3L] = VecNormalize(T3L);
            C1L = drdr*Wgt*Tm3L;
            C1R = C1L;
            
            % Copy into "current" configuration placeholders
            tu3L = Tu3L;
            tu3R = -Tu3L;
            c1L = C1L;
            c1R = C1L;

            %normal vectors
            nLx = tu3L(1);
            nLy = tu3L(2);
            nRx = tu3R(1);
            nRy = tu3R(2);   
            NLx = Tu3L(1);
            NLy = Tu3L(2);         
            nvectL1= [nLx 0   nLy   
                      0  nLy  nLx ];
            nvectR1 = [nRx 0  nRy  
                      0  nRy  nRx ];                   
            nvecL = [nLx; nLy];
            nvecR = [nRx; nRy];

            
            %%%%%%%%%%% Start Left side

            % Method/formulation type branch
            if ielmatL == 51
                
                % Shape function arrays
                for ie = 1:nelL
                  NmatL(1,(ie-1)*ndf+1) = shlL(ie);
                  NmatL(2,(ie-1)*ndf+2) = shlL(ie);

                  BmatL(Bcol1,(ie-1)*ndf+1) = QxyL(ie,col1);
                  BmatL(Bcol2,(ie-1)*ndf+2) = QxyL(ie,col2);
                end
                
            %% Material Models
            % Extract history for left element along interface segment
            if plasversionL == 1
                % Compute input for Radial Return
                eps2d = BmatL*rhspulL; %2D strain
                pn_1 = bulkL*One'*eps2d(1:3);
                eps3d = [eps2d(1); eps2d(2); 0; eps2d(3); 0; 0];
                ephr = nh1-1+14+(l-1)*7; %pointer for plastic strain at pt l
                betahr = nh1-1+14+(l-1)*7+3;
                ahr = nh1-1+14+l*7;
                ep_n = [hr(ephr+1); hr(ephr+2); 0; hr(ephr+3); 0; 0];
                ep_n(3) = -(ep_n(1) + ep_n(2));
                beta_n = [hr(betahr+1); hr(betahr+2); 0; hr(betahr+3); 0; 0];
                beta_n(3) = -(beta_n(1) + beta_n(2));
                a_n = hr(ahr);
                
                [sig_n1,C_n1,ep_n1,beta_n1,a_n1] = J2RadialReturn0(eps3d,ep_n,beta_n,a_n,muL,bulkL,KhardL,HhardL,sigyL);
                 
                
                % Convert output from 3D to 2D
                sigma2L = sig_n1(ind3to2);
                
            else
                % Compute input for Radial Return
                eps2d = BmatL*rhspulL; %2D strain
                pn_1 = bulkL*One'*eps2d(1:3);
                deps2d = BmatL*rhspdulL; %3D enhanced strain
                deps_n1 = [deps2d(1); deps2d(2); deps2d(3); 0];
                ephr = nh1-1+14+(l-1)*7; %pointer for plastic strain at pt l
%                 betahr = nh1-1+(l-1)*7+3;
                ahr = nh1-1+14+l*7;
                ee_n = [hr(ephr+1); hr(ephr+2); hr(ephr+3); hr(ephr+4)];
%                 beta_n = [hr(betahr+1); hr(betahr+2); 0; hr(betahr+3); 0; 0];
%                 beta_n(3) = -(beta_n(1) + beta_n(2));
                a_n = hr(ahr);
                eps_n1 = ee_n + deps_n1;
                
                if plasmodelL == 1
                [DGAMA,LALGVA,RSTAVA,STRES] = SUVM(2,2,[ElemYML,ElemvL],[0,100;sigyL,KhardL*100+sigyL],[ee_n; a_n],eps_n1);
                C_n1 = CTVM(DGAMA,LALGVA(1),2,2,[ElemYML,ElemvL],[0,100;sigyL,KhardL*100+sigyL],RSTAVA,STRES);
                elseif plasmodelL == 2
                [DGAM,LALGVA,RSTAVA,STRES] = SUTR(2,2,[ElemYML,ElemvL],[0,100;sigy,KhardL*100+sigyL],[ee_n; a_n],eps_n1);
                C_n1 = CTTR(LALGVA(1),2,LALGVA,2,[ElemYML,ElemvL],[0,100;sigyL,KhardL*100+sigyL],RSTAVA,eps_n1,STRES);
                end
                
                % Convert output from 3D to 2D
                sig_n1 = [STRES(1); STRES(2); STRES(4); STRES(3); 0; 0];
                sigma2L = sig_n1(ind3to2);
            end
          
            
            % Cauchy stress tensors for both sides
            SmatL1=[sigma2L(1), sigma2L(3)
                    sigma2L(3), sigma2L(2)];
            
            else
                
                % Shape function arrays
                for mm = 1:nelL    
                NmatL(:,3*mm-2:3*mm) = [shlL(mm,1)     0    0       
                                        0        shlL(mm,1) 0]; 

                BmatL(:,3*mm-2:3*mm) = [QxyL(mm,1) 0         0
                                        0         QxyL(mm,2) 0        
                                        QxyL(mm,2) QxyL(mm,1) 0
                                        0 0 shlpL(mm)];    
                end 
           
            % Material Models
            % Extract history for left element along interface segment
            if plasversionL == 1
                % Compute input for Radial Return
                du = BmatL*rhspulL; %2D strain
                pn_1L = du(4);
                eps2d = du(1:3); %2D enhanced strain
                eps3d = [eps2d(1); eps2d(2); 0; eps2d(3); 0; 0];
                ephr = nh1-1+14+(l-1)*7; %pointer for plastic strain at pt l
                betahr = nh1-1+14+(l-1)*7+3;
                ahr = nh1-1+14+l*7;
                ep_n = [hr(ephr+1); hr(ephr+2); 0; hr(ephr+3); 0; 0];
                ep_n(3) = -(ep_n(1) + ep_n(2));
                beta_n = [hr(betahr+1); hr(betahr+2); 0; hr(betahr+3); 0; 0];
                beta_n(3) = -(beta_n(1) + beta_n(2));
                a_n = hr(ahr);
                
                [sig_n1,C_n1,ep_n1,beta_n1,a_n1,sdev3,Cdev_n1] = J2RadialReturn0(eps3d,ep_n,beta_n,a_n,muL,bulkL,KhardL,HhardL,sigyL);
                 
                
                % Convert output from 3D to 2D
                sigma2L = sdev3(ind3to2);
%                 sigma2L = sig_n1(ind3to2);
%                 cmatL = C_n1(ind3to2,ind3to2);
                
                
            else
                % Compute input for Radial Return
                du = BmatL*rhspulL; %2D strain
                pn_1L = du(4);
                eps2d = du(1:3); %2D enhanced strain
                deps2d = BmatL*rhspdulL; %3D enhanced strain
                deps_n1 = [deps2d(1); deps2d(2); deps2d(3); 0];
                ephr = nh1-1+14+(l-1)*7; %pointer for plastic strain at pt l
%                 betahr = nh1-1+(l-1)*7+3;
                ahr = nh1-1+14+l*7;
                ee_n = [hr(ephr+1); hr(ephr+2); hr(ephr+3); hr(ephr+4)];
%                 beta_n = [hr(betahr+1); hr(betahr+2); 0; hr(betahr+3); 0; 0];
%                 beta_n(3) = -(beta_n(1) + beta_n(2));
                a_n = hr(ahr);
                eps_n1 = ee_n + deps_n1;
                
                if plasmodelL == 1
                [DGAMA,LALGVA,RSTAVA,STRES] = SUVM(2,2,[ElemYML,ElemvL],[0,100;sigyL,KhardL*100+sigyL],[ee_n; a_n],eps_n1);
                elseif plasmodelL == 2
                [DGAM,LALGVA,RSTAVA,STRES] = SUTR(2,2,[ElemYML,ElemvL],[0,100;sigy,KhardL*100+sigyL],[ee_n; a_n],eps_n1);
                end
                
                
                % Convert output from 3D to 2D
                sig_n1 = [STRES(1); STRES(2); STRES(4); STRES(3); 0; 0];
                sdev3 = sig_n1 - 1/3*OneOne*sig_n1;
                sigma2L = sdev3(ind3to2);
                
            end
          
            
            % Cauchy stress tensors for both sides
            SmatL1=[sigma2L(1), sigma2L(3)
                    sigma2L(3), sigma2L(2)]+pn_1L*eye(2);
            
            end % method type
            
            %%%%%%%%%%% End Left side
            
            
            
            %%%%%%%%%%% Start Right side

            % Method/formulation type branch
            if ielmatR == 51
                
                % Shape function arrays
                for ie = 1:nelR
                  NmatR(1,(ie-1)*ndf+1) = shlR(ie);
                  NmatR(2,(ie-1)*ndf+2) = shlR(ie);

                  BmatR(Bcol1,(ie-1)*ndf+1) = QxyR(ie,col1);
                  BmatR(Bcol2,(ie-1)*ndf+2) = QxyR(ie,col2);
                end 
                
            % Extract history for right element along interface segment
            if plasversionR == 1
                % Compute input for Radial Return
                eps2d = BmatR*rhspulR; %2D strain
                pn_1 = bulkR*One'*eps2d(1:3);
                eps3d = [eps2d(1); eps2d(2); 0; eps2d(3); 0; 0];
                ephr = nh1-1+35+(l-1)*7; %pointer for plastic strain at pt l
                betahr = nh1-1+35+(l-1)*7+3;
                ahr = nh1-1+35+l*7;
                ep_n = [hr(ephr+1); hr(ephr+2); 0; hr(ephr+3); 0; 0];
                ep_n(3) = -(ep_n(1) + ep_n(2));
                beta_n = [hr(betahr+1); hr(betahr+2); 0; hr(betahr+3); 0; 0];
                beta_n(3) = -(beta_n(1) + beta_n(2));
                a_n = hr(ahr);
                
                [sig_n1,C_n1,ep_n1,beta_n1,a_n1] = J2RadialReturn0(eps3d,ep_n,beta_n,a_n,muR,bulkR,KhardR,HhardR,sigyR);
                
                
                % Convert output from 3D to 2D
                sigma2R = sig_n1(ind3to2);
                
            else
                % Compute input for Radial Return
                eps2d = BmatR*rhspulR; %2D strain
                pn_1 = bulkR*One'*eps2d(1:3);
                deps2d = BmatR*rhspdulR; %3D enhanced strain
                deps_n1 = [deps2d(1); deps2d(2); deps2d(3); 0];
                ephr = nh1-1+35+(l-1)*7; %pointer for plastic strain at pt l
%                 betahr = nh1-1+(l-1)*7+3;
                ahr = nh1-1+35+l*7;
                ee_n = [hr(ephr+1); hr(ephr+2); hr(ephr+3); hr(ephr+4)];
%                 beta_n = [hr(betahr+1); hr(betahr+2); 0; hr(betahr+3); 0; 0];
%                 beta_n(3) = -(beta_n(1) + beta_n(2));
                a_n = hr(ahr);
                eps_n1 = ee_n + deps_n1;
                
                if plasmodelR == 1
                [DGAMA,LALGVA,RSTAVA,STRES] = SUVM(2,2,[ElemYMR,ElemvR],[0,100;sigyR,KhardR*100+sigyR],[ee_n; a_n],eps_n1);
                elseif plasmodelR == 2
                [DGAM,LALGVA,RSTAVA,STRES] = SUTR(2,2,[ElemYMR,ElemvR],[0,100;sigy,KhardR*100+sigyR],[ee_n; a_n],eps_n1);
                end
                 
                
                % Convert output from 3D to 2D
                sig_n1 = [STRES(1); STRES(2); STRES(4); STRES(3); 0; 0];
                sigma2R = sig_n1(ind3to2);
            end
            %% End Material Models
          
            
            % Cauchy stress tensors for both sides
            SmatR1=[sigma2R(1), sigma2R(3)
                    sigma2R(3), sigma2R(2)];
                
            else
                
                % Shape function arrays
                for mm = 1:nelR    
                NmatR(:,3*mm-2:3*mm) = [shlR(mm,1)     0     0      
                                         0        shlR(mm,1) 0]; 

                BmatR(:,3*mm-2:3*mm) = [QxyR(mm,1) 0          0    
                                        0         QxyR(mm,2)  0
                                        QxyR(mm,2) QxyR(mm,1) 0
                                        0 0 shlpR(mm)];                     
                end  
                
            % Extract history for right element along interface segment
            if plasversionR == 1
                % Compute input for Radial Return
                du = BmatR*rhspulR; %2D strain
                pn_1R = du(4);
                eps2d = du(1:3); %2D enhanced strain
                eps3d = [eps2d(1); eps2d(2); 0; eps2d(3); 0; 0];
                ephr = nh1-1+35+(l-1)*7; %pointer for plastic strain at pt l
                betahr = nh1-1+35+(l-1)*7+3;
                ahr = nh1-1+35+l*7;
                ep_n = [hr(ephr+1); hr(ephr+2); 0; hr(ephr+3); 0; 0];
                ep_n(3) = -(ep_n(1) + ep_n(2));
                beta_n = [hr(betahr+1); hr(betahr+2); 0; hr(betahr+3); 0; 0];
                beta_n(3) = -(beta_n(1) + beta_n(2));
                a_n = hr(ahr);
                
                [sig_n1,C_n1,ep_n1,beta_n1,a_n1,sdev3,Cdev_n1] = J2RadialReturn0(eps3d,ep_n,beta_n,a_n,muR,bulkR,KhardR,HhardR,sigyR);
                
                
                % Convert output from 3D to 2D
                sigma2R = sdev3(ind3to2);
%                 sigma2R = sig_n1(ind3to2);
%                 cmatR = C_n1(ind3to2,ind3to2);
                
            else
                % Compute input for Radial Return
                du = BmatR*rhspulR; %2D strain
                pn_1R = du(4);
                eps2d = du(1:3); %2D enhanced strain
                deps2d = BmatR*rhspdulR; %3D enhanced strain
                deps_n1 = [deps2d(1); deps2d(2); deps2d(3); 0];
                ephr = nh1-1+35+(l-1)*7; %pointer for plastic strain at pt l
%                 betahr = nh1-1+(l-1)*7+3;
                ahr = nh1-1+35+l*7;
                ee_n = [hr(ephr+1); hr(ephr+2); hr(ephr+3); hr(ephr+4)];
%                 beta_n = [hr(betahr+1); hr(betahr+2); 0; hr(betahr+3); 0; 0];
%                 beta_n(3) = -(beta_n(1) + beta_n(2));
                a_n = hr(ahr);
                eps_n1 = ee_n + deps_n1;
                
                if plasmodelR == 1
                [DGAMA,LALGVA,RSTAVA,STRES] = SUVM(2,2,[ElemYMR,ElemvR],[0,100;sigyR,KhardR*100+sigyR],[ee_n; a_n],eps_n1);
                elseif plasmodelR == 2
                [DGAM,LALGVA,RSTAVA,STRES] = SUTR(2,2,[ElemYMR,ElemvR],[0,100;sigy,KhardR*100+sigyR],[ee_n; a_n],eps_n1);
                end
                 
                
                % Convert output from 3D to 2D
                sig_n1 = [STRES(1); STRES(2); STRES(4); STRES(3); 0; 0];
                sdev3 = sig_n1 - 1/3*OneOne*sig_n1;
                sigma2R = sdev3(ind3to2);
                
            end
            % End Material Models
          
            
            % Cauchy stress tensors for both sides
            SmatR1=[sigma2R(1), sigma2R(3)
                    sigma2R(3), sigma2R(2)]+pn_1R*eye(2);
            
            end % method type
            
            %%%%%%%%%%% End Right side
            
            
            tvtr = gamL*SmatL1*nvecL-gamR*SmatR1*nvecR; % weighted traction
            jumpu = ulL*shlL - ulR*shlR; % displacement jump - points from Left side to Right side when uL > uR
            jumpu = jumpu(1:2);
            
            % numerical flux - traction vector pointing out from Right side
            jumps(1:2,l) = jumpu;
            tractions_int(1:4,l) = [-tvtr; -tvtr+ep*jumpu];
            
        end
        
        % reorder integration point values from tensor-product to FEM-
        % counter-clockwise so that normal shape functions can be used
%         if nelL == 4 || nelL == 10
%         tractions_int = tractions_int(:,[1 2 3]); % VERIFY THIS
%         else
%         tractions_int = tractions_int(:,[1 2 4 3]);
%         end
        
        % extrapolate values to segment nodes
%         if nelseg == 3
%             plist = [-1/3 5/3 -1/3
%                      -1/3 -1/3 5/3];
%         else % nelseg == 4
            sqr6 = 1/sqrt(0.6);
            plist = [-sqr6 0 sqr6];
%         end
        
        for ll = 1:nelseg
            
            r = plist(1,ll);
            s = -1;
            shpS = shl1d(r,1,2);
            
            tractions(:,ll) = [jumps; tractions_int]*shpS(1:3)';
            
        end
        
        
        % assemble into interface output quantities
        InterQuant(1:6,ElemFlagI,step) = tractions;
        
end