% Tim Truster
% 05/10/2015
%
% Nonlinear DG element for small deformation von Mises plasticity
% Mixed method
% Implementation copied from NL_Elem53_3d

if isw ~= 1
CGtoDGarrays
nelLP = nelL;
nelRP = nelR;

inter = elem - (numel - numSI);
nodeAR = SurfacesI(inter,1);
nodeBR = SurfacesI(inter,2);
nodeAL = SurfacesI(inter,3);
nodeBL = SurfacesI(inter,4);
if length(matepropL)>6
    plasmodel = matepropL(7);
    if length(matepropL)>7
        plasversion = matepropL(8);
    else
        plasversion = 1;2;
    end
else
    plasmodel = 1;
    plasversion = 1;2;
end
end

nitvms = 1;2;
if nitvms == 1 %VMS parameter for the stability tensor rp
pencoeff = 40;200;6;20;4;3;1;1;
elseif nitvms == 2 %Nitsche
pencoeff = 1;4;2;
else %RFB
pencoeff = 1;4;2;
end

if ~exist('modifDG','var') % allows the value to be changed in input file
modifDG = 3;1;4;
end
forceelast = 1;0; % force elastic value of tau


switch isw %Task Switch
%%
    case 1
        
        if ndf > 4
            
            for i = 5:ndf
                lie(i,1) = 0;
            end
            
        end
        
        % History is stored as follows:
        % 1-7 is for left element, interior (tau calculation)
        % 8-14 is for right element, interior
        % 15-35 is for left element, interface/DG integration points
        % 36-56 is for right element, interface
        nh1 = 13*2 + 13*13*2;
        % Store the values of penalty, weighting tensors
        nh3 = 9*3 + 6;
        iste = 4+4+4+3+6+6;
        
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
        One = [1.0d0; 1.0d0; 1.0d0; 0.0d0; 0.0d0; 0.0d0];
        OneOne = ([1.0d0; 1.0d0; 1.0d0; 0.0d0; 0.0d0; 0.0d0]*[1.0d0; 1.0d0; 1.0d0; 0.0d0; 0.0d0; 0.0d0]');
       
        NmatL = zeros(3,nstL);
        BmatL = zeros(7,nstL);
        NmatR = zeros(3,nstR);
        BmatR = zeros(7,nstR);
%         dmatL = zeros(9,3);
%         dmatR = zeros(9,3);
        
        ul_nL = ul_n(1:ndf,1:nelL);
        ul_nR = ul_n(1:ndf,nelL+1:nelL+nelR);
       

        % Determin bounds of integration segment
        xlintL = xlL;
        xlintR = xlR;
        
        % Set Gauss points
        if nelL == 4
            lint = 3;
        elseif nelL == 8
            lint = 4;
        elseif nelL == 10
            lint = 13;
        else
            lint = 9;
        end
        der = 0;
        bf = 0;
        ib = 5;   
        ll = 0; % Counter for history variables

        
        %Form tau matrix: compute on first iteration and store; this gives
        %quadratic convergence when tau isn't recalculated (no
        %linearization of tau required).
        if nitvms == 1
            
        iterset = 1;
        if iter  <=iterset % == 0 %
            
        
        % Extract history from left element interior
        history = hr(nh1:nh1-1+13);
        % Compute stability tensor for left side
        [tauL,history1,intb,CnL] = Tau3d54(xlintL,xlL,ulL,ul_nL,history,matepropL,plasmodel,plasversion,nelL,nen,ndf,forceelast);
        % Store history from left element interior
        hr(nh2:nh2-1+13) = history1;
        % Extract history from right element interior
        history = hr(nh1+13:nh1-1+13*2);
        % Compute stability tensor for right side
        [tauR,history1,intb,CnR] = Tau3d54(xlintR,xlR,ulR,ul_nR,history,matepropR,plasmodel,plasversion,nelR,nen,ndf,forceelast);
        % Store history from right element interior
        hr(nh2+13:nh2-1+13*2) = history1;

        % Integrate edge bubble functions over interface segment
        ebL = 0;
        ebR = 0;
        intedge = 0;
        
        for ie = 1:lint            
% For separate bubble types on T and Q
           if nelL == 4 || nelL == 10
                [w,ss] =  int3d_t(ie,lint,ib);
                ebeL = facebubbleT(ss,nelL);
                [shlL,shldL,shls,be] = shltt(ss,nelL,nelL,0,0);
           else
                [Wgt,ss] = intpntb(ie,lint,ib);
                [shlL,shldL,shls,be] = shlb(ss,nelL,nelL,0,0);
                ebeL = facebubbleQ(ss,nelL);
           end
            
                %Physical location of int pt
%                 xint = XlL(1,:)*shlL;
%                 yint = XlL(2,:)*shlL;
%                 zint = XlL(3,:)*shlL;
% 
%                 xi = POU_Coord3(xint,yint,zint,XlR,1,nelR);   %Get the kesi eta in the right hand side
%                 rR = xi(1);
%                 sR = xi(2);
%                 tR = xi(3);
                
           if nelR == 4 || nelR == 10
                [w,ss] =  int3d_t(ie,lint,ib);
                ebeR = facebubbleT(ss,nelR);
           else
                [Wgt,ss] = intpntb(ie,lint,ib);
                ebeR = facebubbleQ(ss,nelR);
           end           
            
           if nelL == 4 || nelL == 10
                [w,ss] =  int3d_t(ie,lint,ib);
                [shlL,shldL,shls,be] = shltt(ss,nelL,nelL,0,0);
                [PxyL,shgsL,JdetL,bubbleL,XsL] = shgtt(xlL(:,1:nelL),nelL,shldL,shls,nen,0,0,be);               
           else
                [Wgt,ss] = intpntb(ie,lint,ib);
                [shlL,shldL,shls,be] = shlb(ss,nelL,nelL,0,0);
                [PxyL,shgsL,JdetL,bubbleL,XsL] = shgb(xlL(:,1:nelL),nelL,shldL,shls,nen,0,0,be); 
           end
                    
            %Evaluate tangent and normal vectors
                T1L = XsL(:,1);
                [Tm1L, Tu1L] = VecNormalize(T1L);
                T2L = XsL(:,2);
                [Tm2L, Tu2L] = VecNormalize(T2L);
                T3L = VecCrossProd(T1L,T2L);
                [Tm3L, Tu3L] = VecNormalize(T3L);
                
                C1L = Wgt*Tm3L;
            
            if nelL == 10 && ie == 13
                g = [Tu1L; Tu3L; Tu2L];
            Srot = mm10_RT2RVE(g);
            CnR3 = Srot*CnR*Srot';
            CnRrot = CnR3;
            CnL3 = Srot*CnL*Srot';
            CnLrot = CnL3;
            end
                
            ebL = ebL + C1L*ebeL;
            ebR = ebR + C1L*ebeR;
            intedge = intedge + C1L;

        end   
        
        edgeK = tauL*ebL^2 + tauR*ebR^2;
        gamL = ebL^2*(edgeK\tauL);
        gamR = ebR^2*(edgeK\tauR);
        hr(nh3:nh3+8) = [gamL(1,1),gamL(1,2),gamL(1,3),gamL(2,1),gamL(2,2),gamL(2,3),gamL(3,1),gamL(3,2),gamL(3,3)];
        hr(nh3+9:nh3+17) = [gamR(1,1),gamR(1,2),gamR(1,3),gamR(2,1),gamR(2,2),gamR(2,3),gamR(3,1),gamR(3,2),gamR(3,3)];
%         gamL_list(iter+1,:,inter) = [gamL(1,1),gamL(1,2),gamL(1,3),gamL(2,1),gamL(2,2),gamL(2,3),gamL(3,1),gamL(3,2),gamL(3,3)];
%         gamR_list(iter+1,:,inter) = [gamR(1,1),gamR(1,2),gamR(1,3),gamR(2,1),gamR(2,2),gamR(2,3),gamR(3,1),gamR(3,2),gamR(3,3)];
        ep = pencoeff*intedge*inv(edgeK); 
        hr(nh3+18:nh3+26) = [ep(1,1),ep(1,2),ep(1,3),ep(2,1),ep(2,2),ep(2,3),ep(3,1),ep(3,2),ep(3,3)];
        
        else % load taus from previous iterations
        gamL = [hr(nh3) hr(nh3+1) hr(nh3+2)
                hr(nh3+3) hr(nh3+4) hr(nh3+5)
                hr(nh3+6) hr(nh3+7) hr(nh3+8)];
        gamR = [hr(nh3+9) hr(nh3+10) hr(nh3+11)
                hr(nh3+12) hr(nh3+13) hr(nh3+14)
                hr(nh3+15) hr(nh3+16) hr(nh3+17)];
        ep =   [hr(nh3+18) hr(nh3+19) hr(nh3+20)
                hr(nh3+21) hr(nh3+22) hr(nh3+23)
                hr(nh3+24) hr(nh3+25) hr(nh3+26)];
        end %iterset
        
        else % Nitsche method, hard-coded penalty terms
            
        gamL = 0.5*eye(3);
        gamR = 0.5*eye(3);     
        ep = 20000*pencoeff*eye(3);
        
        end %nitvms
        

        rhspulL = reshape(ulL,ndf*nelL,1);
        rhspulR = reshape(ulR,ndf*nelR,1);
        rhspdulL = reshape(ulL-ul_nL,ndf*nelL,1);
        rhspdulR = reshape(ulR-ul_nR,ndf*nelR,1);
        % Main integration loop for DG residual/stiffness terms
        ll=0;           
        for l = 1:lint

                ll = ll + 1;

                % Evaluate  basis functions at integration points
                if nelL == 4 || nelL == 10
                  [Wgt,ss] = int3d_t(l,lint,ib);
                  [shlL,shldL,shls,be] = shltt(ss,nelL,nelL,0,0);             
                  [QxyL,shgsL,JdetL,bubbleL,xsL] = shgtt(xlL(:,1:nelL),nelL,shldL,shls,nen,0,0,be);
                  shlpL = shltt(ss,4,nelL,0,0);
                else
                  [Wgt,ss] = intpntb(l,lint,ib);
                  [shlL,shldL,shls,be] = shlb(ss,nelL,nelL,0,0);
                  [QxyL,shgsL,JdetL,bubbleL,xsL] = shgb(xlL(:,1:nelL),nelL,shldL,shls,nen,0,0,be);
                  shlpL = shlb(ss,8,nelL,0,0);
                end

                %Physical location of int pt
                xint = xlL(1,:)*shlL;
                yint = xlL(2,:)*shlL;
                zint = xlL(3,:)*shlL;

                xi = POU_Coord3(xint,yint,zint,xlR,1,nelR);   %Get the kesi eta in the right hand side
                rR = xi(1);
                sR = xi(2);
                tR = xi(3);

                % Evaluate  basis functions at integration points
                if nelR == 4 || nelR == 10
                  [shlR,shldR,shls,be] = shltt([rR sR tR],nelR,nelR,0,0);               
                  [QxyR,shgsR,JdetR,bubbleR,xsR] = shgtt(xlR(:,1:nelR),nelR,shldR,shls,nen,0,0,be);
                  shlpR = shltt(ss,4,nelR,0,0);
                else
                  [shlR,shldR,shls,be] = shlb([rR sR tR],nelR,nelR,0,0);    
                  [QxyR,shgsR,JdetR,bubbleR,xsR] = shgb(xlR(:,1:nelR),nelR,shldR,shls,nen,0,0,be);
                  shlpR = shlb(ss,4,nelR,0,0);
                end

 % initial stress sigma and material bulk term c_ijkl
           for mm = 1:nelL 
 NmatL(:,4*mm-3:4*mm) = [shlL(mm,1)     0          0 0
                           0        shlL(mm,1)     0 0
                           0            0       shlL(mm,1) 0];
 BmatL(:,4*mm-3:4*mm) = [QxyL(mm,1) 0         0         0
                        0         QxyL(mm,2) 0         0
                        0         0         QxyL(mm,3) 0
                        QxyL(mm,2) QxyL(mm,1) 0         0
                        0         QxyL(mm,3) QxyL(mm,2) 0
                        QxyL(mm,3) 0         QxyL(mm,1) 0
                        0 0 0 shlpL(mm)];  
           end
            
           for mm = 1:nelR    
 NmatR(:,4*mm-3:4*mm) = [shlR(mm,1)     0          0 0
                           0        shlR(mm,1)     0 0
                           0            0       shlR(mm,1) 0];
 BmatR(:,4*mm-3:4*mm) = [QxyR(mm,1) 0         0         0
                        0         QxyR(mm,2) 0         0
                        0         0         QxyR(mm,3) 0
                        QxyR(mm,2) QxyR(mm,1) 0         0
                        0         QxyR(mm,3) QxyR(mm,2) 0
                        QxyR(mm,3) 0         QxyR(mm,1) 0
                        0 0 0 shlpR(mm)];                
           end  
           
                
            %% Material Models
            % Extract history for left element along interface segment
            
                % Compute input for Radial Return
                du = BmatL*rhspulL; %2D strain
                pn_1L = du(7);
                eps3d = du(1:6); %3D enhanced strain
                ephr = nh1-1+26+(l-1)*13; %pointer for plastic strain at pt l
                betahr = nh1-1+26+(l-1)*13+6;
                ahr = nh1-1+26+l*13;
                ep_n = [hr(ephr+1); hr(ephr+2); hr(ephr+3); hr(ephr+4); hr(ephr+5); hr(ephr+6)];
                beta_n = [hr(betahr+1); hr(betahr+2); hr(betahr+3); hr(betahr+4); hr(betahr+5); hr(betahr+6)];
                a_n = hr(ahr);
                
                [sig_n1,C_n1,ep_n1,beta_n1,a_n1,sdev3,Cdev_n1] = J2RadialReturn0(eps3d,ep_n,beta_n,a_n,muL,bulkL,KhardL,HhardL,sigyL);
                 
                if modifDG == 0
                    if a_n1 > 0 && initia == 0
%                     [dmatL1]=dmat3_ep_numer(du,ep_n,beta_n,a_n,ElemYM,Elemv,Khard,Hhard,sigy,Cdev_n1,plasversion,plasmodel);
                    dmatL1 = J2Radial6thOrder(eps3d,ep_n,beta_n,a_n,mu,bulk,Khard,Hhard,sigy); %Analytical 6th order tensor
                    else
                        dmatL1 = zeros(36,6);
                    end
                end
                
                % Convert output from 3D to 2D
                sigma2L = sdev3;
                cmatL = Cdev_n1;
                
                % Store history variables
                ephr = nh2-1+26+(l-1)*13; %pointer for plastic strain at pt l
                betahr = nh2-1+26+(l-1)*13+6;
                ahr = nh2-1+26+l*13;
                hr(ephr+1:ephr+6) = ep_n1(1:6);
                hr(betahr+1:betahr+6) = beta_n1(1:6);
                hr(ahr) = a_n1;
            
                
            % Extract history for right element along interface segment

                % Compute input for Radial Return
                du = BmatR*rhspulR; %2D strain
                pn_1R = du(7);
                eps3d = du(1:6); %3D enhanced strain
                ephr = nh1-1+39+(l-1)*13; %pointer for plastic strain at pt l
                betahr = nh1-1+39+(l-1)*13+6;
                ahr = nh1-1+39+l*13;
                ep_n = [hr(ephr+1); hr(ephr+2); hr(ephr+3); hr(ephr+4); hr(ephr+5); hr(ephr+6)];
                beta_n = [hr(betahr+1); hr(betahr+2); hr(betahr+3); hr(betahr+4); hr(betahr+5); hr(betahr+6)];
                a_n = hr(ahr);
                
                [sig_n1,C_n1,ep_n1,beta_n1,a_n1,sdev3,Cdev_n1] = J2RadialReturn0(eps3d,ep_n,beta_n,a_n,muR,bulkR,KhardR,HhardR,sigyR);
                
                if modifDG == 0 
                    if a_n1 > 0 && initia == 0
%                     [dmatR1]=dmat2_ep_numer(du,ep_n,beta_n,a_n,ElemYM,Elemv,Khard,Hhard,sigy,Cdev_n1,plasversion,plasmodel);
                    dmatR1 = J2Radial6thOrder(eps3d,ep_n,beta_n,a_n,mu,bulk,Khard,Hhard,sigy);
%                     dmatR1 = D22([1 2 4 7 8 10 19 20 22],[1 2 4]);
                    else
                        dmatR1 = zeros(36,6);
                    end
                end
                
                % Convert output from 3D to 2D
                sigma2R = sdev3;
                cmatR = Cdev_n1;
                
                % Store history variables
                ephr = nh2-1+39+(l-1)*13; %pointer for plastic strain at pt l
                betahr = nh2-1+39+(l-1)*13+6;
                ahr = nh2-1+39+l*13;
                hr(ephr+1:ephr+6) = ep_n1(1:6);
                hr(betahr+1:betahr+6) = beta_n1(1:6);
                hr(ahr) = a_n1;
                
            %% End Material Models
          
            
            % Cauchy stress tensors for both sides
            SmatL1=[sigma2L(1), sigma2L(4), sigma2L(6)
                    sigma2L(4), sigma2L(2), sigma2L(5)
                    sigma2L(6), sigma2L(5), sigma2L(3)]+pn_1L*eye(3);
            SmatR1=[sigma2R(1), sigma2R(4), sigma2R(6)
                    sigma2R(4), sigma2R(2), sigma2R(5)
                    sigma2R(6), sigma2R(5), sigma2R(3)]+pn_1R*eye(3);
             
            % Evaluate tangent and normal vectors
            T1L = XsL(:,1);
            [Tm1L, Tu1L] = VecNormalize(T1L);
            T2L = XsL(:,2);
            [Tm2L, Tu2L] = VecNormalize(T2L);
            T3L = VecCrossProd(T1L,T2L);
            [Tm3L, Tu3L] = VecNormalize(T3L);
            C1L = Wgt*Tm3L;
            C1R = C1L;
            
            % Copy into "current" configuration placeholders
            tu3L = Tu3L;
            tu3R = -Tu3L;
            c1L = C1L;
            c1R = C1L;

            %normal vectors
                nLx = -tu3L(1);
                nLy = -tu3L(2);
                nLz = -tu3L(3);
                nRx = -tu3R(1);
                nRy = -tu3R(2);
                nRz = -tu3R(3);          
            nvectL1= [nLx 0  0  nLy  0  nLz 
                      0  nLy 0  nLx nLz  0 
                      0   0 nLz  0  nLy nLx];
            nvectR1 = [nRx 0  0  nRy  0  nRz 
                      0  nRy 0  nRx  nRz 0 
                      0   0 nRz  0  nRy nRx];   
            nvecL = [nLx; nLy; nLz];
            nvecR = [nRx; nRy; nRz];


            % stiffness matrix terms for weighted stress and material moduli
            term18L=(gamL*nvectL1*[cmatL One])';
            term18R=(gamR*nvectR1*[cmatR One])';
            
            
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
             dmatL2 = [reshape(dmatL2,6,6) zeros(6,1); zeros(1,7)];
             dmatR2 = dmatR1*nvectR1'*gamR'*jumpu;
             dmatR2 = [reshape(dmatR2,6,6) zeros(6,1); zeros(1,7)];
             
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
%%
    case 26 % Element stresses: output stability tensors
        
        ElemS = zeros(nestr,1);
        
        gamL = [hr(nh3) hr(nh3+1) hr(nh3+2)
                hr(nh3+3) hr(nh3+4) hr(nh3+5)
                hr(nh3+6) hr(nh3+7) hr(nh3+8)];
        gamR = [hr(nh3+9) hr(nh3+10) hr(nh3+11)
                hr(nh3+12) hr(nh3+13) hr(nh3+14)
                hr(nh3+15) hr(nh3+16) hr(nh3+17)];
        ep =   [hr(nh3+18) hr(nh3+19) hr(nh3+20)
                hr(nh3+21) hr(nh3+22) hr(nh3+23)
                hr(nh3+24) hr(nh3+25) hr(nh3+26)];
        eigCnL = hr(nh3+27:nh3+32);
        eigCnR = hr(nh3+33:nh3+38);
        
        if exist('iprob','var') && iprob == 1
        % Make the last eigenvalue always be the one that is not zero and
        % not 18.666667, for the elastic and the plastic versions
        if abs(eigCnL(1)) > 1e-9
            eigCnL([3 2]) = eigCnL([2 3]);
        else
            if eigCnL(2) > 18.67 || eigCnL(2) < 18.66
                eigCnL([3 2]) = eigCnL([2 3]);
            end
        end
        if abs(eigCnR(1)) > 1e-9
            eigCnR([3 2]) = eigCnR([2 3]);
        else
            if eigCnR(2) > 18.67 || eigCnR(2) < 18.66
                eigCnR([3 2]) = eigCnR([2 3]);
            end
        end
        end
        
            
            
        % Determin unit normal at segment midpoint
            
        % Determin bounds of integration segment
        InterConn2D2 % InterConn2DT %
        
        % Set jacobian for integration space
        drdr = (eL2 - eL1)/drL;
        drdrR = (eR2 - eR1)/drR; % This causes errors for partial interface segments
        m = (eR2-eR1)/(eL1-eL2);
        
        if nelL == 3 || nelL == 6
            [Wgt,litr,lits] = intpntt(1,1,1);
        elseif nelL == 4 || nelL == 9
            [Wgt,litr,lits] = intpntq(1,1,1);
        end
            
        rL = drdr*(litr-roL)+eL1;
       if nelL == 3 || nelL == 6
       [shlL,shldL,shlsL,be] = shlt(rL,lits,nelL,nelL,0,0);
       [PxyL,shgsL,JdetL,bubbleL,XsL] = shgt(xlL(:,1:nelL),nelL,shldL,shlsL,nen,0,0,be);
       else
       [shlL,shldL,shlsL,be] = shlq(rL,lits,nelL,nelL,0,0);
       [PxyL,shgsL,JdetL,bubbleL,XsL] = shgq(xlL(:,1:nelL),nelL,shldL,shlsL,nen,0,0,be);      
       end 
        rR = m*(rL-eL2) + eR1;

        if nelR == 3 || nelR == 6
            sR = 0;
        else %if nelR == 4
            sR = -1;
        end
       if nelR == 3 || nelR == 6            
        [shlR,shldR,shlsR,be] = shlt(rR,sR,nelR,nelR,0,0);
        [PxyR,shgsR,JdetR,bubbleR,XsR] = shgt(xlR(:,1:nelR),nelR,shldR,shlsR,nen,0,0,be);     
       else
        [shlR,shldR,shlsR,be] = shlq(rR,sR,nelR,nelR,0,0);
        [PxyR,shgsR,JdetR,bubbleR,XsR] = shgq(xlR(:,1:nelR),nelR,shldR,shlsR,nen,0,0,be);  
       end

         %Evaluate tangent and normal vectors
        T1L = [XsL(:,1); 0];
        [Tm1L, Tu1L] = VecNormalize(T1L);
        T2L = [0; 0; 1];
        Tm2L = 1;
        Tu2L = T2L';
        T3L = VecCrossProd(T1L,T2L);
        [tm3L, tu3L] = VecNormalize(T3L);
        tm3R = tm3L;
        tu3R = -tu3L;

        %normal vectors
        nLx = tu3L(1);
        nLy = tu3L(2);                             
        nvecL = [nLx; nLy];
        tvecL = [-nLy; nLx];
        nRx = tu3R(1);
        nRy = tu3R(2);                             
        nvecR = [nRx; nRy];
        tvecR = [-nRy; nRx];
        
        % Penalty tensor
        ep_nn = nvecL'*ep*nvecL;
        ep_tt = tvecL'*ep*tvecL;
        ep_nt = tvecL'*ep*nvecL;
        
        tauL = pencoeff*inv(ep)*gamL;
        tauR = pencoeff*inv(ep)*gamR;
        
        % Stability tensor left side
        tL_nn = nvecL'*tauL*nvecL;
        tL_tt = tvecL'*tauL*tvecL;
        tL_nt = tvecL'*tauL*nvecL;
        
        % Stability tensor right side
        tR_nn = nvecR'*tauR*nvecR;
        tR_tt = tvecR'*tauR*tvecR;
        tR_nt = tvecR'*tauR*nvecR;
        
        ElemS(1:4+4+4+3+6+6) = [reshape(ep,4,1); reshape(gamL,4,1); reshape(gamR,4,1); ep_nn; ep_tt; ep_nt; ...
            tL_nn; tL_tt; tL_nt; tR_nn; tR_tt; tR_nt; eigCnL; eigCnR];
            
        
    case 61 % form data structure for interface segments
        
        % The results for DCB3C9U3DGIc at step 1 or 10 do not seem correct
        % or meaningful. Thus, something needs to be fixed.
        % For this portion to work, FormIData -> nenseg needs to be set to
        % 3.

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
              shlpL = shlq(rL,lits,4,nelL,0,0);
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
              shlpR = shlq(rR,sR,4,nelR,0,0);        
          
            
             
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
            if plasversion == 1
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
                
                if plasmodel == 1
                [DGAMA,LALGVA,RSTAVA,STRES] = SUVM(2,2,[ElemYML,ElemvL],[0,100;sigyL,KhardL*100+sigyL],[ee_n; a_n],eps_n1);
                elseif plasmodel == 2
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
            
            %%%%%%%%%%% End Left side
            
            
            
            %%%%%%%%%%% Start Right side
                
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
            if plasversion == 1
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
                
                if plasmodel == 1
                [DGAMA,LALGVA,RSTAVA,STRES] = SUVM(2,2,[ElemYMR,ElemvR],[0,100;sigyR,KhardR*100+sigyR],[ee_n; a_n],eps_n1);
                elseif plasmodel == 2
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