%% DG implementation of large deformaiton 
% Mixed method 2D case
% Tim Truster
% 10/2011
% modified by Pinlei Chen
% 10/13/2013
% for 2D large deformation with interface in it
% have tau and delta in it
% have d_ijklmn in it 
% verified for the noncomforming mesh of patch test
% for body force problem, change the lint
% UIUC
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


nitvms = 1;
if nitvms == 1 %VMS parameter for the stability tensor rp
pencoeff = 1;1;
elseif nitvms == 2 %Nitsche
pencoeff = 1;4;2;
else %RFB
pencoeff = 1;4;2;
end


switch isw %Task Switch
%%
    case 1
        
        if ndf > 3  % 3 for mixed method in 2D
            
            for i = 4:ndf
                lie(i,1) = 0;
            end
            
        end
        
    case 3 %interface stiffness
       % generate the initial parameter 
        ElemKRR = zeros(nstR,nstR);
        ElemKRL = zeros(nstR,nstL);
        ElemKLR = zeros(nstL,nstR);
        ElemKLL = zeros(nstL,nstL);
        ElemFR = zeros(nstR,1);
        ElemFL = zeros(nstL,1);
        
        t1 = zeros(3,1);
        t2 = t1;
        t3 = t1;           
        dtol = 1e-11;
        itchat = 100;12;8; %flag for iteration number to suppress chatter
        
        I1 = [1; 1; 0; 0];
        I2 = eye(2);
        I6 = eye(4);
        P1 = [1 0 0 0
              0 1 0 0
              0 0 1 0];
        P2 = [1 0 0 0 
              0 0 1/2 1/2
              0 0 1/2 -1/2
              0 1 0 0];
        P3 = [1 0 0 0 
              0 0 1/2 -1/2
              0 0 1/2 1/2
              0 1 0 0]; 
        Z2 = zeros(2);
        spvec0 = I1;
        spmat0 = I2; % s_cup
        cpmat1 = I1*I1';
        cpmat2 = diag([-2,-2,-1,0]);
        cpmat0 = cpmat1 + cpmat2;  %c_cup
%         dpmat1 = [cpmat1; cpmat1; zeros(4,4)];
%         dpmat2 =-[6 2 0 0
%                   2 2 0 0
%                   0 0 1 0
%                   0 0 0 0
%                   2 2 0 0
%                   2 6 0 0
%                   0 0 1 0
%                   0 0 0 0
%                   0 0 1 0
%                   0 0 1 0
%                   1 1 0 0
%                   0 0 0 0];
%         dpmat3 = [8 0 0 0
%                   0 0 0 0
%                   0 0 2 0
%                   0 0 0 0
%                   0 0 0 0
%                   0 8 0 0
%                   0 0 2 0
%                   0 0 0 0
%                   0 0 2 0
%                   0 0 2 0
%                   2 2 0 0
%                   0 0 0 0];
        I1d = [1; 1; 0];
        cpmat1d = I1d*I1d';
        dpmat1 = [cpmat1d; cpmat1d; zeros(3,3)];
        dpmat2 =-[6 2 0 
                  2 2 0 
                  0 0 1 
                  2 2 0 
                  2 6 0 
                  0 0 1 
                  0 0 1 
                  0 0 1 
                  1 1 0 ];
        dpmat3 = [8 0 0 
                  0 0 0 
                  0 0 2 
                  0 0 0 
                  0 8 0 
                  0 0 2 
                  0 0 2 
                  0 0 2 
                  2 2 0 ];
        dpmat0 = dpmat1 + dpmat2 + dpmat3;
        ElemEL = matepropL(1);
        ElemvL = matepropL(2);
        ElemER = matepropR(1);
        ElemvR = matepropR(2);
        lam = getlam(matepropL); %For both material is the same at each side
       
        NmatL = zeros(3,nstL);
        BmatL = zeros(5,nstL);
%         BmatL1 = zeros(3,nstL);        
%         BmatL2 = zeros(4,nstL);
%         BmatL3 = zeros(4,nstL);       
        NmatR = zeros(3,nstR);
        BmatR = zeros(5,nstR);
%         BmatR1 = zeros(3,nstR);   
%         BmatR2 = zeros(4,nstR);
%         BmatR3 = zeros(4,nstL);        
        dmatL = zeros(9,3);
        dmatR = zeros(9,3);
     
        % Determin bounds of integration segment
        InterConn2D2 % InterConn2DT % 
        
        % Set jacobian for integration space
        drdr = (eL2 - eL1)/drL;
        drdrR = (eR2 - eR1)/drR;
        m = (eR2-eR1)/(eL1-eL2);
        if nelL == 3 || nelL == 6  
        lint = 3;10;2;3;
        else
         if exist('iprob','var') == 1 && iprob == 6
             lint = 10;4;
         else
             lint = 3;
         end    
         %10 for body force problem; 4 for other problem 
        end
        ideriv = 0;
     


%Form tau matrix: compute on first iteration and store; this gives
        %quadratic convergence when tau isn't recalculated (no
        %linearization of tau required).
   iterset = 1;4;
    if iter  <=iterset % == 0 %
        [tauL,intb] = TauS2_M(xlL,ulL,matepropL,nelL,nelLP,nen,lam,roL,eL1,drdr); %[Y^(-1)]
%        [tauL,intb] = TauS2_1(xlL,ulL,mateprop,nelL,nel2L,nen,lam,roL,eL1,drdr); %[Y^(-1)]        
%       TauListL(elem,:) = [tauL(1),tauL(2),tauL(3),tauL(4),tauL(5),tauL(6),tauL(7),tauL(8),tauL(9)];
        [tauR,intb] = TauS2_M(xlR,ulR,matepropR,nelR,nelRP,nen,lam,roL,eL1,drdrR);
%        [tauR,intb] = TauS2_1(xlR,ulR,mateprop,nelR,nel2R,nen,lam,roL,eL1,drdrR);
%          TauListR(elem,:) = [tauR(1),tauR(2),tauR(3),tauR(4),tauR(5),tauR(6),tauR(7),tauR(8),tauR(9)];

%         tau = tauL + tauR;

        ebL = 0;
        ebR = 0;
        intedge = 0;
        
        for ie = 1:lint            
% For separate bubble types on T and Q
           if nelL == 3 || nelL == 6
                [Wgt,litr,lits] = intpntt(ie,lint,1); 
                 ebeL = edgebubble(litr,lits,nelL);  %edgebubble is for T3 element
%                 rL = drdr*(litr-roL)+eL1;
%                 [shlL,shldL,shlsL,be] = shlt(rL,lits,nelL,nel2L,0,0);
           else
                [Wgt,litr,lits] = intpntq(ie,lint,1);
%                 rL = drdr*(litr-roL)+eL1;
%                 [shlL,shldL,shlsL,be] = shlq(rL,lits,nelL,nel2L,0,0);
                 ebeL = edgebubbleQ(litr,lits,nelL);
           end
                            
           if nelR == 3 || nelR == 6
                [Wgt,litr,lits] = intpntt(ie,lint,1); 
                 ebeR = edgebubble(litr,lits,nelL);
%                 rR = drdr*(litr-roL)+eR1;
%                 [shlR,shldR,shlsR,be] = shlt(rR,lits,nelR,nel2R,0,0);
           else
                [Wgt,litr,lits] = intpntq(ie,lint,1);
%                 rR = drdr*(litr-roL)+eR1;
%                 [shlR,shldR,shlsR,be] = shlq(rR,lits,nelR,nel2R,0,0);
                 ebeR = edgebubbleQ(litr,lits,nelR);
           end         
                    
%            b = edgebubble(litr,0);
           if nelL == 3 || nelL == 6
               rL = drdr*(litr-roL)+eL1;
               [shlL,shldL,shlsL,be] = shlt(rL,lits,nelL,nelL,0,0);
               [PxyL,shgsL,JdetL,bubbleL,XsL] = shgt(xlL(:,1:nelL),nelL,shldL,shlsL,nen,0,0,be); 
           else    
               rL = drdr*(litr-roL)+eL1;
               [shlL,shldL,shlsL,be] = shlq(rL,lits,nelL,nelL,0,0);
               [PxyL,shgsL,JdetL,bubbleL,XsL] = shgq(xlL(:,1:nelL),nelL,shldL,shlsL,nen,0,0,be); 
           end
%            [PxyL,shgsL,JdetL,bubbleL,xsL] = shgq(xlL(:,1:nelL)+ulL(:,1:nelL),nelL,shldL,shlsL,nen,0,0,be);
%             QxyL = PxyL;  
            
            
%             rR = m*(rL-eL2) + eR1;
%             
%             if nelR == 3 || nelR == 6
%                 sR = 0;
%             else %if nelR == 4
%                 sR = -1;
%             end
%             
%             [shlR,shldR,shlsR,be] = shlq(rR,sR,nelR,nel2R,0,0);
%             [PxyR,shgsR,JdetR,bubbleR,XsR] = shgq(xlR(:,1:nelR),nelR,shldR,shlsR,nen,0,0,be);     
% %             [PxyR,shgsR,JdetR,bubbleR,xsR] = shgq(xlR(:,1:nelR)+ulR(:,1:nelR),nelR,shldR,shlsR,nen,0,0,be);           
% %             QxyR = PxyR;
             
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
        
        edgeK = tauL*ebL^2 + tauR*ebR^2;
        gamL = ebL^2*(edgeK\tauL);
        gamL_list(iter+1,:,inter) = [gamL(1,1),gamL(1,2),gamL(2,1),gamL(2,2)];
        gamR = ebR^2*(edgeK\tauR);
        gamR_list(iter+1,:,inter) = [gamR(1,1),gamR(1,2),gamR(2,1),gamR(2,2)];
%         gamL_list(elem,:) = [gamL(1),gamL(2),gamL(3),gamL(4),gamL(5),gamL(6),gamL(7),gamL(8),gamL(9)];
%         gamR_list(elem,:) = [gamR(1),gamR(2),gamR(3),gamR(4),gamR(5),gamR(6),gamR(7),gamR(8),gamR(9)];
        ep = pencoeff*intedge*inv(edgeK); 
        ep_list(iter+1,:,inter) = [ep(1,1),ep(1,2),ep(2,1),ep(2,2)];
%         ep_List(elem,:) = [ep(1),ep(2),ep(3),ep(4),ep(5),ep(6),ep(7),ep(8),ep(9)];
%       ep_List_R(elem,:) = [tauR(1),tauR(2),tauR(3),tauR(4),tauR(5),tauR(6),tauR(7),tauR(8),tauR(9)];         
   else
        gamL = [gamL_list(iterset+1,1,inter) gamL_list(iterset+1,2,inter)
                gamL_list(iterset+1,3,inter) gamL_list(iterset+1,4,inter)];
        gamR = [gamR_list(iterset+1,1,inter) gamR_list(iterset+1,2,inter)
                gamR_list(iterset+1,3,inter) gamR_list(iterset+1,4,inter)];
        ep = [ep_list(iterset+1,1,inter) ep_list(iterset+1,2,inter)
              ep_list(iterset+1,3,inter) ep_list(iterset+1,4,inter)];
   end
%        gamL = [0.5 0
%                 0 0.5];
%        gamR = [0.5 0
%                 0 0.5];     
%        ep = 20*eye(2);
   
%        s = -1;
        ll=0;           
     for l = 1:lint

            ll = ll + 1;
            if nelL == 3 || nelL == 6
                [Wgt,litr,lits] = intpntt(ll,lint,1);
            elseif nelL == 4 || nelL == 9
                [Wgt,litr,lits] = intpntq(ll,lint,1);
            end
                    
%            b = edgebubble(litr,0);
            
            rL = drdr*(litr-roL)+eL1;
           if nelL == 3 || nelL == 6
           [shlL,shldL,shlsL,be] = shlt(rL,lits,nelL,nelL,0,0);
           [QxyL,shgsL,JdetL,bubbleL,XsL] = shgt(xlL(:,1:nelL),nelL,shldL,shlsL,nen,0,0,be); 
           [QxyL,shgsL,JdetL,bubbleL,xsL] = shgt(xlL(:,1:nelL)+ulL(1:2,1:nelL),nelL,shldL,shlsL,nen,0,0,be);
           else
           [shlL,shldL,shlsL,be] = shlq(rL,lits,nelL,nelL,0,0);
           [QxyL,shgsL,JdetL,bubbleL,XsL] = shgq(xlL(:,1:nelL),nelL,shldL,shlsL,nen,0,0,be); 
           [QxyL,shgsL,JdetL,bubbleL,xsL] = shgq(xlL(:,1:nelL)+ulL(1:2,1:nelL),nelL,shldL,shlsL,nen,0,0,be);      
           end
% pressure field
             if nelLP == 3 || nelLP == 6
              [shlpL,shldL,shlsL] = shlt(rL,lits,nelLP,nelL,0,0); %shlpL:shape fun
              [PxyL, shgpL] = shgt(xlL(:,1:nelL)+ulL(1:2,1:nelL),nelLP,shldL,shlsL,nelLP,0,0,be); %PxyL:deriv of shape fun for pressure field
            else
              [shlpL,shldL,shlsL] = shlq(rL,lits,nelLP,nelL,0,0);
              [PxyL, shgpL] = shgq(xlL(:,1:nelL)+ulL(1:2,1:nelL),nelLP,shldL,shlsL,nelLP,0,0,be);
            end           
            
            rR = m*(rL-eL2) + eR1;
            
            if nelR == 3 || nelR == 6
                sR = 0;
            else %if nelR == 4
                sR = -1;
            end
           if nelR == 3 || nelR == 6            
            [shlR,shldR,shlsR,be] = shlt(rR,sR,nelR,nelR,0,0);
            [QxyR,shgsR,JdetR,bubbleR,XsR] = shgt(xlR(:,1:nelR),nelR,shldR,shlsR,nen,0,0,be);     
            [QxyR,shgsR,JdetR,bubbleR,xsR] = shgt(xlR(:,1:nelR)+ulR(1:2,1:nelR),nelR,shldR,shlsR,nen,0,0,be);  
           else
            [shlR,shldR,shlsR,be] = shlq(rR,sR,nelR,nelR,0,0);
            [QxyR,shgsR,JdetR,bubbleR,XsR] = shgq(xlR(:,1:nelR),nelR,shldR,shlsR,nen,0,0,be);     
            [QxyR,shgsR,JdetR,bubbleR,xsR] = shgq(xlR(:,1:nelR)+ulR(1:2,1:nelR),nelR,shldR,shlsR,nen,0,0,be);  
           end

             if nelRP == 3 || nelRP == 6
              [shlpR,shldR,shlsR] = shlt(rR,sR,nelRP,nelR,0,0); %shlpL:shape fun
              [PxyR, shgpR] = shgt(xlR(:,1:nelR)+ulR(1:2,1:nelR),nelRP,shldR,shlsR,nelRP,0,0,be); %PxyL:deriv of shape fun for pressure field
            else
              [shlpR,shldR,shlsR] = shlq(rR,sR,nelRP,nelR,0,0);
              [PxyR, shgpR] = shgq(xlR(:,1:nelR)+ulR(1:2,1:nelR),nelRP,shldR,shlsR,nelRP,0,0,be);
            end  
            
 % initial stress sigma and material bulk term c_ijkl
           for mm = 1:nelL    
%            NmatL(:,2*mm-1:2*mm) = [shlL(mm,1)     0          
%                                     0        shlL(mm,1)];                
           NmatL(:,3*mm-2:3*mm) = [shlL(mm,1)     0       0         
                                    0        shlL(mm,1)   0
                                    0             0     shlpL(mm,1)]; 
%            BmatL(:,2*mm-1:2*mm) = [QxyL(mm,1) 0             
%                                     0         QxyL(mm,2)         
%                                     QxyL(mm,2) QxyL(mm,1)        
%                                      QxyL(mm,2) -QxyL(mm,1)]; 
%            PmatL(:,mm) = shlpL(mm,1) ;  %shape function of pressure
           BmatL(:,3*mm-2:3*mm) = [QxyL(mm,1) 0             0        
                                    0         QxyL(mm,2)    0      
                                    QxyL(mm,2) QxyL(mm,1)   0      
                                     QxyL(mm,2) -QxyL(mm,1) 0
                                     0          0           shlpL(mm,1)]; %similar to Bmat [derivative of shape disp + shape of pressure]
%            BBmatL(:,3*mm-2:3*mm) = [shgsL(1,1) 0         0
%                                     shgsL(1,2) 0         0
%                                     shgsL(1,3) 0         0
%                                      0         shgsL(1,1) 0
%                                      0         shgsL(1,2) 0
%                                      0         shgsL(1,3) 0 
%                                      0         0         PxyL(1,1)
%                                      0         0         PxyL(1,2)]; %%similar to BBmat [second derivative of shape disp + derivative of shape of pressure]
%            BmatL1(:,3*mm-2:3*mm) = [QxyL(mm,1) 0                  
%                                      0         QxyL(mm,2)          
%                                      QxyL(mm,2) QxyL(mm,1)];
%             BmatL2(:,3*mm-2:3*mm)=[QxyL(mm,1)  0  0                   
%                                   QxyL(mm,2)  0   0                                   
%                                     0     QxyL(mm,1) 0                
%                                     0     QxyL(mm,2) 0
%                                     0          0     shlpL(mm,1)  ];   
%            BmatL3(:,3*mm-2:3*mm)=[QxyL(mm,1)  0                     
%                                    0     QxyL(mm,1)                                      
%                                    QxyL(mm,2) 0                     
%                                    0     QxyL(mm,2)];     
           end 
            
           for mm = 1:nelR    
%            NmatR(:,2*mm-1:2*mm) = [shlR(mm,1)     0          
%                                      0        shlR(mm,1)];                
           NmatR(:,3*mm-2:3*mm) = [shlR(mm,1)     0       0        
                                     0        shlR(mm,1)  0
                                     0            0       shlpR(mm,1)]; 
%            BmatR(:,2*mm-1:2*mm) = [QxyR(mm,1) 0             
%                                     0         QxyR(mm,2)         
%                                     QxyR(mm,2) QxyR(mm,1)        
%                                      QxyR(mm,2) -QxyR(mm,1)];   
%            PmatR(:,mm) = shlpR(mm,1) ;  %shape function of pressure                                 
           BmatR(:,3*mm-2:3*mm) = [QxyR(mm,1) 0              0  
                                    0         QxyR(mm,2)     0    
                                    QxyR(mm,2) QxyR(mm,1)    0    
                                     QxyR(mm,2) -QxyR(mm,1)  0
                                     0           0           shlpR(mm,1)];
                               
%            BmatR1(:,3*mm-2:3*mm) = [QxyR(mm,1) 0                  
%                                      0         QxyR(mm,2)          
%                                      QxyR(mm,2) QxyR(mm,1)];
%            BmatR2(:,3*mm-2:3*mm)=[QxyR(mm,1)  0                     
%                                  QxyR(mm,2)  0                                      
%                                    0     QxyR(mm,1)                
%                                    0     QxyR(mm,2)];   
%            BmatR3(:,3*mm-2:3*mm)=[QxyR(mm,1)  0                     
%                                    0     QxyR(mm,1)                                      
%                                    QxyR(mm,2) 0                     
%                                    0     QxyR(mm,2)];                       
           end  
    %plug in the -ul to get the inverse of deformation gradient fi as 'F' in this function             
            [fiL,JxXL,FL] = kine2d(PxyL,-ulL,nelL,0); %this is equivalent to ikine2d
            JxXL = 1/JxXL; %this is equivalent to ikine2d
    %         [F,JxX,fi] = ikine2d(Qxy,ul,nel); %this is equivalent to above two lines
            JdetL = JdetL/JxXL;            
            [theta1L,theta2L,theta3L] = ThetaNS(JxXL,matepropL);
            
            [fiR,JxXR,FR] = kine2d(PxyR,-ulR,nelR,0); %this is equivalent to ikine2d
            JxXR = 1/JxXR; %this is equivalent to ikine2d
    %         [F,JxX,fi] = ikine2d(Qxy,ul,nel); %this is equivalent to above two lines
            JdetR = JdetR/JxXR;               
            [theta1R,theta2R,theta3R] = ThetaNS(JxXR,matepropR);   
            
            spvecL = JxXL*spvec0;
%             spmatL = JxXL*spmat0;
            cpmatL = JxXL*cpmat0;
%             cpmat32L = (theta3*JxXL^2 + theta2*JxXL)*cpmat1 + theta2*JxXL*cpmat2;
             dpmatL = JxXL*dpmat0;            
          
         [sigmaiL, cmatiL] = SigmaCmatNSCST2i(FL,JxXL,matepropL); 
           [dmatiL] = dmat2i(matepropL);                     
            pressL = ulL(3,:)*shlpL;
            sigmapL = pressL*spvecL;
            cmatpL = pressL*cpmatL;
            dmatpL = pressL*dpmatL;
            sigma2L = (sigmaiL + sigmapL)./JxXL;  % [sigma11 sigma22 sigma12 0]
            cmatL = (cmatiL + cmatpL)./JxXL;      % [last one is zero]      
            dmatL1 = (dmatpL - dmatiL)./JxXL;
            SmatL1=[sigma2L(1), sigma2L(3)
                    sigma2L(3), sigma2L(2)];           
            
%            [sigma2R, cmatR] = SigmaCmat2i(FR,JxXR,matepropR,lam);
            spvecR = JxXR*spvec0;
%             spmatR = JxXR*spmat0;
            cpmatR = JxXR*cpmat0;
%             cpmat32R = (theta3*JxXR^2 + theta2*JxXR)*cpmat1 + theta2*JxXR*cpmat2;
             dpmatR = JxXR*dpmat0; 
            
         [sigmaiR, cmatiR] = SigmaCmatNSCST2i(FR,JxXR,matepropR); 
          [dmatiR] = dmat2i(matepropR);          
            pressR = ulR(3,:)*shlpR;
            sigmapR = pressR*spvecR;
            cmatpR = pressR*cpmatR;
            dmatpR = pressR*dpmatR;
            sigma2R = (sigmaiR + sigmapR)./JxXR;
            cmatR = (cmatiR + cmatpR)./JxXR;  
            dmatR1 = (dmatpR - dmatiR)./JxXR;           
            SmatR1=[sigma2R(1), sigma2R(3)
                    sigma2R(3), sigma2R(2)];     

             %Evaluate tangent and normal vectors
            t1L = [xsL(:,1); 0];
            [tm1L, tu1L] = VecNormalize(t1L);
            t2L = [0; 0; 1];
            tm2L = 1;
            tu2L = t2L';
            t3L = VecCrossProd(t1L,t2L);
            [tm3L, tu3L] = VecNormalize(t3L);

                %Evaluate tangent and normal vectors
             t1R = [xsR(:,1); 0];
            [tm1R, tu1R] = VecNormalize(t1R);
            t2R = [0; 0; 1];
            tm2R = 1;
            tu2R = t2R';
            t3R = VecCrossProd(t1R,t2R);
            [tm3R, tu3R] = VecNormalize(t3R);
             c1L = drdr*Wgt*tm3L;
             c1R = drdrR*Wgt*tm3R; 
%              c1L = c1R;
             
             %Evaluate tangent and normal vectors
            T1L = [XsL(:,1); 0];
            [Tm1L, Tu1L] = VecNormalize(T1L);
            T2L = [0; 0; 1];
            Tm2L = 1;
            Tu2L = T2L';
            T3L = VecCrossProd(T1L,T2L);
            [Tm3L, Tu3L] = VecNormalize(T3L);

                %Evaluate tangent and normal vectors
             T1R = [XsR(:,1); 0];
            [Tm1R, Tu1R] = VecNormalize(T1R);
            T2R = [0; 0; 1];
            Tm2R = 1;
            Tu2R = T2R';
            T3R = VecCrossProd(T1R,T2R);
            [Tm3R, Tu3R] = VecNormalize(T3R);
             C1L = drdr*Wgt*Tm3L;
             C1R = drdrR*Wgt*Tm3R;  
%              C1L = C1R;
                %normal vectors
                nLx = tu3L(1);
                nLy = tu3L(2);
                nRx = tu3R(1);
                nRy = tu3R(2);
                %tagent vectors
                tLx = tu1L(1);
                tLy = tu1L(2);
                tRx = tu1R(1);
                tRy = tu1R(2);             
                nvectL1= [nLx 0   nLy   
                          0  nLy  nLx ];
                nvectR1 = [nRx 0  nRy  
                          0  nRy  nRx ]; 
                nvectL2 = [eye(3,3)*nLx zeros(3,3) eye(3,3)*nLy
                           zeros(3,3)    eye(3,3)*nLy eye(3,3)*nLx];
                nvectR2 = [eye(3,3)*nRx zeros(3,3) eye(3,3)*nRy
                           zeros(3,3)    eye(3,3)*nRy eye(3,3)*nRx];                       
                nvecL = [nLx; nLy];
                nvecR = [nRx; nRy];


             rhspulL = reshape(ulL,ndf*nelL,1);
             rhspulR = reshape(ulR,ndf*nelR,1);
                    
            jumpu = NmatL*rhspulL - NmatR*rhspulR;  %jumpu=uL-uR
%% right hand side force terms           
           
           SmatnL=[SmatL1*nvecL zeros(ndm,1)
                    zeros(ndm,1) SmatL1*nvecL];
            
            SmatnR=[SmatR1*nvecR zeros(ndm,1)
                    zeros(ndm,1) SmatR1*nvecR];
           cmatnL=(nvectL1*cmatL(1:3,1:3));
           cmatnR=(nvectR1*cmatR(1:3,1:3));
           cmatnBL=BmatL(1:4,:)'*P2'*[cmatnL zeros(2,3)
                            zeros(2,3)  cmatnL ];               
           cmatnBR=BmatR(1:4,:)'*P2'*[cmatnR    zeros(2,3)
                            zeros(2,3)  cmatnR ];

            term17L=P2'*SmatnL*gamL;
            term17R=P2'*SmatnR*gamR;

            term18L=P1'*(gamL*nvectL1*cmatL(1:3,1:3))';
            term18R=P1'*(gamR*nvectR1*cmatR(1:3,1:3))';

            term28L=NmatL(1:2,:)'*(c1L*gamL*SmatL1*nvecL-c1R*gamR*SmatR1*nvecR);   %average stress term
            term28R=NmatR(1:2,:)'*(c1L*gamL*SmatL1*nvecL-c1R*gamR*SmatR1*nvecR);

% additional terms related to the pressure:
            Sigma2UL = theta2L*spmat0; % sigma^U=U''(J)*I  
            Sigma2UR = theta2R*spmat0; % sigma^U=U''(J)*I 
            SmatUnL =  [Sigma2UL*nvecL zeros(ndm,1)
                        zeros(ndm,1) Sigma2UL*nvecL];     
            SmatUnR =  [Sigma2UR*nvecR zeros(ndm,1)
                        zeros(ndm,1) Sigma2UR*nvecR];                     
            term19L = (Sigma2UL*gamL*nvecL)';
            term19R = (Sigma2UR*gamR*nvecR)';
%            term20L = (theta1L-pressL/lam); % interier part
%            term20R = (theta1R-pressR/lam); % interier part
            
%% stiffness terms

            
             gamajumpuL = gamL*jumpu(1:2);
             gamajumpuR = gamR*jumpu(1:2);             
             term5L=cmatnBL*[eye(3,3)*gamajumpuL(1) eye(3,3)*gamajumpuL(2)]'*(P1*BmatL(1:4,:));
             term5R=cmatnBR*[eye(3,3)*gamajumpuR(1) eye(3,3)*gamajumpuR(2)]'*(P1*BmatR(1:4,:));
            
             sig8L2 = [gamajumpuL(1) gamajumpuL(2)]*cmatnL;

%             sig8L3= [sig8L2(1) 0  sig8L2(3)/2 sig8L2(3)/2
%                     0 sig8L2(2)  sig8L2(3)/2 -sig8L2(3)/2
%                     sig8L2(3)/2  sig8L2(3)/2 (sig8L2(2)+sig8L2(1))/4 (sig8L2(2)-sig8L2(1))/4
%                     sig8L2(3)/2 -sig8L2(3)/2 (sig8L2(2)-sig8L2(1))/4 (sig8L2(2)+sig8L2(1))/4];   

          sig8L3 = [sig8L2(1) 0 sig8L2(3) 0;
                    sig8L2(3) 0 sig8L2(2) 0;
                    0 sig8L2(1) 0 sig8L2(3);
                    0 sig8L2(3) 0 sig8L2(2)];
               
           term8L = BmatL(1:4,:)'*P2'*sig8L3*(P3*BmatL(1:4,:));
                                   
%             sig8R = [cmatnR(1,:) zeros(1,3) 
%                      zeros(1,3)  cmatnR(2,:) ];  
%             sig8R1 = sig8R*[eye(3,3)*jumpu(1) eye(3,3)*jumpu(2)]';  
% 
%              sig8R2 = [sum(sig8R1(:,1)) sum(sig8R1(:,2)) sum(sig8R1(:,3))];
          
             sig8R2 = [gamajumpuR(1) gamajumpuR(2)]*cmatnR;
             
%             sig8R3= [sig8R2(1) 0  sig8R2(3)/2 sig8R2(3)/2
%                     0 sig8R2(2)  sig8R2(3)/2 -sig8R2(3)/2
%                     sig8R2(3)/2  sig8R2(3)/2 (sig8R2(2)+sig8R2(1))/4 (sig8R2(2)-sig8R2(1))/4
%                     sig8R2(3)/2 -sig8R2(3)/2 (sig8R2(2)-sig8R2(1))/4 (sig8R2(2)+sig8R2(1))/4];          
% 
%              term8L=BmatL'*sig8L3*BmatL;   
%              term8R=BmatR'*sig8R3*BmatR;

            sig8R3 = [sig8R2(1) 0 sig8R2(3) 0;
                     sig8R2(3) 0 sig8R2(2) 0;
                     0 sig8R2(1) 0 sig8R2(3);
                     0 sig8R2(3) 0 sig8R2(2)];
               
              term8R =BmatR(1:4,:)'*P2'*sig8R3*(P3*BmatR(1:4,:));
             
              term30L=NmatL(1:2,:)'*ep*jumpu(1:2);
              term30R=NmatR(1:2,:)'*ep*jumpu(1:2);

%              [dmatL]=dmat2(JxXL,matepropL,lam);
%              [dmatR]=dmat2(JxXR,matepropR,lam);        
             dmatL2 = [eye(3,3)*jumpu(1) eye(3,3)*jumpu(2)]*[eye(3,3)*gamL(1,1) eye(3,3)*gamL(1,2)
                                                             eye(3,3)*gamL(2,1) eye(3,3)*gamL(2,2)]*nvectL2*dmatL1;  %missing J here
             dmatR2 = [eye(3,3)*jumpu(1) eye(3,3)*jumpu(2)]*[eye(3,3)*gamR(1,1) eye(3,3)*gamR(1,2)
                                                             eye(3,3)*gamR(2,1) eye(3,3)*gamR(2,2)]*nvectR2*dmatR1;  %missing J here
             
             term7L = BmatL(1:4,:)'*P1'*dmatL2*(P1*BmatL(1:4,:));
             term7R = BmatR(1:4,:)'*P1'*dmatR2*(P1*BmatR(1:4,:));
% perssure stiffness terms
             
             term21L = P2'*(SmatUnL*gamL./theta2L)*jumpu(1:2);  %similar to term17
             term21R = P2'*(SmatUnR*gamR./theta2R)*jumpu(1:2);
             term22L = P1'*(gamL*nvectL1*cpmat0(1:3,1:3))'*jumpu(1:2);  %similar to term18
             term22R = P1'*(gamR*nvectR1*cpmat0(1:3,1:3))'*jumpu(1:2);
             term23L = gamL*spmat0*nvecL; %similar to term28, but only left part
             term23R = gamR*spmat0*nvecR;
             
%% right hand side  

              ElemFL = ElemFL - (-c1L*BmatL(1:4,:)'*(term17L+term18L)*jumpu(1:2) - term28L + C1L*term30L - c1L*BmatL(5,:)'*term19L*jumpu(1:2));  % 
              ElemFR = ElemFR - (c1R* BmatR(1:4,:)'*(term17R+term18R)*jumpu(1:2) + term28R - C1R*term30R + c1R*BmatR(5,:)'*term19R*jumpu(1:2));  %;
          

%% left hand side
%relate to eta & Dd
             ElemKLL = ElemKLL - c1L*BmatL(1:4,:)'*(term17L+term18L)*NmatL(1:2,:) - c1L*NmatL(1:2,:)'*(term17L'+term18L')*BmatL(1:4,:);  
             ElemKLR = ElemKLR + c1L*BmatL(1:4,:)'*(term17L+term18L)*NmatR(1:2,:) + c1R*NmatL(1:2,:)'*(term17R'+term18R')*BmatR(1:4,:);
             ElemKRL = ElemKRL + c1R*BmatR(1:4,:)'*(term17R+term18R)*NmatL(1:2,:) + c1L*NmatR(1:2,:)'*(term17L'+term18L')*BmatL(1:4,:);
             ElemKRR = ElemKRR - c1R*BmatR(1:4,:)'*(term17R+term18R)*NmatR(1:2,:) - c1R*NmatR(1:2,:)'*(term17R'+term18R')*BmatR(1:4,:);
            %rp penalty terms  
             ElemKLL = ElemKLL + C1L*(NmatL(1:2,:)'*ep*NmatL(1:2,:));
             ElemKLR = ElemKLR - C1R*(NmatL(1:2,:)'*ep*NmatR(1:2,:));
             ElemKRL = ElemKRL - C1L*(NmatR(1:2,:)'*ep*NmatL(1:2,:));
             ElemKRR = ElemKRR + C1R*(NmatR(1:2,:)'*ep*NmatR(1:2,:));              
            %relate to jumpu terms additional terms
             ElemKLL = ElemKLL - c1L*(term5L+term5L'+term7L+term8L);  %                     
             ElemKRR = ElemKRR + c1R*(term5R+term5R'+term7R+term8R);  %
%relate to the eta & Dp
             ElemKLL = ElemKLL - c1L*BmatL(1:4,:)'*(term21L+term22L)*NmatL(3,:) - c1L*NmatL(1:2,:)'*term23L*NmatL(3,:);         
             ElemKLR = ElemKLR                                                  + c1R*NmatL(1:2,:)'*term23R*NmatR(3,:);
             ElemKRL = ElemKRL                                                  + c1L*NmatR(1:2,:)'*term23L*NmatL(3,:);
             ElemKRR = ElemKRR + c1R*BmatR(1:4,:)'*(term21R+term22R)*NmatR(3,:) - c1R*NmatR(1:2,:)'*term23R*NmatR(3,:);  
%relate to the q & Du
             ElemKLL = ElemKLL - c1L*NmatL(3,:)'*(theta2*(term21L+term22L))'*BmatL(1:4,:) - c1L*NmatL(3,:)'*(theta2*term23L)'*NmatL(1:2,:); 
             ElemKLR = ElemKLR                                                            + c1L*NmatL(3,:)'*(theta2*term23L)'*NmatR(1:2,:);  
             ElemKRL = ElemKRL                                                            + c1R*NmatR(3,:)'*(theta2*term23R)'*NmatL(1:2,:);  
             ElemKRR = ElemKRR + c1R*NmatR(3,:)'*(theta2*(term21R+term22R))'*BmatR(1:4,:) - c1R*NmatR(3,:)'*(theta2*term23R)'*NmatR(1:2,:); 
%relate to the q & Dp
%              ElemKLL22
%              ElemKLR22
%              ElemKRL22
%              ElemKRR22 
          
    end %lint        
% ElemKLL
            ElemK = [ElemKLL ElemKLR; ElemKRL ElemKRR];
            ElemF = [ElemFL; ElemFR];
            
%%
    case 6
        
       % generate the initial parameter 
        ElemKRR = zeros(nstR,nstR);
        ElemKRL = zeros(nstR,nstL);
        ElemKLR = zeros(nstL,nstR);
        ElemKLL = zeros(nstL,nstL);
        ElemFR = zeros(nstR,1);
        ElemFL = zeros(nstL,1);
        
        t1 = zeros(3,1);
        t2 = t1;
        t3 = t1;           
        dtol = 1e-11;
        itchat = 100;12;8; %flag for iteration number to suppress chatter
        
        I1 = [1; 1; 0; 0];
        I2 = eye(2);
        I6 = eye(4);
        P1 = [1 0 0 0
              0 1 0 0
              0 0 1 0];
        P2 = [1 0 0 0 
              0 0 1/2 1/2
              0 0 1/2 -1/2
              0 1 0 0];
        P3 = [1 0 0 0 
              0 0 1/2 -1/2
              0 0 1/2 1/2
              0 1 0 0]; 
        Z2 = zeros(2);
        spvec0 = I1;
        spmat0 = I2; % s_cup
        cpmat1 = I1*I1';
        cpmat2 = diag([-2,-2,-1,0]);
        cpmat0 = cpmat1 + cpmat2;  %c_cup

        ElemEL = matepropL(1);
        ElemvL = matepropL(2);
        ElemER = matepropR(1);
        ElemvR = matepropR(2);
        lam = getlam(matepropL); %For both material is the same at each side
       
        NmatL = zeros(3,nstL);
        BmatL = zeros(5,nstL);
%         BmatL1 = zeros(3,nstL);        
%         BmatL2 = zeros(4,nstL);
%         BmatL3 = zeros(4,nstL);       
        NmatR = zeros(3,nstR);
        BmatR = zeros(5,nstR);
%         BmatR1 = zeros(3,nstR);   
%         BmatR2 = zeros(4,nstR);
%         BmatR3 = zeros(4,nstL);        
        dmatL = zeros(9,3);
        dmatR = zeros(9,3);
     
        % Determin bounds of integration segment
        InterConn2D2 % InterConn2DT % 
        
        % Set jacobian for integration space
        drdr = (eL2 - eL1)/drL;
        drdrR = (eR2 - eR1)/drR;
        m = (eR2-eR1)/(eL1-eL2);
        if nelL == 3 || nelL == 6  
        lint = 3;10;2;3;
        else
        lint = 3;10;4;10;2;3; %10 for body force problem; 4 for other problem 
        end
        ideriv = 0;
     


%Form tau matrix: compute on first iteration and store; this gives
        %quadratic convergence when tau isn't recalculated (no
        %linearization of tau required).
   iterset =1; 4;
   if iter  <=iterset % == 0 %
        [tauL,intb] = TauS2_M(xlL,ulL,matepropL,nelL,nelLP,nen,lam,roL,eL1,drdr); %[Y^(-1)]
%        [tauL,intb] = TauS2_1(xlL,ulL,mateprop,nelL,nel2L,nen,lam,roL,eL1,drdr); %[Y^(-1)]        
%       TauListL(elem,:) = [tauL(1),tauL(2),tauL(3),tauL(4),tauL(5),tauL(6),tauL(7),tauL(8),tauL(9)];
        [tauR,intb] = TauS2_M(xlR,ulR,matepropR,nelR,nelRP,nen,lam,roL,eL1,drdrR);
%        [tauR,intb] = TauS2_1(xlR,ulR,mateprop,nelR,nel2R,nen,lam,roL,eL1,drdrR);
%          TauListR(elem,:) = [tauR(1),tauR(2),tauR(3),tauR(4),tauR(5),tauR(6),tauR(7),tauR(8),tauR(9)];

%         tau = tauL + tauR;

        ebL = 0;
        ebR = 0;
        intedge = 0;
        
        for ie = 1:lint            
% For separate bubble types on T and Q
           if nelL == 3 || nelL == 6
                [Wgt,litr,lits] = intpntt(ie,lint,1); 
                 ebeL = edgebubble(litr,lits,nelL);  %edgebubble is for T3 element
%                 rL = drdr*(litr-roL)+eL1;
%                 [shlL,shldL,shlsL,be] = shlt(rL,lits,nelL,nel2L,0,0);
           else
                [Wgt,litr,lits] = intpntq(ie,lint,1);
%                 rL = drdr*(litr-roL)+eL1;
%                 [shlL,shldL,shlsL,be] = shlq(rL,lits,nelL,nel2L,0,0);
                 ebeL = edgebubbleQ(litr,lits,nelL);
           end
                            
           if nelR == 3 || nelR == 6
                [Wgt,litr,lits] = intpntt(ie,lint,1); 
                 ebeR = edgebubble(litr,lits,nelL);
%                 rR = drdr*(litr-roL)+eR1;
%                 [shlR,shldR,shlsR,be] = shlt(rR,lits,nelR,nel2R,0,0);
           else
                [Wgt,litr,lits] = intpntq(ie,lint,1);
%                 rR = drdr*(litr-roL)+eR1;
%                 [shlR,shldR,shlsR,be] = shlq(rR,lits,nelR,nel2R,0,0);
                 ebeR = edgebubbleQ(litr,lits,nelR);
           end         
                    
%            b = edgebubble(litr,0);
           if nelL == 3 || nelL == 6
               rL = drdr*(litr-roL)+eL1;
               [shlL,shldL,shlsL,be] = shlt(rL,lits,nelL,nelL,0,0);
               [PxyL,shgsL,JdetL,bubbleL,XsL] = shgt(xlL(:,1:nelL),nelL,shldL,shlsL,nen,0,0,be); 
           else    
               rL = drdr*(litr-roL)+eL1;
               [shlL,shldL,shlsL,be] = shlq(rL,lits,nelL,nelL,0,0);
               [PxyL,shgsL,JdetL,bubbleL,XsL] = shgq(xlL(:,1:nelL),nelL,shldL,shlsL,nen,0,0,be); 
           end
%            [PxyL,shgsL,JdetL,bubbleL,xsL] = shgq(xlL(:,1:nelL)+ulL(:,1:nelL),nelL,shldL,shlsL,nen,0,0,be);
%             QxyL = PxyL;  
            
            
%             rR = m*(rL-eL2) + eR1;
%             
%             if nelR == 3 || nelR == 6
%                 sR = 0;
%             else %if nelR == 4
%                 sR = -1;
%             end
%             
%             [shlR,shldR,shlsR,be] = shlq(rR,sR,nelR,nel2R,0,0);
%             [PxyR,shgsR,JdetR,bubbleR,XsR] = shgq(xlR(:,1:nelR),nelR,shldR,shlsR,nen,0,0,be);     
% %             [PxyR,shgsR,JdetR,bubbleR,xsR] = shgq(xlR(:,1:nelR)+ulR(:,1:nelR),nelR,shldR,shlsR,nen,0,0,be);           
% %             QxyR = PxyR;
             
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
        
        edgeK = tauL*ebL^2 + tauR*ebR^2;
        gamL = ebL^2*(edgeK\tauL);
        gamL_list(iter+1,:,inter) = [gamL(1,1),gamL(1,2),gamL(2,1),gamL(2,2)];
        gamR = ebR^2*(edgeK\tauR);
        gamR_list(iter+1,:,inter) = [gamR(1,1),gamR(1,2),gamR(2,1),gamR(2,2)];
%         gamL_list(elem,:) = [gamL(1),gamL(2),gamL(3),gamL(4),gamL(5),gamL(6),gamL(7),gamL(8),gamL(9)];
%         gamR_list(elem,:) = [gamR(1),gamR(2),gamR(3),gamR(4),gamR(5),gamR(6),gamR(7),gamR(8),gamR(9)];
        ep = pencoeff*intedge*inv(edgeK); 
        ep_list(iter+1,:,inter) = [ep(1,1),ep(1,2),ep(2,1),ep(2,2)];
%         ep_List(elem,:) = [ep(1),ep(2),ep(3),ep(4),ep(5),ep(6),ep(7),ep(8),ep(9)];
%       ep_List_R(elem,:) = [tauR(1),tauR(2),tauR(3),tauR(4),tauR(5),tauR(6),tauR(7),tauR(8),tauR(9)];         
   else
        gamL = [gamL_list(iterset+1,1,inter) gamL_list(iterset+1,2,inter)
                gamL_list(iterset+1,3,inter) gamL_list(iterset+1,4,inter)];
        gamR = [gamR_list(iterset+1,1,inter) gamR_list(iterset+1,2,inter)
                gamR_list(iterset+1,3,inter) gamR_list(iterset+1,4,inter)];
        ep = [ep_list(iterset+1,1,inter) ep_list(iterset+1,2,inter)
              ep_list(iterset+1,3,inter) ep_list(iterset+1,4,inter)];
   end
%        gamL = [0.5 0
%                 0 0.5];
%        gamR = [0.5 0
%                 0 0.5];     
%        ep = 20*eye(2);
   
%        s = -1;
        ll=0;           
     for l = 1:lint

            ll = ll + 1;
            if nelL == 3 || nelL == 6
                [Wgt,litr,lits] = intpntt(ll,lint,1);
            elseif nelL == 4 || nelL == 9
                [Wgt,litr,lits] = intpntq(ll,lint,1);
            end
                    
%            b = edgebubble(litr,0);
            
            rL = drdr*(litr-roL)+eL1;
           if nelL == 3 || nelL == 6
           [shlL,shldL,shlsL,be] = shlt(rL,lits,nelL,nelL,0,0);
           [QxyL,shgsL,JdetL,bubbleL,XsL] = shgt(xlL(:,1:nelL),nelL,shldL,shlsL,nen,0,0,be); 
           [QxyL,shgsL,JdetL,bubbleL,xsL] = shgt(xlL(:,1:nelL)+ulL(1:2,1:nelL),nelL,shldL,shlsL,nen,0,0,be);
           else
           [shlL,shldL,shlsL,be] = shlq(rL,lits,nelL,nelL,0,0);
           [QxyL,shgsL,JdetL,bubbleL,XsL] = shgq(xlL(:,1:nelL),nelL,shldL,shlsL,nen,0,0,be); 
           [QxyL,shgsL,JdetL,bubbleL,xsL] = shgq(xlL(:,1:nelL)+ulL(1:2,1:nelL),nelL,shldL,shlsL,nen,0,0,be);      
           end
% pressure field
             if nelLP == 3 || nelLP == 6
              [shlpL,shldL,shlsL] = shlt(rL,lits,nelLP,nelL,0,0); %shlpL:shape fun
              [PxyL, shgpL] = shgt(xlL(:,1:nelL)+ulL(1:2,1:nelL),nelLP,shldL,shlsL,nelLP,0,0,be); %PxyL:deriv of shape fun for pressure field
            else
              [shlpL,shldL,shlsL] = shlq(rL,lits,nelLP,nelL,0,0);
              [PxyL, shgpL] = shgq(xlL(:,1:nelL)+ulL(1:2,1:nelL),nelLP,shldL,shlsL,nelLP,0,0,be);
            end           
            
            rR = m*(rL-eL2) + eR1;
            
            if nelR == 3 || nelR == 6
                sR = 0;
            else %if nelR == 4
                sR = -1;
            end
           if nelR == 3 || nelR == 6            
            [shlR,shldR,shlsR,be] = shlt(rR,sR,nelR,nelR,0,0);
            [QxyR,shgsR,JdetR,bubbleR,XsR] = shgt(xlR(:,1:nelR),nelR,shldR,shlsR,nen,0,0,be);     
            [QxyR,shgsR,JdetR,bubbleR,xsR] = shgt(xlR(:,1:nelR)+ulR(1:2,1:nelR),nelR,shldR,shlsR,nen,0,0,be);  
           else
            [shlR,shldR,shlsR,be] = shlq(rR,sR,nelR,nelR,0,0);
            [QxyR,shgsR,JdetR,bubbleR,XsR] = shgq(xlR(:,1:nelR),nelR,shldR,shlsR,nen,0,0,be);     
            [QxyR,shgsR,JdetR,bubbleR,xsR] = shgq(xlR(:,1:nelR)+ulR(1:2,1:nelR),nelR,shldR,shlsR,nen,0,0,be);  
           end

             if nelRP == 3 || nelRP == 6
              [shlpR,shldR,shlsR] = shlt(rR,sR,nelRP,nelR,0,0); %shlpL:shape fun
              [PxyR, shgpR] = shgt(xlR(:,1:nelR)+ulR(1:2,1:nelR),nelRP,shldR,shlsR,nelRP,0,0,be); %PxyL:deriv of shape fun for pressure field
            else
              [shlpR,shldR,shlsR] = shlq(rR,sR,nelRP,nelR,0,0);
              [PxyR, shgpR] = shgq(xlR(:,1:nelR)+ulR(1:2,1:nelR),nelRP,shldR,shlsR,nelRP,0,0,be);
            end  
            
 % initial stress sigma and material bulk term c_ijkl
           for mm = 1:nelL    
%            NmatL(:,2*mm-1:2*mm) = [shlL(mm,1)     0          
%                                     0        shlL(mm,1)];                
           NmatL(:,3*mm-2:3*mm) = [shlL(mm,1)     0       0         
                                    0        shlL(mm,1)   0
                                    0             0     shlpL(mm,1)]; 
%            BmatL(:,2*mm-1:2*mm) = [QxyL(mm,1) 0             
%                                     0         QxyL(mm,2)         
%                                     QxyL(mm,2) QxyL(mm,1)        
%                                      QxyL(mm,2) -QxyL(mm,1)]; 
%            PmatL(:,mm) = shlpL(mm,1) ;  %shape function of pressure
           BmatL(:,3*mm-2:3*mm) = [QxyL(mm,1) 0             0        
                                    0         QxyL(mm,2)    0      
                                    QxyL(mm,2) QxyL(mm,1)   0      
                                     QxyL(mm,2) -QxyL(mm,1) 0
                                     0          0           shlpL(mm,1)]; %similar to Bmat [derivative of shape disp + shape of pressure]
%            BBmatL(:,3*mm-2:3*mm) = [shgsL(1,1) 0         0
%                                     shgsL(1,2) 0         0
%                                     shgsL(1,3) 0         0
%                                      0         shgsL(1,1) 0
%                                      0         shgsL(1,2) 0
%                                      0         shgsL(1,3) 0 
%                                      0         0         PxyL(1,1)
%                                      0         0         PxyL(1,2)]; %%similar to BBmat [second derivative of shape disp + derivative of shape of pressure]
%            BmatL1(:,3*mm-2:3*mm) = [QxyL(mm,1) 0                  
%                                      0         QxyL(mm,2)          
%                                      QxyL(mm,2) QxyL(mm,1)];
%             BmatL2(:,3*mm-2:3*mm)=[QxyL(mm,1)  0  0                   
%                                   QxyL(mm,2)  0   0                                   
%                                     0     QxyL(mm,1) 0                
%                                     0     QxyL(mm,2) 0
%                                     0          0     shlpL(mm,1)  ];   
%            BmatL3(:,3*mm-2:3*mm)=[QxyL(mm,1)  0                     
%                                    0     QxyL(mm,1)                                      
%                                    QxyL(mm,2) 0                     
%                                    0     QxyL(mm,2)];     
           end 
            
           for mm = 1:nelR    
%            NmatR(:,2*mm-1:2*mm) = [shlR(mm,1)     0          
%                                      0        shlR(mm,1)];                
           NmatR(:,3*mm-2:3*mm) = [shlR(mm,1)     0       0        
                                     0        shlR(mm,1)  0
                                     0            0       shlpR(mm,1)]; 
%            BmatR(:,2*mm-1:2*mm) = [QxyR(mm,1) 0             
%                                     0         QxyR(mm,2)         
%                                     QxyR(mm,2) QxyR(mm,1)        
%                                      QxyR(mm,2) -QxyR(mm,1)];   
%            PmatR(:,mm) = shlpR(mm,1) ;  %shape function of pressure                                 
           BmatR(:,3*mm-2:3*mm) = [QxyR(mm,1) 0              0  
                                    0         QxyR(mm,2)     0    
                                    QxyR(mm,2) QxyR(mm,1)    0    
                                     QxyR(mm,2) -QxyR(mm,1)  0
                                     0           0           shlpR(mm,1)];
                               
%            BmatR1(:,3*mm-2:3*mm) = [QxyR(mm,1) 0                  
%                                      0         QxyR(mm,2)          
%                                      QxyR(mm,2) QxyR(mm,1)];
%            BmatR2(:,3*mm-2:3*mm)=[QxyR(mm,1)  0                     
%                                  QxyR(mm,2)  0                                      
%                                    0     QxyR(mm,1)                
%                                    0     QxyR(mm,2)];   
%            BmatR3(:,3*mm-2:3*mm)=[QxyR(mm,1)  0                     
%                                    0     QxyR(mm,1)                                      
%                                    QxyR(mm,2) 0                     
%                                    0     QxyR(mm,2)];                       
           end  
    %plug in the -ul to get the inverse of deformation gradient fi as 'F' in this function             
            [fiL,JxXL,FL] = kine2d(PxyL,-ulL,nelL,0); %this is equivalent to ikine2d
            JxXL = 1/JxXL; %this is equivalent to ikine2d
    %         [F,JxX,fi] = ikine2d(Qxy,ul,nel); %this is equivalent to above two lines
            JdetL = JdetL/JxXL;            
            [theta1L,theta2L,theta3L] = ThetaNS(JxXL,matepropL);
            
            [fiR,JxXR,FR] = kine2d(PxyR,-ulR,nelR,0); %this is equivalent to ikine2d
            JxXR = 1/JxXR; %this is equivalent to ikine2d
    %         [F,JxX,fi] = ikine2d(Qxy,ul,nel); %this is equivalent to above two lines
            JdetR = JdetR/JxXR;               
            [theta1R,theta2R,theta3R] = ThetaNS(JxXR,matepropR);   
            
            spvecL = JxXL*spvec0;
%             spmatL = JxXL*spmat0;
            cpmatL = JxXL*cpmat0;
%             cpmat32L = (theta3*JxXL^2 + theta2*JxXL)*cpmat1 + theta2*JxXL*cpmat2;
%             dpmatL = JxXL*dpmat0;            
          
         [sigmaiL, cmatiL] = SigmaCmatNSCST2i(FL,JxXL,matepropL); 
                     
            pressL = ulL(3,:)*shlpL;
            sigmapL = pressL*spvecL;
            cmatpL = pressL*cpmatL;
            sigma2L = (sigmaiL + sigmapL)./JxXL;  % [sigma11 sigma22 sigma12 0]
            cmatL = (cmatiL + cmatpL)./JxXL;      % [last one is zero]      
            
            SmatL1=[sigma2L(1), sigma2L(3)
                    sigma2L(3), sigma2L(2)];           
            
%            [sigma2R, cmatR] = SigmaCmat2i(FR,JxXR,matepropR,lam);
            spvecR = JxXR*spvec0;
%             spmatR = JxXR*spmat0;
            cpmatR = JxXR*cpmat0;
%             cpmat32R = (theta3*JxXR^2 + theta2*JxXR)*cpmat1 + theta2*JxXR*cpmat2;
%             dpmatR = JxXR*dpmat0; 
            
         [sigmaiR, cmatiR] = SigmaCmatNSCST2i(FR,JxXR,matepropR); 
                    
            pressR = ulR(3,:)*shlpR;
            sigmapR = pressR*spvecR;
            cmatpR = pressR*cpmatR;
            sigma2R = (sigmaiR + sigmapR)./JxXR;
            cmatR = (cmatiR + cmatpR)./JxXR;  
            
            SmatR1=[sigma2R(1), sigma2R(3)
                    sigma2R(3), sigma2R(2)];     

             %Evaluate tangent and normal vectors
            t1L = [xsL(:,1); 0];
            [tm1L, tu1L] = VecNormalize(t1L);
            t2L = [0; 0; 1];
            tm2L = 1;
            tu2L = t2L';
            t3L = VecCrossProd(t1L,t2L);
            [tm3L, tu3L] = VecNormalize(t3L);

                %Evaluate tangent and normal vectors
             t1R = [xsR(:,1); 0];
            [tm1R, tu1R] = VecNormalize(t1R);
            t2R = [0; 0; 1];
            tm2R = 1;
            tu2R = t2R';
            t3R = VecCrossProd(t1R,t2R);
            [tm3R, tu3R] = VecNormalize(t3R);
             c1L = drdr*Wgt*tm3L;
             c1R = drdrR*Wgt*tm3R; 
%              c1L = c1R;
             
             %Evaluate tangent and normal vectors
            T1L = [XsL(:,1); 0];
            [Tm1L, Tu1L] = VecNormalize(T1L);
            T2L = [0; 0; 1];
            Tm2L = 1;
            Tu2L = T2L';
            T3L = VecCrossProd(T1L,T2L);
            [Tm3L, Tu3L] = VecNormalize(T3L);

                %Evaluate tangent and normal vectors
             T1R = [XsR(:,1); 0];
            [Tm1R, Tu1R] = VecNormalize(T1R);
            T2R = [0; 0; 1];
            Tm2R = 1;
            Tu2R = T2R';
            T3R = VecCrossProd(T1R,T2R);
            [Tm3R, Tu3R] = VecNormalize(T3R);
             C1L = drdr*Wgt*Tm3L;
             C1R = drdrR*Wgt*Tm3R;  
%              C1L = C1R;
                %normal vectors
                nLx = tu3L(1);
                nLy = tu3L(2);
                nRx = tu3R(1);
                nRy = tu3R(2);
                %tagent vectors
                tLx = tu1L(1);
                tLy = tu1L(2);
                tRx = tu1R(1);
                tRy = tu1R(2);             
                nvectL1= [nLx 0   nLy   
                          0  nLy  nLx ];
                nvectR1 = [nRx 0  nRy  
                          0  nRy  nRx ]; 
                nvectL2 = [eye(3,3)*nLx zeros(3,3) eye(3,3)*nLy
                           zeros(3,3)    eye(3,3)*nLy eye(3,3)*nLx];
                nvectR2 = [eye(3,3)*nRx zeros(3,3) eye(3,3)*nRy
                           zeros(3,3)    eye(3,3)*nRy eye(3,3)*nRx];                       
                nvecL = [nLx; nLy];
                nvecR = [nRx; nRy];


             rhspulL = reshape(ulL,ndf*nelL,1);
             rhspulR = reshape(ulR,ndf*nelR,1);
                    
            jumpu = NmatL*rhspulL - NmatR*rhspulR;  %jumpu=uL-uR
%% right hand side force terms           
           
           SmatnL=[SmatL1*nvecL zeros(ndm,1)
                    zeros(ndm,1) SmatL1*nvecL];
            
            SmatnR=[SmatR1*nvecR zeros(ndm,1)
                    zeros(ndm,1) SmatR1*nvecR];
           cmatnL=(nvectL1*cmatL(1:3,1:3));
           cmatnR=(nvectR1*cmatR(1:3,1:3));
           cmatnBL=BmatL(1:4,:)'*P2'*[cmatnL zeros(2,3)
                            zeros(2,3)  cmatnL ];               
           cmatnBR=BmatR(1:4,:)'*P2'*[cmatnR    zeros(2,3)
                            zeros(2,3)  cmatnR ];

            term17L=P2'*SmatnL*gamL;
            term17R=P2'*SmatnR*gamR;

            term18L=P1'*(gamL*nvectL1*cmatL(1:3,1:3))';
            term18R=P1'*(gamR*nvectR1*cmatR(1:3,1:3))';

            term28L=NmatL(1:2,:)'*(c1L*gamL*SmatL1*nvecL-c1R*gamR*SmatR1*nvecR);   %average stress term
            term28R=NmatR(1:2,:)'*(c1L*gamL*SmatL1*nvecL-c1R*gamR*SmatR1*nvecR);

% additional terms related to the pressure:
            Sigma2UL = theta2L*spmat0; % sigma^U=U''(J)*I  
            Sigma2UR = theta2R*spmat0; % sigma^U=U''(J)*I 
            SmatUnL =  [Sigma2UL*nvecL zeros(ndm,1)
                        zeros(ndm,1) Sigma2UL*nvecL];     
            SmatUnR =  [Sigma2UR*nvecR zeros(ndm,1)
                        zeros(ndm,1) Sigma2UR*nvecR];                     
            term19L = (Sigma2UL*gamL*nvecL)';
            term19R = (Sigma2UR*gamR*nvecR)';
%            term20L = (theta1L-pressL/lam); % interier part
%            term20R = (theta1R-pressR/lam); % interier part
            
            term30L=NmatL(1:2,:)'*ep*jumpu(1:2);
            term30R=NmatR(1:2,:)'*ep*jumpu(1:2);

             
%% right hand side  

              ElemFL = ElemFL - (-c1L*BmatL(1:4,:)'*(term17L+term18L)*jumpu(1:2) - term28L + C1L*term30L - c1L*BmatL(5,:)'*term19L*jumpu(1:2));  % 
              ElemFR = ElemFR - (c1R* BmatR(1:4,:)'*(term17R+term18R)*jumpu(1:2) + term28R - C1R*term30R + c1R*BmatR(5,:)'*term19R*jumpu(1:2));  %;
          
          
    end %lint  
%%

    case -1 % boundary tractions
        
        ElemF = zeros(nst,1);

        % Determine bounds of integration
        
        if nel == 4 || nel == 9
            
        dr = 2;
        ro = -1;
        
        % Upper Limit
        if nodeA == ElemFlag(2)
            eR2 = 1;
        else %
            eR2 = 0;
        end
        % Lower Limit
        if nodeB == ElemFlag(1)
            eR1 = -1;
        else %nodeA == ElemFlag(5)
            eR1 = 0;
        end
        
        elseif nel == 3 || nel == 6
            
        dr = 1;
        ro = 0;
            
        % Upper Limit
        if nodeA == ElemFlag(2)
            eR2 = 1;
        else %nodeA == ElemFlagR(5)
            eR2 = 1/2;
        end
        % Lower Limit
        if nodeB == ElemFlag(1)
            eR1 = 0;
        else %nodeA == ElemFlag(5)
            eR1 = 1/2;
        end
        
        end
        
        % Set jacobian for integration space
        drdr = (eR2 - eR1)/dr;
        
        lint = 10;4; %10; % Use 10 for body force BF2U4M0.m
%         % Load Gauss Points for quadrature
%         if enrich == 1
%             [rlist, rWgts, rnum] = GaussPoints(pr+2);
%             slist = -1;
%             sWgts = 1;
%             snum = 1;
%         else
%             [rlist, rWgts, rnum] = GaussPoints(pr+1);
%             slist = -1;
%             sWgts = 1;
%             snum = 1;
%         end

%         lamda = Elemv*ElemE/((1+Elemv)*(1-2*Elemv));
%         mu = ElemE/(2*(1+Elemv));
%         thick = 1;
        ideriv = 0;
        der = 0;
        bf = 0;
        thick = 1;
        
        for je = 1:lint

            if nel == 3 || nel == 6
                [Wgt,r,s] = intpntt(je,lint,1);
            elseif nel == 4 || nel == 9
                [Wgt,r,s] = intpntq(je,lint,1);
            end
            r = drdr*(r-ro)+eR1;

            if nel == 3 || nel == 6
              [shl,shld,shls,be] = shlt(r,s,nel,nel,der,bf);
              [QXY, shgs, Jdet, be, xs] = shgt(xl,nel,shld,shls,nel,bf,der,be);
            else
              [shl,shld,shls,be] = shlq(r,s,nel,nel,der,bf);
              [QXY, shgs, Jdet, be, xs] = shgq(xl,nel,shld,shls,nel,bf,der,be);
            end
                    
            %Evaluate tangent and normal vectors
            t1 = [xs(:,1); 0];
            [tm1, tu1] = VecNormalize(t1);
            t2 = [0; 0; 1];
            tm2 = 1;
            tu2 = t2';
            t3 = VecCrossProd(t1,t2);
            [tm3, tu3] = VecNormalize(t3);
            
            if exist('iprob','var')
                if iprob == 6
                    mu = 40;
                    lam = 40;
                    X = xl(1,:)*shl;
                    Y = xl(2,:)*shl;
                    P = [ (101*X*lam*(101*X + 100))/10000, (101*Y*mu)/(101*X + 100) - (10201*X*Y*lam)/10000, 0
                         (101*Y*mu)/100, mu - (100*mu)/(101*X + 100) + X*((101*lam)/100 + (101*mu)/100),     0
                              0,                                            0, (101*X*lam*(101*X + 100))/10000];
                    Traction = P*tu3';
                else
                    Traction = traction;
                end
            else
                Traction = traction;
             end
            
            c1 = Wgt*tm3*drdr*thick;
            for o=1:nel
                don = shl(o);
                F = don*Traction';

    %                 Fmtl = F'*t; %Magnitudes of F dot tunit(l=1:3)
    % %                 for l = 1:3
    % %                     for m = 1:3
    %                         Ftl = t*Fmtl'; %Sum of Vectors {F.t(l)}t(m)
    % %                     end
    % %                 end  t*t'*F = eye*F = F

                ElemF(ndf*o-1) = ElemF(ndf*o-1) + F(1)*c1;

                ElemF(ndf*o-0)   = ElemF(ndf*o-0)   + F(2)*c1;

            end %o

        end %ie
        ElemF;
%%


        
%         end %intt
ep;
gamL;
ElemFL;
norm(ElemKLL-ElemKLL');
    case 22
        
        ElemM = zeros(nst);
        ElemF = zeros(nst,1);

        lam = getlam(mateprop);
        
        I1 = [1; 1; 1; 0; 0; 0; 0; 0; 0];
        stresID = [1 2 4 0 1 3];
        spvec0 = I1;
        
        thick = 1;
        
        % Load Guass Integration Points

        if nel == 3
            lint = 7; %minimum of 7 for all integrals in deformed state
            der = 0;
        elseif nel == 4
%             lint = 4;
            lint = 16;
            der = 1;
        elseif nel == 6
            lint = 13; %minimum of 13 for all integrals in deformed state
            der = 1;
        else
%             lint = 9;
            lint = 25;
            der = 0;
        end

        %Integration Loop
        for ll = 1:lint

            %Evaluate first derivatives of basis functions at int. point
            if nel == 3 || nel == 6
              [Wgt,litr,lits] =  intpntt(ll,lint,ib);
              [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
              [Qxy, shgs, Jdet,be] = shgt(xl+ul(1:2,:),nel,shld,shls,nel,bf,der,be);
            else
              [Wgt,litr,lits] =  intpntq(ll,lint,ib);
              [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
              [Qxy, shgs, Jdet,be] = shgq(xl+ul(1:2,:),nel,shld,shls,nel,bf,der,be);
            end

            if nelS == 3 || nelS == 6
              [shlS,shld,shls] = shlt(litr,lits,nelS,nel,0,0);
            else
              [shlS,shld,shls] = shlq(litr,lits,nelS,nel,0,0);
            end
            
            [fi,JxX,F] = kine2d(Qxy,-ul,nel,0); %this is equivalent to ikine2d
            JxX = 1/JxX; %this is equivalent to ikine2d
    %         [F,JxX,fi] = ikine2d(Qxy,ul,nel); %this is equivalent to above two lines
            Jdet = Jdet/JxX;

            sigma = SigmaCmat3([F zeros(2,1); zeros(1,2) 1],JxX,mateprop,lam);
            sigma = sigma/JxX;
            
            % Form B matrix
            Nmat = shlS';

            w = Wgt*Jdet*thick;
            
            if stres <= 3 % stress components
                sigmas = sigma(stresID(stres));
            elseif stres >= 5
                if stres <= 6 % principal stresses
                sigma = [sigma(1) sigma(4) sigma(6); sigma(4) sigma(2) sigma(5); sigma(6) sigma(5) sigma(3)];
                psig = eig(sigma);
                sigmas = psig(stresID(stres));
                else % hydrostatic stress
                sigmas = 1/3*sigma'*I1;
                end
            else % von Mises stress
                trs = sigma'*I1;
                dsig = sigma - 1/3*trs*I1;
                sigmas = sqrt(3/2*(dsig'*dsig+dsig(4)^2));
            end
            
            ElemF = ElemF + w*Nmat'*sigmas;
            
            ElemM = ElemM + w*(Nmat'*Nmat);

        end %je
%%        
    case 25 %Stress Projection2

        ElemS = zeros(nel,npstr+1);
        ElemS2 = zeros(nel,npstr);

        lam = getlam(mateprop);
        
        I1 = [1; 1; 1; 0; 0; 0; 0; 0; 0];
        stresID = [1 2 4 0 1 3];
        spvec0 = I1;
        
        thick = 1;
        
        % Load Guass Integration Points

        if nel == 3
            lint = 1;
            nint = 1;
        elseif nel == 4
%             lint = 4;
            lint = 4;
            nint = 1;
        elseif nel == 6
            lint = 7;
            nint = 3;
        else
            lint = 9;
            nint = 4;
        end
        
        der = 0;
        bf = 1;
        ib = 0;

        %Stress Loop
        for ll = 1:nint

            %Evaluate first derivatives of basis functions at int. point
            if nel == 3 || nel == 6
              [Wgt,litr,lits] =  intpntt(ll,nint,ib);
              [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
              [Qxy, shgs, Jdet,be] = shgt(xl+ul(1:2,:),nel,shld,shls,nel,bf,der,be);
            else
              [Wgt,litr,lits] =  intpntq(ll,nint,ib);
              [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
              [Qxy, shgs, Jdet,be] = shgq(xl+ul(1:2,:),nel,shld,shls,nel,bf,der,be);
            end
            
            [fi,JxX,F] = kine2d(Qxy,-ul,nel,0); %this is equivalent to ikine2d
            JxX = 1/JxX; %this is equivalent to ikine2d
    %         [F,JxX,fi] = ikine2d(Qxy,ul,nel); %this is equivalent to above two lines
            Jdet = Jdet/JxX;
            
            sigma = SigmaCmat3([F zeros(2,1); zeros(1,2) 1],JxX,mateprop,lam);
            sigma = sigma/JxX;
            
            for stres = 1:npstr
            
            if stres <= 3 % stress components
                sigmas = sigma(stresID(stres));
            elseif stres >= 5
                if stres <= 6 % principal stresses
                sigma2 = [sigma(1) sigma(4) sigma(6); sigma(4) sigma(2) sigma(5); sigma(6) sigma(5) sigma(3)];
                psig = eig(sigma2);
                sigmas = psig(stresID(stres));
                else % hydrostatic stress
                sigmas = 1/3*sigma'*I1;
                end
            else % von Mises stress
                trs = sigma'*I1;
                dsig = sigma - 1/3*trs*I1;
                sigmas = sqrt(3/2*(dsig'*dsig));
            end
            
            ElemS2(ll,stres) = sigmas;
            
            end

        end %je
        
        % interpolate stress at nodes
        if nel == 3
            plist = [0 1 0
                     0 0 1];
        elseif nel == 4
            plist = [-1 1 1 -1
                     -1 -1 1 1];
        elseif nel == 6
            plist = [-1/3 5/3 -1/3 2/3 2/3 -1/3
                     -1/3 -1/3 5/3 -1/3 2/3 2/3];
        else
            sqr3 = sqrt(3);
            plist = [-sqr3 sqr3 sqr3 -sqr3 0 sqr3 0 -sqr3 0
                     -sqr3 -sqr3 sqr3 sqr3 -sqr3 0 sqr3 0 0];
        end
        
        for ll = 1:nelS
            
            r = plist(1,ll);
            s = plist(2,ll);
            shpS = sshp2d(r,s,nint);
            
            for stres = 1:npstr
                
                sigmas = ElemS2(1:nint,stres)'*shpS;
                ElemS(ll,stres) = sigmas;
                
            end
            
        end
        
        %Integration Loop
        Vol = 0;
        for ll = 1:lint

            %Evaluate first derivatives of basis functions at int. point
            if nel == 3 || nel == 6
              [Wgt,litr,lits] =  intpntt(ll,lint,ib);
              [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
              [Qxy, shgs, Jdet,be] = shgt(xl+ul(1:2,:),nel,shld,shls,nel,bf,der,be);
            else
              [Wgt,litr,lits] =  intpntq(ll,lint,ib);
              [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
              [Qxy, shgs, Jdet,be] = shgq(xl+ul(1:2,:),nel,shld,shls,nel,bf,der,be);
            end
            
            [fi,JxX,F] = kine2d(Qxy,-ul,nel,0); %this is equivalent to ikine2d
            JxX = 1/JxX; %this is equivalent to ikine2d
    %         [F,JxX,fi] = ikine2d(Qxy,ul,nel); %this is equivalent to above two lines
            Jdet = Jdet/JxX;

            w = Wgt*Jdet*thick;
            
            Vol = Vol + w;

        end %je

        for i = 1:nel
        ElemS(i,npstr+1) = Vol; % use this for weighted average 1; % use this for simple average 
        end
end