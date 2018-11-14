%% DG implementation of large deformaiton 
% pure displacement
%04/17/2013 Pinlei Chen
%% have tau and delta in it
%% have body force in it
%% with d_ijklmn in it
%% lint for interface part =100,interior =1000, traction =100
%% for demage part in it, use the theory in the paper.
PatchE = mateprop(3);
Patchv = mateprop(4);

Bcol1 = [1; 4; 6];
Bcol2 = [2; 4; 5];
Bcol3 = [3; 5; 6];
col1 = [1; 2; 3];
col2 = [2; 1; 3];
col3 = [3; 2; 1];

switch isw %Task Switch

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

        dtol = 1e-11;
        itchat = 100;12;8; %flag for iteration number to suppress chatter
        
%          beta = 1;
%          beta2 = beta^2;
%          sigmax = 0.02e-3;
%          dc = 40;0.1;0.2;
          muf =  0.52;
       beta = 1;0.707;
       beta2 = beta^2;
        sigmax = 22;100;50; %0.01e-3;
       dc = 0.5;40;
%       change sigmaxc and dc to make it a perfect plasticity (make Hc to be a small number)        
        Hc = sigmax/dc;
        rp = 1000*Hc;
%         sigmax = 100;
        ElemEL = matepropL(3);
        ElemvL = matepropL(4);
        ElemER = matepropR(3);
        ElemvR = matepropR(4);
        lamL = getlam(matepropL);
        lamR = getlam(matepropR);
        
        BmatL1 = zeros(6,nstL);        
        BmatR1 = zeros(6,nstR);
        BmatL2 = zeros(9,nstL);
        BmatR2 = zeros(9,nstR);
        term17L = zeros(3,nstL);
        term18L = zeros(3,nstL);
        term17R = zeros(3,nstR);
        term18R = zeros(3,nstR);  
        
        nvectL1 = zeros(3,6);
        nvectR1 = zeros(3,6);        
        NmatL = zeros(3,nstL);
        NmatR = zeros(3,nstR);
        BmatL = zeros(9,nstL);
        BmatR = zeros(9,nstR);
        bnAdN2 = zeros(6,nstR);
        
        % Load Guass Integration Points

        if nel == 4
            lint = 16;11;5;16;
        elseif nel == 8
            lint = 4;16;100;100;8;
        elseif nel == 10
            lint = 14;
%             lint = 27;
        else
            lint = 27;
        end
        der = 0;
        bf = 0;
        ib = 5;   
        ll = 0; % Counter for history variables

        
                %Form tau matrix: compute on first iteration and store; this gives
        %quadratic convergence when tau isn't recalculated (no
        %linearization of tau required).
   iterset = 4;
   if iter  <=iterset % == 0 %
         [tauL,intb] = TauS3(xlL,ulL,matepropL,nel,lint,lamL); %[Y^(-1)]       
%         [tauL,intb] = TauS3_1(xlL,ulL,mateprop,nel,lint,lam); %[Y^(-1)]
%          TauListL(elem,:) = [tauL(1),tauL(2),tauL(3),tauL(4),tauL(5),tauL(6),tauL(7),tauL(8),tauL(9)];
         [tauR,intb] = TauS3(xlR,ulR,matepropR,nel,lint,lamR);
%         [tauR,intb] = TauS3_1(xlR,ulR,mateprop,nel,lint,lam);
%          TauListR(elem,:) = [tauR(1),tauR(2),tauR(3),tauR(4),tauR(5),tauR(6),tauR(7),tauR(8),tauR(9)];

%         tau = tauL + tauR;

        ebL = 0;
        ebR = 0;
        intedge = 0;
        
        for ie = 1:lint            
% For separate bubble types on T and Q
           if nelL == 4 || nelL == 10
                [w,ss] =  int3d_t(l,lint,ib);
                ebeL = facebubbleT(ss);
                [shlL,shldL,shls,be] = shltt([rL sL tL],nelL,nel2L,0,0);
           else
                [Wgt,ss] = intpntb(ie,lint,ib);
                [shlL,shldL,shls,be] = shlb(ss,nelL,nel2L,0,0);
                ebeL = facebubbleQ(ss);
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
                [w,ss] =  int3d_t(l,lint,ib);
                ebeR = facebubbleT(ss);
           else
                [Wgt,ss] = intpntb(ie,lint,ib);
                ebeR = facebubbleQ(ss);
           end           
            
           if nelL == 4 || nelL == 10
                [w,ss] =  int3d_t(l,lint,ib);
                [shlL,shldL,shls,be] = shltt(ss,nelL,nel2L,0,0);
                [PxyL,shgsL,JdetL,bubbleL,XsL] = shgtt(xlL(:,1:nelL),nelL,shldL,shls,nen,0,0,be);                 
                [PxyL,shgsL,JdetL,bubbleL,xsL] = shgtt(xlL(:,1:nelL)+ulL(:,1:nelL),nelL,shldL,shls,nen,0,0,be);
           else
                [Wgt,ss] = intpntb(l,lint,ib);
                [shlL,shldL,shls,be] = shlb(ss,nelL,nel2L,0,0);
                [PxyL,shgsL,JdetL,bubbleL,XsL] = shgb(xlL(:,1:nelL),nelL,shldL,shls,nen,0,0,be); 
                [PxyL,shgsL,JdetL,bubbleL,xsL] = shgb(xlL(:,1:nelL)+ulL(:,1:nelL),nelL,shldL,shls,nen,0,0,be);
           end
                    
            %Evaluate tangent and normal vectors
                T1L = XsL(:,1);
                [Tm1L, Tu1L] = VecNormalize(T1L);
                T2L = XsL(:,2);
                [Tm2L, Tu2L] = VecNormalize(T2L);
                T3L = VecCrossProd(T1L,T2L);
                [Tm3L, Tu3L] = VecNormalize(T3L);
                
                C1L = Wgt*Tm3L;

                t1L = xsL(:,1);
                [tm1L, tu1L] = VecNormalize(t1L);
                t2L = xsL(:,2);
                [tm2L, tu2L] = VecNormalize(t2L);
                t3L = VecCrossProd(t1L,t2L);
                [tm3L, tu3L] = VecNormalize(t3L);

                c1L = Wgt*tm3L;
                
            ebL = ebL + C1L*ebeL;
            ebR = ebR + C1L*ebeR;
            intedge = intedge + C1L;

        end   
        
        pencoeff = 1;  %VMS
        edgeK = tauL*ebL^2 + tauR*ebR^2;
        gamL = ebL^2*(edgeK\tauL);
        gamR = ebR^2*(edgeK\tauR);
        gamL_list(iter+1,:,inter) = [gamL(1,1),gamL(1,2),gamL(1,3),gamL(2,1),gamL(2,2),gamL(2,3),gamL(3,1),gamL(3,2),gamL(3,3)];
        gamR_list(iter+1,:,inter) = [gamR(1,1),gamR(1,2),gamR(1,3),gamR(2,1),gamR(2,2),gamR(2,3),gamR(3,1),gamR(3,2),gamR(3,3)];
        ep = pencoeff*intedge*inv(edgeK); 
        ep_list(iter+1,:,inter) = [ep(1,1),ep(1,2),ep(1,3),ep(2,1),ep(2,2),ep(2,3),ep(3,1),ep(3,2),ep(3,3)];

   else
        gamL = [gamL_list(iterset+1,1,inter) gamL_list(iterset+1,2,inter) gamL_list(iterset+1,3,inter)
                gamL_list(iterset+1,4,inter) gamL_list(iterset+1,5,inter) gamL_list(iterset+1,6,inter)
                gamL_list(iterset+1,7,inter) gamL_list(iterset+1,8,inter) gamL_list(iterset+1,9,inter)];
        gamR = [gamR_list(iterset+1,1,inter) gamR_list(iterset+1,2,inter) gamR_list(iterset+1,3,inter)
                gamR_list(iterset+1,4,inter) gamR_list(iterset+1,5,inter) gamR_list(iterset+1,6,inter)
                gamR_list(iterset+1,7,inter) gamR_list(iterset+1,8,inter) gamR_list(iterset+1,9,inter)];
        ep = [ep_list(iterset+1,1,inter) ep_list(iterset+1,2,inter) ep_list(iterset+1,3,inter)
              ep_list(iterset+1,4,inter) ep_list(iterset+1,5,inter) ep_list(iterset+1,6,inter)
              ep_list(iterset+1,7,inter) ep_list(iterset+1,8,inter) ep_list(iterset+1,9,inter)];
   end
       gamL =0.5*eye(3,3) ;
       gamR =0.5*eye(3,3) ; 
       ep =  rp*eye(3,3) ;   
   

         for l = 1:lint

                ll = ll + 1;

                % Evaluate  basis functions at integration points
                if nelL == 4 || nelL == 10
                  [shlL,shldL,shls,be] = shltt([rL sL tL],nelL,nel2L,0,0);
                  [PxyL,shgsL,JdetL,bubbleL,XsL] = shgtt(xlL(:,1:nelL),nelL,shldL,shls,nen,0,0,be);                 
                  [QxyL,shgsL,JdetL,bubbleL,xsL] = shgtt(xlL(:,1:nelL)+ulL(:,1:nelL),nelL,shldL,shls,nen,0,0,be);
                else
                  [Wgt,ss] = intpntb(l,lint,ib);
                  [shlL,shldL,shls,be] = shlb(ss,nelL,nel2L,0,0);
                  [PxyL,shgsL,JdetL,bubbleL,XsL] = shgb(xlL(:,1:nelL),nelL,shldL,shls,nen,0,0,be); 
                  [QxyL,shgsL,JdetL,bubbleL,xsL] = shgb(xlL(:,1:nelL)+ulL(:,1:nelL),nelL,shldL,shls,nen,0,0,be);
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
                  [shlR,shldR,shls,be] = shltt([rR sR tR],nelR,nel2R,0,0);
                  [PxyR,shgsR,JdetR,bubbleR,XsR] = shgtt(xlR(:,1:nelR),nelR,shldR,shls,nen,0,0,be);                  
                  [QxyR,shgsR,JdetR,bubbleR,xsR] = shgtt(xlR(:,1:nelR)+ulR(:,1:nelR),nelR,shldR,shls,nen,0,0,be);
                  
                else
                  [shlR,shldR,shls,be] = shlb([rR sR tR],nelR,nel2R,0,0);
                  [PxyR,shgsR,JdetR,bubbleR,XsR] = shgb(xlR(:,1:nelR),nelR,shldR,shls,nen,0,0,be);     
                  [QxyR,shgsR,JdetR,bubbleR,xsR] = shgb(xlR(:,1:nelR)+ulR(:,1:nelR),nelR,shldR,shls,nen,0,0,be);
                end


 % initial stress sigma and material bulk term c_ijkl
           for mm = 1:nelL 
 NmatL(:,3*mm-2:3*mm) = [shlL(mm,1)     0          0
                           0        shlL(mm,1)     0
                           0            0       shlL(mm,1) ];
 BmatL(:,3*mm-2:3*mm) = [QxyL(mm,1) 0         0         
                        0         QxyL(mm,2) 0         
                        0         0         QxyL(mm,3) 
                        QxyL(mm,2) QxyL(mm,1) 0         
                        0         QxyL(mm,3) QxyL(mm,2) 
                        QxyL(mm,3) 0         QxyL(mm,1) 
                        QxyL(mm,2) -QxyL(mm,1) 0         
                        0         QxyL(mm,3) -QxyL(mm,2) 
                        -QxyL(mm,3) 0         QxyL(mm,1) ];                
 BmatL1(:,3*mm-2:3*mm) = [QxyL(mm,1) 0         0         
                        0         QxyL(mm,2) 0         
                        0         0         QxyL(mm,3) 
                        QxyL(mm,2) QxyL(mm,1) 0         
                        0         QxyL(mm,3) QxyL(mm,2) 
                        QxyL(mm,3) 0         QxyL(mm,1) ];
BmatL2(:,3*mm-2:3*mm)=[QxyL(mm,1)  0             0         
                       QxyL(mm,2)  0             0         
                       QxyL(mm,3)  0             0          
                        0     QxyL(mm,1)         0         
                        0     QxyL(mm,2)         0         
                        0     QxyL(mm,3)         0
                        0         0         QxyL(mm,1)         
                        0         0         QxyL(mm,2)         
                        0         0         QxyL(mm,3)];
 BmatL4(:,3*mm-2:3*mm) = [PxyL(mm,1) 0         0         
                        0         PxyL(mm,2) 0         
                        0         0         PxyL(mm,3) 
                        PxyL(mm,2) PxyL(mm,1) 0         
                        0         PxyL(mm,3) PxyL(mm,2) 
                        PxyL(mm,3) 0         PxyL(mm,1) 
                        PxyL(mm,2) -PxyL(mm,1) 0         
                        0         PxyL(mm,3) -PxyL(mm,2) 
                        -PxyL(mm,3) 0         PxyL(mm,1) ];   
BmatL5(:,3*mm-2:3*mm)=[PxyL(mm,1)  0             0         
                       PxyL(mm,2)  0             0         
                       PxyL(mm,3)  0             0          
                        0     PxyL(mm,1)         0         
                        0     PxyL(mm,2)         0         
                        0     PxyL(mm,3)         0
                        0         0         PxyL(mm,1)         
                        0         0         PxyL(mm,2)         
                        0         0         PxyL(mm,3)];                    
           end
            
           for mm = 1:nelR    
 NmatR(:,3*mm-2:3*mm) = [shlR(mm,1)     0          0
                           0        shlR(mm,1)     0
                           0            0       shlR(mm,1) ];
 BmatR(:,3*mm-2:3*mm) = [QxyR(mm,1) 0         0         
                        0         QxyR(mm,2) 0         
                        0         0         QxyR(mm,3) 
                        QxyR(mm,2) QxyR(mm,1) 0         
                        0         QxyR(mm,3) QxyR(mm,2) 
                        QxyR(mm,3) 0         QxyR(mm,1) 
                        QxyR(mm,2) -QxyR(mm,1) 0         
                        0         QxyR(mm,3) -QxyR(mm,2) 
                        -QxyR(mm,3) 0         QxyR(mm,1) ];
 BmatR1(:,3*mm-2:3*mm) = [QxyR(mm,1) 0         0         
                        0         QxyR(mm,2) 0         
                        0         0         QxyR(mm,3) 
                        QxyR(mm,2) QxyR(mm,1) 0         
                        0         QxyR(mm,3) QxyR(mm,2) 
                        QxyR(mm,3) 0         QxyR(mm,1) ];
 BmatR2(:,3*mm-2:3*mm)=[QxyR(mm,1)  0             0         
                       QxyR(mm,2)  0             0         
                       QxyR(mm,3)  0             0          
                        0     QxyR(mm,1)         0         
                        0     QxyR(mm,2)         0         
                        0     QxyR(mm,3)         0
                        0         0         QxyR(mm,1)         
                        0         0         QxyR(mm,2)         
                        0         0         QxyR(mm,3)];  
BmatR4(:,3*mm-2:3*mm) = [PxyR(mm,1) 0         0         
                        0         PxyR(mm,2) 0         
                        0         0         PxyR(mm,3) 
                        PxyR(mm,2) PxyR(mm,1) 0         
                        0         PxyR(mm,3) PxyR(mm,2) 
                        PxyR(mm,3) 0         PxyR(mm,1) 
                        PxyR(mm,2) -PxyR(mm,1) 0         
                        0         PxyR(mm,3) -PxyR(mm,2) 
                        -PxyR(mm,3) 0         PxyR(mm,1) ];                   
 BmatR5(:,3*mm-2:3*mm)=[PxyR(mm,1)  0             0         
                       PxyR(mm,2)  0             0         
                       PxyR(mm,3)  0             0          
                        0     PxyR(mm,1)         0         
                        0     PxyR(mm,2)         0         
                        0     PxyR(mm,3)         0
                        0         0         PxyR(mm,1)         
                        0         0         PxyR(mm,2)         
                        0         0         PxyR(mm,3)];  
           end  
    %plug in the -ul to get the inverse of deformation gradient fi as 'F' in this function             
            [fiL,JxXL,FL] = kine3d(QxyL,-ulL,nelL,0); %this is equivalent to ikine2d 
            JxXL = 1/JxXL; %this is equivalent to ikine2d 
    %         [F,JxX,fi] = ikine2d(Qxy,ul,nel); %this is equivalent to above two lines
            JdetxL = JdetL/JxXL;               
%             C1L = Wgt*JdetL;
            [fiR,JxXR,FR] = kine3d(QxyR,-ulR,nelR,0); %this is equivalent to ikine2d
            JxXR = 1/JxXR; %this is equivalent to ikine2d
    %         [F,JxX,fi] = ikine2d(Qxy,ul,nel); %this is equivalent to above two lines
            JdetxR = JdetR/JxXL;                          
%             C1R = Wgt*JdetR;         this is for volume integral  
          
            [sigma2L, cmatL] = SigmaCmat3i(FL,JxXL,matepropL,lamL);
            P = [fiL(1,1)^2 fiL(1,2)^2 fiL(1,3)^2 2*fiL(1,1)*fiL(1,2) 2*fiL(1,2)*fiL(1,3) 2*fiL(1,3)*fiL(1,1)  zeros(1,3)
                 fiL(2,1)^2 fiL(2,2)^2 fiL(2,3)^2 2*fiL(2,1)*fiL(2,2) 2*fiL(2,2)*fiL(2,3) 2*fiL(2,3)*fiL(2,1)  zeros(1,3)
                 fiL(3,1)^2 fiL(3,2)^2 fiL(3,3)^2 2*fiL(3,1)*fiL(3,2) 2*fiL(3,2)*fiL(3,3) 2*fiL(3,3)*fiL(3,1)  zeros(1,3)                 
                 fiL(1,1)*fiL(2,1) fiL(1,2)*fiL(2,2) fiL(1,3)*fiL(2,3) fiL(1,1)*fiL(2,2)+fiL(1,2)*fiL(2,1) fiL(1,2)*fiL(2,3)+fiL(1,3)*fiL(2,2) fiL(1,1)*fiL(2,3)+fiL(1,3)*fiL(2,1) zeros(1,3)
                 fiL(2,1)*fiL(3,1) fiL(2,2)*fiL(3,2) fiL(2,3)*fiL(3,3) fiL(2,1)*fiL(3,2)+fiL(2,2)*fiL(3,1) fiL(2,2)*fiL(3,3)+fiL(2,3)*fiL(3,2) fiL(2,1)*fiL(3,3)+fiL(2,3)*fiL(3,1) zeros(1,3)                 
                 fiL(3,1)*fiL(1,1) fiL(3,2)*fiL(1,2) fiL(3,3)*fiL(1,3) fiL(3,1)*fiL(1,2)+fiL(3,2)*fiL(1,1) fiL(3,2)*fiL(1,3)+fiL(3,3)*fiL(1,2) fiL(3,1)*fiL(1,3)+fiL(3,3)*fiL(1,1) zeros(1,3)                  
                 zeros(3,9)];
            SPKL = P*sigma2L*JxXL;
            CmatL = P*cmatL*JxXL*P';            
            SmatL1=[SPKL(1), SPKL(4), SPKL(6)
                    SPKL(4), SPKL(2), SPKL(5)
                    SPKL(6), SPKL(5), SPKL(3)];        
            
            SmatL = ...
[    sigma2L(1),        0,        0,           sigma2L(4)/2,                 0,           sigma2L(6)/2,           sigma2L(4)/2,                 0,          -sigma2L(6)/2
         0,    sigma2L(2),        0,           sigma2L(4)/2,           sigma2L(5)/2,                 0,          -sigma2L(4)/2,           sigma2L(5)/2,                 0
         0,        0,    sigma2L(3),                 0,           sigma2L(5)/2,           sigma2L(6)/2,                 0,          -sigma2L(5)/2,           sigma2L(6)/2
   sigma2L(4)/2,  sigma2L(4)/2,        0, sigma2L(1)/4 + sigma2L(2)/4,           sigma2L(6)/4,           sigma2L(5)/4, sigma2L(2)/4 - sigma2L(1)/4,           sigma2L(6)/4,          -sigma2L(5)/4
         0,  sigma2L(5)/2,  sigma2L(5)/2,           sigma2L(6)/4, sigma2L(2)/4 + sigma2L(3)/4,           sigma2L(4)/4,          -sigma2L(6)/4, sigma2L(3)/4 - sigma2L(2)/4,           sigma2L(4)/4
   sigma2L(6)/2,        0,  sigma2L(6)/2,           sigma2L(5)/4,           sigma2L(4)/4, sigma2L(1)/4 + sigma2L(3)/4,           sigma2L(5)/4,          -sigma2L(4)/4, sigma2L(1)/4 - sigma2L(3)/4
   sigma2L(4)/2, -sigma2L(4)/2,        0, sigma2L(2)/4 - sigma2L(1)/4,          -sigma2L(6)/4,           sigma2L(5)/4, sigma2L(1)/4 + sigma2L(2)/4,          -sigma2L(6)/4,          -sigma2L(5)/4
         0,  sigma2L(5)/2, -sigma2L(5)/2,           sigma2L(6)/4, sigma2L(3)/4 - sigma2L(2)/4,          -sigma2L(4)/4,          -sigma2L(6)/4, sigma2L(2)/4 + sigma2L(3)/4,          -sigma2L(4)/4
  -sigma2L(6)/2,        0,  sigma2L(6)/2,          -sigma2L(5)/4,           sigma2L(4)/4, sigma2L(1)/4 - sigma2L(3)/4,          -sigma2L(5)/4,          -sigma2L(4)/4, sigma2L(1)/4 + sigma2L(3)/4];
 
            smatL1=[sigma2L(1), sigma2L(4), sigma2L(6)
                    sigma2L(4), sigma2L(2), sigma2L(5)
                    sigma2L(6), sigma2L(5), sigma2L(3)];      
            
            [sigma2R, cmatR] = SigmaCmat3i(FR,JxXR,matepropR,lamR);

            P = [fiR(1,1)^2 fiR(1,2)^2 fiR(1,3)^2 2*fiR(1,1)*fiR(1,2) 2*fiR(1,2)*fiR(1,3) 2*fiR(1,3)*fiR(1,1)  zeros(1,3)
                 fiR(2,1)^2 fiR(2,2)^2 fiR(2,3)^2 2*fiR(2,1)*fiR(2,2) 2*fiR(2,2)*fiR(2,3) 2*fiR(2,3)*fiR(2,1)  zeros(1,3)
                 fiR(3,1)^2 fiR(3,2)^2 fiR(3,3)^2 2*fiR(3,1)*fiR(3,2) 2*fiR(3,2)*fiR(3,3) 2*fiR(3,3)*fiR(3,1)  zeros(1,3)                 
                 fiR(1,1)*fiR(2,1) fiR(1,2)*fiR(2,2) fiR(1,3)*fiR(2,3) fiR(1,1)*fiR(2,2)+fiR(1,2)*fiR(2,1) fiR(1,2)*fiR(2,3)+fiR(1,3)*fiR(2,2) fiR(1,1)*fiR(2,3)+fiR(1,3)*fiR(2,1) zeros(1,3)
                 fiR(2,1)*fiR(3,1) fiR(2,2)*fiR(3,2) fiR(2,3)*fiR(3,3) fiR(2,1)*fiR(3,2)+fiR(2,2)*fiR(3,1) fiR(2,2)*fiR(3,3)+fiR(2,3)*fiR(3,2) fiR(2,1)*fiR(3,3)+fiR(2,3)*fiR(3,1) zeros(1,3)                 
                 fiR(3,1)*fiR(1,1) fiR(3,2)*fiR(1,2) fiR(3,3)*fiR(1,3) fiR(3,1)*fiR(1,2)+fiR(3,2)*fiR(1,1) fiR(3,2)*fiR(1,3)+fiR(3,3)*fiR(1,2) fiR(3,1)*fiR(1,3)+fiR(3,3)*fiR(1,1) zeros(1,3)                  
                 zeros(3,9)];
            SPKR = P*sigma2R*JxXR;
            CmatR = P*cmatR*JxXR*P';
            SmatR1=[SPKR(1), SPKR(4), SPKR(6)
                    SPKR(4), SPKR(2), SPKR(5)
                    SPKR(6), SPKR(5), SPKR(3)];   
                
            SmatR = ...
[    sigma2R(1),        0,        0,           sigma2R(4)/2,                 0,           sigma2R(6)/2,           sigma2R(4)/2,                 0,          -sigma2R(6)/2
         0,    sigma2R(2),        0,           sigma2R(4)/2,           sigma2R(5)/2,                 0,          -sigma2R(4)/2,           sigma2R(5)/2,                 0
         0,        0,    sigma2R(3),                 0,           sigma2R(5)/2,           sigma2R(6)/2,                 0,          -sigma2R(5)/2,           sigma2R(6)/2
   sigma2R(4)/2,  sigma2R(4)/2,        0, sigma2R(1)/4 + sigma2R(2)/4,           sigma2R(6)/4,           sigma2R(5)/4, sigma2R(2)/4 - sigma2R(1)/4,           sigma2R(6)/4,          -sigma2R(5)/4
         0,  sigma2R(5)/2,  sigma2R(5)/2,           sigma2R(6)/4, sigma2R(2)/4 + sigma2R(3)/4,           sigma2R(4)/4,          -sigma2R(6)/4, sigma2R(3)/4 - sigma2R(2)/4,           sigma2R(4)/4
   sigma2R(6)/2,        0,  sigma2R(6)/2,           sigma2R(5)/4,           sigma2R(4)/4, sigma2R(1)/4 + sigma2R(3)/4,           sigma2R(5)/4,          -sigma2R(4)/4, sigma2R(1)/4 - sigma2R(3)/4
   sigma2R(4)/2, -sigma2R(4)/2,        0, sigma2R(2)/4 - sigma2R(1)/4,          -sigma2R(6)/4,           sigma2R(5)/4, sigma2R(1)/4 + sigma2R(2)/4,          -sigma2R(6)/4,          -sigma2R(5)/4
         0,  sigma2R(5)/2, -sigma2R(5)/2,           sigma2R(6)/4, sigma2R(3)/4 - sigma2R(2)/4,          -sigma2R(4)/4,          -sigma2R(6)/4, sigma2R(2)/4 + sigma2R(3)/4,          -sigma2R(4)/4
  -sigma2R(6)/2,        0,  sigma2R(6)/2,          -sigma2R(5)/4,           sigma2R(4)/4, sigma2R(1)/4 - sigma2R(3)/4,          -sigma2R(5)/4,          -sigma2R(4)/4, sigma2R(1)/4 + sigma2R(3)/4];
            smatR1=[sigma2R(1), sigma2R(4), sigma2R(6)
                    sigma2R(4), sigma2R(2), sigma2R(5)
                    sigma2R(6), sigma2R(5), sigma2R(3)]; 
       

                %Evaluate tangent and normal vectors
                % jecobian j=dx/dkesi for surface integral for current
                % config
                t1L = xsL(:,1);
                [tm1L, tu1L] = VecNormalize(t1L);
                t2L = xsL(:,2);
                [tm2L, tu2L] = VecNormalize(t2L);
                t3L = VecCrossProd(t1L,t2L);
                [tm3L, tu3L] = VecNormalize(t3L);
                t1R = xsR(:,1);
                [tm1R, tu1R] = VecNormalize(t1R);
                t2R = xsR(:,2);
                [tm2R, tu2R] = VecNormalize(t2R);
                t3R = VecCrossProd(t1R,t2R);
                [tm3R, tu3R] = VecNormalize(t3R);
                c1L = Wgt*tm3L;
                c1R = Wgt*tm3R;
                % jecobian j=dx/dkesi for surface integral for material
                % config
                T1L = XsL(:,1);
                [Tm1L, Tu1L] = VecNormalize(T1L);
                T2L = XsL(:,2);
                [Tm2L, Tu2L] = VecNormalize(T2L);
                T3L = VecCrossProd(T1L,T2L);
                [Tm3L, Tu3L] = VecNormalize(T3L);
                T1R = XsR(:,1);
                [Tm1R, Tu1R] = VecNormalize(T1R);
                T2R = XsR(:,2);
                [Tm2R, Tu2R] = VecNormalize(T2R);
                T3R = VecCrossProd(T1R,T2R);
                [Tm3R, Tu3R] = VecNormalize(T3R);
                C1L = Wgt*Tm3L;
                C1R = Wgt*Tm3R;                
                %normal vectors
                nLx = -tu3L(1);
                nLy = -tu3L(2);
                nLz = -tu3L(3);
                nRx = -tu3R(1);
                nRy = -tu3R(2);
                nRz = -tu3R(3);
                NLx = -Tu3L(1);
                NLy = -Tu3L(2);
                NLz = -Tu3L(3);
                NRx = -Tu3R(1);
                NRy = -Tu3R(2);
                NRz = -Tu3R(3);                
                %tagent vectors
                tLx = tu1L(1);
                tLy = tu1L(2);
                tLz = tu1L(3);
                tRx = tu1R(1);
                tRy = tu1R(2);
                tRz = tu1R(3);                
                nvectL = [nLx 0  0  nLy  0  nLz nLy 0  -nLz
                          0  nLy 0  nLx nLz  0 -nLx nLz  0
                          0   0 nLz  0  nLy nLx  0 -nLy nLx];
                nvectL1= [nLx 0  0  nLy  0  nLz 
                          0  nLy 0  nLx nLz  0 
                          0   0 nLz  0  nLy nLx];
                nvectR = [nRx 0  0  nRy  0  nRz nRy 0  -nRz
                          0  nRy 0  nRx  nRz 0 -nRx nRz  0
                          0   0 nRz  0  nRy nRx  0 -nRy nRx]; 
                nvectR1 = [nRx 0  0  nRy  0  nRz 
                          0  nRy 0  nRx  nRz 0 
                          0   0 nRz  0  nRy nRx];  
                nvectL2 = [eye(6,6)*nLx zeros(6,6)    zeros(6,6)    eye(6,6)*nLy  zeros(6,6)    eye(6,6)*nLz
                           zeros(6,6)    eye(6,6)*nLy zeros(6,6)    eye(6,6)*nLx  eye(6,6)*nLz  zeros(6,6)
                           zeros(6,6)    zeros(6,6)   eye(6,6)*nLz   zeros(6,6)   eye(6,6)*nLy  eye(6,6)*nLx];
                nvectR2 = [eye(6,6)*nRx zeros(6,6)    zeros(6,6)    eye(6,6)*nRy  zeros(6,6)    eye(6,6)*nRz
                           zeros(6,6)    eye(6,6)*nRy zeros(6,6)    eye(6,6)*nRx  eye(6,6)*nRz  zeros(6,6)
                           zeros(6,6)    zeros(6,6)   eye(6,6)*nRz   zeros(6,6)   eye(6,6)*nRy  eye(6,6)*nRx];
                nvecL = [nLx; nLy; nLz];
                nvecR = [nRx; nRy; nRz];
                NvecL = [NLx; NLy; NLz];
                NvecR = [NRx; NRy; NRz];
                NvectL1= [NLx 0  0  NLy  0  NLz 
                          0  NLy 0  NLx NLz  0 
                          0   0 NLz  0  NLy NLx];     
                NvectR1 = [NRx 0  0  NRy  0  NRz 
                          0  NRy 0  NRx  NRz 0 
                          0   0 NRz  0  NRy NRx];                        
%additional stiffness term

            smatnL=[smatL1*nvecL zeros(ndm,2)
                    zeros(ndm,1) smatL1*nvecL zeros(ndm,1)
                    zeros(ndm,2) smatL1*nvecL];
            
            smatnR=[smatR1*nvecR zeros(ndm,2)
                    zeros(ndm,1) smatR1*nvecR zeros(ndm,1)
                    zeros(ndm,2) smatR1*nvecR];

           cmatnL=gamL*(nvectL1*cmatL(1:6,1:6));
           cmatnR=gamR*(nvectR1*cmatR(1:6,1:6));
           CmatnL=gamL*(NvectL1*CmatL(1:6,1:6));
           CmatnR=gamR*(NvectR1*CmatR(1:6,1:6)); 
           
           cmatnBL=BmatL2'*[cmatnL zeros(3,12)
                    zeros(3,6)   cmatnL   zeros(3,6)
                    zeros(3,12)  cmatnL ]; 
                
           cmatnBR=BmatR2'*[cmatnR    zeros(3,12)
                    zeros(3,6)   cmatnR   zeros(3,6)
                    zeros(3,12)  cmatnR ];

           CmatnBL=BmatL2'*[CmatnL zeros(3,12)
                    zeros(3,6)   CmatnL   zeros(3,6)
                    zeros(3,12)  CmatnL ]; 
                
           CmatnBR=BmatR2'*[CmatnR    zeros(3,12)
                    zeros(3,6)   CmatnR   zeros(3,6)
                    zeros(3,12)  CmatnR ];
                
            term17L=gamL*smatnL'*BmatL2;
            term17R=gamR*smatnR'*BmatR2;

        I2 = eye(3);   
        P2 = [1 0 0  0   0   0    0   0    0
              0 0 0 1/2  0   0  -1/2  0    0
              0 0 0  0   0  1/2   0   0   1/2
              0 0 0 1/2  0   0   1/2  0    0
              0 1 0  0   0   0    0   0    0
              0 0 0  0  1/2  0    0  -1/2  0
              0 0 0  0   0  1/2   0   0   -1/2
              0 0 0  0  1/2  0    0   1/2  0   
              0 0 1  0   0   0    0   0    0] ;
            term32 = smatL1*nvecL;
            term32M = [term32(1)*I2 term32(2)*I2 term32(3)*I2];
            term33L = gamL*term32M*P2*BmatL;
            term32 = smatR1*nvecR;
            term32M = [term32(1)*I2 term32(2)*I2 term32(3)*I2];
            term33R = gamR*term32M*P2*BmatR;
            Term32 = SmatL1*NvecL;
            Term32M = [Term32(1)*I2 Term32(2)*I2 Term32(3)*I2];
            Term33L = gamL*Term32M*P2*BmatL4;
            Term32 = SmatR1*NvecR;
            Term32M = [Term32(1)*I2 Term32(2)*I2 Term32(3)*I2];
            Term33R = gamR*Term32M*P2*BmatR4;
            
            term18L=gamL*nvectL1*cmatL(1:6,1:6)*BmatL1;
            term18R=gamR*nvectR1*cmatR(1:6,1:6)*BmatR1;

            bigFL = [FL(1,1)  0        0      FL(2,1)  0       0      FL(3,1)   0        0  
                     0       FL(1,2)   0          0   FL(2,2)  0        0     FL(3,2)    0
                     0        0      FL(1,3)      0    0      FL(2,3)   0       0      FL(3,3)                     
                     FL(1,2) FL(1,1)   0      FL(2,2) FL(2,1)  0      FL(3,2) FL(3,1)    0
                     0       FL(1,3) FL(1,2)      0   FL(2,3) FL(2,2)   0     FL(3,3)  FL(3,2)
                     FL(1,3)  0      FL(1,1)  FL(2,3)  0      FL(2,1) FL(3,3)   0      FL(3,1)];
            Term18L=gamL*FL*NvectL1*CmatL(1:6,1:6)*bigFL*BmatL5;
            bigFR = [FR(1,1)  0        0      FR(2,1)  0       0      FR(3,1)   0        0  
                     0       FR(1,2)   0          0   FR(2,2)  0        0     FR(3,2)    0
                     0        0      FR(1,3)      0    0      FR(2,3)   0       0      FR(3,3)                     
                     FR(1,2) FR(1,1)   0      FR(2,2) FR(2,1)  0      FR(3,2) FR(3,1)    0
                     0       FR(1,3) FR(1,2)      0   FR(2,3) FR(2,2)   0     FR(3,3)  FR(3,2)
                     FR(1,3)  0      FR(1,1)  FR(2,3)  0      FR(2,1) FR(3,3)   0      FR(3,1)];
            Term18R=gamR*FR*NvectR1*CmatR(1:6,1:6)*bigFR*BmatR5;
            
            term28L=(c1R*gamR*smatR1*nvecR-c1L*gamL*smatL1*nvecL);   %average stress term
            term28R=(c1R*gamR*smatR1*nvecR-c1L*gamL*smatL1*nvecL);  
            
            Term28L=(gamR*FR*SmatR1*NvecR-gamL*FL*SmatL1*NvecL);   %average stress term
            Term28R=(gamR*FR*SmatR1*NvecR-gamL*FL*SmatL1*NvecL);   %average stress term

             rhspulL = reshape(ulL,ndf*nelL,1);
             rhspulR = reshape(ulR,ndf*nelR,1);
                   
            jumpu = NmatR*rhspulR - NmatL*rhspulL;  %jumpu=uR-uL


%% damage part
            % Load history
            damhr = nh1-1+(ll-1)*7;
            dmaxhr = nh1-1+(ll-1)*7+4;
            chat1hr = nh2-1+(ll-1)*7+6;
            chat2hr = nh2-1+(ll-1)*7+7;
            dvec = [hr(damhr+1); hr(damhr+2); hr(damhr+3)];
            dmax = hr(dmaxhr); 
            bonddbond1 = hr(chat1hr);
            chatter = hr(chat2hr);
                
            NBmatL = (Term33L+Term18L) - rp*NmatL;
            NBmatR = (Term33R+Term18R) - rp*NmatR;
            rhspulL = reshape(ulL,ndf*nelL,1);
            rhspulR = reshape(ulR,ndf*nelR,1);
            
%             if iprob == 3
%                 cvec = -0.0022*[1; 1; 1; 0; 0; 0];
%             end
            
            tvtr = -Term28L; %average stress
            normsig = NvecL'*tvtr;
            dn = jumpu'*NvecL;
            tn = normsig+rp*dn;
            
            tss = tvtr + rp*jumpu - tn*NvecL;              %shear traction v
            tt = sqrt(tss'*tss);                          %shear traction
            
            if tn >= 0 % tension
                
                tvec = tvtr + rp*jumpu;
                tb = sqrt(tn^2+beta2*tt^2);
                
            else % compression
                
                tvec = tss;
                tb = beta*tt;
                
            end
            
            normtvec = sqrt(tvec'*tvec);
            if dmax >= dc
                psik = 0;
            else
                psik = sigmax - Hc*dmax;
            end

%%            
            if tn >= 0 % tension
              
              if dmax == 0 && sqrt(tn^2+tt^2/beta2) <= sigmax % perfect adhesion
                
                dvec = zeros(3,1);
                dddt = zeros(3,3);
                dmax = 0;
                bonddbond2 = 0;
                
              elseif dmax == 0 && sqrt(tn^2+tt^2/beta2) > sigmax % damage from initial adhesion
                
                [dvec,dtdd,deq] = damageNR(tvec,NvecL,dvec,normtvec,sigmax,dc,rp,beta2);
                if deq >= dc % opening/sliding
                    dvec = tvec/rp;
                    dddt = 1/rp*eye(3);
                    dn = dvec'*	NvecL;
                    dss = dvec - dn*NvecL;
                    dt = sqrt(dss'*dss);
                    deq = sqrt(dn^2+beta2*dt^2);
                else
                    dddt = inv(dtdd);
                end
                dmax = max(dmax,deq);
                bonddbond2 = -1;
                
              elseif dmax > 0 && tb < rp*dmax % crumpling
                
                dvec = tvec/rp;
                dddt = 1/rp*eye(3);
                dn = dvec'*NvecL;
                dss = dvec - dn*NvecL;
                dt = sqrt(dss'*dss);
                deq = sqrt(dn^2+beta2*dt^2);
                dmax = max(dmax,deq);
                bonddbond2 = -1;
                
              elseif dmax > 0 && rp*dmax <= tb && (tn^2/(rp*dmax + psik)^2 + beta2*tt^2/(rp*dmax + beta2*psik)^2 - 1) <= dtol %unloading
                
                [dvec,dddt,deq] = unloadingNR(tn,tt,tss,NvecL,dmax,sigmax,dc,rp,beta2);
                bonddbond2 = -1;
                
              elseif (tn^2/(rp*dmax + psik)^2 + beta2*tt^2/(rp*dmax + beta2*psik) - 1)^2 > dtol % damage
                
                [dvec,dtdd,deq] = damageNR(tvec,NvecL,dvec,normtvec,sigmax,dc,rp,beta2);
                if deq >= dc % opening/sliding
                    dvec = tvec/rp;
                    dddt = 1/rp*eye(3);
                    dn = dvec'*NvecL;
                    dss = dvec - dn*NvecL;
                    dt = sqrt(dss'*dss);
                    deq = sqrt(dn^2+beta2*dt^2);
                else
                    dddt = inv(dtdd);
                end
                dmax = max(dmax,deq);
                bonddbond2 = -1;
                
              end
              
            else % compression
                 
              if dmax == 0 && tt/beta <= sigmax % perfect adhesion
                
                dvec = zeros(3,1);
                dddt = zeros(3,3);
                dmax = 0;
                bonddbond2 = 0;
                
              elseif dmax > 0 && tb < rp*dmax % crumpling
                
                dvec = tvec/rp;
                dddt = 1/rp*(eye(3) - NvecL*NvecL');
                deq = beta*sqrt(dvec'*dvec);
                dmax = max(dmax,deq);
                bonddbond2 = -1;
                
              elseif dmax > 0 && rp*dmax <= tb && (tb - (rp*dmax + beta2*psik)) <= dtol %unloading
                
                ntilda = tss/tt;
                dvec = (dmax/beta)*ntilda;
                dddt = (dmax/beta)/tt*(eye(3) - ntilda*ntilda' - NvecL*NvecL');
%                 deq = sqrt(dvec'*dvec);
%                 dmax = max(dmax,deq);
                bonddbond2 = -1;
                
              elseif dtol < (tb - (rp*dmax + beta2*psik)) % damage
                
                deq = beta*(tt - beta*sigmax)/(rp - beta2*Hc); % assume only damage
                ntilda = tss/tt;
                if deq >= dc % opening/sliding
                    dvec = tvec/rp;
                    deq = beta*sqrt(dvec'*dvec);
                    dddt = 1/rp*(eye(3) - NvecL*NvecL');
                else % only damage
                    dvec = (deq/beta)*ntilda;
                    dddt = 1/(1*(rp - beta2*Hc))*(ntilda*ntilda') + (deq/beta)/tt*(eye(3) - ntilda*ntilda' - NvecL*NvecL');
                end
                dmax = deq;
                bonddbond2 = -1;
                
              end
                
            end
            
%%          additional stiffness terms
              u_k = jumpu - dvec;
              term5L=cmatnBL*[eye(6,6)*(u_k(1)) eye(6,6)*(u_k(2)) eye(6,6)*(u_k(3))]'*BmatL1;
              term5R=cmatnBR*[eye(6,6)*(u_k(1)) eye(6,6)*(u_k(2)) eye(6,6)*(u_k(3))]'*BmatR1;

             [dmatL]=dmat3(JxXL,mateprop,lam);
             [dmatR]=dmat3(JxXR,mateprop,lam);   
             dmatL2 = [eye(6,6)*(u_k(1)) eye(6,6)*(u_k(2)) eye(6,6)*(u_k(3))]*[eye(6,6)*gamL(1,1) eye(6,6)*gamL(1,2) eye(6,6)*gamL(1,3)
                                                                                                             eye(6,6)*gamL(2,1) eye(6,6)*gamL(2,2) eye(6,6)*gamL(2,3)
                                                                                                             eye(6,6)*gamL(3,1) eye(6,6)*gamL(3,2) eye(6,6)*gamL(3,3)]*nvectL2*dmatL/JxXL;
             dmatR2 = [eye(6,6)*(u_k(1)) eye(6,6)*(u_k(2)) eye(6,6)*(u_k(3))]*[eye(6,6)*gamR(1,1) eye(6,6)*gamR(1,2) eye(6,6)*gamR(1,3)
                                                                                                             eye(6,6)*gamR(2,1) eye(6,6)*gamR(2,2) eye(6,6)*gamR(2,3)
                                                                                                             eye(6,6)*gamR(3,1) eye(6,6)*gamR(3,2) eye(6,6)*gamR(3,3)]*nvectR2*dmatR/JxXR;
             
             term7L = BmatL1'*dmatL2*BmatL1;
             term7R = BmatR1'*dmatR2*BmatR1;             

            sig8L = [cmatnL(1,:) zeros(1,12) 
                     zeros(1,6) cmatnL(2,:) zeros(1,6) 
                     zeros(1,12)             cmatnL(3,:) ];  
            sig8L1 = sig8L*[ones(6,6)*(u_k(1)) eye(6,6)*(u_k(2)) eye(6,6)*(u_k(3))]';  

            sig8L2 = [sum(sig8L1(:,1)) sum(sig8L1(:,2)) sum(sig8L1(:,3)) sum(sig8L1(:,4)) sum(sig8L1(:,5)) sum(sig8L1(:,6))];
                                    
            sig8L3= ...
[    sig8L2(1),        0,        0,           sig8L2(4)/2,                 0,           sig8L2(6)/2,           sig8L2(4)/2,                 0,          -sig8L2(6)/2
         0,    sig8L2(2),        0,           sig8L2(4)/2,           sig8L2(5)/2,                 0,          -sig8L2(4)/2,           sig8L2(5)/2,                 0
         0,        0,    sig8L2(3),                 0,           sig8L2(5)/2,           sig8L2(6)/2,                 0,          -sig8L2(5)/2,           sig8L2(6)/2
   sig8L2(4)/2,  sig8L2(4)/2,        0, sig8L2(1)/4 + sig8L2(2)/4,           sig8L2(6)/4,           sig8L2(5)/4, sig8L2(2)/4 - sig8L2(1)/4,           sig8L2(6)/4,          -sig8L2(5)/4
         0,  sig8L2(5)/2,  sig8L2(5)/2,           sig8L2(6)/4, sig8L2(2)/4 + sig8L2(3)/4,           sig8L2(4)/4,          -sig8L2(6)/4, sig8L2(3)/4 - sig8L2(2)/4,           sig8L2(4)/4
   sig8L2(6)/2,        0,  sig8L2(6)/2,           sig8L2(5)/4,           sig8L2(4)/4, sig8L2(1)/4 + sig8L2(3)/4,           sig8L2(5)/4,          -sig8L2(4)/4, sig8L2(1)/4 - sig8L2(3)/4
   sig8L2(4)/2, -sig8L2(4)/2,        0, sig8L2(2)/4 - sig8L2(1)/4,          -sig8L2(6)/4,           sig8L2(5)/4, sig8L2(1)/4 + sig8L2(2)/4,          -sig8L2(6)/4,          -sig8L2(5)/4
         0,  sig8L2(5)/2, -sig8L2(5)/2,           sig8L2(6)/4, sig8L2(3)/4 - sig8L2(2)/4,          -sig8L2(4)/4,          -sig8L2(6)/4, sig8L2(2)/4 + sig8L2(3)/4,          -sig8L2(4)/4
  -sig8L2(6)/2,        0,  sig8L2(6)/2,          -sig8L2(5)/4,           sig8L2(4)/4, sig8L2(1)/4 - sig8L2(3)/4,          -sig8L2(5)/4,          -sig8L2(4)/4, sig8L2(1)/4 + sig8L2(3)/4];                   
                                   
            sig8R = [cmatnR(1,:) zeros(1,12) 
                     zeros(1,6) cmatnR(2,:) zeros(1,6) 
                     zeros(1,12)             cmatnR(3,:) ];  
            sig8R1 = sig8R*[ones(6,6)*(u_k(1)) eye(6,6)*(u_k(2)) eye(6,6)*(u_k(3))]';  

            sig8R2 = [sum(sig8R1(:,1)) sum(sig8R1(:,2)) sum(sig8R1(:,3)) sum(sig8R1(:,4)) sum(sig8R1(:,5)) sum(sig8R1(:,6))];
                                    
            sig8R3= ...
[    sig8R2(1),        0,        0,           sig8R2(4)/2,                 0,           sig8R2(6)/2,           sig8R2(4)/2,                 0,          -sig8R2(6)/2
         0,    sig8R2(2),        0,           sig8R2(4)/2,           sig8R2(5)/2,                 0,          -sig8R2(4)/2,           sig8R2(5)/2,                 0
         0,        0,    sig8R2(3),                 0,           sig8R2(5)/2,           sig8R2(6)/2,                 0,          -sig8R2(5)/2,           sig8R2(6)/2
   sig8R2(4)/2,  sig8R2(4)/2,        0, sig8R2(1)/4 + sig8R2(2)/4,           sig8R2(6)/4,           sig8R2(5)/4, sig8R2(2)/4 - sig8R2(1)/4,           sig8R2(6)/4,          -sig8R2(5)/4
         0,  sig8R2(5)/2,  sig8R2(5)/2,           sig8R2(6)/4, sig8R2(2)/4 + sig8R2(3)/4,           sig8R2(4)/4,          -sig8R2(6)/4, sig8R2(3)/4 - sig8R2(2)/4,           sig8R2(4)/4
   sig8R2(6)/2,        0,  sig8R2(6)/2,           sig8R2(5)/4,           sig8R2(4)/4, sig8R2(1)/4 + sig8R2(3)/4,           sig8R2(5)/4,          -sig8R2(4)/4, sig8R2(1)/4 - sig8R2(3)/4
   sig8R2(4)/2, -sig8R2(4)/2,        0, sig8R2(2)/4 - sig8R2(1)/4,          -sig8R2(6)/4,           sig8R2(5)/4, sig8R2(1)/4 + sig8R2(2)/4,          -sig8R2(6)/4,          -sig8R2(5)/4
         0,  sig8R2(5)/2, -sig8R2(5)/2,           sig8R2(6)/4, sig8R2(3)/4 - sig8R2(2)/4,          -sig8R2(4)/4,          -sig8R2(6)/4, sig8R2(2)/4 + sig8R2(3)/4,          -sig8R2(4)/4
  -sig8R2(6)/2,        0,  sig8R2(6)/2,          -sig8R2(5)/4,           sig8R2(4)/4, sig8R2(1)/4 - sig8R2(3)/4,          -sig8R2(5)/4,          -sig8R2(4)/4, sig8R2(1)/4 + sig8R2(3)/4];              

             term8L=BmatL'*sig8L3*BmatL;   
             term8R=BmatR'*sig8R3*BmatR;
             
             
            ElemFL = ElemFL - ( - C1L*NmatL'*(tvtr + rp*(jumpu - dvec)) + C1L*(Term33L'+Term18L')*(jumpu - dvec));
            ElemFR = ElemFR - ( + C1L*NmatR'*(tvtr + rp*(jumpu - dvec)) + -C1R*(Term33R'+Term18R')*(jumpu - dvec));

             ElemKLL = ElemKLL - c1L*(term17L'+term18L')*NmatL  -  c1L*NmatL'*(term17L+term18L);
             ElemKLR = ElemKLR + c1L*(term17L'+term18L')*NmatR  +  c1R*NmatL'*(term17R+term18R);
             ElemKRL = ElemKRL + c1R*(term17R'+term18R')*NmatL  +  c1L*NmatR'*(term17L+term18L);
             ElemKRR = ElemKRR - c1R*(term18R'+term17R')*NmatR  -  c1R*NmatR'*(term17R+term18R);
             
             ElemKLL = ElemKLL + C1L*(NmatL'*rp*NmatL);
             ElemKLR = ElemKLR - C1R*(NmatL'*rp*NmatR);
             ElemKRL = ElemKRL - C1L*(NmatR'*rp*NmatL);
             ElemKRR = ElemKRR + C1R*(NmatR'*rp*NmatR); 

%              ElemKLL = ElemKLL - C1L*NmatL'*(Term33L+Term18L) - C1L*NBmatL'*NmatL;
%              ElemKLR = ElemKLR + C1L*NmatR'*(Term33L+Term18L) + C1R*NBmatR'*NmatL;
%              ElemKRL = ElemKRL + C1R*NmatL'*(Term33R+Term18R) + C1L*NBmatL'*NmatR;
%              ElemKRR = ElemKRR - C1R*NmatR'*(Term33R+Term18R) - C1R*NBmatR'*NmatR;

            %relate to jumpu terms &damage terms most complex terms
                 ElemKLL = ElemKLL + c1L*(term5L+term5L'+term7L+term8L);%                    
                 ElemKRR = ElemKRR - c1R*(term5R+term5R'+term7R+term8R);%
                 
                % Damage terms

                ElemKLL = ElemKLL - C1L*NBmatL'*dddt*NBmatL;
                ElemKLR = ElemKLR + C1L*NBmatL'*dddt*NBmatR;
                ElemKRL = ElemKRL + C1R*NBmatR'*dddt*NBmatL;
                ElemKRR = ElemKRR - C1R*NBmatR'*dddt*NBmatR;              
            
            % Store history
            if dmax > 0 && tn >= 0
                if dmax > dc
                    dbond = 2;
                else
                    dbond = 1;
                end
            elseif dmax > 0 && tn < 0
                if dmax > dc
                    dbond = -2;
                else
                    dbond = -1;
                end
            else
                dbond = 0;
            end
            if iter <= itchat
              if sign(bonddbond1) ~= sign(bonddbond2)
                chatter = 1;
              else
                chatter = 0;
              end
            end
            damhr = nh2-1+(ll-1)*7;
            dmaxhr = nh2-1+(ll-1)*7+4;
            dbonhr = nh2-1+(ll-1)*7+5;
            chat1hr = nh2-1+(ll-1)*7+6;
            chat2hr = nh2-1+(ll-1)*7+7;
            hr(damhr+1) = dvec(1);
            hr(damhr+2) = dvec(2);
            hr(damhr+3) = dvec(3);
            hr(dmaxhr) = dmax;
            hr(dbonhr) = dbond;
            hr(chat1hr) = bonddbond2;
            hr(chat2hr) = chatter;
            numdbond = numdbond + abs(dbond);
             

 % end of damage part
         
            
            end %lint       
% ElemKLL
ElemK = [ElemKLL ElemKLR; ElemKRL ElemKRR];
ElemF = [ElemFL; ElemFR];
norm(ElemKLL-ElemKLL',1);
bonddbond2


    case 11  %error estimation 
        
        ElemE = zeros(numEn,1);
        %Set integration number
        lint = IntPoint3(nel);
        ib = 0;
        bf = 0;
        der = 0;

        lam = getlam(mateprop);
        mu = PatchE/(2*(1+Patchv));
        thick = 1;
        iderswit = 2;
        fbx = 0;
        fby = 0;
        fbz = 0;
        hx = 0;
        hy = 0;
        hz = 0;
        el2el = zeros(4,1);
        eprixel = zeros(4,1);
        epriyel = zeros(4,1);
        eprizel = zeros(4,1);
        el2fine = zeros(13,1);
        
        for l = 1:lint

                %Evaluate first derivatives of basis functions at int. point
                if nel == 4 || nel == 10
                  [Wgt,ss] =  int3d_t(ll,lint,ib);
                  [shl,shld,shls,be] = shltt(ss,nel,nen,der,bf);
                  [Qxy, shgs, Jdet] = shgtt(xl+ul,nel,shld,shls,nel,bf,der,be);
                else
                  [Wgt,ss] =  intpntb(ll,lint,ib);
                  [shl,shld,shls,be] = shlb(ss,nel,nen,der,bf);
                  [Qxy, shgs, Jdet] = shgb(xl+ul,nel,shld,shls,nel,bf,der,be);
                end
                
            for mm = 1:nel  
             Nmat(:,3*mm-2:3*mm) = [shl(mm,1)     0          0
                                      0        shl(mm,1)     0
                                      0            0       shl(mm,1) ];  % for body force part
             Bmat(:,3*mm-2:3*mm) = [Qxy(mm,1) 0         0         
                                      0         Qxy(mm,2) 0         
                                      0         0         Qxy(mm,3) 
                                    Qxy(mm,2) Qxy(mm,1) 0         
                                      0         Qxy(mm,3) Qxy(mm,2) 
                                   Qxy(mm,3) 0         Qxy(mm,1) 
                                   Qxy(mm,2) -Qxy(mm,1) 0         
                                      0         Qxy(mm,3) -Qxy(mm,2) 
                                   -Qxy(mm,3) 0         Qxy(mm,1) ];
            end
%                 
            [fi,JxX,F] = kine3d(Qxy,-ul,nel,0); %this is equivalent to ikine2d
            JxX = 1/JxX; %this is equivalent to ikine2d
    %         [F,JxX,fi] = ikine2d(Qxy,ul,nel); %this is equivalent to above two lines
            Jdet = Jdet/JxX;

           W_int = (lam/2*log(JxX)-mu)*log(JxX)+mu/2*(trace(F'*F)-3);
            c1 = Jdet*w*thick;

            xint = xl(1,:)*shl;
            yint = xl(2,:)*shl;
            zint = xl(3,:)*shl;
            dux = ul*Qxy(:,1);
            duy = ul*Qxy(:,2);
            duz = ul*Qxy(:,3);
            u = ul*shl;

            %Compute value of exact fields at int. point
            if iprob == 1
                [ue,duex,duey,duez] = uexact_bar3(xint,yint,zint,PatchE,Patchv,rho);  % exact
%                 solution for bar            
%             elseif iprob == 2
%                 [ue,duex,duey,duez] = uexact_fluid(xint,yint,zint);
            elseif iprob ==3
               [ue,duex,duey,duez] = uexact_cantilever3(xint,yint,zint,PatchE,Patchv);  % exact  
             else
                ue = zeros(4,1);duex=ue;duey=ue;duez=ue;
             end

            %Add standard int. point error to element standard error
            for in = 1:4
                un   = c1 * ( (u(in)-ue(in))^2 );
                upnx   = c1 * ( (dux(in)-duex(in))^2 );
                upny   = c1 * ( (duy(in)-duey(in))^2 );
                upnz   = c1 * ( (duz(in)-duez(in))^2 );
                el2el(in)   = el2el(in)   + un;
                eprixel(in) = eprixel(in) + upnx;
                epriyel(in) = epriyel(in) + upny;
                eprizel(in) = eprizel(in) + upnz;
            end


        end %l

        for in= 1:3
            ElemE(in) = el2el(in);
            ElemE(in+4) = eprixel(in);
            ElemE(in+8) = epriyel(in);
            ElemE(in+12) = eprizel(in);
        end

    case 12  %internal energy 
        
        %Set integration number
        lint = IntPoint3(nel);
        ib = 0;
        bf = 0;
        der = 0;
        Bmat = zeros(9,3*nel);
        lam = getlam(mateprop);
        mu = PatchE/(2*(1+Patchv));
        thick = 1;
        
        for l = 1:lint

                %Evaluate first derivatives of basis functions at int. point
                if nel == 4 || nel == 10
                  [Wgt,ss] =  int3d_t(ll,lint,ib);
                  [shl,shld,shls,be] = shltt(ss,nel,nen,der,bf);
                  [Qxy, shgs, Jdet] = shgtt(xl+ul,nel,shld,shls,nel,bf,der,be);
                else
                  [Wgt,ss] =  intpntb(ll,lint,ib);
                  [shl,shld,shls,be] = shlb(ss,nel,nen,der,bf);
                  [Qxy, shgs, Jdet] = shgb(xl+ul,nel,shld,shls,nel,bf,der,be);
                end
                
            for mm = 1:nel  
             Nmat(:,3*mm-2:3*mm) = [shl(mm,1)     0          0
                                      0        shl(mm,1)     0
                                      0            0       shl(mm,1) ];  % for body force part
             Bmat(:,3*mm-2:3*mm) = [Qxy(mm,1) 0         0         
                                      0         Qxy(mm,2) 0         
                                      0         0         Qxy(mm,3) 
                                    Qxy(mm,2) Qxy(mm,1) 0         
                                      0         Qxy(mm,3) Qxy(mm,2) 
                                   Qxy(mm,3) 0         Qxy(mm,1) 
                                   Qxy(mm,2) -Qxy(mm,1) 0         
                                      0         Qxy(mm,3) -Qxy(mm,2) 
                                   -Qxy(mm,3) 0         Qxy(mm,1) ];
            end
%                 
            [fi,JxX,F] = kine3d(Qxy,-ul,nel,0); %this is equivalent to ikine2d
            JxX = 1/JxX; %this is equivalent to ikine2d
    %         [F,JxX,fi] = ikine2d(Qxy,ul,nel); %this is equivalent to above two lines
            Jdet = Jdet/JxX;

           Energy_int = (lam/2*log(JxX)-mu)*log(JxX)+mu/2*(trace(F'*F)-3);
           
        end


    case 21
 
        ElemK = zeros(nst);

        %Set integration number
        lint = IntPoint3(nel);
        ib = 0;
        bf = 0;
        der = 0;

%         lam = Patchv*PatchE/((1+Patchv)*(1-2*Patchv));
%         mu = PatchE/(2*(1+Patchv));
        thick = 1;
        Dmat = mu*diag([2 2 2 1 1 1]) + lam*[1; 1; 1; 0; 0; 0]*[1 1 1 0 0 0];
        Bmat = zeros(6,3*nel);

        for l = 1:lint                    
                    
            % Evaluate 1-D basis functions at integration points
            if nel == 4 || nel == 10
              [w,lit] =  int3d_t(l,lint,ib);
              [shl,shld,shls,be] = shltt(lit,nel,nel,der,bf);
              [shg, shgs, Jdet] = shgtt(xl,nel,shld,shls,nen,bf,der,be);
%                   shg = [shg'; shl'];
%                   shgs = shgs';
            else
              [w,lit] =  intpntb(l,lint,ib);
              [shl,shld,shls,be] = shlb(lit,nel,nel,der,bf);
              [shg, shgs, Jdet, be] = shgb(xl,nel,shld,shls,nen,bf,der,be);
%                   shg = [shg'; shl'];
%                   shgs = shgs';
            end
            
            if Jdet < 0
                elem
                Jdet %#ok<NOPTS>
            end
            c1 = Jdet*w*thick;
            
            % Form B matrix
            for ie = 1:nel
                
            Bmat(:,(ie-1)*4+1:4*ie) = [shg(ie,1) 0         0         0
                                       0         shg(ie,2) 0         0
                                       0         0         shg(ie,3) 0
                                       shg(ie,2) shg(ie,1) 0         0
                                       0         shg(ie,3) shg(ie,2) 0
                                       shg(ie,3) 0         shg(ie,1) 0
                                       0         0         0         shl(ie)];
                                    
            end
    
            ElemK = ElemK + c1*(Bmat'*Dmat*Bmat);

        end %je
ElemK;

    case -1
        
        ElemF = zeros(nst,1);
        t1 = zeros(3,1);
        t2 = t1;
        t3 = t1;
    
        if nel == 4 || nel == 10
            lint = 13;
        else
            lint = 100; %16; % Use 100 for body force patchtest4_BF &patchtest5_BF  
        end
        der = 0;
        bf = 0;
        ib = 5;
        
        % Integration Loop
        for je = 1:lint

            if nel == 4 || nel == 10
              [Wgt,ss] =  int3d_t(je,lint,ib);
              [shl,shld,shls,be] = shltt(ss,nel,nel,der,bf);
              [QXY,shgs,Jdet,be,sx] = shgtt(xl,nel,shld,shls,nel,bf,der,be);
            else
              [Wgt,ss] =  intpntb(je,lint,ib);
              [shl,shld,shls,be] = shlb(ss,nel,nel,der,bf);
              [QXY,shgs,Jdet,be,sx] = shgb(xl,nel,shld,shls,nel,bf,der,be);
            end
          
            %Evaluate tangent and normal vectors
            t1 = sx(:,2);
            [tm1, tu1] = VecNormalize(t1);
            t2 = sx(:,1);
            [tm2, tu2] = VecNormalize(t2);
            t3 = VecCrossProd(t1,t2);
            [tm3, tu3] = VecNormalize(t3);
            
            if exist('iprob','var')
                if iprob == 6
                    mu = 40;
                    lam = getlam(mateprop);
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

            %Force components are positive in positive coord. direction
            c1 = Wgt*tm3;
            for o=1:nel
                don = shl(o);
                F = don*Traction';  %traction for the one without body force

                ElemF(ndf*o-2) = ElemF(ndf*o-2) + F(1)*c1;

                ElemF(ndf*o-1) = ElemF(ndf*o-1) + F(2)*c1;

                ElemF(ndf*o-0) = ElemF(ndf*o-0) + F(3)*c1;

            end %o
            
        end %je
        ElemF;
 %%
    case -2
        
        ElemF = zeros(nst,1);
        
        lint = 13;

        nil = surfacesi(1);

        for intt = 1:nil %integrate on left domain
                
            trinum = surfacesi(intt+1);
            xit = zeros(ndm,3);
            for j = 1:3
                node = ixt(trinum,j);
                for i = 1:ndm
                    xit(i,j) = Coordinates(node,i);
                end
            end
        
            for l = 1:lint


                %Integration point, weight, jacobian
                [Wgt,litr,lits] =  intpntt(l,lint,0);
                [shl,shld,shls,be] = shlt(litr,lits,3,3,0,0);
                [shg, shgs, Jdet, be, xs] = shgt(xit,3,shld,shls,3,0,0,be);

                %Physical location of int pt
                xint = xit(1,:)*shl;
                yint = xit(2,:)*shl;

                xi = POU_Coord(xint,yint,xl,1,4);
                r = xi(1);
                s = xi(2);
                t = 1;
                ss = [r s t];
                % FIX THIS FOR TETS

                % Evaluate  basis functions at integration points
                if nel == 4 || nel == 10
                  shl = shltt(ss,nel,nel,0,0);
                else
                  shl = shlb(ss,nel,nel,0,0);
                end

                %Evaluate tangent and normal vectors
                t1 = [xs(:,1); 0];
                [tm1, tu1] = VecNormalize(t1);
                t2 = [xs(:,2); 0];
                [tm2, tu2] = VecNormalize(t2);
                t3 = VecCrossProd(t1,t2);
                [tm3, tu3] = VecNormalize(t3);
    %             t = [tu1' tu2' tu3'];
            
                if iprob == 5

                else
                    Traction = traction;
                end

                c1 = Wgt*tm3;

                for o=1:nel

                    don = shl(o);
                    F = don*Traction';

                    ElemF(ndf*o-2) = ElemF(ndf*o-2) + F(1)*c1;

                    ElemF(ndf*o-1) = ElemF(ndf*o-1) + F(2)*c1;

                    ElemF(ndf*o-0) = ElemF(ndf*o-0) + F(3)*c1;

                end %o

            end %ie
            
        end %intt
        ElemF;
%%        


    case 51 % Volume stress/strain homogenization
        
        ElemSS = zeros(17,1);

%         %Set integration number
%         lint = IntPoint3(nel);
%         ib = 0;
%         bf = 0;
%         der = 0;
% 
% %         lam = Patchv*PatchE/((1+Patchv)*(1-2*Patchv));
% %         mu = PatchE/(2*(1+Patchv));
%         thick = 1;
%         fbx = 0;
%         fby = 0;
%         fbz = 0;
%         Dmat = mu*diag([2 2 2 1 1 1]) + lam*[1; 1; 1; 0; 0; 0]*[1 1 1 0 0 0];
%         Nmat = zeros(3,3*nel);
%         Bmat = zeros(6,3*nel);

        if nel == 4
            lint = 16;11;5;16;
        elseif nel == 8
            lint = 8; 16;1000;%1000 for body force
        elseif nel == 10
            lint = 14;
%             lint = 27;
        else
            lint = 27;
        end
        der = 0;
        bf = 0;
        ib = 0;  
        ll = 0;
        for l = 1:lint                    

                ll = ll + 1;
                
               %Evaluate first derivatives of basis functions at int. point
                if nel == 4 || nel == 10
                  [Wgt,ss] =  int3d_t(ll,lint,ib);
                  [shl,shld,shls,be] = shltt(ss,nel,nen,der,bf);
                  [Pxy, Shgs, JdetC] = shgtt(xl,nel,shld,shls,nel,bf,der,be);                 
                  [Qxy, shgs, Jdet] = shgtt(xl+ul,nel,shld,shls,nel,bf,der,be);
                else
                  [Wgt,ss] =  intpntb(ll,lint,ib);
                  [shl,shld,shls,be] = shlb(ss,nel,nen,der,bf);
                  [Pxy, Shgs, JdetC] = shgb(xl,nel,shld,shls,nel,bf,der,be);                  
                  [Qxy, shgs, Jdet] = shgb(xl+ul,nel,shld,shls,nel,bf,der,be);
                end




 % initial stress sigma and material bulk term c_ijkl
%            for mm = 1:nelL 
% %  NmatL(:,3*mm-2:3*mm) = [shlL(mm,1)     0          0
% %                            0        shlL(mm,1)     0
% %                            0            0       shlL(mm,1) ];
% %  BmatL(:,3*mm-2:3*mm) = [QxyL(mm,1) 0         0         
% %                         0         QxyL(mm,2) 0         
% %                         0         0         QxyL(mm,3) 
% %                         QxyL(mm,2) QxyL(mm,1) 0         
% %                         0         QxyL(mm,3) QxyL(mm,2) 
% %                         QxyL(mm,3) 0         QxyL(mm,1) 
% %                         QxyL(mm,2) -QxyL(mm,1) 0         
% %                         0         QxyL(mm,3) -QxyL(mm,2) 
% %                         -QxyL(mm,3) 0         QxyL(mm,1) ];                
% %  BmatL1(:,3*mm-2:3*mm) = [QxyL(mm,1) 0         0         
% %                         0         QxyL(mm,2) 0         
% %                         0         0         QxyL(mm,3) 
% %                         QxyL(mm,2) QxyL(mm,1) 0         
% %                         0         QxyL(mm,3) QxyL(mm,2) 
% %                         QxyL(mm,3) 0         QxyL(mm,1) ];
% % BmatL2(:,3*mm-2:3*mm)=[QxyL(mm,1)  0             0         
% %                        QxyL(mm,2)  0             0         
% %                        QxyL(mm,3)  0             0          
% %                         0     QxyL(mm,1)         0         
% %                         0     QxyL(mm,2)         0         
% %                         0     QxyL(mm,3)         0
% %                         0         0         QxyL(mm,1)         
% %                         0         0         QxyL(mm,2)         
% %                         0         0         QxyL(mm,3)];
%  BmatL4(:,3*mm-2:3*mm) = [PxyL(mm,1) 0         0         
%                         0         PxyL(mm,2) 0         
%                         0         0         PxyL(mm,3) 
%                         PxyL(mm,2) PxyL(mm,1) 0         
%                         0         PxyL(mm,3) PxyL(mm,2) 
%                         PxyL(mm,3) 0         PxyL(mm,1) 
%                         PxyL(mm,2) -PxyL(mm,1) 0         
%                         0         PxyL(mm,3) -PxyL(mm,2) 
%                         -PxyL(mm,3) 0         PxyL(mm,1) ];   
% BmatL5(:,3*mm-2:3*mm)=[PxyL(mm,1)  0             0         
%                        PxyL(mm,2)  0             0         
%                        PxyL(mm,3)  0             0          
%                         0     PxyL(mm,1)         0         
%                         0     PxyL(mm,2)         0         
%                         0     PxyL(mm,3)         0
%                         0         0         PxyL(mm,1)         
%                         0         0         PxyL(mm,2)         
%                         0         0         PxyL(mm,3)];                    
%            end
%             
%            for mm = 1:nelR    
% %  NmatR(:,3*mm-2:3*mm) = [shlR(mm,1)     0          0
% %                            0        shlR(mm,1)     0
% %                            0            0       shlR(mm,1) ];
% %  BmatR(:,3*mm-2:3*mm) = [QxyR(mm,1) 0         0         
% %                         0         QxyR(mm,2) 0         
% %                         0         0         QxyR(mm,3) 
% %                         QxyR(mm,2) QxyR(mm,1) 0         
% %                         0         QxyR(mm,3) QxyR(mm,2) 
% %                         QxyR(mm,3) 0         QxyR(mm,1) 
% %                         QxyR(mm,2) -QxyR(mm,1) 0         
% %                         0         QxyR(mm,3) -QxyR(mm,2) 
% %                         -QxyR(mm,3) 0         QxyR(mm,1) ];
% %  BmatR1(:,3*mm-2:3*mm) = [QxyR(mm,1) 0         0         
% %                         0         QxyR(mm,2) 0         
% %                         0         0         QxyR(mm,3) 
% %                         QxyR(mm,2) QxyR(mm,1) 0         
% %                         0         QxyR(mm,3) QxyR(mm,2) 
% %                         QxyR(mm,3) 0         QxyR(mm,1) ];
% %  BmatR2(:,3*mm-2:3*mm)=[QxyR(mm,1)  0             0         
% %                        QxyR(mm,2)  0             0         
% %                        QxyR(mm,3)  0             0          
% %                         0     QxyR(mm,1)         0         
% %                         0     QxyR(mm,2)         0         
% %                         0     QxyR(mm,3)         0
% %                         0         0         QxyR(mm,1)         
% %                         0         0         QxyR(mm,2)         
% %                         0         0         QxyR(mm,3)];  
% BmatR4(:,3*mm-2:3*mm) = [PxyR(mm,1) 0         0         
%                         0         PxyR(mm,2) 0         
%                         0         0         PxyR(mm,3) 
%                         PxyR(mm,2) PxyR(mm,1) 0         
%                         0         PxyR(mm,3) PxyR(mm,2) 
%                         PxyR(mm,3) 0         PxyR(mm,1) 
%                         PxyR(mm,2) -PxyR(mm,1) 0         
%                         0         PxyR(mm,3) -PxyR(mm,2) 
%                         -PxyR(mm,3) 0         PxyR(mm,1) ];                   
%  BmatR5(:,3*mm-2:3*mm)=[PxyR(mm,1)  0             0         
%                        PxyR(mm,2)  0             0         
%                        PxyR(mm,3)  0             0          
%                         0     PxyR(mm,1)         0         
%                         0     PxyR(mm,2)         0         
%                         0     PxyR(mm,3)         0
%                         0         0         PxyR(mm,1)         
%                         0         0         PxyR(mm,2)         
%                         0         0         PxyR(mm,3)];  
%            end  
           
            [fi,JxX,F] = kine3d(Qxy,-ul,nel,0); %this is equivalent to ikine2d
            JxX = 1/JxX; %this is equivalent to ikine2d
    %         [F,JxX,fi] = ikine2d(Qxy,ul,nel); %this is equivalent to above two lines
            Jdet = JdetC/JxX;
                
            c1 = Wgt*Jdet;
          
            [sigma2, cmat] = SigmaCmat3(F,JxX,mateprop,lam);
            P = [fi(1,1)^2 fi(1,2)^2 fi(1,3)^2 2*fi(1,1)*fi(1,2) 2*fi(1,2)*fi(1,3) 2*fi(1,3)*fi(1,1)  zeros(1,3)
                 fi(2,1)^2 fi(2,2)^2 fi(2,3)^2 2*fi(2,1)*fi(2,2) 2*fi(2,2)*fi(2,3) 2*fi(2,3)*fi(2,1)  zeros(1,3)
                 fi(3,1)^2 fi(3,2)^2 fi(3,3)^2 2*fi(3,1)*fi(3,2) 2*fi(3,2)*fi(3,3) 2*fi(3,3)*fi(3,1)  zeros(1,3)                 
                 fi(1,1)*fi(2,1) fi(1,2)*fi(2,2) fi(1,3)*fi(2,3) fi(1,1)*fi(2,2)+fi(1,2)*fi(2,1) fi(1,2)*fi(2,3)+fi(1,3)*fi(2,2) fi(1,1)*fi(2,3)+fi(1,3)*fi(2,1) zeros(1,3)
                 fi(2,1)*fi(3,1) fi(2,2)*fi(3,2) fi(2,3)*fi(3,3) fi(2,1)*fi(3,2)+fi(2,2)*fi(3,1) fi(2,2)*fi(3,3)+fi(2,3)*fi(3,2) fi(2,1)*fi(3,3)+fi(2,3)*fi(3,1) zeros(1,3)                 
                 fi(3,1)*fi(1,1) fi(3,2)*fi(1,2) fi(3,3)*fi(1,3) fi(3,1)*fi(1,2)+fi(3,2)*fi(1,1) fi(3,2)*fi(1,3)+fi(3,3)*fi(1,2) fi(3,1)*fi(1,3)+fi(3,3)*fi(1,1) zeros(1,3)                  
                 zeros(3,9)];
            SPK = P*sigma2;   
            Smat1=[SPK(1), SPK(4), SPK(6)
                    SPK(4), SPK(2), SPK(5)
                    SPK(6), SPK(5), SPK(3)];  
                
            FPK = F*Smat1;  
%             C = [F(1,1)^2 F(1,2)^2 F(1,3)^2 2*F(1,1)*F(1,2) 2*F(1,2)*F(1,3) 2*F(1,3)*F(1,1)  
%                  F(2,1)^2 F(2,2)^2 F(2,3)^2 2*F(2,1)*F(2,2) 2*F(2,2)*F(2,3) 2*F(2,3)*F(2,1)  
%                  F(3,1)^2 F(3,2)^2 F(3,3)^2 2*F(3,1)*F(3,2) 2*F(3,2)*F(3,3) 2*F(3,3)*F(3,1)             
%                  F(1,1)*F(2,1) F(1,2)*F(2,2) F(1,3)*F(2,3) F(1,1)*F(2,2)+F(1,2)*F(2,1) F(1,2)*F(2,3)+F(1,3)*F(2,2) F(1,1)*F(2,3)+F(1,3)*F(2,1) 
%                  F(2,1)*F(3,1) F(2,2)*F(3,2) F(2,3)*F(3,3) F(2,1)*F(3,2)+F(2,2)*F(3,1) F(2,2)*F(3,3)+F(2,3)*F(3,2) F(2,1)*F(3,3)+F(2,3)*F(3,1)              
%                  F(3,1)*F(1,1) F(3,2)*F(1,2) F(3,3)*F(1,3) F(3,1)*F(1,2)+F(3,2)*F(1,1) F(3,2)*F(1,3)+F(3,3)*F(1,2) F(3,1)*F(1,3)+F(3,3)*F(1,1) ]';
            E_strain = 1/2*(F'*F-eye(3,3));
         
%             epsil = Bmat*reshape(ul,ndf*nel,1);
%             stres = Dmat*epsil;
            volum = c1;
            
            ElemSS(1:6) = ElemSS(1:6) + c1*[E_strain(1,1);E_strain(2,2);E_strain(3,3);E_strain(1,2);E_strain(2,3);E_strain(3,1)];
            ElemSS(7:12) = ElemSS(7:12) + c1*[FPK(1,1);FPK(2,2);FPK(3,3);FPK(1,2);FPK(2,3);FPK(3,1)];
            ElemSS(13) = ElemSS(13) + volum;

        end %je
            lamda_step = (Coordinates(1,1)*2+DispList(1,1,step+1))/(Coordinates(1,1)*2);
            ElemSS(14) = log(lamda_step);
            ElemSS(15) = ElemSS(7)./((Coordinates(1,1)*2+DispList(1,1,step+1))*Coordinates(1,1)*2*20*lamda_step);
            ElemSS(16) = ElemSS(8)./((Coordinates(1,1)*2+DispList(1,1,step+1))*Coordinates(1,1)*2*20);
            ElemSS(17) = ElemSS(9)./((Coordinates(1,1)*2+DispList(1,1,step+1))*Coordinates(1,1)*2*20);            
    case 52 % Surface strain homogenization
        
        ElemSS = zeros(13,1);
        t1 = zeros(3,1);
        t2 = t1;
        t3 = t1;
        
        ElemEL = matepropL(1);
        ElemvL = matepropL(2);
        ElemER = matepropR(1);
        ElemvR = matepropR(2);
        lamdaR = ElemvR*ElemER/((1+ElemvR)*(1-2*ElemvR));
        muR = ElemER/(2*(1+ElemvR));
        lamdaL = ElemvL*ElemEL/((1+ElemvL)*(1-2*ElemvL));
        muL = ElemEL/(2*(1+ElemvL));
        DmatR = muR*diag([2 2 2 1 1 1]) + lamdaR*[1; 1; 1; 0; 0; 0]*[1 1 1 0 0 0];
        DmatL = muL*diag([2 2 2 1 1 1]) + lamdaL*[1; 1; 1; 0; 0; 0]*[1 1 1 0 0 0];
        
        NmatL = zeros(3,nstL);
        BmatL = zeros(6,nstL);
        bnAdN1 = zeros(6,nstL);
        N1 = zeros(3,nstL);
        NmatR = zeros(3,nstR);
        BmatR = zeros(6,nstR);
        bnAdN2 = zeros(6,nstR);
        N2 = zeros(3,nstR);
        
        lint = 4;
        ll = 0; % Counter for history variables

%         nil = surfacesi(1);
% 
%         for intt = 1:nil %integrate on left domain
%                 
%             trinum = surfacesi(intt+1);
%             xit = zeros(ndm,3);
%             for j = 1:3
%                 node = ixt(trinum,j);
%                 for i = 1:ndm
%                     xit(i,j) = Coordinates(node,i);
%                 end
%             end
%         
            for l = 1:lint


               ll = ll + 1;

                % Evaluate  basis functions at integration points
                if nelL == 4 || nelL == 10
                  [shlL,shldL,shls,be] = shltt([rL sL tL],nelL,nel2L,0,0);
                  [PxyL,shgsL,JdetL,bubbleL,XsL] = shgtt(xlL(:,1:nelL),nelL,shldL,shls,nen,0,0,be);                 
                  [QxyL,shgsL,JdetL,bubbleL,xsL] = shgtt(xlL(:,1:nelL)+ulL(:,1:nelL),nelL,shldL,shls,nen,0,0,be);
                else
                  [Wgt,ss] = intpntb(l,lint,ib);
                  [shlL,shldL,shls,be] = shlb(ss,nelL,nel2L,0,0);
                  [PxyL,shgsL,JdetL,bubbleL,XsL] = shgb(xlL(:,1:nelL),nelL,shldL,shls,nen,0,0,be); 
                  [QxyL,shgsL,JdetL,bubbleL,xsL] = shgb(xlL(:,1:nelL)+ulL(:,1:nelL),nelL,shldL,shls,nen,0,0,be);
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
                  [shlR,shldR,shls,be] = shltt([rR sR tR],nelR,nel2R,0,0);
                  [PxyR,shgsR,JdetR,bubbleR,XsR] = shgtt(xlR(:,1:nelR),nelR,shldR,shls,nen,0,0,be);                  
                  [QxyR,shgsR,JdetR,bubbleR,xsR] = shgtt(xlR(:,1:nelR)+ulR(:,1:nelR),nelR,shldR,shls,nen,0,0,be);
                  
                else
                  [shlR,shldR,shls,be] = shlb([rR sR tR],nelR,nel2R,0,0);
                  [PxyR,shgsR,JdetR,bubbleR,XsR] = shgb(xlR(:,1:nelR),nelR,shldR,shls,nen,0,0,be);     
                  [QxyR,shgsR,JdetR,bubbleR,xsR] = shgb(xlR(:,1:nelR)+ulR(:,1:nelR),nelR,shldR,shls,nen,0,0,be);
                end


 % initial stress sigma and material bulk term c_ijkl
           for mm = 1:nelL 
 NmatL(:,3*mm-2:3*mm) = [shlL(mm,1)     0          0
                           0        shlL(mm,1)     0
                           0            0       shlL(mm,1) ];
%  BmatL(:,3*mm-2:3*mm) = [QxyL(mm,1) 0         0         
%                         0         QxyL(mm,2) 0         
%                         0         0         QxyL(mm,3) 
%                         QxyL(mm,2) QxyL(mm,1) 0         
%                         0         QxyL(mm,3) QxyL(mm,2) 
%                         QxyL(mm,3) 0         QxyL(mm,1) 
%                         QxyL(mm,2) -QxyL(mm,1) 0         
%                         0         QxyL(mm,3) -QxyL(mm,2) 
%                         -QxyL(mm,3) 0         QxyL(mm,1) ];                
%  BmatL1(:,3*mm-2:3*mm) = [QxyL(mm,1) 0         0         
%                         0         QxyL(mm,2) 0         
%                         0         0         QxyL(mm,3) 
%                         QxyL(mm,2) QxyL(mm,1) 0         
%                         0         QxyL(mm,3) QxyL(mm,2) 
%                         QxyL(mm,3) 0         QxyL(mm,1) ];
% BmatL2(:,3*mm-2:3*mm)=[QxyL(mm,1)  0             0         
%                        QxyL(mm,2)  0             0         
%                        QxyL(mm,3)  0             0          
%                         0     QxyL(mm,1)         0         
%                         0     QxyL(mm,2)         0         
%                         0     QxyL(mm,3)         0
%                         0         0         QxyL(mm,1)         
%                         0         0         QxyL(mm,2)         
%                         0         0         QxyL(mm,3)];
 BmatL4(:,3*mm-2:3*mm) = [PxyL(mm,1) 0         0         
                        0         PxyL(mm,2) 0         
                        0         0         PxyL(mm,3) 
                        PxyL(mm,2) PxyL(mm,1) 0         
                        0         PxyL(mm,3) PxyL(mm,2) 
                        PxyL(mm,3) 0         PxyL(mm,1) 
                        PxyL(mm,2) -PxyL(mm,1) 0         
                        0         PxyL(mm,3) -PxyL(mm,2) 
                        -PxyL(mm,3) 0         PxyL(mm,1) ];   
BmatL5(:,3*mm-2:3*mm)=[PxyL(mm,1)  0             0         
                       PxyL(mm,2)  0             0         
                       PxyL(mm,3)  0             0          
                        0     PxyL(mm,1)         0         
                        0     PxyL(mm,2)         0         
                        0     PxyL(mm,3)         0
                        0         0         PxyL(mm,1)         
                        0         0         PxyL(mm,2)         
                        0         0         PxyL(mm,3)];                    
           end
            
           for mm = 1:nelR    
 NmatR(:,3*mm-2:3*mm) = [shlR(mm,1)     0          0
                           0        shlR(mm,1)     0
                           0            0       shlR(mm,1) ];
%  BmatR(:,3*mm-2:3*mm) = [QxyR(mm,1) 0         0         
%                         0         QxyR(mm,2) 0         
%                         0         0         QxyR(mm,3) 
%                         QxyR(mm,2) QxyR(mm,1) 0         
%                         0         QxyR(mm,3) QxyR(mm,2) 
%                         QxyR(mm,3) 0         QxyR(mm,1) 
%                         QxyR(mm,2) -QxyR(mm,1) 0         
%                         0         QxyR(mm,3) -QxyR(mm,2) 
%                         -QxyR(mm,3) 0         QxyR(mm,1) ];
%  BmatR1(:,3*mm-2:3*mm) = [QxyR(mm,1) 0         0         
%                         0         QxyR(mm,2) 0         
%                         0         0         QxyR(mm,3) 
%                         QxyR(mm,2) QxyR(mm,1) 0         
%                         0         QxyR(mm,3) QxyR(mm,2) 
%                         QxyR(mm,3) 0         QxyR(mm,1) ];
%  BmatR2(:,3*mm-2:3*mm)=[QxyR(mm,1)  0             0         
%                        QxyR(mm,2)  0             0         
%                        QxyR(mm,3)  0             0          
%                         0     QxyR(mm,1)         0         
%                         0     QxyR(mm,2)         0         
%                         0     QxyR(mm,3)         0
%                         0         0         QxyR(mm,1)         
%                         0         0         QxyR(mm,2)         
%                         0         0         QxyR(mm,3)];  
BmatR4(:,3*mm-2:3*mm) = [PxyR(mm,1) 0         0         
                        0         PxyR(mm,2) 0         
                        0         0         PxyR(mm,3) 
                        PxyR(mm,2) PxyR(mm,1) 0         
                        0         PxyR(mm,3) PxyR(mm,2) 
                        PxyR(mm,3) 0         PxyR(mm,1) 
                        PxyR(mm,2) -PxyR(mm,1) 0         
                        0         PxyR(mm,3) -PxyR(mm,2) 
                        -PxyR(mm,3) 0         PxyR(mm,1) ];                   
 BmatR5(:,3*mm-2:3*mm)=[PxyR(mm,1)  0             0         
                       PxyR(mm,2)  0             0         
                       PxyR(mm,3)  0             0          
                        0     PxyR(mm,1)         0         
                        0     PxyR(mm,2)         0         
                        0     PxyR(mm,3)         0
                        0         0         PxyR(mm,1)         
                        0         0         PxyR(mm,2)         
                        0         0         PxyR(mm,3)];  
           end  
    %plug in the -ul to get the inverse of deformation gradient fi as 'F' in this function             
            [fi,JxXL,FL] = kine3d(QxyL,-ulL,nelL,0); %this is equivalent to ikine2d 
            JxXL = 1/JxXL; %this is equivalent to ikine2d 
    %         [F,JxX,fi] = ikine2d(Qxy,ul,nel); %this is equivalent to above two lines
            JdetxL = JdetL/JxXL;               
%             C1L = Wgt*JdetL;
            [fiR,JxXR,FR] = kine3d(QxyR,-ulR,nelR,0); %this is equivalent to ikine2d
            JxXR = 1/JxXR; %this is equivalent to ikine2d
    %         [F,JxX,fi] = ikine2d(Qxy,ul,nel); %this is equivalent to above two lines
            JdetxR = JdetR/JxXL;                          
%             C1R = Wgt*JdetR;         this is for volume integral  
          
            [sigma2L, cmatL] = SigmaCmat3i(FL,JxXL,mateprop,lam);
            P = [fiL(1,1)^2 fiL(1,2)^2 fiL(1,3)^2 2*fiL(1,1)*fiL(1,2) 2*fiL(1,2)*fiL(1,3) 2*fiL(1,3)*fiL(1,1)  zeros(1,3)
                 fiL(2,1)^2 fiL(2,2)^2 fiL(2,3)^2 2*fiL(2,1)*fiL(2,2) 2*fiL(2,2)*fiL(2,3) 2*fiL(2,3)*fiL(2,1)  zeros(1,3)
                 fiL(3,1)^2 fiL(3,2)^2 fiL(3,3)^2 2*fiL(3,1)*fiL(3,2) 2*fiL(3,2)*fiL(3,3) 2*fiL(3,3)*fiL(3,1)  zeros(1,3)                 
                 fiL(1,1)*fiL(2,1) fiL(1,2)*fiL(2,2) fiL(1,3)*fiL(2,3) fiL(1,1)*fiL(2,2)+fiL(1,2)*fiL(2,1) fiL(1,2)*fiL(2,3)+fiL(1,3)*fiL(2,2) fiL(1,1)*fiL(2,3)+fiL(1,3)*fiL(2,1) zeros(1,3)
                 fiL(2,1)*fiL(3,1) fiL(2,2)*fiL(3,2) fiL(2,3)*fiL(3,3) fiL(2,1)*fiL(3,2)+fiL(2,2)*fiL(3,1) fiL(2,2)*fiL(3,3)+fiL(2,3)*fiL(3,2) fiL(2,1)*fiL(3,3)+fiL(2,3)*fiL(3,1) zeros(1,3)                 
                 fiL(3,1)*fiL(1,1) fiL(3,2)*fiL(1,2) fiL(3,3)*fiL(1,3) fiL(3,1)*fiL(1,2)+fiL(3,2)*fiL(1,1) fiL(3,2)*fiL(1,3)+fiL(3,3)*fiL(1,2) fiL(3,1)*fiL(1,3)+fiL(3,3)*fiL(1,1) zeros(1,3)                  
                 zeros(3,9)];
            SPKL = P*sigma2L*JxXL;
            CmatL = P*cmatL*JxXL*P';            
            SmatL1=[SPKL(1), SPKL(4), SPKL(6)
                    SPKL(4), SPKL(2), SPKL(5)
                    SPKL(6), SPKL(5), SPKL(3)];        
                        
          
  
            
            [sigma2R, cmatR] = SigmaCmat3i(FR,JxXR,mateprop,lam);

            P = [fiR(1,1)^2 fiR(1,2)^2 fiR(1,3)^2 2*fiR(1,1)*fiR(1,2) 2*fiR(1,2)*fiR(1,3) 2*fiR(1,3)*fiR(1,1)  zeros(1,3)
                 fiR(2,1)^2 fiR(2,2)^2 fiR(2,3)^2 2*fiR(2,1)*fiR(2,2) 2*fiR(2,2)*fiR(2,3) 2*fiR(2,3)*fiR(2,1)  zeros(1,3)
                 fiR(3,1)^2 fiR(3,2)^2 fiR(3,3)^2 2*fiR(3,1)*fiR(3,2) 2*fiR(3,2)*fiR(3,3) 2*fiR(3,3)*fiR(3,1)  zeros(1,3)                 
                 fiR(1,1)*fiR(2,1) fiR(1,2)*fiR(2,2) fiR(1,3)*fiR(2,3) fiR(1,1)*fiR(2,2)+fiR(1,2)*fiR(2,1) fiR(1,2)*fiR(2,3)+fiR(1,3)*fiR(2,2) fiR(1,1)*fiR(2,3)+fiR(1,3)*fiR(2,1) zeros(1,3)
                 fiR(2,1)*fiR(3,1) fiR(2,2)*fiR(3,2) fiR(2,3)*fiR(3,3) fiR(2,1)*fiR(3,2)+fiR(2,2)*fiR(3,1) fiR(2,2)*fiR(3,3)+fiR(2,3)*fiR(3,2) fiR(2,1)*fiR(3,3)+fiR(2,3)*fiR(3,1) zeros(1,3)                 
                 fiR(3,1)*fiR(1,1) fiR(3,2)*fiR(1,2) fiR(3,3)*fiR(1,3) fiR(3,1)*fiR(1,2)+fiR(3,2)*fiR(1,1) fiR(3,2)*fiR(1,3)+fiR(3,3)*fiR(1,2) fiR(3,1)*fiR(1,3)+fiR(3,3)*fiR(1,1) zeros(1,3)                  
                 zeros(3,9)];
            SPKR = P*sigma2R*JxXR;
            CmatR = P*cmatR*JxXR*P';
            SmatR1=[SPKR(1), SPKR(4), SPKR(6)
                    SPKR(4), SPKR(2), SPKR(5)
                    SPKR(6), SPKR(5), SPKR(3)];   
                
                ElemSS(1:6) = ElemSS(1:6) + c1*epsil;
            
            end %lint
%         
%         end %intt
%%        
    case 25 %Stress Projection2

        ElemS = zeros(nel,npstr+1);
        ElemS2 = zeros(nel,npstr);

         lam = Patchv*PatchE/((1+Patchv)*(1-2*Patchv));
         mu = PatchE/(2*(1+Patchv));
        thick = 1;
        fbx = 0;
        fby = 0;
        fbz = 0;
        Dmat = mu*diag([2 2 2 1 1 1]) + lam*[1; 1; 1; 0; 0; 0]*[1 1 1 0 0 0];
        Nmat = zeros(3,3*nel);
        Bmat = zeros(6,3*nel);
        I1 = [1; 1; 1; 0; 0; 0];
        
        % Load Guass Integration Points

        if nel == 4
            lint = 1;
            nint = 1;
        elseif nel == 8
%             lint = 4;
            lint = 8;
            nint = 1;
        elseif nel == 10
            lint = 11;
            nint = 4;
        else
            lint = 27;
            nint = 8;
        end
        
        der = 0;
        bf = 0;
        ib = 0;

        %Stress Loop
        for ll = 1:nint

            %Evaluate first derivatives of basis functions at int. point
            if nel == 4 || nel == 10
              [Wgt,ss] =  int3d_t(ll,nint,ib);
              [shl,shld,shls,be] = shltt(ss,nel,nel,der,bf);
              [shg,shgs,Jdet,be] = shgtt(xl,nel,shld,shls,nel,bf,der,be);
            else
              [Wgt,ss] =  intpntb(ll,nint,ib);
              [shl,shld,shls,be] = shlb(ss,nel,nel,der,bf);
              [shg,shgs,Jdet,be] = shgb(xl,nel,shld,shls,nel,bf,der,be);
            end
            
            % Form B matrix
            for ie = 1:nel
              
              Bmat(Bcol1,(ie-1)*3+1) = shg(ie,col1);
              Bmat(Bcol2,(ie-1)*3+2) = shg(ie,col2);
              Bmat(Bcol3,3*ie      ) = shg(ie,col3);
                 
            end
            
            epsil = Bmat*reshape(ul,ndf*nel,1);
            sigma = Dmat*epsil;
            
            for stres = 1:npstr
            
            if stres <= 6 % stress components
                sigmas = sigma(stres);
            elseif stres >= 8
                if stres <= 10 % principal stresses
                sigma2 = [sigma(1) sigma(4) sigma(6); sigma(4) sigma(2) sigma(5); sigma(6) sigma(5) sigma(3)];
                psig = eig(sigma2);
                sigmas = psig(stres-7);
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
        if nel == 4
            plist = [1 0 0 0
                     0 1 0 0
                     0 0 0 1];
        elseif nel == 8
            plist = [-1 1 1 -1 -1 1 1 -1
                     -1 -1 1 1 -1 -1 1 1
                     -1 -1 -1 -1 1 1 1 1];
        elseif nel == 10
            plist = [ 1.927050983124842 -0.309016994374947 -0.309016994374947 -0.309016994374947  0.809016994374947 -0.309016994374947  0.809016994374947  0.809016994374947 -0.309016994374947 -0.309016994374947
                     -0.309016994374947  1.927050983124842 -0.309016994374947 -0.309016994374947  0.809016994374947  0.809016994374947 -0.309016994374947 -0.309016994374947  0.809016994374947 -0.309016994374947
                     -0.309016994374947 -0.309016994374947 -0.309016994374947  1.927050983124842 -0.309016994374947 -0.309016994374947 -0.309016994374947  0.809016994374947  0.809016994374947  0.809016994374947];
        else
            sqr3 = sqrt(3);
            plist = [-sqr3 sqr3 sqr3 -sqr3 -sqr3 sqr3 sqr3 -sqr3 0 sqr3 0 -sqr3 0 sqr3 0 -sqr3 -sqr3 sqr3 sqr3 -sqr3 0 0 -sqr3 sqr3 0 0 0
                     -sqr3 -sqr3 sqr3 sqr3 -sqr3 -sqr3 sqr3 sqr3 -sqr3 0 sqr3 0 -sqr3 0 sqr3 0 -sqr3 -sqr3 sqr3 sqr3 0 0 0 0 -sqr3 sqr3 0
                     -sqr3 -sqr3 -sqr3 -sqr3 sqr3 sqr3 sqr3 sqr3 -sqr3 -sqr3 -sqr3 -sqr3 sqr3 sqr3 sqr3 sqr3 0 0 0 0 -sqr3 sqr3 0 0 0 0 0];
        end
        
        for ll = 1:nelS
            
            r = plist(1,ll);
            s = plist(2,ll);
            t = plist(3,ll);
            shpS = sshp3d(r,s,t,nint);
            
            for stres = 1:npstr
                
                sigmas = ElemS2(1:nint,stres)'*shpS;
                ElemS(ll,stres) = sigmas;
                
            end
            
        end
        
        %Integration Loop
        Vol = 0;
        for ll = 1:lint

            %Evaluate first derivatives of basis functions at int. point
            if nel == 4 || nel == 10
              [Wgt,ss] =  int3d_t(ll,lint,ib);
              [shl,shld,shls,be] = shltt(ss,nel,nel,der,bf);
              [Qxy,shgs,Jdet,be] = shgtt(xl,nel,shld,shls,nel,bf,der,be);
            else
              [Wgt,ss] =  intpntb(ll,lint,ib);
              [shl,shld,shls,be] = shlb(ss,nel,nel,der,bf);
              [Qxy,shgs,Jdet,be] = shgb(xl,nel,shld,shls,nel,bf,der,be);
            end

            w = Wgt*Jdet*thick;
            
            Vol = Vol + w;

        end %je

        for i = 1:nel
        ElemS(i,npstr+1) = 1; % use this for simple average Vol; % use this for weighted average 
        end
        
    case 60
        
        numhr = 4;
        ElemI = zeros(15,numhr);
        t1 = zeros(3,1);
        t2 = t1;
        t3 = t1;
        
        dtol = 1e-11;
        
%         beta = 1;
%         beta2 = beta^2;
%         sigmax = 100;
%         dc = 0.2;
%         beta = 0.707;
%         beta2 = beta^2;
%         sigmax = 0.01e-3;
%         dc = 20;
%         
%         Hc = sigmax/dc;
%         rp = 100*Hc;
        ElemEL = matepropL(1);
        ElemvL = matepropL(2);
        ElemER = matepropR(1);
        ElemvR = matepropR(2);
        lamdaR = ElemvR*ElemER/((1+ElemvR)*(1-2*ElemvR));
        muR = ElemER/(2*(1+ElemvR));
        lamdaL = ElemvL*ElemEL/((1+ElemvL)*(1-2*ElemvL));
        muL = ElemEL/(2*(1+ElemvL));
        DmatR = muR*diag([2 2 2 1 1 1]) + lamdaR*[1; 1; 1; 0; 0; 0]*[1 1 1 0 0 0];
        DmatL = muL*diag([2 2 2 1 1 1]) + lamdaL*[1; 1; 1; 0; 0; 0]*[1 1 1 0 0 0];
        
        NmatL = zeros(3,nstL);
        BmatL = zeros(6,nstL);
        bnAdN1 = zeros(6,nstL);
        N1 = zeros(3,nstL);
        NmatR = zeros(3,nstR);
        BmatR = zeros(6,nstR);
        bnAdN2 = zeros(6,nstR);
        N2 = zeros(3,nstR);
        
        lint = 4;3;
        ll = 0; % Counter for history variables
        
            for l = 1:lint

                ll = ll + 1;

                % Evaluate  basis functions at integration points
                if nelL == 4 || nelL == 10
                  [shlL,shldL,shls,be] = shltt([rL sL tL],nelL,nel2L,0,0);
                  PxyL = shgtt(xlL(:,1:nelL),nelL,shldL,shls,nen,0,0,be);
                else
                  [Wgt,ss] = intpntb(l,lint,5);
                  [shlL,shldL,shls,be] = shlb(ss,nelL,nel2L,0,0);
                  [PxyL,shgs,Jdet,bubble,xs] = shgb(xlL(:,1:nelL),nelL,shldL,shls,nen,0,0,be);
                end
                QxyL = PxyL;

                %Physical location of int pt
                xint = xlL(1,:)*shlL;
                yint = xlL(2,:)*shlL;
                zint = xlL(3,:)*shlL;

                xi = POU_Coord3(xint,yint,zint,xlR,1,nelR);
                rR = xi(1);
                sR = xi(2);
                tR = xi(3);

                % Evaluate  basis functions at integration points
                if nelR == 4 || nelR == 10
                  [shlR,shldR,shls,be] = shltt([rR sR tR],nelR,nel2R,0,0);
                  PxyR = shgtt(xlR(:,1:nelR),nelR,shldR,shls,nen,0,0,be);
                else
                  [shlR,shldR,shls,be] = shlb([rR sR tR],nelR,nel2R,0,0);
                  PxyR = shgb(xlR(:,1:nelR),nelR,shldR,shls,nen,0,0,be);
                end
                QxyR = PxyR;

                %Evaluate tangent and normal vectors
                t1 = xs(:,1);
                [tm1, tu1] = VecNormalize(t1);
                t2 = xs(:,2);
                [tm2, tu2] = VecNormalize(t2);
                t3 = VecCrossProd(t1,t2);
                [tm3, tu3] = VecNormalize(t3);
                nLx = -tu3(1);
                nLy = -tu3(2);
                nLz = -tu3(3);
                nRx = tu3(1);
                nRy = tu3(2);
                nRz = tu3(3);
                tLx = tu1(1);
                tLy = tu1(2);
                tLz = tu1(3);
                nvect = [nLx 0 0 nLy 0 nLz
                         0 nLy 0 nLx nLz 0
                         0 0 nLz 0 nLy nLx]; %- ?
                nvec = [nLx; nLy; nLz];

                c1 = Wgt*tm3;
                
            for i = 1:nelL
%                 NmatL(:,3*i-2:3*i) = shlL(i)*eye(3);
%                 BmatL(:,3*i-2:3*i) = [QxyL(i,1) 0 0 
%                                       0 QxyL(i,2) 0 
%                                       0 0 QxyL(i,3) 
%                                       QxyL(i,2) QxyL(i,1) 0 
%                                       0 QxyL(i,3) QxyL(i,2) 
%                                       QxyL(i,3) 0 QxyL(i,1) ];
                NmatL(1,(i-1)*3+1) = shlL(i);
                NmatL(2,(i-1)*3+2) = shlL(i);
                NmatL(3,3*i      ) = shlL(i);
                BmatL(Bcol1,(i-1)*3+1) = QxyL(i,col1);
                BmatL(Bcol2,(i-1)*3+2) = QxyL(i,col2);
                BmatL(Bcol3,3*i      ) = QxyL(i,col3);
            end
            
            for i = 1:nelR
%                 NmatR(:,3*i-2:3*i) = shlR(i)*eye(3);
%                 BmatR(:,3*i-2:3*i) = [QxyR(i,1) 0 0 
%                                       0 QxyR(i,2) 0 
%                                       0 0 QxyR(i,3) 
%                                       QxyR(i,2) QxyR(i,1) 0 
%                                       0 QxyR(i,3) QxyR(i,2) 
%                                       QxyR(i,3) 0 QxyR(i,1)];
                NmatR(1,(i-1)*3+1) = shlR(i);
                NmatR(2,(i-1)*3+2) = shlR(i);
                NmatR(3,3*i      ) = shlR(i);
                BmatR(Bcol1,(i-1)*3+1) = QxyR(i,col1);
                BmatR(Bcol2,(i-1)*3+2) = QxyR(i,col2);
                BmatR(Bcol3,3*i      ) = QxyR(i,col3);
            end
            
            % Load history
            damhr = nh1-1+(ll-1)*7;
            dmaxhr = nh1-1+(ll-1)*7+4;
            dvec = [hr(damhr+1); hr(damhr+2); hr(damhr+3)];
            dmax = hr(dmaxhr);
                
            bnAdN1 = nvect*DmatL*BmatL;
            bnAdN2 = nvect*DmatR*BmatR;
            rhspulL = reshape(ulL,ndf*nelL,1);
            rhspulR = reshape(ulR,ndf*nelR,1);
            
            tvtr = 1/2*(bnAdN1*rhspulL + bnAdN2*rhspulR);
            normsig = nvec'*tvtr;
            jumpu = NmatR*rhspulR - NmatL*rhspulL;
            dn = jumpu'*nvec;
            tn = normsig+rp*dn;
            
            if tn >= 0 % tension
                
                tvec = tvtr + rp*jumpu;
                
            else % compression
                
                tvec = tvtr + rp*jumpu - tn*nvec;
                
            end
            
            normtvec = sqrt(tvec'*tvec);
            if dmax >= dc
                psik = 0;
            else
                psik = sigmax - Hc*dmax;
            end
            
            if xint > 0
                if yint > 0
                    theta = atan(yint/xint);
                else
                    theta = 2*pi + atan(yint/xint);
                end
            else
                if yint > 0
                    theta = pi + atan(yint/xint);
                else
                    theta = pi + atan(yint/xint);
                end
            end
            mvec = [-nvec(2); nvec(1); 0];
            
%             ElemI(:,ll) = [theta
%                            dn
%                            jumpu'*mvec
%                            normsig
%                            tvtr'*mvec
%                            tn
%                            (tvtr + rp*jumpu)'*mvec
%                            normtvec
%                            dvec'*nvec
%                            dvec'*mvec];
%             ElemI(:,ll) = [theta
%                            dn
%                            jumpu'*mvec
%                            normsig
%                            tvtr'*mvec
%                            (tvtr + rp*(jumpu-dvec))'*nvec
%                            (tvtr + rp*(jumpu-dvec))'*mvec
%                            normtvec
%                            dvec'*nvec
%                            dvec'*mvec];
           ElemI(:,ll) = [xint
                          yint
                          zint
                          normtvec
                          dn
                          jumpu(1)
                          jumpu(2)
                          jumpu(3)
                          tvtr(1)
                          tvtr(2)
                          tvtr(3)
                          dvec(1)
                          dvec(2)
                          dvec(3)
                          bonddbond2];
            
            end %lint
        
%         end %intt
        
end %Task Switch
