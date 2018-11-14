%% DG implementation of large deformaiton 
% pure displacement
%04/17/2013 Pinlei Chen
%% have tau and delta in it
%% have body force in it
%% verified for patchtest4_BF and patchtest5_BF but not quadratic converge for DG
%% with d_ijklmn in it
%% lint for interface part =100,interior =1000, traction =100
% Modified 01/01/2014 by Tim, verified for quadratic convergence under all
% conditions, using TT8U2DG2.m
% 3/21/2014 - TJT - Fixed all terms to verify convergence rate for general
% gamL/R tensors. Debugged by comparing with NL_Elem21_2d_2

Bcol1 = [1; 4; 6];
Bcol2 = [2; 4; 5];
Bcol3 = [3; 5; 6];
col1 = [1; 2; 3];
col2 = [2; 1; 3];
col3 = [3; 2; 1];

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

if ~exist('modifDG','var') % allows the value to be changed in input file
modifDG = 0;1;2;3;
end
% modifDG = 0 for full, 1 for penalty only, 2 for incomplete interior
% penalty (IIGP), 3 for using initial Cijkl for nonsymm interior penalty, 4
% for using Cijkl_n-1 for nonsymm interior penalty

nitvms = 1;
if nitvms == 1 %VMS parameter for the stability tensor rp
if ~exist('pencoeff','var') % allows the value to be changed in input file
pencoeff = 4;20;6;3;1;1;
end
elseif nitvms == 2 %Nitsche
pencoeff = 1;4;2;
else %RFB
pencoeff = 1;4;2;
end
        
        P1 = [eye(6,6) zeros(6,3)];
        P3 = [1 0 0 0 0 0 0 0 0
              0 0 0 1/2 0 0 -1/2 0 0
              0 0 0 0 0 1/2 0 0 1/2
              0 0 0 1/2 0 0 1/2 0 0
              0 1 0 0 0 0 0 0 0
              0 0 0 0 1/2 0 0 -1/2 0
              0 0 0 0 0 1/2 0 0 -1/2
              0 0 0 0 1/2 0 0 1/2 0
              0 0 1 0 0 0 0 0 0];
        P2 = [1 0 0 0 0 0 0 0 0
              0 0 0 1/2 0 0 1/2 0 0
              0 0 0 0 0 1/2 0 0 -1/2
              0 0 0 1/2 0 0 -1/2 0 0
              0 1 0 0 0 0 0 0 0
              0 0 0 0 1/2 0 0 1/2 0
              0 0 0 0 0 1/2 0 0 1/2
              0 0 0 0 1/2 0 0 -1/2 0
              0 0 1 0 0 0 0 0 0];

switch isw %Task Switch
    case 1
        
        if ndf > 3
            
            for i = 4:ndf
                lie(i,1) = 0;
            end
            
        end
        
        nh3 = 27;
          
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

        ElemEL = matepropL(1);
        ElemvL = matepropL(2);
        ElemER = matepropR(1);
        ElemvR = matepropR(2);
        lam = getlam(matepropL);
       
        
        BmatL1 = zeros(6,nstL);        
        BmatR1 = zeros(6,nstR);
        BmatLi = zeros(9,nstL);
        BmatRi = zeros(9,nstR);
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
        
        % Load Guass Integration Points

        if nelL == 4
            lint = 3;11;5;16;
        elseif nelL == 8
            lint = 4;100;4;8;
        elseif nelL == 10
            lint = 13;14;
%             lint = 27;
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
   iterset = 1;3;4;
   if iter  <=iterset % == 0 %
         [tauL,intb] = TauS3_PC(xlL,ulL,matepropL,nelL,lint,lam); %[Y^(-1)]       
%         [tauL,intb] = TauS3_1(xlL,ulL,mateprop,nel,lint,lam); %[Y^(-1)]
%          TauListL(elem,:) = [tauL(1),tauL(2),tauL(3),tauL(4),tauL(5),tauL(6),tauL(7),tauL(8),tauL(9)];
         [tauR,intb] = TauS3_PC(xlR,ulR,matepropR,nelR,lint,lam);
%         [tauR,intb] = TauS3_1(xlR,ulR,mateprop,nel,lint,lam);
%          TauListR(elem,:) = [tauR(1),tauR(2),tauR(3),tauR(4),tauR(5),tauR(6),tauR(7),tauR(8),tauR(9)];


        % Modifications to tau
        if exist('diagt','var') && diagt == 1 % make tauL and tauR diagonal
            if exist('scalt','var') && scalt == 1 % use a scalar of tau times the identity
                if exist('minmaxt','var') && minmaxt == -1 % minimum entry
                    tau = diag(inv(tauL));
                    tau = min(tau)*ones(2,1);
                    tauL = inv(diag(tau));
                    tau = diag(inv(tauR));
                    tau = min(tau)*ones(2,1);
                    tauR = inv(diag(tau));
                else %if minmaxt == 1 % maximum entry
                    tau = diag(inv(tauL));
                    tau = max(tau)*ones(2,1);
                    tauL = inv(diag(tau));
                    tau = diag(inv(tauR));
                    tau = max(tau)*ones(2,1);
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
                    tau = min(tau)*ones(2,1);
                    tauL = inv(diag(tau));
                    tau = diag(inv(tauR));
                    tau = min(tau)*ones(2,1);
                    tauR = inv(diag(tau));
                else %if minmaxt == 1 % maximum entry
                    tau = diag(inv(tauL));
                    tau = max(tau)*ones(2,1);
                    tauL = inv(diag(tau));
                    tau = diag(inv(tauR));
                    tau = max(tau)*ones(2,1);
                    tauR = inv(diag(tau));
                end
            end
%             else % equal only
                tau = inv(tauL + tauR);
                tauL = 1/2*(tauL*tau*tauR + tauR*tau*tauL);
                tauR = tauL;
%             end
        end

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
%         ep_list(iter+1,:,inter) = [ep(1,1),ep(1,2),ep(1,3),ep(2,1),ep(2,2),ep(2,3),ep(3,1),ep(3,2),ep(3,3)];

   else
        gamL = [hr(nh3) hr(nh3+1) hr(nh3+2)
                hr(nh3+3) hr(nh3+4) hr(nh3+5)
                hr(nh3+6) hr(nh3+7) hr(nh3+8)];
        gamR = [hr(nh3+9) hr(nh3+10) hr(nh3+11)
                hr(nh3+12) hr(nh3+13) hr(nh3+14)
                hr(nh3+15) hr(nh3+16) hr(nh3+17)];
        ep =   [hr(nh3+18) hr(nh3+19) hr(nh3+20)
                hr(nh3+21) hr(nh3+22) hr(nh3+23)
                hr(nh3+24) hr(nh3+25) hr(nh3+26)];
%         gamL = [gamL_list(iterset+1,1,inter) gamL_list(iterset+1,2,inter) gamL_list(iterset+1,3,inter)
%                 gamL_list(iterset+1,4,inter) gamL_list(iterset+1,5,inter) gamL_list(iterset+1,6,inter)
%                 gamL_list(iterset+1,7,inter) gamL_list(iterset+1,8,inter) gamL_list(iterset+1,9,inter)];
%         gamR = [gamR_list(iterset+1,1,inter) gamR_list(iterset+1,2,inter) gamR_list(iterset+1,3,inter)
%                 gamR_list(iterset+1,4,inter) gamR_list(iterset+1,5,inter) gamR_list(iterset+1,6,inter)
%                 gamR_list(iterset+1,7,inter) gamR_list(iterset+1,8,inter) gamR_list(iterset+1,9,inter)];
%         ep = [ep_list(iterset+1,1,inter) ep_list(iterset+1,2,inter) ep_list(iterset+1,3,inter)
%               ep_list(iterset+1,4,inter) ep_list(iterset+1,5,inter) ep_list(iterset+1,6,inter)
%               ep_list(iterset+1,7,inter) ep_list(iterset+1,8,inter) ep_list(iterset+1,9,inter)];
   end
%        gamL =0.5*eye(3,3) ;
%        gamR =0.5*eye(3,3) ; 
%        ep = 2000*eye(3,3) ;   
   

         for l = 1:lint

                ll = ll + 1;

                % Evaluate  basis functions at integration points
                if nelL == 4 || nelL == 10
                  [Wgt,ss] = int3d_t(l,lint,ib);
                  [shlL,shldL,shls,be] = shltt(ss,nelL,nelL,0,0);
                  [PxyL,shgsL,JdetL,bubbleL,XsL] = shgtt(xlL(:,1:nelL),nelL,shldL,shls,nen,0,0,be);                 
                  [QxyL,shgsL,JdetL,bubbleL,xsL] = shgtt(xlL(:,1:nelL)+ulL(:,1:nelL),nelL,shldL,shls,nen,0,0,be);
                else
                  [Wgt,ss] = intpntb(l,lint,ib);
                  [shlL,shldL,shls,be] = shlb(ss,nelL,nelL,0,0);
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
                  [shlR,shldR,shls,be] = shltt([rR sR tR],nelR,nelR,0,0);
                  [PxyR,shgsR,JdetR,bubbleR,XsR] = shgtt(xlR(:,1:nelR),nelR,shldR,shls,nen,0,0,be);                  
                  [QxyR,shgsR,JdetR,bubbleR,xsR] = shgtt(xlR(:,1:nelR)+ulR(:,1:nelR),nelR,shldR,shls,nen,0,0,be);
                  
                else
                  [shlR,shldR,shls,be] = shlb([rR sR tR],nelR,nelR,0,0);
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
           end  
           
           if modifDG == 3
           for mm = 1:nelL  
BmatLi(:,3*mm-2:3*mm)=[PxyL(mm,1) 0         0         
                        0         PxyL(mm,2) 0         
                        0         0         PxyL(mm,3) 
                        PxyL(mm,2) PxyL(mm,1) 0         
                        0         PxyL(mm,3) PxyL(mm,2) 
                        PxyL(mm,3) 0         PxyL(mm,1) 
                        PxyL(mm,2) -PxyL(mm,1) 0         
                        0         PxyL(mm,3) -PxyL(mm,2) 
                        -PxyL(mm,3) 0         PxyL(mm,1) ];
           end
            
           for mm = 1:nelR    
 BmatRi(:,3*mm-2:3*mm)=[PxyR(mm,1) 0         0         
                        0         PxyR(mm,2) 0         
                        0         0         PxyR(mm,3) 
                        PxyR(mm,2) PxyR(mm,1) 0         
                        0         PxyR(mm,3) PxyR(mm,2) 
                        PxyR(mm,3) 0         PxyR(mm,1) 
                        PxyR(mm,2) -PxyR(mm,1) 0         
                        0         PxyR(mm,3) -PxyR(mm,2) 
                        -PxyR(mm,3) 0         PxyR(mm,1) ];                  
           end  
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
          
            [sigma2L, cmatL] = SigmaCmat3i(FL,JxXL,matepropL,lam);
            [sigma2Li, cmatLi] = SigmaCmat3i(eye(3),1,matepropL,lam);
            
            SmatL1=[sigma2L(1), sigma2L(4), sigma2L(6)
                    sigma2L(4), sigma2L(2), sigma2L(5)
                    sigma2L(6), sigma2L(5), sigma2L(3)];        
            SmatL1i=[sigma2Li(1), sigma2Li(4), sigma2Li(6)
                    sigma2Li(4), sigma2Li(2), sigma2Li(5)
                    sigma2Li(6), sigma2Li(5), sigma2Li(3)];
            
            [sigma2R, cmatR] = SigmaCmat3i(FR,JxXR,matepropR,lam);
            [sigma2Ri, cmatRi] = SigmaCmat3i(eye(3),1,matepropR,lam);
            
            SmatR1=[sigma2R(1), sigma2R(4), sigma2R(6)
                    sigma2R(4), sigma2R(2), sigma2R(5)
                    sigma2R(6), sigma2R(5), sigma2R(3)]; 
            SmatR1i=[sigma2Ri(1), sigma2Ri(4), sigma2Ri(6)
                    sigma2Ri(4), sigma2Ri(2), sigma2Ri(5)
                    sigma2Ri(6), sigma2Ri(5), sigma2Ri(3)];       

                %Evaluate tangent and normal vectors
                % jacobian j=dx/dkesi for surface integral for current
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
                % jacobian j=dx/dkesi for surface integral for material
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
                %Normal vectors
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
%                 nvectL = [nLx 0  0  nLy  0  nLz nLy 0  -nLz
%                           0  nLy 0  nLx nLz  0 -nLx nLz  0
%                           0   0 nLz  0  nLy nLx  0 -nLy nLx];
                nvectL1= [nLx 0  0  nLy  0  nLz 
                          0  nLy 0  nLx nLz  0 
                          0   0 nLz  0  nLy nLx];
%                 nvectR = [nRx 0  0  nRy  0  nRz nRy 0  -nRz
%                           0  nRy 0  nRx  nRz 0 -nRx nRz  0
%                           0   0 nRz  0  nRy nRx  0 -nRy nRx]; 
                nvectR1 = [nRx 0  0  nRy  0  nRz 
                          0  nRy 0  nRx  nRz 0 
                          0   0 nRz  0  nRy nRx];  
%                 nvectL2 = [eye(6,6)*nLx zeros(6,6)    zeros(6,6)    eye(6,6)*nLy  zeros(6,6)    eye(6,6)*nLz
%                            zeros(6,6)    eye(6,6)*nLy zeros(6,6)    eye(6,6)*nLx  eye(6,6)*nLz  zeros(6,6)
%                            zeros(6,6)    zeros(6,6)   eye(6,6)*nLz   zeros(6,6)   eye(6,6)*nLy  eye(6,6)*nLx];
%                 nvectR2 = [eye(6,6)*nRx zeros(6,6)    zeros(6,6)    eye(6,6)*nRy  zeros(6,6)    eye(6,6)*nRz
%                            zeros(6,6)    eye(6,6)*nRy zeros(6,6)    eye(6,6)*nRx  eye(6,6)*nRz  zeros(6,6)
%                            zeros(6,6)    zeros(6,6)   eye(6,6)*nRz   zeros(6,6)   eye(6,6)*nRy  eye(6,6)*nRx];
                nvecL = [nLx; nLy; nLz];
                nvecR = [nRx; nRy; nRz];
                NvectL1= [NLx 0  0  NLy  0  NLz 
                          0  NLy 0  NLx NLz  0 
                          0   0 NLz  0  NLy NLx];
                NvectR1 = [NRx 0  0  NRy  0  NRz 
                          0  NRy 0  NRx  NRz 0 
                          0   0 NRz  0  NRy NRx]; 


            SmatnL=[SmatL1*nvecL zeros(ndm,2)
                    zeros(ndm,1) SmatL1*nvecL zeros(ndm,1)
                    zeros(ndm,2) SmatL1*nvecL];
                %%%%%NOTE: I will need to compute N at u_n-1 to make these
                %%%%%terms for the constant version
            
            SmatnR=[SmatR1*nvecR zeros(ndm,2)
                    zeros(ndm,1) SmatR1*nvecR zeros(ndm,1)
                    zeros(ndm,2) SmatR1*nvecR];

           cmatnL=(nvectL1*cmatL(1:6,1:6));
           cmatnR=(nvectR1*cmatR(1:6,1:6));
           cmatnBL=BmatL'*P2'*[cmatnL zeros(3,12)
                    zeros(3,6)   cmatnL   zeros(3,6)
                    zeros(3,12)  cmatnL ]; 
                
           cmatnBR=BmatR'*P2'*[cmatnR    zeros(3,12)
                    zeros(3,6)   cmatnR   zeros(3,6)
                    zeros(3,12)  cmatnR ];

            term17L=P2'*SmatnL*gamL';
            term17R=P2'*SmatnR*gamR';
            term17Li = 0*term17L;
            term17Ri = 0*term17R;

            term18L=P1'*(gamL*nvectL1*cmatL(1:6,1:6))';
            term18R=P1'*(gamR*nvectR1*cmatR(1:6,1:6))';
            term18Li=P1'*(gamL*NvectL1*cmatLi(1:6,1:6))';
            term18Ri=P1'*(gamR*NvectR1*cmatRi(1:6,1:6))';

            term28L=NmatL'*(c1L*gamL*SmatL1*nvecL-c1R*gamR*SmatR1*nvecR);   %average stress term
            term28R=NmatR'*(c1L*gamL*SmatL1*nvecL-c1R*gamR*SmatR1*nvecR);           

             rhspulL = reshape(ulL,ndf*nelL,1);
             rhspulR = reshape(ulR,ndf*nelR,1);
                   
            jumpu = NmatL*rhspulL - NmatR*rhspulR;  %jumpu=uL-uR
             
             term30L=NmatL'*ep*jumpu;
             term30R=NmatR'*ep*jumpu;
            

            % Unexpected terms, involves jumpu and 6th order tensor
             gamajumpuL = gamL'*jumpu;
             gamajumpuR = gamR'*jumpu;             
             term5L=cmatnBL*[eye(6,6)*gamajumpuL(1) eye(6,6)*gamajumpuL(2) eye(6,6)*gamajumpuL(3)]'*(P1*BmatL);
             term5R=cmatnBR*[eye(6,6)*gamajumpuR(1) eye(6,6)*gamajumpuR(2) eye(6,6)*gamajumpuR(3)]'*(P1*BmatR);

             sig8L2 = [gamajumpuL(1) gamajumpuL(2) gamajumpuL(3)]*cmatnL;
             
          sig8L3 = [sig8L2(1) 0 0 sig8L2(4) 0 0 sig8L2(6) 0 0 ;
                    sig8L2(4) 0 0 sig8L2(2) 0 0 sig8L2(5) 0 0 ;
                    sig8L2(6) 0 0 sig8L2(5) 0 0 sig8L2(3) 0 0 ;
                    0 sig8L2(1) 0 0 sig8L2(4) 0 0 sig8L2(6) 0;
                    0 sig8L2(4) 0 0 sig8L2(2) 0 0 sig8L2(5) 0;
                    0 sig8L2(6) 0 0 sig8L2(5) 0 0 sig8L2(3) 0;
                    0 0 sig8L2(1) 0 0 sig8L2(4) 0 0 sig8L2(6);
                    0 0 sig8L2(4) 0 0 sig8L2(2) 0 0 sig8L2(5);
                    0 0 sig8L2(6) 0 0 sig8L2(5) 0 0 sig8L2(3)];
               
           term8L = BmatL'*P2'*sig8L3*(P3*BmatL);
           
             sig8R2 = [gamajumpuR(1) gamajumpuR(2) gamajumpuR(3)]*cmatnR;
             
          sig8R3 = [sig8R2(1) 0 0 sig8R2(4) 0 0 sig8R2(6) 0 0 ;
                    sig8R2(4) 0 0 sig8R2(2) 0 0 sig8R2(5) 0 0 ;
                    sig8R2(6) 0 0 sig8R2(5) 0 0 sig8R2(3) 0 0 ;
                    0 sig8R2(1) 0 0 sig8R2(4) 0 0 sig8R2(6) 0;
                    0 sig8R2(4) 0 0 sig8R2(2) 0 0 sig8R2(5) 0;
                    0 sig8R2(6) 0 0 sig8R2(5) 0 0 sig8R2(3) 0;
                    0 0 sig8R2(1) 0 0 sig8R2(4) 0 0 sig8R2(6);
                    0 0 sig8R2(4) 0 0 sig8R2(2) 0 0 sig8R2(5);
                    0 0 sig8R2(6) 0 0 sig8R2(5) 0 0 sig8R2(3)];

              term8R =BmatR'*P2'*sig8R3*(P3*BmatR);
              
             [dmatL]=dmat3(JxXL,matepropL,lam);
             [dmatR]=dmat3(JxXR,matepropR,lam);   
%              dmatL2 = [eye(6,6)*jumpu(1) eye(6,6)*jumpu(2) eye(6,6)*jumpu(3)]*[eye(6,6)*gamL(1,1) eye(6,6)*gamL(1,2) eye(6,6)*gamL(1,3)
%                                                                                eye(6,6)*gamL(2,1) eye(6,6)*gamL(2,2) eye(6,6)*gamL(2,3)
%                                                                                eye(6,6)*gamL(3,1) eye(6,6)*gamL(3,2) eye(6,6)*gamL(3,3)]*nvectL2*dmatL/JxXL;
%              dmatR2 = [eye(6,6)*jumpu(1) eye(6,6)*jumpu(2) eye(6,6)*jumpu(3)]*[eye(6,6)*gamR(1,1) eye(6,6)*gamR(1,2) eye(6,6)*gamR(1,3)
%                                                                                eye(6,6)*gamR(2,1) eye(6,6)*gamR(2,2) eye(6,6)*gamR(2,3)
%                                                                                eye(6,6)*gamR(3,1) eye(6,6)*gamR(3,2) eye(6,6)*gamR(3,3)]*nvectR2*dmatR/JxXR;
             % Revised by TJT to use reshape function and smaller
             % multiplications
             dmatL2 = dmatL/JxXL*nvectL1'*gamL'*jumpu;
             dmatL2 = reshape(dmatL2,6,6);
             dmatR2 = dmatR/JxXR*nvectR1'*gamR'*jumpu;
             dmatR2 = reshape(dmatR2,6,6);
             
             term7L = BmatL'*P1'*dmatL2*(P1*BmatL);
             term7R = BmatR'*P1'*dmatR2*(P1*BmatR);   
             
             
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
             
             if modifDG >= 1 % IVMDG
                 
             % average stress terms
                 
             ElemFL = ElemFL-(-term28L);
             ElemFR = ElemFR-(+term28R);
% 
             ElemKLL = ElemKLL  - c1L*NmatL'*(term17L'+term18L')*BmatL;
             ElemKLR = ElemKLR  + c1R*NmatL'*(term17R'+term18R')*BmatR;
             ElemKRL = ElemKRL  + c1L*NmatR'*(term17L'+term18L')*BmatL;
             ElemKRR = ElemKRR  - c1R*NmatR'*(term17R'+term18R')*BmatR;
             
             if modifDG == 2 %RVMDG
                 
             % constant nonsymmetric terms

             ElemFL = ElemFL-(-C1L*BmatLi'*(term17Li+term18Li)*jumpu);
             ElemFR = ElemFR-(C1R* BmatRi'*(term17Ri+term18Ri)*jumpu);

             ElemKLL = ElemKLL - C1L*BmatLi'*(term17Li+term18Li)*NmatL;
             ElemKLR = ElemKLR + C1L*BmatLi'*(term17Li+term18Li)*NmatR;
             ElemKRL = ElemKRL + C1R*BmatRi'*(term17Ri+term18Ri)*NmatL;
             ElemKRR = ElemKRR - C1R*BmatRi'*(term17Ri+term18Ri)*NmatR;
             
             elseif modifDG == 3 %VMDGs
                 
             % nonsymmetric terms

             ElemFL = ElemFL-(-c1L*BmatL'*(term17L+term18L)*jumpu);
             ElemFR = ElemFR-(c1R* BmatR'*(term17R+term18R)*jumpu);

             ElemKLL = ElemKLL - c1L*BmatL'*(term17L+term18L)*NmatL;
             ElemKLR = ElemKLR + c1L*BmatL'*(term17L+term18L)*NmatR;
             ElemKRL = ElemKRL + c1R*BmatR'*(term17R+term18R)*NmatL;
             ElemKRR = ElemKRR - c1R*BmatR'*(term17R+term18R)*NmatR;
             
             elseif modifDG == 4 %IVMDGs
                 
             % nonsymmetric terms

             ElemKLL = ElemKLL - c1L*BmatL'*(term17L+term18L)*NmatL;
             ElemKLR = ElemKLR + c1L*BmatL'*(term17L+term18L)*NmatR;
             ElemKRL = ElemKRL + c1R*BmatR'*(term17R+term18R)*NmatL;
             ElemKRR = ElemKRR - c1R*BmatR'*(term17R+term18R)*NmatR;
             
             elseif modifDG == 5 %RVMDGs
                 
             % constant nonsymmetric terms

             ElemFL = ElemFL-(-C1L*BmatLi'*(term17Li+term18Li)*jumpu);
             ElemFR = ElemFR-(C1R* BmatRi'*(term17Ri+term18Ri)*jumpu);
                 
             % nonsymmetric terms

             ElemKLL = ElemKLL - c1L*BmatL'*(term17L+term18L)*NmatL;
             ElemKLR = ElemKLR + c1L*BmatL'*(term17L+term18L)*NmatR;
             ElemKRL = ElemKRL + c1R*BmatR'*(term17R+term18R)*NmatL;
             ElemKRR = ElemKRR - c1R*BmatR'*(term17R+term18R)*NmatR;
             
             end
             end

             else % full method VMDG
             
                ElemFL = ElemFL-(-term28L+C1L*term30L);
                ElemFR = ElemFR-(+term28R-C1R*term30R);  
                
                ElemFL = ElemFL-(-c1L*BmatL'*(term17L+term18L)*jumpu);
                ElemFR = ElemFR-(c1R* BmatR'*(term17R+term18R)*jumpu); 
 
                ElemKLL = ElemKLL + C1L*(NmatL'*ep*NmatL);
                ElemKLR = ElemKLR - C1R*(NmatL'*ep*NmatR);
                ElemKRL = ElemKRL - C1L*(NmatR'*ep*NmatL);
                ElemKRR = ElemKRR + C1R*(NmatR'*ep*NmatR);   

             ElemKLL = ElemKLL  - c1L*NmatL'*(term17L'+term18L')*BmatL;
             ElemKLR = ElemKLR  + c1R*NmatL'*(term17R'+term18R')*BmatR;
             ElemKRL = ElemKRL  + c1L*NmatR'*(term17L'+term18L')*BmatL;
             ElemKRR = ElemKRR  - c1R*NmatR'*(term17R'+term18R')*BmatR;             

             ElemKLL = ElemKLL - c1L*BmatL'*(term17L+term18L)*NmatL;
             ElemKLR = ElemKLR + c1L*BmatL'*(term17L+term18L)*NmatR;
             ElemKRL = ElemKRL + c1R*BmatR'*(term17R+term18R)*NmatL;
             ElemKRR = ElemKRR - c1R*BmatR'*(term17R+term18R)*NmatR;           

            %relate to jumpu terms additional terms
                ElemKLL = ElemKLL - c1L*(term5L+term5L'+term7L+term8L);  %                     
                ElemKRR = ElemKRR + c1R*(term5R+term5R'+term7R+term8R);  %

             end
%             if l == 1 && elem == 87 && step > 0
%                % find the first PK stress and Cijkl in refrence config. for
%                % left side
%                P = tranr4(fiL,fiL); % transformation tensor F_iI*F_jJ
%                P = [P' zeros(6,3); zeros(3,9)];
%                Se = P*(sigma2L*JxXL); % second PK stress
%                Se_m = [Se(1) Se(4) Se(6); Se(4) Se(2) Se(5); Se(6) Se(5) Se(3)];
%                Cmat_e = P*(cmatL*JxXL)*P';
%                Pe = FL*Se_m; % first PK stress 
%                
%                Mate_tau(:,:,step) = ep; 
%                Mate_tauL(:,:,step) = tauL; 
%                Mate_tauR(:,:,step) = tauR; 
%                Mate_smat(:,:,step) = Se_m;
%                Mate_cmat(:,:,step) = Cmat_e;
%                Mate_F(:,:,step) = FL;
%             end
%             if jumpu>jumpu_max
%                 jumpu_max = jumpu;
%             end
            end %lint
%          jumpu_max = [jumpu_max; jumpu];
%         end %intt
% ElemKLL
            ElemK = [ElemKLL ElemKLR; ElemKRL ElemKRR];
            ElemF = [ElemFL; ElemFR];
 
    case 6  % post processing assemble internal force
        
        ElemFR = zeros(nstR,1);
        ElemFL = zeros(nstL,1);
        t1 = zeros(3,1);
        t2 = t1;
        t3 = t1;

        ElemEL = matepropL(1);
        ElemvL = matepropL(2);
        ElemER = matepropR(1);
        ElemvR = matepropR(2);
        lam = getlam(matepropL);
       
        
        BmatL1 = zeros(6,nstL);        
        BmatR1 = zeros(6,nstR);
        BmatLi = zeros(9,nstL);
        BmatRi = zeros(9,nstR);
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
        
        % Load Guass Integration Points

        if nelL == 4
            lint = 3;11;5;16;
        elseif nelL == 8
            lint = 4;100;4;8;
        elseif nelL == 10
            lint = 14;
%             lint = 27;
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
   iterset = 1;3;4;
   if iter  <=iterset % == 0 %
         [tauL,intb] = TauS3_PC(xlL,ulL,matepropL,nelL,lint,lam); %[Y^(-1)]       
%         [tauL,intb] = TauS3_1(xlL,ulL,mateprop,nel,lint,lam); %[Y^(-1)]
%          TauListL(elem,:) = [tauL(1),tauL(2),tauL(3),tauL(4),tauL(5),tauL(6),tauL(7),tauL(8),tauL(9)];
         [tauR,intb] = TauS3_PC(xlR,ulR,matepropR,nelR,lint,lam);
%         [tauR,intb] = TauS3_1(xlR,ulR,mateprop,nel,lint,lam);
%          TauListR(elem,:) = [tauR(1),tauR(2),tauR(3),tauR(4),tauR(5),tauR(6),tauR(7),tauR(8),tauR(9)];


        % Modifications to tau
        if exist('diagt','var') && diagt == 1 % make tauL and tauR diagonal
            if exist('scalt','var') && scalt == 1 % use a scalar of tau times the identity
                if exist('minmaxt','var') && minmaxt == -1 % minimum entry
                    tau = diag(inv(tauL));
                    tau = min(tau)*ones(2,1);
                    tauL = inv(diag(tau));
                    tau = diag(inv(tauR));
                    tau = min(tau)*ones(2,1);
                    tauR = inv(diag(tau));
                else %if minmaxt == 1 % maximum entry
                    tau = diag(inv(tauL));
                    tau = max(tau)*ones(2,1);
                    tauL = inv(diag(tau));
                    tau = diag(inv(tauR));
                    tau = max(tau)*ones(2,1);
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
                    tau = min(tau)*ones(2,1);
                    tauL = inv(diag(tau));
                    tau = diag(inv(tauR));
                    tau = min(tau)*ones(2,1);
                    tauR = inv(diag(tau));
                else %if minmaxt == 1 % maximum entry
                    tau = diag(inv(tauL));
                    tau = max(tau)*ones(2,1);
                    tauL = inv(diag(tau));
                    tau = diag(inv(tauR));
                    tau = max(tau)*ones(2,1);
                    tauR = inv(diag(tau));
                end
            end
%             else % equal only
                tau = inv(tauL + tauR);
                tauL = 1/2*(tauL*tau*tauR + tauR*tau*tauL);
                tauR = tauL;
%             end
        end

        ebL = 0;
        ebR = 0;
        intedge = 0;
        
        for ie = 1:lint            
% For separate bubble types on T and Q
           if nelL == 4 || nelL == 10
                [w,ss] =  int3d_t(ie,lint,ib);
                ebeL = facebubbleT(ss);
                [shlL,shldL,shls,be] = shltt(ss,nelL,nelL,0,0);
           else
                [Wgt,ss] = intpntb(ie,lint,ib);
                [shlL,shldL,shls,be] = shlb(ss,nelL,nelL,0,0);
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
                [w,ss] =  int3d_t(ie,lint,ib);
                ebeR = facebubbleT(ss);
           else
                [Wgt,ss] = intpntb(ie,lint,ib);
                ebeR = facebubbleQ(ss);
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
%         ep_list(iter+1,:,inter) = [ep(1,1),ep(1,2),ep(1,3),ep(2,1),ep(2,2),ep(2,3),ep(3,1),ep(3,2),ep(3,3)];

   else
        gamL = [hr(nh3) hr(nh3+1) hr(nh3+2)
                hr(nh3+3) hr(nh3+4) hr(nh3+5)
                hr(nh3+6) hr(nh3+7) hr(nh3+8)];
        gamR = [hr(nh3+9) hr(nh3+10) hr(nh3+11)
                hr(nh3+12) hr(nh3+13) hr(nh3+14)
                hr(nh3+15) hr(nh3+16) hr(nh3+17)];
        ep =   [hr(nh3+18) hr(nh3+19) hr(nh3+20)
                hr(nh3+21) hr(nh3+22) hr(nh3+23)
                hr(nh3+24) hr(nh3+25) hr(nh3+26)];
%         gamL = [gamL_list(iterset+1,1,inter) gamL_list(iterset+1,2,inter) gamL_list(iterset+1,3,inter)
%                 gamL_list(iterset+1,4,inter) gamL_list(iterset+1,5,inter) gamL_list(iterset+1,6,inter)
%                 gamL_list(iterset+1,7,inter) gamL_list(iterset+1,8,inter) gamL_list(iterset+1,9,inter)];
%         gamR = [gamR_list(iterset+1,1,inter) gamR_list(iterset+1,2,inter) gamR_list(iterset+1,3,inter)
%                 gamR_list(iterset+1,4,inter) gamR_list(iterset+1,5,inter) gamR_list(iterset+1,6,inter)
%                 gamR_list(iterset+1,7,inter) gamR_list(iterset+1,8,inter) gamR_list(iterset+1,9,inter)];
%         ep = [ep_list(iterset+1,1,inter) ep_list(iterset+1,2,inter) ep_list(iterset+1,3,inter)
%               ep_list(iterset+1,4,inter) ep_list(iterset+1,5,inter) ep_list(iterset+1,6,inter)
%               ep_list(iterset+1,7,inter) ep_list(iterset+1,8,inter) ep_list(iterset+1,9,inter)];
   end
%        gamL =0.5*eye(3,3) ;
%        gamR =0.5*eye(3,3) ; 
%        ep = 2000*eye(3,3) ;   
   

         for l = 1:lint

                ll = ll + 1;

                % Evaluate  basis functions at integration points
                if nelL == 4 || nelL == 10
                  [Wgt,ss] = int3d_t(l,lint,ib);
                  [shlL,shldL,shls,be] = shltt(ss,nelL,nelL,0,0);
                  [PxyL,shgsL,JdetL,bubbleL,XsL] = shgtt(xlL(:,1:nelL),nelL,shldL,shls,nen,0,0,be);                 
                  [QxyL,shgsL,JdetL,bubbleL,xsL] = shgtt(xlL(:,1:nelL)+ulL(:,1:nelL),nelL,shldL,shls,nen,0,0,be);
                else
                  [Wgt,ss] = intpntb(l,lint,ib);
                  [shlL,shldL,shls,be] = shlb(ss,nelL,nelL,0,0);
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
                  [shlR,shldR,shls,be] = shltt([rR sR tR],nelR,nelR,0,0);
                  [PxyR,shgsR,JdetR,bubbleR,XsR] = shgtt(xlR(:,1:nelR),nelR,shldR,shls,nen,0,0,be);                  
                  [QxyR,shgsR,JdetR,bubbleR,xsR] = shgtt(xlR(:,1:nelR)+ulR(:,1:nelR),nelR,shldR,shls,nen,0,0,be);
                  
                else
                  [shlR,shldR,shls,be] = shlb([rR sR tR],nelR,nelR,0,0);
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
           end  
           
           if modifDG == 3
           for mm = 1:nelL  
BmatLi(:,3*mm-2:3*mm)=[PxyL(mm,1) 0         0         
                        0         PxyL(mm,2) 0         
                        0         0         PxyL(mm,3) 
                        PxyL(mm,2) PxyL(mm,1) 0         
                        0         PxyL(mm,3) PxyL(mm,2) 
                        PxyL(mm,3) 0         PxyL(mm,1) 
                        PxyL(mm,2) -PxyL(mm,1) 0         
                        0         PxyL(mm,3) -PxyL(mm,2) 
                        -PxyL(mm,3) 0         PxyL(mm,1) ];
           end
            
           for mm = 1:nelR    
 BmatRi(:,3*mm-2:3*mm)=[PxyR(mm,1) 0         0         
                        0         PxyR(mm,2) 0         
                        0         0         PxyR(mm,3) 
                        PxyR(mm,2) PxyR(mm,1) 0         
                        0         PxyR(mm,3) PxyR(mm,2) 
                        PxyR(mm,3) 0         PxyR(mm,1) 
                        PxyR(mm,2) -PxyR(mm,1) 0         
                        0         PxyR(mm,3) -PxyR(mm,2) 
                        -PxyR(mm,3) 0         PxyR(mm,1) ];                  
           end  
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
          
            [sigma2L, cmatL] = SigmaCmat3i(FL,JxXL,matepropL,lam);
            [sigma2Li, cmatLi] = SigmaCmat3i(eye(3),1,matepropL,lam);
            
            SmatL1=[sigma2L(1), sigma2L(4), sigma2L(6)
                    sigma2L(4), sigma2L(2), sigma2L(5)
                    sigma2L(6), sigma2L(5), sigma2L(3)];        
            SmatL1i=[sigma2Li(1), sigma2Li(4), sigma2Li(6)
                    sigma2Li(4), sigma2Li(2), sigma2Li(5)
                    sigma2Li(6), sigma2Li(5), sigma2Li(3)];
            
            [sigma2R, cmatR] = SigmaCmat3i(FR,JxXR,matepropR,lam);
            [sigma2Ri, cmatRi] = SigmaCmat3i(eye(3),1,matepropR,lam);
            
            SmatR1=[sigma2R(1), sigma2R(4), sigma2R(6)
                    sigma2R(4), sigma2R(2), sigma2R(5)
                    sigma2R(6), sigma2R(5), sigma2R(3)]; 
            SmatR1i=[sigma2Ri(1), sigma2Ri(4), sigma2Ri(6)
                    sigma2Ri(4), sigma2Ri(2), sigma2Ri(5)
                    sigma2Ri(6), sigma2Ri(5), sigma2Ri(3)];       

                %Evaluate tangent and normal vectors
                % jacobian j=dx/dkesi for surface integral for current
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
                % jacobian j=dx/dkesi for surface integral for material
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
                %Normal vectors
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
%                 nvectL2 = [eye(6,6)*nLx zeros(6,6)    zeros(6,6)    eye(6,6)*nLy  zeros(6,6)    eye(6,6)*nLz
%                            zeros(6,6)    eye(6,6)*nLy zeros(6,6)    eye(6,6)*nLx  eye(6,6)*nLz  zeros(6,6)
%                            zeros(6,6)    zeros(6,6)   eye(6,6)*nLz   zeros(6,6)   eye(6,6)*nLy  eye(6,6)*nLx];
%                 nvectR2 = [eye(6,6)*nRx zeros(6,6)    zeros(6,6)    eye(6,6)*nRy  zeros(6,6)    eye(6,6)*nRz
%                            zeros(6,6)    eye(6,6)*nRy zeros(6,6)    eye(6,6)*nRx  eye(6,6)*nRz  zeros(6,6)
%                            zeros(6,6)    zeros(6,6)   eye(6,6)*nRz   zeros(6,6)   eye(6,6)*nRy  eye(6,6)*nRx];
                nvecL = [nLx; nLy; nLz];
                nvecR = [nRx; nRy; nRz];
                NvectL1= [NLx 0  0  NLy  0  NLz 
                          0  NLy 0  NLx NLz  0 
                          0   0 NLz  0  NLy NLx];
                NvectR1 = [NRx 0  0  NRy  0  NRz 
                          0  NRy 0  NRx  NRz 0 
                          0   0 NRz  0  NRy NRx]; 


            SmatnL=[SmatL1*nvecL zeros(ndm,2)
                    zeros(ndm,1) SmatL1*nvecL zeros(ndm,1)
                    zeros(ndm,2) SmatL1*nvecL];
                %%%%%NOTE: I will need to compute N at u_n-1 to make these
                %%%%%terms for the constant version
            
            SmatnR=[SmatR1*nvecR zeros(ndm,2)
                    zeros(ndm,1) SmatR1*nvecR zeros(ndm,1)
                    zeros(ndm,2) SmatR1*nvecR];

           cmatnL=(nvectL1*cmatL(1:6,1:6));
           cmatnR=(nvectR1*cmatR(1:6,1:6));
           cmatnBL=BmatL'*P2'*[cmatnL zeros(3,12)
                    zeros(3,6)   cmatnL   zeros(3,6)
                    zeros(3,12)  cmatnL ]; 
                
           cmatnBR=BmatR'*P2'*[cmatnR    zeros(3,12)
                    zeros(3,6)   cmatnR   zeros(3,6)
                    zeros(3,12)  cmatnR ];

            term17L=P2'*SmatnL*gamL';
            term17R=P2'*SmatnR*gamR';
            term17Li = 0*term17L;
            term17Ri = 0*term17R;

            term18L=P1'*(gamL*nvectL1*cmatL(1:6,1:6))';
            term18R=P1'*(gamR*nvectR1*cmatR(1:6,1:6))';
            term18Li=P1'*(gamL*NvectL1*cmatLi(1:6,1:6))';
            term18Ri=P1'*(gamR*NvectR1*cmatRi(1:6,1:6))';

            term28L=NmatL'*(c1L*gamL*SmatL1*nvecL-c1R*gamR*SmatR1*nvecR);   %average stress term
            term28R=NmatR'*(c1L*gamL*SmatL1*nvecL-c1R*gamR*SmatR1*nvecR);           

             rhspulL = reshape(ulL,ndf*nelL,1);
             rhspulR = reshape(ulR,ndf*nelR,1);
                   
            jumpu = NmatL*rhspulL - NmatR*rhspulR;  %jumpu=uL-uR
             
             term30L=NmatL'*ep*jumpu;
             term30R=NmatR'*ep*jumpu;
             
             
            % Combine contributions into element force vector and stiffness
            % matrix
             if exist('modifDG','var') && modifDG > 0 % modify the DG terms used in the formulation
             
             % Penalty terms
                 
             ElemFL = ElemFL-(+C1L*term30L);
             ElemFR = ElemFR-(-C1R*term30R);
             
             if modifDG > 1 
                 
             % average stress terms
                 
             ElemFL = ElemFL-(-term28L);
             ElemFR = ElemFR-(+term28R);
             
             if modifDG == 3 
                 
             % constant nonsymmetric terms

             ElemFL = ElemFL-(-C1L*BmatLi'*(term17Li+term18Li)*jumpu);
             ElemFR = ElemFR-(C1R* BmatRi'*(term17Ri+term18Ri)*jumpu);
             
             elseif modifDG > 3 
                 
             % nonsymmetric terms

             ElemFL = ElemFL-(-c1L*BmatL'*(term17L+term18L)*jumpu);
             ElemFR = ElemFR-(c1R* BmatR'*(term17R+term18R)*jumpu);
             
             end
             end

             else % full method
             
                ElemFL = ElemFL-(-term28L+C1L*term30L);
                ElemFR = ElemFR-(+term28R-C1R*term30R);  
                
                ElemFL = ElemFL-(-c1L*BmatL'*(term17L+term18L)*jumpu);
                ElemFR = ElemFR-(c1R* BmatR'*(term17R+term18R)*jumpu);

             end
%             if l == 1 && elem == 87 && step > 0
%                % find the first PK stress and Cijkl in refrence config. for
%                % left side
%                P = tranr4(fiL,fiL); % transformation tensor F_iI*F_jJ
%                P = [P' zeros(6,3); zeros(3,9)];
%                Se = P*(sigma2L*JxXL); % second PK stress
%                Se_m = [Se(1) Se(4) Se(6); Se(4) Se(2) Se(5); Se(6) Se(5) Se(3)];
%                Cmat_e = P*(cmatL*JxXL)*P';
%                Pe = FL*Se_m; % first PK stress 
%                
%                Mate_tau(:,:,step) = ep; 
%                Mate_tauL(:,:,step) = tauL; 
%                Mate_tauR(:,:,step) = tauR; 
%                Mate_smat(:,:,step) = Se_m;
%                Mate_cmat(:,:,step) = Cmat_e;
%                Mate_F(:,:,step) = FL;
%             end
%             if jumpu>jumpu_max
%                 jumpu_max = jumpu;
%             end
            end %lint
%          jumpu_max = [jumpu_max; jumpu];
%         end %intt
% ElemKLL
            ElemF = [ElemFL; ElemFR];

    case -1
        
        ElemF = zeros(nst,1);
 %%
    case -2
        
        ElemF = zeros(nst,1);
%%        


    case 51 % Volume stress/strain homogenization
        
        ElemSS = zeros(13,1);

        %Set integration number
        lint = IntPoint3(nel);
        ib = 0;
        bf = 0;
        der = 0;

%         lam = Patchv*PatchE/((1+Patchv)*(1-2*Patchv));
%         mu = PatchE/(2*(1+Patchv));
        thick = 1;
        fbx = 0;
        fby = 0;
        fbz = 0;
        Dmat = mu*diag([2 2 2 1 1 1]) + lam*[1; 1; 1; 0; 0; 0]*[1 1 1 0 0 0];
        Nmat = zeros(3,3*nel);
        Bmat = zeros(6,3*nel);

        for l = 1:lint                    

            % Evaluate 1-D basis functions at integration points
            if nel == 4 || nel == 10
              [w,ss] =  int3d_t(l,lint,ib);
              [shl,shld,shls,be] = shltt(ss,nel,nel,der,bf);
              [shg,shgs,Jdet,be] = shgtt(xl,nel,shld,shls,nel,bf,der,be);
            else
              [w,ss] =  intpntb(l,lint,ib);
              [shl,shld,shls,be] = shlb(ss,nel,nel,der,bf);
              [shg,shgs,Jdet,be] = shgb(xl,nel,shld,shls,nel,bf,der,be);
            end
            
            if Jdet < 0
                elem
                Jdet %#ok<NOPTS>
            end
            c1 = Jdet*w*thick;
            
            % Form B matrix
            for ie = 1:nel
              
              Bmat(Bcol1,(ie-1)*3+1) = shg(ie,col1);
              Bmat(Bcol2,(ie-1)*3+2) = shg(ie,col2);
              Bmat(Bcol3,3*ie      ) = shg(ie,col3);
                 
            end
            
            epsil = Bmat*reshape(ul,ndf*nel,1);
            stres = Dmat*epsil;
            volum = c1;
            
            ElemSS(1:6) = ElemSS(1:6) + c1*[1; 1; 1; 1/2; 1/2; 1/2].*epsil;
            ElemSS(7:12) = ElemSS(7:12) + c1*stres;
            ElemSS(13) = ElemSS(13) + volum;

        end %je

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
        
        lint = 3;
        ll = 0; % Counter for history variables

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

                ll = ll + 1;

                %Integration point, weight, jacobian
                [Wgt,litr,lits] =  intpntt(l,lint,0);
                [shl,shld,shls,be] = shlt(litr,lits,3,3,0,0);
                [shg, shgs, Jdet, be, xs] = shgt(xit,3,shld,shls,3,0,0,be);
            
                %Physical location of int pt
                xint = xit(1,:)*shl;
                yint = xit(2,:)*shl;
                zint = xit(3,:)*shl;

                xi = POU_Coord3(xint,yint,zint,xlL,1,nelL);
                rL = xi(1);
                sL = xi(2);
                tL = xi(3);

                % Evaluate  basis functions at integration points
                if nelL == 4 || nelL == 10
                  [shlL,shldL,shls,be] = shltt([rL sL tL],nelL,nel2L,0,0);
                  PxyL = shgtt(xlL(:,1:nelL),nelL,shldL,shls,nen,0,0,be);
                else
                  [shlL,shldL,shls,be] = shlb([rL sL tL],nelL,nel2L,0,0);
                  PxyL = shgb(xlL(:,1:nelL),nelL,shldL,shls,nen,0,0,be);
                end
                QxyL = PxyL;

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
                nLx = tu3(1);
                nLy = tu3(2);
                nLz = tu3(3);
                nRx = -tu3(1);
                nRy = -tu3(2);
                nRz = -tu3(3);
                tLx = tu1(1);
                tLy = tu1(2);
                tLz = tu1(3);
                nvect = [nLx 0 0 nLy 0 nLz
                         0 nLy 0 nLx nLz 0
                         0 0 nLz 0 nLy nLx]; %- ?

                c1 = Wgt*tm3;
                
                for i = 1:nelL
                    NmatL(1,(i-1)*3+1) = shlL(i);
                    NmatL(2,(i-1)*3+2) = shlL(i);
                    NmatL(3,3*i      ) = shlL(i);
                    BmatL(Bcol1,(i-1)*3+1) = QxyL(i,col1);
                    BmatL(Bcol2,(i-1)*3+2) = QxyL(i,col2);
                    BmatL(Bcol3,3*i      ) = QxyL(i,col3);
                end

                for i = 1:nelR
                    NmatR(1,(i-1)*3+1) = shlR(i);
                    NmatR(2,(i-1)*3+2) = shlR(i);
                    NmatR(3,3*i      ) = shlR(i);
                    BmatR(Bcol1,(i-1)*3+1) = QxyR(i,col1);
                    BmatR(Bcol2,(i-1)*3+2) = QxyR(i,col2);
                    BmatR(Bcol3,3*i      ) = QxyR(i,col3);
                end
                
%                 dam = 1;
                nvec = [nLx; nLy; nLz];
                jumpu = NmatL*reshape(ulL,ndf*nelL,1) - NmatR*reshape(ulR,ndf*nelR,1);
                epsili = -1/2*(nvec*jumpu' + jumpu*nvec');
                epsil = [epsili(1,1); epsili(2,2); epsili(3,3); epsili(1,2); epsili(2,3); epsili(3,1)];
                
                ElemSS(1:6) = ElemSS(1:6) + c1*epsil;
            
            end %lint
        
        end %intt
        
    case 61 % form data structure for interface segments

        % Set the number of interface quantities per node to be stored
        if ~exist('numIQ','var')
            numIQ = 9; % 3 disp-jump, 3 traction, 3 numerical-flux
        else
            numIQ = max(numIQ,9);
        end
        
        segment = segment + 1;
        IElemSeg(1,inter) = segment;
        IElemSeg(2,inter) = segment;
        
        % Always use face of left element for the nodes
        % Triangular or Quadrilateral face
        if nelL == 4 || nelL == 10
            Iix(segment,1:3) = (segnode+1:segnode+3);
            ICoordinates(segnode+1:segnode+3,1:ndm) = xlL(1:ndm,1:3)';
            segnode = segnode + 3;
        else
            Iix(segment,1:4) = (segnode+1:segnode+4);
            ICoordinates(segnode+1:segnode+4,1:ndm) = xlL(1:ndm,1:4)';
            segnode = segnode + 4;
        end
        Iix(segment,nenseg+1) = ma; % good habit to copy material ID too
        
    case 60 % output interface quantities for plotting
        
        % get segment number for DG element (usually segment=inter)
        segment = IElemSeg(1,inter);
    
        % get nodes on the interface segment
        ElemFlagI = Iix(segment,1:nenseg);
        xlI = ICoordinates(ElemFlagI,1:ndm)';
        
        lam = getlam(matepropL);
        
        % compute quantities of interest
        % jump in displacement
        if nelL == 4 || nelL == 10
            nelseg = 3;
            plist = [0 1 0
                     0 0 1];
        else % nelL == 8 || nelL == 27
            nelseg = 4;
            plist = [-1 1 1 -1
                     -1 -1 1 1];
        end
        
        jumps = zeros(3,nelseg);
        for l = 1:nelseg
            
            r = plist(1,l);
            s = plist(2,l);
            if nelseg == 3
              shlS = shlt(r,s,nelseg,nelseg,0,0);
            else
              shlS = shlq(r,s,nelseg,nelseg,0,0);
            end

            %Physical location of int pt
            xint = xlI(1,:)*shlS;
            yint = xlI(2,:)*shlS;
            zint = xlI(3,:)*shlS;

            xi = POU_Coord3(xint,yint,zint,xlL,1,nelL);
            rL = xi(1);
            sL = xi(2);
            tL = xi(3);

            % Evaluate  basis functions at integration points
            if nelL == 4 || nelL == 10
              shlL = shltt([rL sL tL],nelL,nelL,0,0);
            else
              shlL = shlb([rL sL tL],nelL,nelL,0,0);
            end

            xi = POU_Coord3(xint,yint,zint,xlR,1,nelR);   %Get the kesi eta in the right hand side
            rR = xi(1);
            sR = xi(2);
            tR = xi(3);

            % Evaluate  basis functions at integration points
            if nelR == 4 || nelR == 10
              shlR = shltt([rR sR tR],nelR,nelR,0,0);
            else
              shlR = shlb([rR sR tR],nelR,nelR,0,0);
            end
            
            jumps(1:3,l) = ulL*shlL - ulR*shlR;
            
        end
        
        
        % interface traction, extrapolated from integration points
        % pull out weighting tensors and penalty parameter for this segment
        gamL = [hr(nh3) hr(nh3+1) hr(nh3+2)
                hr(nh3+3) hr(nh3+4) hr(nh3+5)
                hr(nh3+6) hr(nh3+7) hr(nh3+8)];
        gamR = [hr(nh3+9) hr(nh3+10) hr(nh3+11)
                hr(nh3+12) hr(nh3+13) hr(nh3+14)
                hr(nh3+15) hr(nh3+16) hr(nh3+17)];
        ep =   [hr(nh3+18) hr(nh3+19) hr(nh3+20)
                hr(nh3+21) hr(nh3+22) hr(nh3+23)
                hr(nh3+24) hr(nh3+25) hr(nh3+26)];
            
        tractions_int = zeros(6,nelseg); % tractions at integration points
        tractions = zeros(6,nelseg); % tractions at segment nodes
        for l = 1:nelseg
            
            % Evaluate  basis functions at integration points
            if nelL == 4 || nelL == 10
              [Wgt,ss] = int3d_t(l,nelseg,4);
              [shlL,shldL,shls,be] = shltt(ss,nelL,nelL,0,0);                
              [QxyL,shgsL,JdetL,bubbleL,xsL] = shgtt(xlL(:,1:nelL)+ulL(:,1:nelL),nelL,shldL,shls,nen,0,0,be);
            else
              [Wgt,ss] = intpntb(l,nelseg,5);
              [shlL,shldL,shls,be] = shlb(ss,nelL,nelL,0,0); 
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
              [shlR,shldR,shls,be] = shltt([rR sR tR],nelR,nelR,0,0);                  
              [QxyR,shgsR,JdetR,bubbleR,xsR] = shgtt(xlR(:,1:nelR)+ulR(:,1:nelR),nelR,shldR,shls,nen,0,0,be);

            else
              [shlR,shldR,shls,be] = shlb([rR sR tR],nelR,nelR,0,0);
              [QxyR,shgsR,JdetR,bubbleR,xsR] = shgb(xlR(:,1:nelR)+ulR(:,1:nelR),nelR,shldR,shls,nen,0,0,be);
            end
            
            % deformation gradient           
            [fiL,JxXL,FL] = kine3d(QxyL,-ulL,nelL,0); %this is equivalent to ikine2d 
            JxXL = 1/JxXL; %this is equivalent to ikine2d 
            
            [fiR,JxXR,FR] = kine3d(QxyR,-ulR,nelR,0); %this is equivalent to ikine2d
            JxXR = 1/JxXR; %this is equivalent to ikine2d
                      
            % Cauchy stresses
            sigma2L = SigmaCmat3i(FL,JxXL,matepropL,lam);
            
            SmatL1=[sigma2L(1), sigma2L(4), sigma2L(6)
                    sigma2L(4), sigma2L(2), sigma2L(5)
                    sigma2L(6), sigma2L(5), sigma2L(3)];
            
            [sigma2R, cmatR] = SigmaCmat3i(FR,JxXR,matepropR,lam);
            
            SmatR1=[sigma2R(1), sigma2R(4), sigma2R(6)
                    sigma2R(4), sigma2R(2), sigma2R(5)
                    sigma2R(6), sigma2R(5), sigma2R(3)];
            
            % unit vectors
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
            %normal vectors
            nLx = -tu3L(1);
            nLy = -tu3L(2);
            nLz = -tu3L(3);
            nRx = -tu3R(1);
            nRy = -tu3R(2);
            nRz = -tu3R(3);
            nvecL = [nLx; nLy; nLz];
            nvecR = [nRx; nRy; nRz];
            
            tvtr = gamL*SmatL1*nvecL-gamR*SmatR1*nvecR; % weighted traction
            jumpu = ulL*shlL - ulR*shlR; % displacement jump - points from Left side to Right side when uL > uR
            
            % numerical flux - traction vector pointing out from Right side
            tractions_int(1:6,l) = [-tvtr; -tvtr+ep*jumpu];
            
        end
        
        % reorder integration point values from tensor-product to FEM-
        % counter-clockwise so that normal shape functions can be used
        if nelL == 4 || nelL == 10
        tractions_int = tractions_int(:,[1 2 3]); % VERIFY THIS
        else
        tractions_int = tractions_int(:,[1 2 4 3]);
        end
        
        % extrapolate values to segment nodes
        if nelseg == 3
            plist = [-1/3 5/3 -1/3
                     -1/3 -1/3 5/3];
        else % nelseg == 4
            sqr3 = sqrt(3);
            plist = [-sqr3 sqr3 sqr3 -sqr3
                     -sqr3 -sqr3 sqr3 sqr3];
        end
        
        for ll = 1:nelseg
            
            r = plist(1,ll);
            s = plist(2,ll);
            shpS = sshp2d(r,s,nelseg);
            
            tractions(:,ll) = tractions_int*shpS;
            
        end
        
        
        % assemble into interface output quantities
        InterQuant(1:9,ElemFlagI,step) = [jumps; tractions];
        
end %Task Switch
