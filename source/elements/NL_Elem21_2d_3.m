%% DG implementation of large deformaiton 
% pure displacement 2D case
% Tim Truster
% 10/2011
% modified by Pinlei Chen
% 04/23/2013
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

if ~exist('modifDG','var')
modifDG = 0;1;2;
end
% modifDG = 0 for full, 1 for penalty only, 2 for incomplete interior
% penalty (IIGP), 3 for using initial Cijkl for nonsymm interior penalty, 4
% for using Cijkl_n-1 for nonsymm interior penalty

nitvms = 1;
if nitvms == 1 %VMS parameter for the stability tensor rp
pencoeff = 3;1;4;20;1;
elseif nitvms == 2 %Nitsche
pencoeff = 1;4;2;
else %RFB
pencoeff = 1;4;2;
end


switch isw %Task Switch
%%
    case 1
        
        if ndf > 2
            
            for i = 3:ndf
                lie(i,1) = 0;
            end
            
        end
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
        

        ElemEL = matepropL(1);
        ElemvL = matepropL(2);
        ElemER = matepropR(1);
        ElemvR = matepropR(2);
        lam = getlam(matepropL);
       
        NmatL = zeros(2,nstL);
        BmatL = zeros(4,nstL);
       
    
        NmatR = zeros(2,nstR);
        BmatR = zeros(4,nstR);
 
       
        dmatL = zeros(9,3);
        dmatR = zeros(9,3);
       

        % Determin bounds of integration segment
        InterConn2D2 % InterConn2DT %
        
        % Set jacobian for integration space
        drdr = (eL2 - eL1)/drL;
        drdrR = (eR2 - eR1)/drR; % This causes errors for partial interface segments
        m = (eR2-eR1)/(eL1-eL2);
        if nelL == 3 || nelL == 6  
        lint = 3;10;2;3;
        else
        lint = 10;3;4;10;2;3; %10 for body force problem; 4 for other problem 
        end
        ideriv = 0;

%Form tau matrix: compute on first iteration and store; this gives
        %quadratic convergence when tau isn't recalculated (no
        %linearization of tau required).
   iterset = 1;
   if iter  <=iterset % == 0 %
       if exist('dropp','var') && dropp == 1
       lamt = 0;
       else
       lamt = lam;
       end
        [tauL,intb] = TauS2(xlintL,xlL,ulL,matepropL,nelL,nen,lamt); %[Y^(-1)]
%        [tauL,intb] = TauS2_1(xlL,ulL,mateprop,nelL,nel2L,nen,lam,roL,eL1,drdr); %[Y^(-1)]        
%       TauListL(elem,:) = [tauL(1),tauL(2),tauL(3),tauL(4),tauL(5),tauL(6),tauL(7),tauL(8),tauL(9)];
        [tauR,intb] = TauS2(xlintR,xlR,ulR,matepropR,nelR,nen,lamt);
%        [tauR,intb] = TauS2_1(xlR,ulR,mateprop,nelR,nel2R,nen,lam,roL,eL1,drdrR);
%          TauListR(elem,:) = [tauR(1),tauR(2),tauR(3),tauR(4),tauR(5),tauR(6),tauR(7),tauR(8),tauR(9)];

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
%         ep = 20*eye(2);
   
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
           [PxyL,shgsL,JdetL,bubbleL,XsL] = shgt(xlL(:,1:nelL),nelL,shldL,shlsL,nen,0,0,be); 
           [PxyL,shgsL,JdetL,bubbleL,xsL] = shgt(xlL(:,1:nelL)+ulL(:,1:nelL),nelL,shldL,shlsL,nen,0,0,be);
           else
           [shlL,shldL,shlsL,be] = shlq(rL,lits,nelL,nelL,0,0);
           [PxyL,shgsL,JdetL,bubbleL,XsL] = shgq(xlL(:,1:nelL),nelL,shldL,shlsL,nen,0,0,be); 
           [PxyL,shgsL,JdetL,bubbleL,xsL] = shgq(xlL(:,1:nelL)+ulL(:,1:nelL),nelL,shldL,shlsL,nen,0,0,be);      
           end
            QxyL = PxyL;  
            
            
            rR = m*(rL-eL2) + eR1;
            
            if nelR == 3 || nelR == 6
                sR = 0;
            else %if nelR == 4
                sR = -1;
            end
           if nelR == 3 || nelR == 6            
            [shlR,shldR,shlsR,be] = shlt(rR,sR,nelR,nelR,0,0);
            [PxyR,shgsR,JdetR,bubbleR,XsR] = shgt(xlR(:,1:nelR),nelR,shldR,shlsR,nen,0,0,be);     
            [PxyR,shgsR,JdetR,bubbleR,xsR] = shgt(xlR(:,1:nelR)+ulR(:,1:nelR),nelR,shldR,shlsR,nen,0,0,be);  
           else
            [shlR,shldR,shlsR,be] = shlq(rR,sR,nelR,nelR,0,0);
            [PxyR,shgsR,JdetR,bubbleR,XsR] = shgq(xlR(:,1:nelR),nelR,shldR,shlsR,nen,0,0,be);     
            [PxyR,shgsR,JdetR,bubbleR,xsR] = shgq(xlR(:,1:nelR)+ulR(:,1:nelR),nelR,shldR,shlsR,nen,0,0,be);  
           end
            QxyR = PxyR;

 % initial stress sigma and material bulk term c_ijkl
           for mm = 1:nelL    
           NmatL(:,2*mm-1:2*mm) = [shlL(mm,1)     0          
                                    0        shlL(mm,1)]; 
                          
           BmatL(:,2*mm-1:2*mm) = [QxyL(mm,1) 0             
                                    0         QxyL(mm,2)         
                                    QxyL(mm,2) QxyL(mm,1)        
                                     QxyL(mm,2) -QxyL(mm,1)];
%            BmatL1(:,2*mm-1:2*mm) = [QxyL(mm,1) 0                  
%                                      0         QxyL(mm,2)          
%                                      QxyL(mm,2) QxyL(mm,1)];
%            BmatL2(:,2*mm-1:2*mm)=[QxyL(mm,1)  0                     
%                                  QxyL(mm,2)  0                                      
%                                    0     QxyL(mm,1)                
%                                    0     QxyL(mm,2) ];   
%            BmatL3(:,2*mm-1:2*mm)=[QxyL(mm,1)  0                     
%                                    0     QxyL(mm,1)                                      
%                                    QxyL(mm,2) 0                     
%                                    0     QxyL(mm,2)];     
           end 
            
           for mm = 1:nelR    
           NmatR(:,2*mm-1:2*mm) = [shlR(mm,1)     0          
                                     0        shlR(mm,1)]; 
                          
           BmatR(:,2*mm-1:2*mm) = [QxyR(mm,1) 0             
                                    0         QxyR(mm,2)         
                                    QxyR(mm,2) QxyR(mm,1)        
                                     QxyR(mm,2) -QxyR(mm,1)];
%            BmatR1(:,2*mm-1:2*mm) = [QxyR(mm,1) 0                  
%                                      0         QxyR(mm,2)          
%                                      QxyR(mm,2) QxyR(mm,1)];
%            BmatR2(:,2*mm-1:2*mm)=[QxyR(mm,1)  0                     
%                                  QxyR(mm,2)  0                                      
%                                    0     QxyR(mm,1)                
%                                    0     QxyR(mm,2)];   
%            BmatR3(:,2*mm-1:2*mm)=[QxyR(mm,1)  0                     
%                                    0     QxyR(mm,1)                                      
%                                    QxyR(mm,2) 0                     
%                                    0     QxyR(mm,2)];                       
           end  
    %plug in the -ul to get the inverse of deformation gradient fi as 'F' in this function             
            [fiL,JxXL,FL] = kine2d(QxyL,-ulL,nelL,0); %this is equivalent to ikine2d
            JxXL = 1/JxXL; %this is equivalent to ikine2d
    %         [F,JxX,fi] = ikine2d(Qxy,ul,nel); %this is equivalent to above two lines
            JdetL = JdetL/JxXL;            

            
            [fiR,JxXR,FR] = kine2d(QxyR,-ulR,nelR,0); %this is equivalent to ikine2d
            JxXR = 1/JxXR; %this is equivalent to ikine2d
    %         [F,JxX,fi] = ikine2d(Qxy,ul,nel); %this is equivalent to above two lines
            JdetR = JdetR/JxXR;               
          
            
          
            [sigma2L, cmatL] = SigmaCmat2i(FL,JxXL,matepropL,lam);
            
%             SmatL = [sigma2L(1) 0  sigma2L(3)/2 sigma2L(3)/2
%                     0 sigma2L(2)  sigma2L(3)/2 -sigma2L(3)/2
%                     sigma2L(3)/2  sigma2L(3)/2 (sigma2L(2)+sigma2L(1))/4 (sigma2L(2)-sigma2L(1))/4
%                     sigma2L(3)/2 -sigma2L(3)/2 (sigma2L(2)-sigma2L(1))/4 (sigma2L(2)+sigma2L(1))/4];            
            
            SmatL1=[sigma2L(1), sigma2L(3)
                    sigma2L(3), sigma2L(2)];           
            
            [sigma2R, cmatR] = SigmaCmat2i(FR,JxXR,matepropR,lam);
                                 
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
%additional stiffness term


           
           
           SmatnL=[SmatL1*nvecL zeros(ndm,1)
                    zeros(ndm,1) SmatL1*nvecL];
            
            SmatnR=[SmatR1*nvecR zeros(ndm,1)
                    zeros(ndm,1) SmatR1*nvecR];
           cmatnL=(nvectL1*cmatL(1:3,1:3));
           cmatnR=(nvectR1*cmatR(1:3,1:3));
           cmatnBL=BmatL'*P2'*[cmatnL zeros(2,3)
                            zeros(2,3)  cmatnL ];               
           cmatnBR=BmatR'*P2'*[cmatnR    zeros(2,3)
                            zeros(2,3)  cmatnR ];

            term17L=P2'*SmatnL*gamL;
            term17R=P2'*SmatnR*gamR;

            term18L=P1'*(gamL*nvectL1*cmatL(1:3,1:3))';
            term18R=P1'*(gamR*nvectR1*cmatR(1:3,1:3))';

            term28L=NmatL'*(c1L*gamL*SmatL1*nvecL-c1R*gamR*SmatR1*nvecR);   %average stress term
            term28R=NmatR'*(c1L*gamL*SmatL1*nvecL-c1R*gamR*SmatR1*nvecR);

             rhspulL = reshape(ulL,ndf*nelL,1);
             rhspulR = reshape(ulR,ndf*nelR,1);
                    
            jumpu = NmatL*rhspulL - NmatR*rhspulR;  %jumpu=uL-uR


%%          assume dealta1=delata2=0.5

%              jumpu1=[jumpu(1) 0    jumpu(2)   
%                         0  jumpu(2)  jumpu(1)];  
             gamajumpuL = gamL*jumpu;
             gamajumpuR = gamR*jumpu;             
             term5L=cmatnBL*[eye(3,3)*gamajumpuL(1) eye(3,3)*gamajumpuL(2)]'*(P1*BmatL);
             term5R=cmatnBR*[eye(3,3)*gamajumpuR(1) eye(3,3)*gamajumpuR(2)]'*(P1*BmatR);
            
             sig8L2 = [gamajumpuL(1) gamajumpuL(2)]*cmatnL;

%             sig8L3= [sig8L2(1) 0  sig8L2(3)/2 sig8L2(3)/2
%                     0 sig8L2(2)  sig8L2(3)/2 -sig8L2(3)/2
%                     sig8L2(3)/2  sig8L2(3)/2 (sig8L2(2)+sig8L2(1))/4 (sig8L2(2)-sig8L2(1))/4
%                     sig8L2(3)/2 -sig8L2(3)/2 (sig8L2(2)-sig8L2(1))/4 (sig8L2(2)+sig8L2(1))/4];   

          sig8L3 = [sig8L2(1) 0 sig8L2(3) 0;
                    sig8L2(3) 0 sig8L2(2) 0;
                    0 sig8L2(1) 0 sig8L2(3);
                    0 sig8L2(3) 0 sig8L2(2)];
               
           term8L = BmatL'*P2'*sig8L3*(P3*BmatL);
                                   
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
               
              term8R =BmatR'*P2'*sig8R3*(P3*BmatR);
             
              term30L=NmatL'*ep*jumpu;
              term30R=NmatR'*ep*jumpu;

             [dmatL1]=dmat2_no_p(JxXL,matepropL,lam);
             [dmatR1]=dmat2_no_p(JxXR,matepropR,lam);        
             dmatL2 = [eye(3,3)*jumpu(1) eye(3,3)*jumpu(2)]*[eye(3,3)*gamL(1,1) eye(3,3)*gamL(1,2)
                                                             eye(3,3)*gamL(2,1) eye(3,3)*gamL(2,2)]*nvectL2*dmatL1/JxXL;  %missing J here
             dmatR2 = [eye(3,3)*jumpu(1) eye(3,3)*jumpu(2)]*[eye(3,3)*gamR(1,1) eye(3,3)*gamR(1,2)
                                                             eye(3,3)*gamR(2,1) eye(3,3)*gamR(2,2)]*nvectR2*dmatR1/JxXR;  %missing J here
             
             term7L = BmatL'*P1'*dmatL2*(P1*BmatL);
             term7R = BmatR'*P1'*dmatR2*(P1*BmatR);

            if exist('modifDG','var') && modifDG > 0 % modify the DG terms used in the formulation
             
             % Penalty terms
                 
             ElemFL = ElemFL-(+C1L*term30L);
             ElemFR = ElemFR-(-C1R*term30R);
              
             ElemKLL = ElemKLL + C1L*(NmatL'*ep*NmatL);
             ElemKLR = ElemKLR - C1R*(NmatL'*ep*NmatR);
             ElemKRL = ElemKRL - C1L*(NmatR'*ep*NmatL);
             ElemKRR = ElemKRR + C1R*(NmatR'*ep*NmatR);  
             
             if modifDG > 1 
                 
             % average stress terms
                 
             ElemFL = ElemFL-(-term28L);
             ElemFR = ElemFR-(+term28R);
% 
             ElemKLL = ElemKLL  - c1L*NmatL'*(term17L'+term18L')*BmatL;
             ElemKLR = ElemKLR  + c1R*NmatL'*(term17R'+term18R')*BmatR;
             ElemKRL = ElemKRL  + c1L*NmatR'*(term17L'+term18L')*BmatL;
             ElemKRR = ElemKRR  - c1R*NmatR'*(term17R'+term18R')*BmatR;
             
             if modifDG == 3 
                 
             % constant nonsymmetric terms

             ElemFL = ElemFL-(-C1L*BmatLi'*(term17Li+term18Li)*jumpu);
             ElemFR = ElemFR-(C1R* BmatRi'*(term17Ri+term18Ri)*jumpu);

             ElemKLL = ElemKLL - C1L*BmatLi'*(term17Li+term18Li)*NmatL;
             ElemKLR = ElemKLR + C1L*BmatLi'*(term17Li+term18Li)*NmatR;
             ElemKRL = ElemKRL + C1R*BmatRi'*(term17Ri+term18Ri)*NmatL;
             ElemKRR = ElemKRR - C1R*BmatRi'*(term17Ri+term18Ri)*NmatR;
             
             elseif modifDG > 3 
                 
             % nonsymmetric terms

             ElemFL = ElemFL-(-c1L*BmatL'*(term17L+term18L)*jumpu);
             ElemFR = ElemFR-(c1R* BmatR'*(term17R+term18R)*jumpu);

             ElemKLL = ElemKLL - c1L*BmatL'*(term17L+term18L)*NmatL;
             ElemKLR = ElemKLR + c1L*BmatL'*(term17L+term18L)*NmatR;
             ElemKRL = ElemKRL + c1R*BmatR'*(term17R+term18R)*NmatL;
             ElemKRR = ElemKRR - c1R*BmatR'*(term17R+term18R)*NmatR;
             
             end
             end 

             else % full method            
             
             
                ElemFL = ElemFL-(-c1L*BmatL'*(term17L+term18L)*jumpu-term28L+C1L*term30L);
                ElemFR = ElemFR-(c1R* BmatR'*(term17R+term18R)*jumpu+term28R-C1R*term30R);

%               ElemFL = ElemFL-(-term28L+C1L*term30L);
%               ElemFR = ElemFR-(+term28R-C1R*term30R);

%              ElemFL = ElemFL-(-c1L*(term17L'+term18L')*jumpu+C1L*term30L);
%              ElemFR = ElemFR-(c1R*(term17R'+term18R')*jumpu-C1R*term30R);

                ElemKLL = ElemKLL - c1L*BmatL'*(term17L+term18L)*NmatL - c1L*NmatL'*(term17L'+term18L')*BmatL;
                ElemKLR = ElemKLR + c1L*BmatL'*(term17L+term18L)*NmatR + c1R*NmatL'*(term17R'+term18R')*BmatR;
                ElemKRL = ElemKRL + c1R*BmatR'*(term17R+term18R)*NmatL + c1L*NmatR'*(term17L'+term18L')*BmatL;
                ElemKRR = ElemKRR - c1R*BmatR'*(term17R+term18R)*NmatR - c1R*NmatR'*(term17R'+term18R')*BmatR;

%              ElemKLL = ElemKLL - c1L*(term17L'+term18L')*NmatL;
%              ElemKLR = ElemKLR + c1L*(term17L'+term18L')*NmatR;
%              ElemKRL = ElemKRL + c1R*(term17R'+term18R')*NmatL;
%              ElemKRR = ElemKRR - c1R*(term17R'+term18R')*NmatR;
% 
%              ElemKLL = ElemKLL  - c1L*NmatL'*(term17L+term18L);
%              ElemKLR = ElemKLR  + c1R*NmatL'*(term17R+term18R);
%              ElemKRL = ElemKRL  + c1L*NmatR'*(term17L+term18L);
%              ElemKRR = ElemKRR  - c1R*NmatR'*(term17R+term18R);
 
                ElemKLL = ElemKLL + C1L*(NmatL'*ep*NmatL);
                ElemKLR = ElemKLR - C1R*(NmatL'*ep*NmatR);
                ElemKRL = ElemKRL - C1L*(NmatR'*ep*NmatL);
                ElemKRR = ElemKRR + C1R*(NmatR'*ep*NmatR);              

            %relate to jumpu terms additional terms
                ElemKLL = ElemKLL - c1L*(term5L+term5L'+term7L+term8L);  %                     
                ElemKRR = ElemKRR + c1R*(term5R+term5R'+term7R+term8R);  %
           end
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
        
%          beta = 1;
%          beta2 = beta^2;
%          sigmax = 100;
%          dc = 0.2;
% %        beta = 0.707;
% %        beta2 = beta^2;
%         sigmaxc = 100; %0.01e-3;
% %        dc = 40;
%         
%         Hc = sigmax/dc;
%         rp = 4000;1*Hc;

        ElemEL = matepropL(1);
        ElemvL = matepropL(2);
        ElemER = matepropR(1);
        ElemvR = matepropR(2);
        lam = getlam(matepropL);
       
        NmatL = zeros(2,nstL);
        BmatL = zeros(4,nstL);
      
        NmatR = zeros(2,nstR);
        BmatR = zeros(4,nstR);
       
        dmatL = zeros(9,3);
        dmatR = zeros(9,3);
%         ib = 1; 
%         % Determine bounds of integration, right
%         
%         if nelR == 4 || nelR == 9
%             
%             t1 = [(xlR(1,2)-xlR(1,1)); (xlR(2,2)-xlR(2,1)); 0];
%             t2 = [(xlR(1,4)-xlR(1,1)); (xlR(2,4)-xlR(2,1)); 0];
%             t3 = VecCrossProd(t1,t2);
%             aR = VecNormalize(t3);
%             hR = sqrt((xlR(1,2)-xlR(1,1))^2+(xlR(2,2)-xlR(2,1))^2)/aR;
%             
%             drR = 2;
%             roR = -1;
% 
%             % Upper Limit
%             if nodeAR == ElemFlagR(2)
%                 eR2 = 1;
%                 xlintR(:,2) = xlR(:,2);
%                 xlintR(:,3) = xlR(:,3);
%             elseif nodeAR == ElemFlagR(nel2R)
%                 eR2 = epR;
%             elseif nodeAR == -1 %no enrichment but DG instead
%                 xy = xlL(:,1);
%                 POUxi = POU_Coord(xy(1),xy(2),xlR,1,nelR);
%                 eR2 = POUxi(1);
%                 xlintR(:,2) = xy;
%             elseif nelR == 9 && nodeAR == ElemFlagR(5)
%                 eR2 = 0;
%                 xlintR(:,2) = xlR(:,5);
%             end
%             % Lower Limit
%             if nodeBR == ElemFlagR(1)
%                 eR1 = -1;
%                 xlintR(:,1) = xlR(:,1);
%                 xlintR(:,3) = xlR(:,4);
%             elseif nodeBR == ElemFlagR(nel2R)
%                 eR1 = epR;
%             elseif nodeBR == -1 %no enrichment but DG instead
%                 xy = xlL(:,2);
%                 POUxi = POU_Coord(xy(1),xy(2),xlR,1,nelR);
%                 eR1 = POUxi(1);
%                 xlintR(:,1) = xy;
%                 xlintR(:,3) = xlR(:,3);
%             elseif nelR == 9 && nodeBR == ElemFlagR(5)
%                 eR1 = 0;
%                 xlintR(:,1) = xlR(:,5);
%                 xlintR(:,3) = xlR(:,3);
%             end
%         
%         elseif nelR == 3 || nelR == 6
%             
%             t1 = [(xlR(1,2)-xlR(1,1)); (xlR(2,2)-xlR(2,1)); 0];
%             t2 = [(xlR(1,3)-xlR(1,1)); (xlR(2,3)-xlR(2,1)); 0];
%             t3 = VecCrossProd(t1,t2);
%             aR = VecNormalize(t3);
%             aR = aR/2;
%             hR = sqrt((xlR(1,2)-xlR(1,1))^2+(xlR(2,2)-xlR(2,1))^2)/aR;
%             
%             drR = 1;
%             roR = 0;
% 
%             % Upper Limit
%             if nodeAR == ElemFlagR(2)
%                 eR2 = 1;
%                 xlintR(:,2) = xlR(:,2);
%             elseif nodeAR == ElemFlagR(nel2R)
%                 eR2 = epR;
%             elseif nodeAR == -1 %no enrichment but DG instead
%                 xy = xlL(:,1);
%                 POUxi = POU_Coord(xy(1),xy(2),xlR,1,nelR);
%                 eR2 = POUxi(1);
%                 xlintR(:,2) = xy;
%             elseif nelR == 6 && nodeAR == ElemFlagR(4)
%                 eR2 = 1/2;
%                 xlintR(:,2) = xlR(:,4);
%             end
%             % Lower Limit
%             if nodeBR == ElemFlagR(1)
%                 eR1 = 0;
%                 xlintR(:,1) = xlR(:,1);
%             elseif nodeBR == ElemFlagR(nel2R)
%                 eR1 = epR;
%             elseif nodeBR == -1 %no enrichment but DG instead
%                 xy = xlL(:,2);
%                 POUxi = POU_Coord(xy(1),xy(2),xlR,1,nelR);
%                 eR1 = POUxi(1);
%                 xlintR(:,1) = xy;
%             elseif nelR == 6 && nodeBR == ElemFlagR(4)
%                 eR1 = 1/2;
%                 xlintR(:,1) = xlR(:,4);
%             end
%             xlintR(:,3) = xlR(:,3);
%         
%         end
%         
%         % Determine bounds of integration, left
%         
%         if nelL == 4 || nelL == 9
%             
%             t1 = [(xlL(1,2)-xlL(1,1)); (xlL(2,2)-xlL(2,1)); 0];
%             t2 = [(xlL(1,4)-xlL(1,1)); (xlL(2,4)-xlL(2,1)); 0];
%             t3 = VecCrossProd(t1,t2);
%             aL = VecNormalize(t3);
%             hL = sqrt((xlL(1,2)-xlL(1,1))^2+(xlL(2,2)-xlL(2,1))^2)/aL;
%             
%             drL = 2;
%             roL = -1;
% 
%             % Upper Limit
%             if nodeAL == ElemFlagL(1)
%                 eL1 = -1;
%                 xlintL(:,1) = xlL(:,1);
%                 xlintL(:,3) = xlL(:,4);
%             elseif nodeAL == ElemFlagL(nel2L)
%                 eL1 = epL;
%             elseif nodeAL == -1 %no enrichment but DG instead
%                 xy = xlR(:,2);
%                 POUxi = POU_Coord(xy(1),xy(2),xlL,1,nelL);
%                 eL1 = POUxi(1);
%                 xlintL(:,1) = xlL(:,1);
%             elseif nelL == 9 && nodeAL == ElemFlagL(5)
%                 eL1 = 0;
%                 xlintL(:,1) = xlL(:,5);
%             end
%             % Lower Limit
%             if nodeBL == ElemFlagL(2)
%                 eL2 = 1;
%                 xlintL(:,2) = xlL(:,2);
%                 xlintL(:,3) = xlL(:,3);
%             elseif nodeBL == ElemFlagL(nel2L)
%                 eL2 = epL;
%             elseif nodeBL == -1 %no enrichment but DG instead
%                 xy = xlR(:,1);
%                 POUxi = POU_Coord(xy(1),xy(2),xlL,1,nelL);
%                 eL2 = POUxi(1);
%                 xlintL(:,2) = xlL(:,2);
%                 xlintL(:,3) = xlL(:,4);
%             elseif nelL == 9 && nodeB == ElemFlagL(5)
%                 eL2 = 0;
%                 xlintL(:,2) = xlL(:,5);
%                 xlintL(:,3) = xlL(:,4);
%             end
%         
%         elseif nelL == 3 || nelL == 6
%             
%             t1 = [(xlL(1,2)-xlL(1,1)); (xlL(2,2)-xlL(2,1)); 0];
%             t2 = [(xlL(1,3)-xlL(1,1)); (xlL(2,3)-xlL(2,1)); 0];
%             t3 = VecCrossProd(t1,t2);
%             aL = VecNormalize(t3);
%             aL = aL/2;
%             hL = sqrt((xlL(1,2)-xlL(1,1))^2+(xlL(2,2)-xlL(2,1))^2)/aL;
%             
%             drL = 1;
%             roL = 0;
% 
%             % Upper Limit
%             if nodeAL == ElemFlagL(1)
%                 eL1 = 0;
%                 xlintL(:,1) = xlL(:,1);                
%             elseif nodeAL == ElemFlagL(nel2L)
%                 eL1 = epL;
%             elseif nodeAL == -1 %no enrichment but DG instead
%                 xy = xlR(:,2);
%                 POUxi = POU_Coord(xy(1),xy(2),xlL,1,nelL);
%                 eL1 = POUxi(1);
%                 xlintL(:,1) = xy;  
%             elseif nelL == 6 && nodeAL == ElemFlagL(4)
%                 eL1 = 1/2;
%                 xlintL(:,1) = xlL(:,4);
%             end
%             % Lower Limit
%             if nodeBL == ElemFlagL(2)
%                 eL2 = 1;
%                 xlintL(:,2) = xlL(:,2);  
%             elseif nodeBL == ElemFlagL(nel2L)
%                 eL2 = epL;
%             elseif nodeBL == -1 %no enrichment but DG instead
%                 xy = xlR(:,1);
%                 POUxi = POU_Coord(xy(1),xy(2),xlL,1,nelL);
%                 eL2 = POUxi(1);
%                 xlintL(:,2) = xy;  
%             elseif nelL == 6 && nodeBL == ElemFlagL(4)
%                 eL2 = 1/2;
%                 xlintL(:,1) = xlL(:,4);
%             end
%             xlintL(:,3) = xlL(:,3);  
%         
%         end
        
%         h = 2/(hR + hL);
% %         h = 0.125;        

        % Determin bounds of integration segment
        InterConn2D2 % InterConn2DT % 
        
        % Set jacobian for integration space
        drdr = (eL2 - eL1)/drL;
        drdrR = (eR2 - eR1)/drR;
        m = (eR2-eR1)/(eL1-eL2);
        if nelL == 3 || nelL == 6  
        lint = 3;10;2;3;
        else
        lint = 10;4;10;2;3; %10 for body force problem; 4 for other problem 
        end
        ideriv = 0;
     
%         etauL = TauEE2d(xlintL,DmatL,lintt6);
%         etauR = TauEE2d(xlintR,DmatR,lintt6);
% %         tau = tauL;
%         etau = etauL + etauR;
% %         ep = tauR/tauL;
%         ep = 10*max(muL,muR)/h;
%         eb = 0;
%         Kinv = 2*eye(2);

%Form tau matrix: compute on first iteration and store; this gives
        %quadratic convergence when tau isn't recalculated (no
        %linearization of tau required).
   iterset = 4;
   if iter  <=iterset % == 0 %
%         [tauL,intb] = TauS2(xlL,ulL,matepropL,nelL,nen,lam,roL,eL1,drdr); %[Y^(-1)]
% %        [tauL,intb] = TauS2_1(xlL,ulL,mateprop,nelL,nel2L,nen,lam,roL,eL1,drdr); %[Y^(-1)]        
% %       TauListL(elem,:) = [tauL(1),tauL(2),tauL(3),tauL(4),tauL(5),tauL(6),tauL(7),tauL(8),tauL(9)];
%         [tauR,intb] = TauS2(xlR,ulR,matepropR,nelR,nen,lam,roL,eL1,drdrR);
%        [tauR,intb] = TauS2_1(xlR,ulR,mateprop,nelR,nel2R,nen,lam,roL,eL1,drdrR);
%          TauListR(elem,:) = [tauR(1),tauR(2),tauR(3),tauR(4),tauR(5),tauR(6),tauR(7),tauR(8),tauR(9)];
        [tauL,intb] = TauS2(xlintL,xlL,ulL,matepropL,nelL,nen,lamt); %[Y^(-1)]
%        [tauL,intb] = TauS2_1(xlL,ulL,mateprop,nelL,nel2L,nen,lam,roL,eL1,drdr); %[Y^(-1)]        
%       TauListL(elem,:) = [tauL(1),tauL(2),tauL(3),tauL(4),tauL(5),tauL(6),tauL(7),tauL(8),tauL(9)];
        [tauR,intb] = TauS2(xlintR,xlR,ulR,matepropR,nelR,nen,lamt);
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
%        ep = 100*eye(2);
   
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
           [PxyL,shgsL,JdetL,bubbleL,XsL] = shgt(xlL(:,1:nelL),nelL,shldL,shlsL,nen,0,0,be); 
           [PxyL,shgsL,JdetL,bubbleL,xsL] = shgt(xlL(:,1:nelL)+ulL(:,1:nelL),nelL,shldL,shlsL,nen,0,0,be);
           else
           [shlL,shldL,shlsL,be] = shlq(rL,lits,nelL,nelL,0,0);
           [PxyL,shgsL,JdetL,bubbleL,XsL] = shgq(xlL(:,1:nelL),nelL,shldL,shlsL,nen,0,0,be); 
           [PxyL,shgsL,JdetL,bubbleL,xsL] = shgq(xlL(:,1:nelL)+ulL(:,1:nelL),nelL,shldL,shlsL,nen,0,0,be);      
           end
            QxyL = PxyL;  
            
            
            rR = m*(rL-eL2) + eR1;
            
            if nelR == 3 || nelR == 6
                sR = 0;
            else %if nelR == 4
                sR = -1;
            end
           if nelR == 3 || nelR == 6            
            [shlR,shldR,shlsR,be] = shlt(rR,sR,nelR,nelR,0,0);
            [PxyR,shgsR,JdetR,bubbleR,XsR] = shgt(xlR(:,1:nelR),nelR,shldR,shlsR,nen,0,0,be);     
            [PxyR,shgsR,JdetR,bubbleR,xsR] = shgt(xlR(:,1:nelR)+ulR(:,1:nelR),nelR,shldR,shlsR,nen,0,0,be);  
           else
            [shlR,shldR,shlsR,be] = shlq(rR,sR,nelR,nelR,0,0);
            [PxyR,shgsR,JdetR,bubbleR,XsR] = shgq(xlR(:,1:nelR),nelR,shldR,shlsR,nen,0,0,be);     
            [PxyR,shgsR,JdetR,bubbleR,xsR] = shgq(xlR(:,1:nelR)+ulR(:,1:nelR),nelR,shldR,shlsR,nen,0,0,be);  
           end
            QxyR = PxyR;

 % initial stress sigma and material bulk term c_ijkl
           for mm = 1:nelL    
           NmatL(:,2*mm-1:2*mm) = [shlL(mm,1)     0          
                                    0        shlL(mm,1)]; 
                          
           BmatL(:,2*mm-1:2*mm) = [QxyL(mm,1) 0             
                                    0         QxyL(mm,2)         
                                    QxyL(mm,2) QxyL(mm,1)        
                                     QxyL(mm,2) -QxyL(mm,1)];
%            BmatL1(:,2*mm-1:2*mm) = [QxyL(mm,1) 0                  
%                                      0         QxyL(mm,2)          
%                                      QxyL(mm,2) QxyL(mm,1)];
%            BmatL2(:,2*mm-1:2*mm)=[QxyL(mm,1)  0                     
%                                  QxyL(mm,2)  0                                      
%                                    0     QxyL(mm,1)                
%                                    0     QxyL(mm,2) ];   
%            BmatL3(:,2*mm-1:2*mm)=[QxyL(mm,1)  0                     
%                                    0     QxyL(mm,1)                                      
%                                    QxyL(mm,2) 0                     
%                                    0     QxyL(mm,2)];     
           end 
            
           for mm = 1:nelR    
           NmatR(:,2*mm-1:2*mm) = [shlR(mm,1)     0          
                                     0        shlR(mm,1)]; 
                          
           BmatR(:,2*mm-1:2*mm) = [QxyR(mm,1) 0             
                                    0         QxyR(mm,2)         
                                    QxyR(mm,2) QxyR(mm,1)        
                                     QxyR(mm,2) -QxyR(mm,1)];
%            BmatR1(:,2*mm-1:2*mm) = [QxyR(mm,1) 0                  
%                                      0         QxyR(mm,2)          
%                                      QxyR(mm,2) QxyR(mm,1)];
%            BmatR2(:,2*mm-1:2*mm)=[QxyR(mm,1)  0                     
%                                  QxyR(mm,2)  0                                      
%                                    0     QxyR(mm,1)                
%                                    0     QxyR(mm,2)];   
%            BmatR3(:,2*mm-1:2*mm)=[QxyR(mm,1)  0                     
%                                    0     QxyR(mm,1)                                      
%                                    QxyR(mm,2) 0                     
%                                    0     QxyR(mm,2)];                       
           end  
    %plug in the -ul to get the inverse of deformation gradient fi as 'F' in this function             
            [fiL,JxXL,FL] = kine2d(QxyL,-ulL,nelL,0); %this is equivalent to ikine2d
            JxXL = 1/JxXL; %this is equivalent to ikine2d
    %         [F,JxX,fi] = ikine2d(Qxy,ul,nel); %this is equivalent to above two lines
            JdetL = JdetL/JxXL;            

            
            [fiR,JxXR,FR] = kine2d(QxyR,-ulR,nelR,0); %this is equivalent to ikine2d
            JxXR = 1/JxXR; %this is equivalent to ikine2d
    %         [F,JxX,fi] = ikine2d(Qxy,ul,nel); %this is equivalent to above two lines
            JdetR = JdetR/JxXR;               
          
            
          
            [sigma2L, cmatL] = SigmaCmat2i(FL,JxXL,matepropL,lam);
            
%             SmatL = [sigma2L(1) 0  sigma2L(3)/2 sigma2L(3)/2
%                     0 sigma2L(2)  sigma2L(3)/2 -sigma2L(3)/2
%                     sigma2L(3)/2  sigma2L(3)/2 (sigma2L(2)+sigma2L(1))/4 (sigma2L(2)-sigma2L(1))/4
%                     sigma2L(3)/2 -sigma2L(3)/2 (sigma2L(2)-sigma2L(1))/4 (sigma2L(2)+sigma2L(1))/4];            
            
            SmatL1=[sigma2L(1), sigma2L(3)
                    sigma2L(3), sigma2L(2)];           
            
            [sigma2R, cmatR] = SigmaCmat2i(FR,JxXR,matepropR,lam);
                                 
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
%additional stiffness term
          
           SmatnL=[SmatL1*nvecL zeros(ndm,1)
                    zeros(ndm,1) SmatL1*nvecL];
            
            SmatnR=[SmatR1*nvecR zeros(ndm,1)
                    zeros(ndm,1) SmatR1*nvecR];
           cmatnL=(nvectL1*cmatL(1:3,1:3));
           cmatnR=(nvectR1*cmatR(1:3,1:3));
           cmatnBL=(P2*BmatL)'*[cmatnL zeros(2,3)
                            zeros(2,3)  cmatnL ];               
           cmatnBR=(P2*BmatR)'*[cmatnR    zeros(2,3)
                            zeros(2,3)  cmatnR ];

            term17L=gamL*SmatnL'*(P2*BmatL);
            term17R=gamR*SmatnR'*(P2*BmatR);

            term18L=gamL*nvectL1*cmatL(1:3,1:3)*(P1*BmatL);
            term18R=gamR*nvectR1*cmatR(1:3,1:3)*(P1*BmatR);

            term28L=NmatL'*(c1L*gamL*SmatL1*nvecL-c1R*gamR*SmatR1*nvecR);   %average stress term
            term28R=NmatR'*(c1L*gamL*SmatL1*nvecL-c1R*gamR*SmatR1*nvecR);

             rhspulL = reshape(ulL,ndf*nelL,1);
             rhspulR = reshape(ulR,ndf*nelR,1);
                    
            jumpu = NmatL*rhspulL - NmatR*rhspulR;  %jumpu=uL-uR
             
              term30L=NmatL'*ep*jumpu;
              term30R=NmatR'*ep*jumpu;
             
              ElemFL = ElemFL-(-c1L*(term17L'+term18L')*jumpu-term28L+C1L*term30L);
              ElemFR = ElemFR-(c1R*(term17R'+term18R')*jumpu+term28R-C1R*term30R);           


     end
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