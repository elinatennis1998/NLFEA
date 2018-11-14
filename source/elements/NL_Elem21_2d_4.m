%% DG implementation of large deformaiton 
% pure displacement 2D case
% Tim Truster
% 10/2011
% modified by Pinlei Chen
% 04/23/2013
% for 2D large deformation with interface in it
% have tau and delta in it
% have d_ijklmn in it 
% this is a term by term version
% for body force problem, change the lint
% NOT suitable for different element type (Q4, Q9 T3 type)
% 03/21/2014 - TJT - verified to give quadratic convergence for general Q4
% deformation with patchtest4_BF_2d, including general and non-symmetric
% weight tensors gamL and gamR. The extra code lines in the file were also
% removed.
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
pencoeff = 3;1; 20;4;1;
elseif nitvms == 2 %Nitsche
pencoeff = 1;4;2;
else %RFB
pencoeff = 1;4;2;
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

switch isw %Task Switch
%%
    case 1
        
        if ndf > 2
            
            for i = 3:ndf
                lie(i,1) = 0;
            end
            
        end
        
        nh3 = 12;
        
    case 3 %interface stiffness

        
        PatchEL = matepropL(4);
        PatchvL = matepropL(5);
        muL = PatchEL/(2*(1+PatchvL));%80.19;
        lamL = PatchvL*PatchEL/((1+PatchvL)*(1-2*PatchvL));%4e5;
        PatchER = matepropR(4);
        PatchvR = matepropR(5);
        muR = PatchER/(2*(1+PatchvR));%80.19;
        lamR = PatchvR*PatchER/((1+PatchvR)*(1-2*PatchvR));%4e5;
        I2 = eye(2);
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
        hr(nh3:nh3+3) = [gamL(1,1),gamL(1,2),gamL(2,1),gamL(2,2)];
%         gamL_list(iter+1,:,inter) = [gamL(1,1),gamL(1,2),gamL(2,1),gamL(2,2)];
        gamR = ebR^2*(edgeK\tauR);
        hr(nh3+4:nh3+7) = [gamR(1,1),gamR(1,2),gamR(2,1),gamR(2,2)];
%         gamR_list(iter+1,:,inter) = [gamR(1,1),gamR(1,2),gamR(2,1),gamR(2,2)];
%         gamL_list(elem,:) = [gamL(1),gamL(2),gamL(3),gamL(4),gamL(5),gamL(6),gamL(7),gamL(8),gamL(9)];
%         gamR_list(elem,:) = [gamR(1),gamR(2),gamR(3),gamR(4),gamR(5),gamR(6),gamR(7),gamR(8),gamR(9)];
        ep = pencoeff*intedge*inv(edgeK); 
        hr(nh3+8:nh3+11) = [ep(1,1),ep(1,2),ep(2,1),ep(2,2)];
%         ep_list(iter+1,:,inter) = [ep(1,1),ep(1,2),ep(2,1),ep(2,2)];
%         ep_List(elem,:) = [ep(1),ep(2),ep(3),ep(4),ep(5),ep(6),ep(7),ep(8),ep(9)];
%       ep_List_R(elem,:) = [tauR(1),tauR(2),tauR(3),tauR(4),tauR(5),tauR(6),tauR(7),tauR(8),tauR(9)];         
   else
        gamL = [hr(nh3) hr(nh3+1)
                hr(nh3+2) hr(nh3+3)];
        gamR = [hr(nh3+4) hr(nh3+5)
                hr(nh3+6) hr(nh3+7)];
        ep =   [hr(nh3+8) hr(nh3+9)
                hr(nh3+10) hr(nh3+11)];
%         gamL = [gamL_list(iterset+1,1,inter) gamL_list(iterset+1,2,inter)
%                 gamL_list(iterset+1,3,inter) gamL_list(iterset+1,4,inter)];
%         gamR = [gamR_list(iterset+1,1,inter) gamR_list(iterset+1,2,inter)
%                 gamR_list(iterset+1,3,inter) gamR_list(iterset+1,4,inter)];
%         ep = [ep_list(iterset+1,1,inter) ep_list(iterset+1,2,inter)
%               ep_list(iterset+1,3,inter) ep_list(iterset+1,4,inter)];
   end
%        gamL = [0.5 0
%                 0 0.5];
%        gamR = [0.5 0
%                 0 0.5];     
%        ep = 200*eye(2);
   
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
%           [PxyL,shgsL,JdetL,bubbleL,xsL] = shgt(xlL(:,1:nelL)+ulL(:,1:nelL),nelL,shldL,shlsL,nen,0,0,be);
           else
           [shlL,shldL,shlsL,be] = shlq(rL,lits,nelL,nelL,0,0);
           [PxyL,shgsL,JdetL,bubbleL,XsL] = shgq(xlL(:,1:nelL),nelL,shldL,shlsL,nen,0,0,be); 
%           [PxyL,shgsL,JdetL,bubbleL,xsL] = shgq(xlL(:,1:nelL)+ulL(:,1:nelL),nelL,shldL,shlsL,nen,0,0,be);      
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
%            [PxyR,shgsR,JdetR,bubbleR,xsR] = shgt(xlR(:,1:nelR)+ulR(:,1:nelR),nelR,shldR,shlsR,nen,0,0,be);  
           else
            [shlR,shldR,shlsR,be] = shlq(rR,sR,nelR,nelR,0,0);
            [PxyR,shgsR,JdetR,bubbleR,XsR] = shgq(xlR(:,1:nelR),nelR,shldR,shlsR,nen,0,0,be);     
%            [PxyR,shgsR,JdetR,bubbleR,xsR] = shgq(xlR(:,1:nelR)+ulR(:,1:nelR),nelR,shldR,shlsR,nen,0,0,be);  
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
%             [fInvL,JxXL,FL] = kine2d(QxyL,-ulL,nelL,0); %this is equivalent to ikine2d
              [FL,JxXL,fInvL] = kine2d(QxyL,ulL,nelL,0);           
%             JxXL = 1/JxXL; %this is equivalent to ikine2d
    %         [F,JxX,fi] = ikine2d(Qxy,ul,nel); %this is equivalent to above two lines
%             JdetL = JdetL/JxXL;            
            PL = muL*(FL - fInvL') + lamL*(JxXL^2-JxXL)*fInvL';
            
%            [fInvR,JxXR,FR] = kine2d(QxyR,-ulR,nelR,0); %this is equivalent to ikine2d
            [FR,JxXR,fInvR] = kine2d(QxyR,ulR,nelR,0);
%             JxXR = 1/JxXR; %this is equivalent to ikine2d
    %         [F,JxX,fi] = ikine2d(Qxy,ul,nel); %this is equivalent to above two lines
%             JdetR = JdetR/JxXR;               
            PR = muR*(FR - fInvR') + lamR*(JxXR^2-JxXR)*fInvR';           
             
             %Evaluate tangent and normal vectors
            T1L = [XsL(:,1); 0];
            [Tm1L, Tu1L] = VecNormalize(T1L);
            T2L = [0; 0; 1];
            Tm2L = 1;
            Tu2L = T2L';
            T3L = VecCrossProd(T1L,T2L);
            [Tm3L, Tu3L] = VecNormalize(T3L);
            %normal vectors
            NLx = Tu3L(1);
            NLy = Tu3L(2);
            Nvec = [NLx; NLy];         
            C1L = drdr*Wgt*Tm3L;
            C1R = C1L;
%             c1 = C1L;
            % Nanson formula
            t3L=JxXL*[fInvL zeros(2,1); zeros(1,2) 1]'*T3L';
            [tm3L, tu3L] = VecNormalize(t3L);
            t3R=JxXR*[fInvR zeros(2,1); zeros(1,2) 1]'*-T3L';
            [tm3R, tu3R] = VecNormalize(t3R);
            c1L = Wgt*drdr*tm3L;
            c1R = Wgt*drdr*tm3R;

                %normal vectors
                nLx = tu3L(1);
                nLy = tu3L(2);
                nRx = tu3R(1);
                nRy = tu3R(2);
                %tagent vectors
%                 tLx = tu1L(1);
%                 tLy = tu1L(2);
%                 tRx = tu1R(1);
%                 tRy = tu1R(2);             
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

             rhspulL = reshape(ulL,ndf*nelL,1);
             rhspulR = reshape(ulR,ndf*nelR,1);
                    
            jumpu = NmatR*rhspulR - NmatL*rhspulL;  %jumpu=uR-uL 
                 
            % make gam's symmetric to help with debugging; they are not
            % always symmetric, so gamL(i,k) /= gamL(k,i)
%              gamL = 1/2*(gamL + gamL');
%              gamR = 1/2*(gamR + gamR');
            tvtr = gamL*PL*Nvec + gamR*PR*Nvec; %average stress
            % copy some terms into other variables that are the same name
            % as in NL_Elem2_3dMDG.m
            nvec = Nvec;
            c1 = C1L;

            % Not needed for linear problems
            ElemFL = ElemFL - C1L*( - NmatL'*(tvtr + ep*jumpu));%  + bnAdN1'*jumpu);% 
            ElemFR = ElemFR - C1R*( + NmatR'*(tvtr + ep*jumpu));%  + bnAdN2'*jumpu);% 
%             ElemFL = ElemFL - C1L*( - NmatL'*(ep*jumpu));%  + bnAdN1'*jumpu);% 
%             ElemFR = ElemFR - C1R*( + NmatR'*(ep*jumpu));%  + bnAdN2'*jumpu);% 

            % Non-symmetric term, enforcing traction continuity
            WLi = 0;
            WRi = 0;
            dULj = 0;
            dURj = 0;
            ElemKLL1 = zeros(nstL,nstL);
            ElemKLR1 = zeros(nstL,nstR);
            ElemKRL1 = zeros(nstR,nstL);
            ElemKRR1 = zeros(nstR,nstR);
            for Aw = 1:4
                WLi = shlL(Aw);
                WRi = shlR(Aw);
                for Bu = 1:4
                    dULj(1) = QxyL(Bu,1);
                    dULj(2) = QxyL(Bu,2);
                    dURj(1) = QxyR(Bu,1);
                    dURj(2) = QxyR(Bu,2);
                    for i = 1:2
                        for j = 1:2
                            termLL = 0;
                            termLR = 0;
                            termRL = 0;
                            termRR = 0;
                            for k = 1:2
                                for I = 1:2
                                    for J = 1:2
                                        ALkIjJ = muL*(I2(k,j)*I2(I,J)+fInvL(I,j)*fInvL(J,k)) ...
                                               + lamL*((2*JxXL^2-JxXL)*fInvL(I,k)*fInvL(J,j) - (JxXL^2-JxXL)*fInvL(I,j)*fInvL(J,k));
                                        ARkIjJ = muR*(I2(k,j)*I2(I,J)+fInvR(I,j)*fInvR(J,k)) ...
                                               + lamR*((2*JxXR^2-JxXR)*fInvR(I,k)*fInvR(J,j) - (JxXR^2-JxXR)*fInvR(I,j)*fInvR(J,k));
                                        termLL = termLL - c1*WLi*gamL(i,k)*nvec(I)*ALkIjJ*dULj(J);
                                        termLR = termLR - c1*WLi*gamR(i,k)*nvec(I)*ARkIjJ*dURj(J);
                                        termRL = termRL + c1*WRi*gamL(i,k)*nvec(I)*ALkIjJ*dULj(J);
                                        termRR = termRR + c1*WRi*gamR(i,k)*nvec(I)*ARkIjJ*dURj(J);
                                    end
                                end
                            end
                            ElemKLL1(ndf*(Aw-1)+i,ndf*(Bu-1)+j) = ElemKLL1(ndf*(Aw-1)+i,ndf*(Bu-1)+j) + termLL;
                            ElemKLR1(ndf*(Aw-1)+i,ndf*(Bu-1)+j) = ElemKLR1(ndf*(Aw-1)+i,ndf*(Bu-1)+j) + termLR;
                            ElemKRL1(ndf*(Aw-1)+i,ndf*(Bu-1)+j) = ElemKRL1(ndf*(Aw-1)+i,ndf*(Bu-1)+j) + termRL;
                            ElemKRR1(ndf*(Aw-1)+i,ndf*(Bu-1)+j) = ElemKRR1(ndf*(Aw-1)+i,ndf*(Bu-1)+j) + termRR;
                        end
                    end
                end
            end

%             JxXL = 1;
%             JxXR = 1;
%             fInvL = I2;
%             fInvR = I2;
            
            % The version below was copied by Tim on 3/21 from
            % NL_Elem2_3dMDG.m and all terms verified for quadratic
            % convergence when gamL = general and gamR = general
            % Symmetrizing term, "stresses" for weighting function
            WLi = 0;
            WRi = 0;
            for Aw = 1:4
                WLi(1) = QxyL(Aw,1);
                WLi(2) = QxyL(Aw,2);
                WRi(1) = QxyR(Aw,1);
                WRi(2) = QxyR(Aw,2);
                for i = 1:2
                    for j = 1:2
                        termL = 0;
                        termR = 0;
                        for k = 1:2
                        for I = 1:2
                            for J = 1:2
                                ALiIkJ = muL*(I2(i,k)*I2(I,J)+fInvL(I,k)*fInvL(J,i)) ...
                                       + lamL*((2*JxXL^2-JxXL)*fInvL(I,i)*fInvL(J,k) - (JxXL^2-JxXL)*fInvL(I,k)*fInvL(J,i));
                                ARiIkJ = muR*(I2(i,k)*I2(I,J)+fInvR(I,k)*fInvR(J,i)) ...
                                       + lamR*((2*JxXR^2-JxXR)*fInvR(I,i)*fInvR(J,k) - (JxXR^2-JxXR)*fInvR(I,k)*fInvR(J,i));
                                termL = termL + c1*jumpu(j)*gamL(j,k)*nvec(J)*ALiIkJ*WLi(I);
                                termR = termR + c1*jumpu(j)*gamR(j,k)*nvec(J)*ARiIkJ*WRi(I);
                            end
                        end
                        end
                        ElemFL(ndf*(Aw-1)+i) = ElemFL(ndf*(Aw-1)+i) - termL;
                        ElemFR(ndf*(Aw-1)+i) = ElemFR(ndf*(Aw-1)+i) - termR;
                    end
                end
            end
% %             ElemFL = ElemFL - c1*( + bnAdN1'*jumpu);% );% 
% %             ElemFR = ElemFR - c1*( + bnAdN2'*jumpu);% );% 

            % Symmetrizing term, "stresses" for weighting function
            WLi = 0;
            WRi = 0;
            dULj = 0;
            dURj = 0;
            ElemKLL2 = zeros(nstL,nstL);
            ElemKLR2 = zeros(nstL,nstR);
            ElemKRL2 = zeros(nstR,nstL);
            ElemKRR2 = zeros(nstR,nstR);
            for Aw = 1:4
                WLi(1) = QxyL(Aw,1);
                WLi(2) = QxyL(Aw,2);
                WRi(1) = QxyR(Aw,1);
                WRi(2) = QxyR(Aw,2);
                for Bu = 1:4
                    dULj = shlL(Bu);
                    dURj = shlR(Bu);
                    for i = 1:2
                        for j = 1:2
                            termLL = 0;
                            termLR = 0;
                            termRL = 0;
                            termRR = 0;
                            for k = 1:2
                            for I = 1:2
                                for J = 1:2
                                    ALiIkJ = muL*(I2(i,k)*I2(I,J)+fInvL(I,k)*fInvL(J,i)) ...
                                           + lamL*((2*JxXL^2-JxXL)*fInvL(I,i)*fInvL(J,k) - (JxXL^2-JxXL)*fInvL(I,k)*fInvL(J,i));
                                    ARiIkJ = muR*(I2(i,k)*I2(I,J)+fInvR(I,k)*fInvR(J,i)) ...
                                           + lamR*((2*JxXR^2-JxXR)*fInvR(I,i)*fInvR(J,k) - (JxXR^2-JxXR)*fInvR(I,k)*fInvR(J,i));
                                    termLL = termLL - c1*dULj*gamL(j,k)*nvec(J)*ALiIkJ*WLi(I);
                                    termRL = termRL - c1*dULj*gamR(j,k)*nvec(J)*ARiIkJ*WRi(I);
                                    termLR = termLR + c1*dURj*gamL(j,k)*nvec(J)*ALiIkJ*WLi(I);
                                    termRR = termRR + c1*dURj*gamR(j,k)*nvec(J)*ARiIkJ*WRi(I);
                                end
                            end
                            end
                            ElemKLL2(ndf*(Aw-1)+i,ndf*(Bu-1)+j) = ElemKLL2(ndf*(Aw-1)+i,ndf*(Bu-1)+j) + termLL;
                            ElemKLR2(ndf*(Aw-1)+i,ndf*(Bu-1)+j) = ElemKLR2(ndf*(Aw-1)+i,ndf*(Bu-1)+j) + termLR;
                            ElemKRL2(ndf*(Aw-1)+i,ndf*(Bu-1)+j) = ElemKRL2(ndf*(Aw-1)+i,ndf*(Bu-1)+j) + termRL;
                            ElemKRR2(ndf*(Aw-1)+i,ndf*(Bu-1)+j) = ElemKRR2(ndf*(Aw-1)+i,ndf*(Bu-1)+j) + termRR;
                        end
                    end
                end
            end

            % Unexpected terms, involves jumpu and 6th order tensor
            WLi = 0;
            WRi = 0;
            dULj = 0;
            dURj = 0;
            ElemKLL3 = zeros(nstL,nstL);
            ElemKRR3 = zeros(nstR,nstR);
            for Aw = 1:4
                WLi(1) = QxyL(Aw,1);
                WLi(2) = QxyL(Aw,2);
                WRi(1) = QxyR(Aw,1);
                WRi(2) = QxyR(Aw,2);
                for Bu = 1:4
                    dULj(1) = QxyL(Bu,1);
                    dULj(2) = QxyL(Bu,2);
                    dURj(1) = QxyR(Bu,1);
                    dURj(2) = QxyR(Bu,2);
                    for i = 1:2
                        for j = 1:2
                            termLL = 0;
                            termRR = 0;
                            for I = 1:2
                                for J = 1:2
                                    termDgLL = 0;
                                    termDgRR = 0;
                                    for kk = 1:2
                                    for k = 1:2
                                        for K = 1:2
DLiIjJkK = - muL*(fInvL(I,k)*fInvL(K,j)*fInvL(J,i)+fInvL(I,j)*fInvL(J,k)*fInvL(K,i)) ...
           + lamL*((4*JxXL^2-JxXL)*fInvL(I,i)*fInvL(J,j)*fInvL(K,k) + (JxXL^2-JxXL)*(fInvL(I,k)*fInvL(K,j)*fInvL(J,i) + fInvL(I,j)*fInvL(J,k)*fInvL(K,i))...
           - (2*JxXL^2-JxXL)*(fInvL(I,j)*fInvL(J,i)*fInvL(K,k) + fInvL(I,k)*fInvL(K,i)*fInvL(J,j) + fInvL(I,i)*fInvL(J,k)*fInvL(K,j)));
DRiIjJkK = - muR*(fInvR(I,k)*fInvR(K,j)*fInvR(J,i)+fInvR(I,j)*fInvR(J,k)*fInvR(K,i)) ...
           + lamR*((4*JxXR^2-JxXR)*fInvR(I,i)*fInvR(J,j)*fInvR(K,k) + (JxXR^2-JxXR)*(fInvR(I,k)*fInvR(K,j)*fInvR(J,i) + fInvR(I,j)*fInvR(J,k)*fInvR(K,i))...
           - (2*JxXR^2-JxXR)*(fInvR(I,j)*fInvR(J,i)*fInvR(K,k) + fInvR(I,k)*fInvR(K,i)*fInvR(J,j) + fInvR(I,i)*fInvR(J,k)*fInvR(K,j)));
                                            termDgLL = termDgLL + gamL(kk,k)*DLiIjJkK*jumpu(kk)*nvec(K);
                                            termDgRR = termDgRR + gamR(kk,k)*DRiIjJkK*jumpu(kk)*nvec(K);
                                        end
                                    end
                                    end
                                    termLL = termLL + c1*dULj(J)*termDgLL*WLi(I);
                                    termRR = termRR + c1*dURj(J)*termDgRR*WRi(I);
                                end
                            end
                            ElemKLL3(ndf*(Aw-1)+i,ndf*(Bu-1)+j) = ElemKLL3(ndf*(Aw-1)+i,ndf*(Bu-1)+j) + termLL;
                            ElemKRR3(ndf*(Aw-1)+i,ndf*(Bu-1)+j) = ElemKRR3(ndf*(Aw-1)+i,ndf*(Bu-1)+j) + termRR;
                        end
                    end
                end
            end

            ElemKLL = ElemKLL + ElemKLL1 + ElemKLL2 + ElemKLL3 + c1*(NmatL'*ep*NmatL);%
            ElemKLR = ElemKLR + ElemKLR1 + ElemKLR2 - c1*(NmatL'*ep*NmatR);
            ElemKRL = ElemKRL + ElemKRL1 + ElemKRL2 - c1*(NmatR'*ep*NmatL);
            ElemKRR = ElemKRR + ElemKRR1 + ElemKRR2 + ElemKRR3 + c1*(NmatR'*ep*NmatR);%    

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
        ElemK = zeros(nst,nst);
        ElemF = zeros(nst,1);
%%

    case -1 % boundary tractions
        
        ElemF = zeros(nst,1);

end