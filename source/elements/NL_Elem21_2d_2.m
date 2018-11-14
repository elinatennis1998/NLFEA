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
% 03/21/2014 - TJT - verified to give quadratic convergence for general Q4
% deformation with patchtest4_BF_2d, including general and non-symmetric
% weight tensors gamL and gamR. The extra code lines in the file were also
% removed. Also verified that the stiffness matrix for full method remains
% symmetric; the transpose of gamL/R is important for keeping that feature.
% Modified formulations with tau and force/stiffness terms also verified.
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
        iste = 4+4+4+3+6;
        
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
        

        ElemEL = matepropL(1);
        ElemvL = matepropL(2);
        ElemER = matepropR(1);
        ElemvR = matepropR(2);
        lam = getlam(matepropL);
       
        NmatL = zeros(2,nstL);
        BmatL = zeros(4,nstL);
        BmatLi = zeros(4,nstL);
        NmatR = zeros(2,nstR);
        BmatR = zeros(4,nstR);
        BmatRi = zeros(4,nstR);
        dmatL = zeros(9,3);
        dmatR = zeros(9,3);
       

        % Determin bounds of integration segment
        InterConn2D2 % InterConn2DT %
        
        % Set jacobian for integration space
        drdr = (eL2 - eL1)/drL;
        drdrR = (eR2 - eR1)/drR; % This causes errors for partial interface segments
        m = (eR2-eR1)/(eL1-eL2);
%         if nelL == 3 || nelL == 6  
%         lint = 3;10;2;3;
%         else
%         lint = 3;10;4;10;2;3; %10 for body force problem; 4 for other problem 
%         end
        if exist('iprob','var') == 1 && iprob == 6
            lintDG = 10;
        else
            if ~exist('lintDG','var') % allows the value to be changed in input file
            lintDG = 3;10;2;3;
            end
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
        
        for ie = 1:lintDG            
% For separate bubble types on T and Q
           if nelL == 3 || nelL == 6
                [Wgt,litr,lits] = intpntt(ie,lintDG,1); 
                 ebeL = edgebubble(litr,lits,nelL);  %edgebubble is for T3 element
%                 rL = drdr*(litr-roL)+eL1;
%                 [shlL,shldL,shlsL,be] = shlt(rL,lits,nelL,nel2L,0,0);
           else
                [Wgt,litr,lits] = intpntq(ie,lintDG,1);
%                 rL = drdr*(litr-roL)+eL1;
%                 [shlL,shldL,shlsL,be] = shlq(rL,lits,nelL,nel2L,0,0);
                 ebeL = edgebubbleQ(litr,lits,nelL);
           end
                            
           if nelR == 3 || nelR == 6
                [WgtR,litr,lits] = intpntt(ie,lintDG,1); 
                 ebeR = edgebubble(litr,lits,nelR);
%                 rR = drdr*(litr-roL)+eR1;
%                 [shlR,shldR,shlsR,be] = shlt(rR,lits,nelR,nel2R,0,0);
           else
                [WgtR,litr,lits] = intpntq(ie,lintDG,1);
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
%        ep = 20*eye(2);
   
%        s = -1;
        ll=0;           
     for l = 1:lintDG

            ll = ll + 1;
            if nelL == 3 || nelL == 6
                [Wgt,litr,lits] = intpntt(ll,lintDG,1);
            elseif nelL == 4 || nelL == 9
                [Wgt,litr,lits] = intpntq(ll,lintDG,1);
            end
                    
%            b = edgebubble(litr,0);
            
            rL = drdr*(litr-roL)+eL1;
           if nelL == 3 || nelL == 6
           [shlL,shldL,shlsL,be] = shlt(rL,lits,nelL,nelL,0,0);
           [PxyL,shgsL,JdetL,bubbleL,XsL] = shgt(xlL(:,1:nelL),nelL,shldL,shlsL,nen,0,0,be); 
           [QxyL,shgsL,JdetL,bubbleL,xsL] = shgt(xlL(:,1:nelL)+ulL(:,1:nelL),nelL,shldL,shlsL,nen,0,0,be);
           else
           [shlL,shldL,shlsL,be] = shlq(rL,lits,nelL,nelL,0,0);
           [PxyL,shgsL,JdetL,bubbleL,XsL] = shgq(xlL(:,1:nelL),nelL,shldL,shlsL,nen,0,0,be); 
           [QxyL,shgsL,JdetL,bubbleL,xsL] = shgq(xlL(:,1:nelL)+ulL(:,1:nelL),nelL,shldL,shlsL,nen,0,0,be);      
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
            [QxyR,shgsR,JdetR,bubbleR,xsR] = shgt(xlR(:,1:nelR)+ulR(:,1:nelR),nelR,shldR,shlsR,nen,0,0,be);  
           else
            [shlR,shldR,shlsR,be] = shlq(rR,sR,nelR,nelR,0,0);
            [PxyR,shgsR,JdetR,bubbleR,XsR] = shgq(xlR(:,1:nelR),nelR,shldR,shlsR,nen,0,0,be);     
            [QxyR,shgsR,JdetR,bubbleR,xsR] = shgq(xlR(:,1:nelR)+ulR(:,1:nelR),nelR,shldR,shlsR,nen,0,0,be);  
           end

 % initial stress sigma and material bulk term c_ijkl
           for mm = 1:nelL    
           NmatL(:,2*mm-1:2*mm) = [shlL(mm,1)     0          
                                    0        shlL(mm,1)]; 
                          
           BmatL(:,2*mm-1:2*mm) = [QxyL(mm,1) 0             
                                    0         QxyL(mm,2)         
                                    QxyL(mm,2) QxyL(mm,1)        
                                     QxyL(mm,2) -QxyL(mm,1)];    
           end 
            
           for mm = 1:nelR    
           NmatR(:,2*mm-1:2*mm) = [shlR(mm,1)     0          
                                     0        shlR(mm,1)]; 
                          
           BmatR(:,2*mm-1:2*mm) = [QxyR(mm,1) 0             
                                    0         QxyR(mm,2)         
                                    QxyR(mm,2) QxyR(mm,1)        
                                     QxyR(mm,2) -QxyR(mm,1)];                     
           end  
           
           if modifDG == 2 || modifDG == 5 % frozen terms using initial configuration
           for mm = 1:nelL    
                          
           BmatLi(:,2*mm-1:2*mm) = [PxyL(mm,1) 0             
                                    0         PxyL(mm,2)         
                                    PxyL(mm,2) PxyL(mm,1)        
                                    PxyL(mm,2) -PxyL(mm,1)];    
           end 
            
           for mm = 1:nelR     
                          
           BmatRi(:,2*mm-1:2*mm) = [PxyR(mm,1) 0             
                                    0         PxyR(mm,2)         
                                    PxyR(mm,2) PxyR(mm,1)        
                                    PxyR(mm,2) -PxyR(mm,1)];                     
           end  
           
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
          
            
            % Cuchy stress tensors and material moduli for both sides
            [sigma2L, cmatL] = SigmaCmat2i(FL,JxXL,matepropL,lam);
            [sigma2Li, cmatLi] = SigmaCmat2i(eye(2),1,matepropL,lam);
            
            SmatL1=[sigma2L(1), sigma2L(3)
                    sigma2L(3), sigma2L(2)];           
            SmatL1i=[sigma2Li(1), sigma2Li(3)
                    sigma2Li(3), sigma2Li(2)];
            
            [sigma2R, cmatR] = SigmaCmat2i(FR,JxXR,matepropR,lam);
            [sigma2Ri, cmatRi] = SigmaCmat2i(eye(2),1,matepropR,lam);
                                 
            SmatR1=[sigma2R(1), sigma2R(3)
                    sigma2R(3), sigma2R(2)];     
            SmatR1i=[sigma2Ri(1), sigma2Ri(3)
                    sigma2Ri(3), sigma2Ri(2)];
             
             %Evaluate tangent and normal vectors
            T1L = [XsL(:,1); 0];
            [Tm1L, Tu1L] = VecNormalize(T1L);
            T2L = [0; 0; 1];
            Tm2L = 1;
            Tu2L = T2L';
            T3L = VecCrossProd(T1L,T2L);
            [Tm3L, Tu3L] = VecNormalize(T3L);
            C1L = drdr*Wgt*Tm3L;
            C1R = C1L;
            
            % Nanson formula
            t3L=JxXL*[fiL zeros(2,1); zeros(1,2) 1]'*T3L';
            [tm3L, tu3L] = VecNormalize(t3L);
            t3R=JxXR*[fiR zeros(2,1); zeros(1,2) 1]'*-T3L';
            [tm3R, tu3R] = VecNormalize(t3R);
            c1L = Wgt*drdr*tm3L;
            c1R = Wgt*drdr*tm3R;

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
            NvectL1= [NLx 0  NLy 
                      0  NLy NLx];
            NvectR1 = -[NLx 0 NLy
                      0  NLy NLx]; 
%             nvectL2 = [eye(3,3)*nLx zeros(3,3) eye(3,3)*nLy
%                        zeros(3,3)    eye(3,3)*nLy eye(3,3)*nLx];
%             nvectR2 = [eye(3,3)*nRx zeros(3,3) eye(3,3)*nRy
%                        zeros(3,3)    eye(3,3)*nRy eye(3,3)*nRx];                       
            nvecL = [nLx; nLy];
            nvecR = [nRx; nRy];


           % stiffness matrix terms for weighted stress and material moduli
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

            term17L=P2'*SmatnL*gamL';
            term17R=P2'*SmatnR*gamR';
            term17Li = 0*term17L;
            term17Ri = 0*term17R;

            term18L=P1'*(gamL*nvectL1*cmatL(1:3,1:3))';
            term18R=P1'*(gamR*nvectR1*cmatR(1:3,1:3))';
            term18Li=P1'*(gamL*NvectL1*cmatLi(1:3,1:3))';
            term18Ri=P1'*(gamR*NvectR1*cmatRi(1:3,1:3))';
            
            
            %average stress term for force vector
            term28L=NmatL'*(c1L*gamL*SmatL1*nvecL-c1R*gamR*SmatR1*nvecR);
            term28R=NmatR'*(c1L*gamL*SmatL1*nvecL-c1R*gamR*SmatR1*nvecR);

             rhspulL = reshape(ulL,ndf*nelL,1);
             rhspulR = reshape(ulR,ndf*nelR,1);
                    
            jumpu = NmatL*rhspulL - NmatR*rhspulR;  %jumpu=uL-uR

              term30L=NmatL'*ep*jumpu;
              term30R=NmatR'*ep*jumpu;

            % Unexpected terms, involves jumpu and 6th order tensor
             gamajumpuL = gamL'*jumpu;
             gamajumpuR = gamR'*jumpu;             
             term5L=cmatnBL*[eye(3,3)*gamajumpuL(1) eye(3,3)*gamajumpuL(2)]'*(P1*BmatL);
             term5R=cmatnBR*[eye(3,3)*gamajumpuR(1) eye(3,3)*gamajumpuR(2)]'*(P1*BmatR);
            
             sig8L2 = [gamajumpuL(1) gamajumpuL(2)]*cmatnL;  

          sig8L3 = [sig8L2(1) 0 sig8L2(3) 0;
                    sig8L2(3) 0 sig8L2(2) 0;
                    0 sig8L2(1) 0 sig8L2(3);
                    0 sig8L2(3) 0 sig8L2(2)];
               
           term8L = BmatL'*P2'*sig8L3*(P3*BmatL);
          
             sig8R2 = [gamajumpuR(1) gamajumpuR(2)]*cmatnR;

            sig8R3 = [sig8R2(1) 0 sig8R2(3) 0;
                     sig8R2(3) 0 sig8R2(2) 0;
                     0 sig8R2(1) 0 sig8R2(3);
                     0 sig8R2(3) 0 sig8R2(2)];
               
              term8R =BmatR'*P2'*sig8R3*(P3*BmatR);

             [dmatL1]=dmat2_no_p(JxXL,matepropL,lam);
             [dmatR1]=dmat2_no_p(JxXR,matepropR,lam);        
%              dmatL2 = [eye(3,3)*jumpu(1) eye(3,3)*jumpu(2)]*[eye(3,3)*gamL(1,1) eye(3,3)*gamL(1,2)
%                                                              eye(3,3)*gamL(2,1) eye(3,3)*gamL(2,2)]*nvectL2*dmatL1/JxXL;  %missing J here
%              dmatR2 = [eye(3,3)*jumpu(1) eye(3,3)*jumpu(2)]*[eye(3,3)*gamR(1,1) eye(3,3)*gamR(1,2)
%                                                              eye(3,3)*gamR(2,1) eye(3,3)*gamR(2,2)]*nvectR2*dmatR1/JxXR;  %missing J here
             % Revised by TJT to use reshape function and smaller
             % multiplications
             dmatL2 = dmatL1/JxXL*nvectL1'*gamL'*jumpu;
             dmatL2 = reshape(dmatL2,3,3);
             dmatR2 = dmatR1/JxXR*nvectR1'*gamR'*jumpu;
             dmatR2 = reshape(dmatR2,3,3);
             
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
    end %lint        
% ElemKLL
            ElemK = [ElemKLL ElemKLR; ElemKRL ElemKRR];
            ElemF = [ElemFL; ElemFR];
            
%%
    case 6

       % generate the initial parameter 
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
       
        NmatL = zeros(2,nstL);
        BmatL = zeros(4,nstL);
        NmatR = zeros(2,nstR);
        BmatR = zeros(4,nstR);
       

        % Determin bounds of integration segment
        InterConn2D2 % InterConn2DT %
        
        % Set jacobian for integration space
        drdr = (eL2 - eL1)/drL;
        drdrR = (eR2 - eR1)/drR; % This causes errors for partial interface segments
        m = (eR2-eR1)/(eL1-eL2);
%         if nelL == 3 || nelL == 6  
%         lint = 3;10;2;3;
%         else
%         lint = 10;3;4;10;2;3; %10 for body force problem; 4 for other problem 
%         end
        if exist('iprob','var') == 1 && iprob == 6
            lintDG = 10;
        else
            if ~exist('lintDG','var') % allows the value to be changed in input file
            lintDG = 3;10;2;3;
            end
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
        
        for ie = 1:lintDG            
% For separate bubble types on T and Q
           if nelL == 3 || nelL == 6
                [Wgt,litr,lits] = intpntt(ie,lintDG,1); 
                 ebeL = edgebubble(litr,lits,nelL);  %edgebubble is for T3 element
%                 rL = drdr*(litr-roL)+eL1;
%                 [shlL,shldL,shlsL,be] = shlt(rL,lits,nelL,nel2L,0,0);
           else
                [Wgt,litr,lits] = intpntq(ie,lintDG,1);
%                 rL = drdr*(litr-roL)+eL1;
%                 [shlL,shldL,shlsL,be] = shlq(rL,lits,nelL,nel2L,0,0);
                 ebeL = edgebubbleQ(litr,lits,nelL);
           end
                            
           if nelR == 3 || nelR == 6
                [WgtR,litr,lits] = intpntt(ie,lintDG,1); 
                 ebeR = edgebubble(litr,lits,nelR);
%                 rR = drdr*(litr-roL)+eR1;
%                 [shlR,shldR,shlsR,be] = shlt(rR,lits,nelR,nel2R,0,0);
           else
                [WgtR,litr,lits] = intpntq(ie,lintDG,1);
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
%        ep = 20*eye(2);
   
%        s = -1;
        ll=0;           
     for l = 1:lintDG

            ll = ll + 1;
            if nelL == 3 || nelL == 6
                [Wgt,litr,lits] = intpntt(ll,lintDG,1);
            elseif nelL == 4 || nelL == 9
                [Wgt,litr,lits] = intpntq(ll,lintDG,1);
            end
                    
%            b = edgebubble(litr,0);
            
            rL = drdr*(litr-roL)+eL1;
           if nelL == 3 || nelL == 6
           [shlL,shldL,shlsL,be] = shlt(rL,lits,nelL,nelL,0,0);
%            [PxyL,shgsL,JdetL,bubbleL,XsL] = shgt(xlL(:,1:nelL),nelL,shldL,shlsL,nen,0,0,be); 
           [PxyL,shgsL,JdetL,bubbleL,xsL] = shgt(xlL(:,1:nelL)+ulL(:,1:nelL),nelL,shldL,shlsL,nen,0,0,be);
           else
           [shlL,shldL,shlsL,be] = shlq(rL,lits,nelL,nelL,0,0);
%            [PxyL,shgsL,JdetL,bubbleL,XsL] = shgq(xlL(:,1:nelL),nelL,shldL,shlsL,nen,0,0,be); 
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
%             [PxyR,shgsR,JdetR,bubbleR,XsR] = shgt(xlR(:,1:nelR),nelR,shldR,shlsR,nen,0,0,be);     
            [PxyR,shgsR,JdetR,bubbleR,xsR] = shgt(xlR(:,1:nelR)+ulR(:,1:nelR),nelR,shldR,shlsR,nen,0,0,be);  
           else
            [shlR,shldR,shlsR,be] = shlq(rR,sR,nelR,nelR,0,0);
%             [PxyR,shgsR,JdetR,bubbleR,XsR] = shgq(xlR(:,1:nelR),nelR,shldR,shlsR,nen,0,0,be);     
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
           end 
            
           for mm = 1:nelR    
           NmatR(:,2*mm-1:2*mm) = [shlR(mm,1)     0          
                                     0        shlR(mm,1)]; 
                          
           BmatR(:,2*mm-1:2*mm) = [QxyR(mm,1) 0             
                                    0         QxyR(mm,2)         
                                    QxyR(mm,2) QxyR(mm,1)        
                                     QxyR(mm,2) -QxyR(mm,1)];                     
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
          
            
            % Cuchy stress tensors and material moduli for both sides
            [sigma2L, cmatL] = SigmaCmat2i(FL,JxXL,matepropL,lam);
            
            SmatL1=[sigma2L(1), sigma2L(3)
                    sigma2L(3), sigma2L(2)];           
            
            [sigma2R, cmatR] = SigmaCmat2i(FR,JxXR,matepropR,lam);
                                 
            SmatR1=[sigma2R(1), sigma2R(3)
                    sigma2R(3), sigma2R(2)];     
             
             %Evaluate tangent and normal vectors
            T1L = [XsL(:,1); 0];
            [Tm1L, Tu1L] = VecNormalize(T1L);
            T2L = [0; 0; 1];
            Tm2L = 1;
            Tu2L = T2L';
            T3L = VecCrossProd(T1L,T2L);
            [Tm3L, Tu3L] = VecNormalize(T3L);
            C1L = drdr*Wgt*Tm3L;
            C1R = C1L;
            
            % Nanson formula
            t3L=JxXL*[fiL zeros(2,1); zeros(1,2) 1]'*T3L';
            [tm3L, tu3L] = VecNormalize(t3L);
            t3R=JxXR*[fiR zeros(2,1); zeros(1,2) 1]'*-T3L';
            [tm3R, tu3R] = VecNormalize(t3R);
            c1L = Wgt*drdr*tm3L;
            c1R = Wgt*drdr*tm3R;

            %normal vectors
            nLx = tu3L(1);
            nLy = tu3L(2);
            nRx = tu3R(1);
            nRy = tu3R(2);            
            nvectL1= [nLx 0   nLy   
                      0  nLy  nLx ];
            nvectR1 = [nRx 0  nRy  
                      0  nRy  nRx ]; 
%             nvectL2 = [eye(3,3)*nLx zeros(3,3) eye(3,3)*nLy
%                        zeros(3,3)    eye(3,3)*nLy eye(3,3)*nLx];
%             nvectR2 = [eye(3,3)*nRx zeros(3,3) eye(3,3)*nRy
%                        zeros(3,3)    eye(3,3)*nRy eye(3,3)*nRx];                       
            nvecL = [nLx; nLy];
            nvecR = [nRx; nRy];


           % stiffness matrix terms for weighted stress and material moduli
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

            term17L=P2'*SmatnL*gamL';
            term17R=P2'*SmatnR*gamR';

            term18L=P1'*(gamL*nvectL1*cmatL(1:3,1:3))';
            term18R=P1'*(gamR*nvectR1*cmatR(1:3,1:3))';
            
            
            %average stress term for force vector
            term28L=NmatL'*(c1L*gamL*SmatL1*nvecL-c1R*gamR*SmatR1*nvecR);
            term28R=NmatR'*(c1L*gamL*SmatL1*nvecL-c1R*gamR*SmatR1*nvecR);

             rhspulL = reshape(ulL,ndf*nelL,1);
             rhspulR = reshape(ulR,ndf*nelR,1);
                    
            jumpu = NmatL*rhspulL - NmatR*rhspulR;  %jumpu=uL-uR

              term30L=NmatL'*ep*jumpu;
              term30R=NmatR'*ep*jumpu;


            % Combine contributions into element force vector
             if exist('modifDG','var') && modifDG > 0 % modify the DG terms used in the formulation
             
             % Penalty terms
                 
             ElemFL = ElemFL-(+C1L*term30L);
             ElemFR = ElemFR-(-C1R*term30R);
             
             if modifDG >= 1 
                 
             % average stress terms
                 
             ElemFL = ElemFL-(-term28L);
             ElemFR = ElemFR-(+term28R);
             
             if modifDG == 2 
                 
             % constant nonsymmetric terms

             ElemFL = ElemFL-(-C1L*BmatLi'*(term17Li+term18Li)*jumpu);
             ElemFR = ElemFR-(C1R* BmatRi'*(term17Ri+term18Ri)*jumpu);
             
             elseif modifDG == 3 
                 
             % nonsymmetric terms

             ElemFL = ElemFL-(-c1L*BmatL'*(term17L+term18L)*jumpu);
             ElemFR = ElemFR-(c1R* BmatR'*(term17R+term18R)*jumpu);
             
             elseif modifDG == 4 
             
             elseif modifDG == 5 
                 
             % constant nonsymmetric terms

             ElemFL = ElemFL-(-C1L*BmatLi'*(term17Li+term18Li)*jumpu);
             ElemFR = ElemFR-(C1R* BmatRi'*(term17Ri+term18Ri)*jumpu);
             
             end
             end

             else % full method
             
                ElemFL = ElemFL-(-term28L+C1L*term30L);
                ElemFR = ElemFR-(+term28R-C1R*term30R);  
                
                ElemFL = ElemFL-(-c1L*BmatL'*(term17L+term18L)*jumpu);
                ElemFR = ElemFR-(c1R* BmatR'*(term17R+term18R)*jumpu);

             end
    end %lint        
% ElemKLL
            ElemF = [ElemFL; ElemFR];
%%

    case -1 % boundary tractions
        
        ElemF = zeros(nst,1);
        
%%
    case 26 % Element stresses: output stability tensors
        
        ElemS = zeros(nestr,1);
        
        gamL = [hr(nh3) hr(nh3+1)
                hr(nh3+2) hr(nh3+3)];
        gamR = [hr(nh3+4) hr(nh3+5)
                hr(nh3+6) hr(nh3+7)];
        ep =   [hr(nh3+8) hr(nh3+9)
                hr(nh3+10) hr(nh3+11)];
            
            
        % Determin unit normal at segment midpoint
            
        % Determin bounds of integration segment
        InterConn2D2 % InterConn2DT %
        
        % Set jacobian for integration space
        drdr = (eL2 - eL1)/drL;
        drdrR = (eR2 - eR1)/drR; % This causes errors for partial interface segments
        m = (eR2-eR1)/(eL1-eL2);
        
        if nelL == 3 || nelL == 6
            [Wgt,litr,lits] = intpntt(ll,1,1);
        elseif nelL == 4 || nelL == 9
            [Wgt,litr,lits] = intpntq(ll,1,1);
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
            
    %plug in the -ul to get the inverse of deformation gradient fi as 'F' in this function
        [FL,JxXL,fiL] = kine2d(PxyL,ulL,nelL,0);
        [FR,JxXR,fiR] = kine2d(PxyR,ulR,nelR,0);

         %Evaluate tangent and normal vectors
        T1L = [XsL(:,1); 0];
        [Tm1L, Tu1L] = VecNormalize(T1L);
        T2L = [0; 0; 1];
        Tm2L = 1;
        Tu2L = T2L';
        T3L = VecCrossProd(T1L,T2L);
        [Tm3L, Tu3L] = VecNormalize(T3L);

        % Nanson formula
        t3L=JxXL*[fiL zeros(2,1); zeros(1,2) 1]'*T3L';
        [tm3L, tu3L] = VecNormalize(t3L);
        t3R=JxXR*[fiR zeros(2,1); zeros(1,2) 1]'*-T3L';
        [tm3R, tu3R] = VecNormalize(t3R);

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
        
        ElemS(1:4+4+4+3+6) = [reshape(ep,4,1); reshape(gamL,4,1); reshape(gamR,4,1); ep_nn; ep_tt; ep_nt; ...
            tL_nn; tL_tt; tL_nt; tR_nn; tR_tt; tR_nt];
            
        
    case 61 % form data structure for interface segments

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
        % Triangular or Quadrilateral face
%         if nelL == 3 || nelL == 6
            Iix(segment,1:2) = (segnode+1:segnode+2);
            ICoordinates(segnode+1:segnode+2,1:ndm) = xlL(1:ndm,1:2)';
            segnode = segnode + 2;
%         else
%             Iix(segment,1:4) = (segnode+1:segnode+4);
%             ICoordinates(segnode+1:segnode+4,1:ndm) = xlL(1:ndm,1:4)';
%             segnode = segnode + 4;
%         end
        Iix(segment,nenseg+1) = ma; % good habit to copy material ID too
        
    case 60 % output interface quantities for plotting
        
        % get segment number for DG element (usually segment=inter)
        segment = IElemSeg(1,inter);
    
        % get nodes on the interface segment
        ElemFlagI = Iix(segment,1:nenseg);
        xlI = ICoordinates(ElemFlagI,1:ndm)';
        
        lam = getlam(matepropL);
        

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
            
            nelseg = 2;
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
                    
%            b = edgebubble(litr,0);
            
            rL = drdr*(litr-roL)+eL1;
           if nelL == 3 || nelL == 6
           [shlL,shldL,shlsL,be] = shlt(rL,lits,nelL,nelL,0,0);
           [PxyL,shgsL,JdetL,bubbleL,XsL] = shgt(xlL(:,1:nelL),nelL,shldL,shlsL,nen,0,0,be); 
           [QxyL,shgsL,JdetL,bubbleL,xsL] = shgt(xlL(:,1:nelL)+ulL(:,1:nelL),nelL,shldL,shlsL,nen,0,0,be);
           else
           [shlL,shldL,shlsL,be] = shlq(rL,lits,nelL,nelL,0,0);
           [PxyL,shgsL,JdetL,bubbleL,XsL] = shgq(xlL(:,1:nelL),nelL,shldL,shlsL,nen,0,0,be); 
           [QxyL,shgsL,JdetL,bubbleL,xsL] = shgq(xlL(:,1:nelL)+ulL(:,1:nelL),nelL,shldL,shlsL,nen,0,0,be);      
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
            [QxyR,shgsR,JdetR,bubbleR,xsR] = shgt(xlR(:,1:nelR)+ulR(:,1:nelR),nelR,shldR,shlsR,nen,0,0,be);  
           else
            [shlR,shldR,shlsR,be] = shlq(rR,sR,nelR,nelR,0,0);
            [PxyR,shgsR,JdetR,bubbleR,XsR] = shgq(xlR(:,1:nelR),nelR,shldR,shlsR,nen,0,0,be);     
            [QxyR,shgsR,JdetR,bubbleR,xsR] = shgq(xlR(:,1:nelR)+ulR(:,1:nelR),nelR,shldR,shlsR,nen,0,0,be);  
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
          
            
            % Cuchy stress tensors and material moduli for both sides
            [sigma2L, cmatL] = SigmaCmat2i(FL,JxXL,matepropL,lam);
            [sigma2Li, cmatLi] = SigmaCmat2i(eye(2),1,matepropL,lam);
            
            SmatL1=[sigma2L(1), sigma2L(3)
                    sigma2L(3), sigma2L(2)];           
            SmatL1i=[sigma2Li(1), sigma2Li(3)
                    sigma2Li(3), sigma2Li(2)];
            
            [sigma2R, cmatR] = SigmaCmat2i(FR,JxXR,matepropR,lam);
            [sigma2Ri, cmatRi] = SigmaCmat2i(eye(2),1,matepropR,lam);
                                 
            SmatR1=[sigma2R(1), sigma2R(3)
                    sigma2R(3), sigma2R(2)];     
            SmatR1i=[sigma2Ri(1), sigma2Ri(3)
                    sigma2Ri(3), sigma2Ri(2)];
             
             %Evaluate tangent and normal vectors
            T1L = [XsL(:,1); 0];
            [Tm1L, Tu1L] = VecNormalize(T1L);
            T2L = [0; 0; 1];
            Tm2L = 1;
            Tu2L = T2L';
            T3L = VecCrossProd(T1L,T2L);
            [Tm3L, Tu3L] = VecNormalize(T3L);
            C1L = drdr*Wgt*Tm3L;
            C1R = C1L;
            
            % Nanson formula
            t3L=JxXL*[fiL zeros(2,1); zeros(1,2) 1]'*T3L';
            [tm3L, tu3L] = VecNormalize(t3L);
            t3R=JxXR*[fiR zeros(2,1); zeros(1,2) 1]'*-T3L';
            [tm3R, tu3R] = VecNormalize(t3R);
            c1L = Wgt*drdr*tm3L;
            c1R = Wgt*drdr*tm3R;  
            %normal vectors
            nLx = tu3L(1);
            nLy = tu3L(2);
            nRx = tu3R(1);
            nRy = tu3R(2);                     
            nvecL = [nLx; nLy];
            nvecR = [nRx; nRy];
            
            tvtr = gamL*SmatL1*nvecL-gamR*SmatR1*nvecR; % weighted traction
            jumpu = ulL*shlL - ulR*shlR; % displacement jump - points from Left side to Right side when uL > uR
            
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
            sqr3 = sqrt(3);
            plist = [-sqr3 sqr3];
%         end
        
        for ll = 1:nelseg
            
            r = plist(1,ll);
            s = -1;
            shpS = sshp2d(r,s,4);
            
            tractions(:,ll) = [jumps; tractions_int]*shpS(1:2);
            
        end
        
        
        % assemble into interface output quantities
        InterQuant(1:6,ElemFlagI,step) = tractions;
        
end