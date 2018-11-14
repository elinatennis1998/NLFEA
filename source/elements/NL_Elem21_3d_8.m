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
% 1/24/2014 - converted to triangular sector integration

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
surfacesi = SurfacesI(inter,9:10);
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
        
        nh3 = 9*3*4; %maximum of 4 triangles per DG element
          
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
        bnAdN2 = zeros(6,nstR);
        
        
        % Load Guass Integration Points
        lint = 3;7;13; %at least need 7 pts for quadratic elements

        nil = surfacesi(2) - surfacesi(1) + 1;

        for intt = 1:nil %integrate on left domain
                
            trinum = surfacesi(1) + intt-1;
            xit = zeros(ndm,3);
            for j = 1:3
                node = ixt(trinum,j);
                for i = 1:ndm
                    xit(i,j) = xintt(node,i);
                end
            end
        
            % Determin bounds of integration segment and sector
            InterConn3D
            
           iterset = 3;1;4;
           if iter  <=iterset % == 0 %   
       
            % Perform integration of various matrices
    % For separate bubble types on T and B
            if nelL == 8 || nelL == 27
         [tauL,intb] = TauS3_PCW(xlintL,xlL,ulL,matepropL,nelL,lint,2,lam); %[Y^(-1)] 
            else
            error('no tetrahedral tau yet')
            end
            if nelR == 8 || nelR == 27
         [tauR,intb] = TauS3_PCW(xlintR,xlR,ulR,matepropR,nelR,lint,2,lam);
            else
            error('no tetrahedral tau yet')
            end


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

            % THIS LOOP COULD BE EVALUATED IN CLOSED-FORM
            for ie = 1:lint


                % For separate bubble types on T and B
                if nelL == 4 || nelL == 10
                    
                elseif nelL == 8 || nelL == 27
                    [Wgt,litr,lits] =  intpntt(ie,lint,0);
                    [shl,shld] = shlt(litr,lits,3,3,0,0);
                    ebeL = facebubbleW([litr lits -1]);
                end

                if nelL == 4 || nelL == 10
                    
                elseif nelL == 8 || nelL == 27
                    [Wgt,litr,lits] =  intpntt(ie,lint,0);
                    ebeR = facebubbleW([litr lits -1]);
                end

                %Evaluate tangent and normal vectors
                t1 = xit*shld(:,1);%[xs(:,1); 0];
                [tm1, tu1] = VecNormalize(t1);
                t2 = xit*shld(:,2);%[xs(:,2); 0];
                [tm2, tu2] = VecNormalize(t2);
                t3 = VecCrossProd(t1,t2);
                [tm3, tu3] = VecNormalize(t3);

                c1 = Wgt*tm3;

                ebL = ebL + c1*ebeL;
                ebR = ebR + c1*ebeR;
                intedge = intedge + c1;

            end
        
            edgeK = tauL*ebL^2 + tauR*ebR^2;
            gamL = ebL^2*(edgeK\tauL);
            gamR = ebR^2*(edgeK\tauR);
            gamLhr = nh3-1+(intt-1)*27; %pointer for gamL at triangle intt
            gamRhr = nh3-1+(intt-1)*27+9;
            ephr = nh3-1+(intt-1)*27+18;
            hr(gamLhr+1:gamLhr+9) = [gamL(1,1),gamL(1,2),gamL(1,3),gamL(2,1),gamL(2,2),gamL(2,3),gamL(3,1),gamL(3,2),gamL(3,3)];
            hr(gamRhr+1:gamRhr+9) = [gamR(1,1),gamR(1,2),gamR(1,3),gamR(2,1),gamR(2,2),gamR(2,3),gamR(3,1),gamR(3,2),gamR(3,3)];
            ep = pencoeff*intedge*inv(edgeK); 
            hr(ephr+1:ephr+9) = [ep(1,1),ep(1,2),ep(1,3),ep(2,1),ep(2,2),ep(2,3),ep(3,1),ep(3,2),ep(3,3)];

           else
            gamLhr = nh3-1+(intt-1)*27; %pointer for gamL at triangle intt
            gamRhr = nh3-1+(intt-1)*27+9;
            ephr = nh3-1+(intt-1)*27+18;
            gamL = [hr(gamLhr+1) hr(gamLhr+2) hr(gamLhr+3)
                    hr(gamLhr+4) hr(gamLhr+5) hr(gamLhr+6)
                    hr(gamLhr+7) hr(gamLhr+8) hr(gamLhr+9)];
            gamR = [hr(gamRhr+1) hr(gamRhr+2) hr(gamRhr+3)
                    hr(gamRhr+4) hr(gamRhr+5) hr(gamRhr+6)
                    hr(gamRhr+7) hr(gamRhr+8) hr(gamRhr+9)];
            ep =   [hr(ephr+1) hr(ephr+2) hr(ephr+3)
                    hr(ephr+4) hr(ephr+5) hr(ephr+6)
                    hr(ephr+7) hr(ephr+8) hr(ephr+9)];
       end
%        gamL =0.5*eye(3,3) ;
%        gamR =0.5*eye(3,3) ; 
%        ep = 200*eye(3,3) ;   
   

         for l = 1:lint

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
                  [shlL,shldL,shls,be] = shltt([rL sL tL],nelL,nelL,0,0);
                  [PxyL,shgsL,JdetL,bubbleL,XsL] = shgtt(xlL(:,1:nelL),nelL,shldL,shls,nen,0,0,be);                 
                  [QxyL,shgsL,JdetL,bubbleL,xsL] = shgtt(xlL(:,1:nelL)+ulL(:,1:nelL),nelL,shldL,shls,nen,0,0,be);
                else
                  [shlL,shldL,shls,be] = shlb([rL sL tL],nelL,nelL,0,0);
                  [PxyL,shgsL,JdetL,bubbleL,XsL] = shgb(xlL(:,1:nelL),nelL,shldL,shls,nen,0,0,be); 
                  [QxyL,shgsL,JdetL,bubbleL,xsL] = shgb(xlL(:,1:nelL)+ulL(:,1:nelL),nelL,shldL,shls,nen,0,0,be);
                end

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
                t1 = xit*shld(:,1);%[xs(:,1); 0];
                [tm1, tu1] = VecNormalize(t1);
                t2 = xit*shld(:,2);%[xs(:,2); 0];
                [tm2, tu2] = VecNormalize(t2);
                t3 = VecCrossProd(t1,t2);
                [tm3, tu3] = VecNormalize(t3);
                Tm3L = tm3;
                Tm3R = tm3;
                Tu3L = tu3;
                Tu3R = -tu3;
                C1L = Wgt*Tm3L;
                C1R = Wgt*Tm3R;  
                % Nanson formula
                t3L=JxXL*fiL'*t3';
                [tm3L, tu3L] = VecNormalize(t3L);
                t3R=JxXR*fiR'*-t3';
                [tm3R, tu3R] = VecNormalize(t3R);
                c1L = Wgt*tm3L;
                c1R = Wgt*tm3R;

                % verify that unit vector points outward from L
                % elements are oriented so that nodes 1-4 are on bottom face of
                % the template element
                xcheck = xlL(:,1) + Tu3L';
                xi = POU_Coord3(xcheck(1),xcheck(2),xcheck(3),xlL,1,nelL);
                if ((nelL == 4 || nelL == 10) && xi(3) < 0) ...
                || ((nelL == 8 || nelL == 27) && xi(3) < -1)
                    nLx = tu3L(1);
                    nLy = tu3L(2);
                    nLz = tu3L(3);
                    NLx = Tu3L(1);
                    NLy = Tu3L(2);
                    NLz = Tu3L(3);
                else
                    nLx = -tu3L(1);
                    nLy = -tu3L(2);
                    nLz = -tu3L(3);
                    NLx = -Tu3L(1);
                    NLy = -Tu3L(2);
                    NLz = -Tu3L(3);
                end
                xcheck = xlR(:,1) + Tu3R';
                xi = POU_Coord3(xcheck(1),xcheck(2),xcheck(3),xlR,1,nelR);
                if ((nelR == 4 || nelR == 10) && xi(3) < 0) ...
                || ((nelR == 8 || nelR == 27) && xi(3) < -1)
                    nRx = tu3R(1);
                    nRy = tu3R(2);
                    nRz = tu3R(3);
                    NRx = Tu3R(1);
                    NRy = Tu3R(2);
                    NRz = Tu3R(3);
                else
                    nRx = -tu3R(1);
                    nRy = -tu3R(2);
                    nRz = -tu3R(3);
                    NRx = -Tu3R(1);
                    NRy = -Tu3R(2);
                    NRz = -Tu3R(3);
                end               
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
%             if jumpu>jumpu_max
%                 jumpu_max = jumpu;
%             end
            end %lint
            
        end %intt
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
        bnAdN2 = zeros(6,nstR);
        
        
        % Load Guass Integration Points
        lint = 7;13;3; %at least need 7 pts for quadratic elements

        nil = surfacesi(2) - surfacesi(1) + 1;

        for intt = 1:nil %integrate on left domain
                
            trinum = surfacesi(1) + intt-1;
            xit = zeros(ndm,3);
            for j = 1:3
                node = ixt(trinum,j);
                for i = 1:ndm
                    xit(i,j) = xintt(node,i);
                end
            end
        
            % Determin bounds of integration segment and sector
            InterConn3D
            
           iterset = 3;1;4;
           if iter  <=iterset % == 0 %   
       
            % Perform integration of various matrices
    % For separate bubble types on T and B
            if nelL == 8 || nelL == 27
         [tauL,intb] = TauS3_PCW(xlintL,xlL,ulL,matepropL,nelL,lint,2,lam); %[Y^(-1)] 
            else
            error('no tetrahedral tau yet')
            end
            if nelR == 8 || nelR == 27
         [tauR,intb] = TauS3_PCW(xlintR,xlR,ulR,matepropR,nelR,lint,2,lam);
            else
            error('no tetrahedral tau yet')
            end


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

            % THIS LOOP COULD BE EVALUATED IN CLOSED-FORM
            for ie = 1:lint


                % For separate bubble types on T and B
                if nelL == 4 || nelL == 10
                    
                elseif nelL == 8 || nelL == 27
                    [Wgt,litr,lits] =  intpntt(ie,lint,0);
                    [shl,shld] = shlt(litr,lits,3,3,0,0);
                    ebeL = facebubbleW([litr lits -1]);
                end

                if nelL == 4 || nelL == 10
                    
                elseif nelL == 8 || nelL == 27
                    [Wgt,litr,lits] =  intpntt(ie,lint,0);
                    ebeR = facebubbleW([litr lits -1]);
                end

                %Evaluate tangent and normal vectors
                t1 = xit*shld(:,1);%[xs(:,1); 0];
                [tm1, tu1] = VecNormalize(t1);
                t2 = xit*shld(:,2);%[xs(:,2); 0];
                [tm2, tu2] = VecNormalize(t2);
                t3 = VecCrossProd(t1,t2);
                [tm3, tu3] = VecNormalize(t3);

                c1 = Wgt*tm3;

                ebL = ebL + c1*ebeL;
                ebR = ebR + c1*ebeR;
                intedge = intedge + c1;

            end
        
            edgeK = tauL*ebL^2 + tauR*ebR^2;
            gamL = ebL^2*(edgeK\tauL);
            gamR = ebR^2*(edgeK\tauR);
            gamLhr = nh3-1+(intt-1)*27; %pointer for gamL at triangle intt
            gamRhr = nh3-1+(intt-1)*27+9;
            ephr = nh3-1+(intt-1)*27+18;
            hr(gamLhr+1:gamLhr+9) = [gamL(1,1),gamL(1,2),gamL(1,3),gamL(2,1),gamL(2,2),gamL(2,3),gamL(3,1),gamL(3,2),gamL(3,3)];
            hr(gamRhr+1:gamRhr+9) = [gamR(1,1),gamR(1,2),gamR(1,3),gamR(2,1),gamR(2,2),gamR(2,3),gamR(3,1),gamR(3,2),gamR(3,3)];
            ep = pencoeff*intedge*inv(edgeK); 
            hr(ephr+1:ephr+9) = [ep(1,1),ep(1,2),ep(1,3),ep(2,1),ep(2,2),ep(2,3),ep(3,1),ep(3,2),ep(3,3)];

           else
            gamLhr = nh3-1+(intt-1)*27; %pointer for gamL at triangle intt
            gamRhr = nh3-1+(intt-1)*27+9;
            ephr = nh3-1+(intt-1)*27+18;
            gamL = [hr(gamLhr+1) hr(gamLhr+2) hr(gamLhr+3)
                    hr(gamLhr+4) hr(gamLhr+5) hr(gamLhr+6)
                    hr(gamLhr+7) hr(gamLhr+8) hr(gamLhr+9)];
            gamR = [hr(gamRhr+1) hr(gamRhr+2) hr(gamRhr+3)
                    hr(gamRhr+4) hr(gamRhr+5) hr(gamRhr+6)
                    hr(gamRhr+7) hr(gamRhr+8) hr(gamRhr+9)];
            ep =   [hr(ephr+1) hr(ephr+2) hr(ephr+3)
                    hr(ephr+4) hr(ephr+5) hr(ephr+6)
                    hr(ephr+7) hr(ephr+8) hr(ephr+9)];
       end
%        gamL =0.5*eye(3,3) ;
%        gamR =0.5*eye(3,3) ; 
%        ep = 200*eye(3,3) ;   
   

         for l = 1:lint

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
                  [shlL,shldL,shls,be] = shltt([rL sL tL],nelL,nelL,0,0);
                  [PxyL,shgsL,JdetL,bubbleL,XsL] = shgtt(xlL(:,1:nelL),nelL,shldL,shls,nen,0,0,be);                 
                  [QxyL,shgsL,JdetL,bubbleL,xsL] = shgtt(xlL(:,1:nelL)+ulL(:,1:nelL),nelL,shldL,shls,nen,0,0,be);
                else
                  [shlL,shldL,shls,be] = shlb([rL sL tL],nelL,nelL,0,0);
                  [PxyL,shgsL,JdetL,bubbleL,XsL] = shgb(xlL(:,1:nelL),nelL,shldL,shls,nen,0,0,be); 
                  [QxyL,shgsL,JdetL,bubbleL,xsL] = shgb(xlL(:,1:nelL)+ulL(:,1:nelL),nelL,shldL,shls,nen,0,0,be);
                end

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
                t1 = xit*shld(:,1);%[xs(:,1); 0];
                [tm1, tu1] = VecNormalize(t1);
                t2 = xit*shld(:,2);%[xs(:,2); 0];
                [tm2, tu2] = VecNormalize(t2);
                t3 = VecCrossProd(t1,t2);
                [tm3, tu3] = VecNormalize(t3);
                Tm3L = tm3;
                Tm3R = tm3;
                Tu3L = tu3;
                Tu3R = -tu3;
                C1L = Wgt*Tm3L;
                C1R = Wgt*Tm3R;  
                % Nanson formula
                t3L=JxXL*fiL'*t3';
                [tm3L, tu3L] = VecNormalize(t3L);
                t3R=JxXR*fiR'*-t3';
                [tm3R, tu3R] = VecNormalize(t3R);
                c1L = Wgt*tm3L;
                c1R = Wgt*tm3R;

                % verify that unit vector points outward from L
                % elements are oriented so that nodes 1-4 are on bottom face of
                % the template element
                xcheck = xlL(:,1) + Tu3L';
                xi = POU_Coord3(xcheck(1),xcheck(2),xcheck(3),xlL,1,nelL);
                if ((nelL == 4 || nelL == 10) && xi(3) < 0) ...
                || ((nelL == 8 || nelL == 27) && xi(3) < -1)
                    nLx = tu3L(1);
                    nLy = tu3L(2);
                    nLz = tu3L(3);
                    NLx = Tu3L(1);
                    NLy = Tu3L(2);
                    NLz = Tu3L(3);
                else
                    nLx = -tu3L(1);
                    nLy = -tu3L(2);
                    nLz = -tu3L(3);
                    NLx = -Tu3L(1);
                    NLy = -Tu3L(2);
                    NLz = -Tu3L(3);
                end
                xcheck = xlR(:,1) + Tu3R';
                xi = POU_Coord3(xcheck(1),xcheck(2),xcheck(3),xlR,1,nelR);
                if ((nelR == 4 || nelR == 10) && xi(3) < 0) ...
                || ((nelR == 8 || nelR == 27) && xi(3) < -1)
                    nRx = tu3R(1);
                    nRy = tu3R(2);
                    nRz = tu3R(3);
                    NRx = Tu3R(1);
                    NRy = Tu3R(2);
                    NRz = Tu3R(3);
                else
                    nRx = -tu3R(1);
                    nRy = -tu3R(2);
                    nRz = -tu3R(3);
                    NRx = -Tu3R(1);
                    NRy = -Tu3R(2);
                    NRz = -Tu3R(3);
                end               
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
%             if jumpu>jumpu_max
%                 jumpu_max = jumpu;
%             end
            end %lint
            
        end %intt
%          jumpu_max = [jumpu_max; jumpu];
%         end %intt
% ElemKLL
            ElemF = [ElemFL; ElemFR];

    case -1
        
        ElemF = zeros(nst,1);
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
        
        ElemSS = zeros(13,1);

    case 52 % Surface strain homogenization
        
        ElemSS = zeros(13,1);
        
    case 61 % form data structure for interface segments

        % Set the number of interface quantities per node to be stored
        if ~exist('numIQ','var')
            numIQ = 9; % 3 disp-jump, 3 traction, 3 numerical-flux
        else
            numIQ = max(numIQ,9);
        end

        nil = surfacesi(2) - surfacesi(1) + 1; % number of triangles in interface segment
        
        IElemSeg(1,inter) = segment + 1;
        IElemSeg(2,inter) = segment + nil;

        for intt = 1:nil %integrate on left domain
                
            trinum = surfacesi(1) + intt-1;
            xit = zeros(ndm,3);
            for j = 1:3
                node = ixt(trinum,j);
                for i = 1:ndm
                    xit(i,j) = xintt(node,i);
                end
            end
            
            segment = segment + 1;
            Iix(segment,1:3) = (segnode+1:segnode+3);
            ICoordinates(segnode+1:segnode+3,1:ndm) = xit';
            segnode = segnode + 3;
            Iix(segment,nenseg+1) = ma; % good habit to copy material ID too
            
        end
        
        
    case 60 % output interface quantities for plotting
        
        lam = getlam(matepropL);
        
        lint = 3; %at least need 7 pts for quadratic elements

        nil = IElemSeg(2,inter) - IElemSeg(1,inter) + 1;

        for intt = 1:nil %integrate on left domain
               
            % get segment number for DG element
            segment = IElemSeg(1,inter) + intt - 1;
    
            % get nodes on the interface segment
            ElemFlagI = Iix(segment,1:nenseg);
            xlI = ICoordinates(ElemFlagI,1:ndm)';
            nelseg = 3;
            
        
            % compute quantities of interest
            % jump in displacement
            % interface traction, extrapolated from integration points
            
            % pull out weighting tensors and penalty parameter for this segment
            gamLhr = nh3-1+(intt-1)*27; %pointer for gamL at triangle intt
            gamRhr = nh3-1+(intt-1)*27+9;
            ephr = nh3-1+(intt-1)*27+18;
            gamL = [hr(gamLhr+1) hr(gamLhr+2) hr(gamLhr+3)
                    hr(gamLhr+4) hr(gamLhr+5) hr(gamLhr+6)
                    hr(gamLhr+7) hr(gamLhr+8) hr(gamLhr+9)];
            gamR = [hr(gamRhr+1) hr(gamRhr+2) hr(gamRhr+3)
                    hr(gamRhr+4) hr(gamRhr+5) hr(gamRhr+6)
                    hr(gamRhr+7) hr(gamRhr+8) hr(gamRhr+9)];
            ep =   [hr(ephr+1) hr(ephr+2) hr(ephr+3)
                    hr(ephr+4) hr(ephr+5) hr(ephr+6)
                    hr(ephr+7) hr(ephr+8) hr(ephr+9)];

            jumps = zeros(3,nelseg);
            tractions_int = zeros(9,nelseg); % tractions at integration points
            tractions = zeros(6,nelseg); % tractions at segment nodes
            for l = 1:nelseg

                %Integration point, weight, jacobian
                [Wgt,litr,lits] =  intpntt(l,lint,0);
                [shl,shld,shls,be] = shlt(litr,lits,3,3,0,0);

                %Physical location of int pt
                xint = xlI(1,:)*shl;
                yint = xlI(2,:)*shl;
                zint = xlI(3,:)*shl;

                xi = POU_Coord3(xint,yint,zint,xlL,1,nelL);
                rL = xi(1);
                sL = xi(2);
                tL = xi(3);

                % Evaluate  basis functions at integration points
                if nelL == 4 || nelL == 10
                  [shlL,shldL,shls,be] = shltt([rL sL tL],nelL,nelL,0,0);               
                  [QxyL,shgsL,JdetL,bubbleL,xsL] = shgtt(xlL(:,1:nelL)+ulL(:,1:nelL),nelL,shldL,shls,nen,0,0,be);
                else
                  [shlL,shldL,shls,be] = shlb([rL sL tL],nelL,nelL,0,0);
                  [QxyL,shgsL,JdetL,bubbleL,xsL] = shgb(xlL(:,1:nelL)+ulL(:,1:nelL),nelL,shldL,shls,nen,0,0,be);
                end

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
                %Evaluate tangent and normal vectors
                t1 = xlI*shld(:,1);%[xs(:,1); 0];
                [tm1, tu1] = VecNormalize(t1);
                t2 = xlI*shld(:,2);%[xs(:,2); 0];
                [tm2, tu2] = VecNormalize(t2);
                t3 = VecCrossProd(t1,t2);
                [tm3, tu3] = VecNormalize(t3);
                Tm3L = tm3;
                Tm3R = tm3;
                Tu3L = tu3;
                Tu3R = -tu3; 
                % Nanson formula
                t3L=JxXL*fiL'*t3';
                [tm3L, tu3L] = VecNormalize(t3L);
                t3R=JxXR*fiR'*-t3';
                [tm3R, tu3R] = VecNormalize(t3R);

                % verify that unit vector points outward from L
                % elements are oriented so that nodes 1-4 are on bottom face of
                % the template element
                xcheck = xlL(:,1) + Tu3L';
                xi = POU_Coord3(xcheck(1),xcheck(2),xcheck(3),xlL,1,nelL);
                if ((nelL == 4 || nelL == 10) && xi(3) < 0) ...
                || ((nelL == 8 || nelL == 27) && xi(3) < -1)
                    nLx = tu3L(1);
                    nLy = tu3L(2);
                    nLz = tu3L(3);
                    NLx = Tu3L(1);
                    NLy = Tu3L(2);
                    NLz = Tu3L(3);
                else
                    nLx = -tu3L(1);
                    nLy = -tu3L(2);
                    nLz = -tu3L(3);
                    NLx = -Tu3L(1);
                    NLy = -Tu3L(2);
                    NLz = -Tu3L(3);
                end
                xcheck = xlR(:,1) + Tu3R';
                xi = POU_Coord3(xcheck(1),xcheck(2),xcheck(3),xlR,1,nelR);
                if ((nelR == 4 || nelR == 10) && xi(3) < 0) ...
                || ((nelR == 8 || nelR == 27) && xi(3) < -1)
                    nRx = tu3R(1);
                    nRy = tu3R(2);
                    nRz = tu3R(3);
                    NRx = Tu3R(1);
                    NRy = Tu3R(2);
                    NRz = Tu3R(3);
                else
                    nRx = -tu3R(1);
                    nRy = -tu3R(2);
                    nRz = -tu3R(3);
                    NRx = -Tu3R(1);
                    NRy = -Tu3R(2);
                    NRz = -Tu3R(3);
                end
                nvecL = [nLx; nLy; nLz];
                nvecR = [nRx; nRy; nRz];

                tvtr = gamL*SmatL1*nvecL-gamR*SmatR1*nvecR; % weighted traction
                jumpu = ulL*shlL - ulR*shlR; % displacement jump - points from Left side to Right side when uL > uR

                % numerical flux - traction vector pointing out from Right side
                tractions_int(1:9,l) = [jumpu; -tvtr; -tvtr+ep*jumpu];

            end

%             % reorder integration point values from tensor-product to FEM-
%             % counter-clockwise so that normal shape functions can be used
%             tractions_int = tractions_int(:,[1 2 3]); % VERIFY THIS

            % extrapolate values to segment nodes
            plist = [-1/3 5/3 -1/3
                     -1/3 -1/3 5/3];

            for ll = 1:nelseg

                r = plist(1,ll);
                s = plist(2,ll);
                shpS = sshp2d(r,s,nelseg);

                jumps(:,ll) = tractions_int(1:3,:)*shpS;
                tractions(:,ll) = tractions_int(4:9,:)*shpS;

            end


            % assemble into interface output quantities
            InterQuant(1:9,ElemFlagI,step) = [jumps; tractions];
            
        end
        
end %Task Switch
