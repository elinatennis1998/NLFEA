
% 01/21/13
% Coulomb friction penalty
% desktop/laptop versions combined 2/7/13 where the interface prop section from the laptop
% was replaced over the outdated version in the desktop; nothing else from the desktop
% version, which was where the simulations were run, was modified.

if isw ~= 1
CGtoDGarrays
end

Bcol1 = [1; 4; 6];
Bcol2 = [2; 4; 5];
Bcol3 = [3; 5; 6];
col1 = [1; 2; 3];
col2 = [2; 1; 3];
col3 = [3; 2; 1];
diagt = 1; % make tauL and tauR diagonal
scalt = 1; % use a scalar of tau times the identity
minmaxt = -1;
pencoeff = 3;
nitvms = 1;

switch isw %Task Switch
    
    case 1
        
        if ndf > 3
            
            for i = 4:ndf
                lie(i,1) = 0;
            end
            
        end
        
        iste = 21;
%%        
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
        itchatn = 8;100;12; %flag for iteration number to suppress chatter
        
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
        
%         eN = 100*muL;%/(0.001^2);
%         eN = 30*muR/((xlR(:,1)-xlR(:,7))'*(xlR(:,1)-xlR(:,7)));
        if length(mateprop)>4
            eT = mateprop(5);
        else
            eT = 1;
        end
        
        NmatL = zeros(3,nstL);
        BmatL = zeros(6,nstL);
        bnAdN1 = zeros(6,nstL);
        N1 = zeros(3,nstL);
        NmatR = zeros(3,nstR);
        BmatR = zeros(6,nstR);
        bnAdN2 = zeros(6,nstR);
        N2 = zeros(3,nstR);
        
        if nelL == 4
            lint = 3;
        elseif nelL == 10
            lint = 13;7;
        else
            lint = 4;
        end
        ib = 5;
        

        % Compute the tau integrals
        [tauL,intbL,volL] = TauE1_3d(xlL,DmatL,nelL);
        [tauR,intbR,volR] = TauE1_3d(xlR,DmatR,nelR);


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
                [Wgt,ss] =  int3d_t(ie,lint,ib);
                ebeL = facebubbleT(ss,nelL);
                [shlL,shldL,shls,be] = shltt(ss,nelL,nelL,0,0);
           else
                [Wgt,ss] = intpntb(ie,lint,ib);
                [shlL,shldL,shls,be] = shlb(ss,nelL,nelL,0,0);
                ebeL = facebubbleQ(ss,nelL);
           end
                
           if nelR == 4 || nelR == 10
                [Wgt,ss] =  int3d_t(ie,lint,ib);
                ebeR = facebubbleT(ss,nelR);
           else
                [Wgt,ss] = intpntb(ie,lint,ib);
                ebeR = facebubbleQ(ss,nelR);
           end           
            
           if nelL == 4 || nelL == 10
                [Wgt,ss] =  int3d_t(ie,lint,ib);
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
        
        if nitvms == 1
        % VMS
        edgeK = tauL*ebL^2 + tauR*ebR^2;
        gamL = ebL^2*(edgeK\tauL);
        gamR = ebR^2*(edgeK\tauR);
        ep = pencoeff*intedge*inv(edgeK); %#ok<MINV>
        elseif nitvms == 2
        % Nitsche
        h = 2/(intedge/volR + intedge/volL);
        gamL = 1/2*eye(3);
        gamR = 1/2*eye(3);
        ep = 1000;pencoeff*eye(3)*max(muL,muR)/h;
        end
        eN = ep(1,1);
        
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
%             % Find node closest to center of triangle to make bubble func
%             if nelL == 4 || nelL == 10
%                 xitL = [xit xlL(:,4)];
%             else
%                 shl = shlt(1/3,1/3,3,3,0,0);
%                 xcp = xit*shl;
%                 maxdis = norm(xlL(:,5)-xcp,2);
%                 maxind = 5;
%                 for i = 6:8
%                     dis = norm(xlL(:,i)-xcp,2);
%                     if maxdis <= dis
%                         maxdis = dis;
%                         maxind = i;
%                     end
%                 end
%                 xitL = [xit xlL(:,maxind)];
%             end
%             
%             if nelR == 4 || nelR == 10
%                 xitR = [xit xlR(:,4)];
%             else
%                 shl = shlt(1/3,1/3,3,3,0,0);
%                 xcp = xit*shl;
%                 maxdis = norm(xlR(:,5)-xcp,2);
%                 maxind = 5;
%                 for i = 6:8
%                     dis = norm(xlR(:,i)-xcp,2);
%                     if maxdis <= dis
%                         maxdis = dis;
%                         maxind = i;
%                     end
%                 end
%                 xitR = [xit xlR(:,maxind)];
%             end
%             
%             % Compute taus
%             etauL = TauEE3d(xitL,DmatL,14);
%             etauR = TauEE3d(xitR,DmatR,14);
%             etau = etauL + etauR;
%     %         ep = tauR/tauL;
%             ep = 1;
%             eb = 0;
%             
%             for l = 1:lint
% 
%                 %Integration point, weight, jacobian
%                 [Wgt,litr,lits] =  intpntt(l,lint,0);
%                 [shl,shld,shls,be] = shlt(litr,lits,3,3,0,0);
%                 [shg, shgs, Jdet, be, xs] = shgt(xit,3,shld,shls,3,0,0,be);
% 
% %                 if Jdet < 0 % the orientation doesn't matter
% %                     inter
% %                     Jdet = -Jdet;
% %                 end
%                 
%                 bee = facebubble([litr,lits,0]);
%             
%                 %Physical location of int pt
%                 xint = xit(1,:)*shl;
%                 yint = xit(2,:)*shl;
%                 zint = xit(3,:)*shl;
% 
%                 xi = POU_Coord3(xint,yint,zint,xlL,1,nelL);
%                 rL = xi(1);
%                 sL = xi(2);
%                 tL = xi(3);
% 
%                 % Evaluate  basis functions at integration points
%                 if nelL == 4 || nelL == 10
%                   [shlL,shldL,shls,be] = shltt([rL sL tL],nelL,nel2L,0,0);
%                   PxyL = shgtt(xlL(:,1:nelL),nelL,shldL,shls,nen,0,0,be);
%                 else
%                   [shlL,shldL,shls,be] = shlb([rL sL tL],nelL,nel2L,0,0);
%                   PxyL = shgb(xlL(:,1:nelL),nelL,shldL,shls,nen,0,0,be);
%                 end
%                 QxyL = PxyL;
% 
%                 %Evaluate tangent and normal vectors
%                 t1 = [xs(:,1); 0];
%                 [tm1, tu1] = VecNormalize(t1);
%                 t2 = [xs(:,2); 0];
%                 [tm2, tu2] = VecNormalize(t2);
%                 t3 = VecCrossProd(t1,t2);
%                 [tm3, tu3] = VecNormalize(t3);
% %                 t1 = xs(1,:)';
% %                 [tm1, tu1] = VecNormalize(t1);
% %                 t2 = xs(2,:)';
% %                 [tm2, tu2] = VecNormalize(t2);
% %                 t3 = VecCrossProd(t1,t2);
% %                 [tm3, tu3] = VecNormalize(t3);
%                 nLx = tu3(1);
%                 nLy = tu3(2);
%                 nLz = tu3(3);
%                 nRx = -tu3(1);
%                 nRy = -tu3(2);
%                 nRz = -tu3(3);
%                 tLx = tu1(1);
%                 tLy = tu1(2);
%                 tLz = tu1(3);
%                 nvec = [nLx; nLy; nLz];
% 
%                 c1 = Wgt*tm3;
%                 
%                 eb = eb + c1*bee;
%                 
%             end
% 
%             Kinv = inv(etau)/eb^2;
%             Km = max([sum(Kinv(1,:)) sum(Kinv(1,:)) sum(Kinv(1,:))]);
% %             rp = max(Km,rp);
%             Ksinv = 1/(nvec'*etau*nvec)/eb^2;
        
            for l = 1:lint

                ll = ll + 1;

%                 %Integration point, weight, jacobian
%                 [Wgt,litr,lits] =  intpntt(l,lint,0);
%                 [shl,shld,shls,be] = shlt(litr,lits,3,3,0,0);
%                 [shg, shgs, Jdet, be, xs] = shgt(xit,3,shld,shls,3,0,0,be);
% 
% %                 if Jdet < 0 % the orientation doesn't matter
% %                     inter
% %                     Jdet = -Jdet;
% %                 end
%                 
%                 bee = facebubble([litr,lits,0]);
%             
%                 %Physical location of int pt
%                 xint = xit(1,:)*shl;
%                 yint = xit(2,:)*shl;
%                 zint = xit(3,:)*shl;
% 
%                 xi = POU_Coord3(xint,yint,zint,xlL,1,nelL);
%                 rL = xi(1);
%                 sL = xi(2);
%                 tL = xi(3);

                % Evaluate  basis functions at integration points
                if nelL == 4 || nelL == 10
                  [Wgt,ss] = int3d_t(l,lint,1);
                  [shlL,shldL,shls,be] = shltt(ss,nelL,nelL,0,0);
                  [PxyL,shgs,Jdet,bubble,xs] = shgtt(xlL(:,1:nelL),nelL,shldL,shls,nen,0,0,be);
                else
                  [Wgt,ss] = intpntb(l,lint,5);
                  [shlL,shldL,shls,be] = shlb(ss,nelL,nelL,0,0);
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
                  [shlR,shldR,shls,be] = shltt([rR sR tR],nelR,nelR,0,0);
                  PxyR = shgtt(xlR(:,1:nelR),nelR,shldR,shls,nen,0,0,be);
                else
                  [shlR,shldR,shls,be] = shlb([rR sR tR],nelR,nelR,0,0);
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
                nmat = nvec*nvec';

                c1 = Wgt*tm3;
                
            for i = 1:nelL
%                 NmatL(:,3*i-2:3*i) = shlL(i)*eye(3);
%                 BmatL(:,3*i-2:3*i) = [QxyL(i,1) 0 0 
%                                       0 QxyL(i,2) 0 
%                                       0 0 QxyL(i,3) 
%                                       QxyL(i,2) QxyL(i,1) 0 
%                                       0 QxyL(i,3) QxyL(i,2) 
%                                       QxyL(i,3) 0 QxyL(i,1) ];
                NmatL(1,(i-1)*ndf+1) = shlL(i);
                NmatL(2,(i-1)*ndf+2) = shlL(i);
                NmatL(3,(i-1)*ndf+3) = shlL(i);
                BmatL(Bcol1,(i-1)*ndf+1) = QxyL(i,col1);
                BmatL(Bcol2,(i-1)*ndf+2) = QxyL(i,col2);
                BmatL(Bcol3,(i-1)*ndf+3) = QxyL(i,col3);
            end
            
            for i = 1:nelR
%                 NmatR(:,3*i-2:3*i) = shlR(i)*eye(3);
%                 BmatR(:,3*i-2:3*i) = [QxyR(i,1) 0 0 
%                                       0 QxyR(i,2) 0 
%                                       0 0 QxyR(i,3) 
%                                       QxyR(i,2) QxyR(i,1) 0 
%                                       0 QxyR(i,3) QxyR(i,2) 
%                                       QxyR(i,3) 0 QxyR(i,1)];
                NmatR(1,(i-1)*ndf+1) = shlR(i);
                NmatR(2,(i-1)*ndf+2) = shlR(i);
                NmatR(3,(i-1)*ndf+3) = shlR(i);
                BmatR(Bcol1,(i-1)*ndf+1) = QxyR(i,col1);
                BmatR(Bcol2,(i-1)*ndf+2) = QxyR(i,col2);
                BmatR(Bcol3,(i-1)*ndf+3) = QxyR(i,col3);
            end
                
            bnAdN1 = nvect*DmatL*BmatL;
            bnAdN2 = nvect*DmatR*BmatR;
            rhspulL = reshape(ulL,ndf*nelL,1);
            rhspulR = reshape(ulR,ndf*nelR,1);
            
            tvtr = 1/2*(bnAdN1*rhspulL + bnAdN2*rhspulR); %average stress
            sn = nvec'*tvtr;                              %normal stress
            tss = tvtr - sn*nvec;                         %shear stress
            jumpu = NmatR*rhspulR - NmatL*rhspulL;        %displacement jump
            gn = jumpu'*nvec;                             %normal gap
            tN_k = sn + eN*gn;                            %normal traction
            gt = jumpu - gn*nvec;                         %tangential gap
            
            tT_tr = eT*gt;      %Elastic trial force
            
                tT_k = tT_tr;

                ElemFL = ElemFL - c1*( - NmatL'*((tN_k)*nvec) + 1/2*bnAdN1'*gn*nvec);
                ElemFR = ElemFR - c1*( + NmatR'*((tN_k)*nvec) + 1/2*bnAdN2'*gn*nvec);
                
                ElemFL = ElemFL - c1*( - NmatL'*(tT_k));
                ElemFR = ElemFR - c1*( + NmatR'*(tT_k));

                ElemKLL = ElemKLL - c1/2*NmatL'*nmat*bnAdN1 - c1/2*bnAdN1'*nmat*NmatL;
                ElemKLR = ElemKLR - c1/2*NmatL'*nmat*bnAdN2 + c1/2*bnAdN1'*nmat*NmatR;
                ElemKRL = ElemKRL + c1/2*NmatR'*nmat*bnAdN1 - c1/2*bnAdN2'*nmat*NmatL;
                ElemKRR = ElemKRR + c1/2*NmatR'*nmat*bnAdN2 + c1/2*bnAdN2'*nmat*NmatR;

                ElemKLL = ElemKLL + c1*eN*(NmatL'*nmat*NmatL);
                ElemKLR = ElemKLR - c1*eN*(NmatL'*nmat*NmatR);
                ElemKRL = ElemKRL - c1*eN*(NmatR'*nmat*NmatL);
                ElemKRR = ElemKRR + c1*eN*(NmatR'*nmat*NmatR);

                    dddt = eye(3) - nmat;

                    ElemKLL = ElemKLL + c1*eT*NmatL'*dddt*NmatL;
                    ElemKLR = ElemKLR - c1*eT*NmatL'*dddt*NmatR;
                    ElemKRL = ElemKRL - c1*eT*NmatR'*dddt*NmatL;
                    ElemKRR = ElemKRR + c1*eT*NmatR'*dddt*NmatR;
            
            
            
            end %lint
        
%         end %intt
ElemFL;
            ElemK = [ElemKLL ElemKLR; ElemKRL ElemKRR];
            ElemF = [ElemFL; ElemFR];
    
    case 51 % Volume stress/strain homogenization
        
        ElemSS = zeros(13,1);

        %Set integration number
        lint = IntPoint3(nel);
        ib = 0;
        bf = 0;
        der = 0;

        lam = Patchv*PatchE/((1+Patchv)*(1-2*Patchv));
        mu = PatchE/(2*(1+Patchv));
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
        
        if nelL == 4
            numhr = 3;
        else
            numhr = 4;
        end
        ElemI = zeros(14,numhr);
        t1 = zeros(3,1);
        t2 = t1;
        t3 = t1;
        
        dtol = 1e-11;
        eN = 30*muL/((xlR(:,1)-xlR(:,7))'*(xlR(:,1)-xlR(:,7)));
        
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
        
        if nelL == 4
            lint = 3;
        else
            lint = 4;
        end
        ll = 0; % Counter for history variables
        
            for l = 1:lint

                ll = ll + 1;

                % Evaluate  basis functions at integration points
                if nelL == 4 || nelL == 10
                  [Wgt,ss] = int3d_t(l,lint,1);
                  [shlL,shldL,shls,be] = shltt(ss,nelL,nel2L,0,0);
                  [PxyL,shgs,Jdet,bubble,xs] = shgtt(xlL(:,1:nelL),nelL,shldL,shls,nen,0,0,be);
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
            gshr = nh1-1+(ll-1)*18; %%%%%% THINK MORE ABOUT THIS
            tthr = nh1-1+(ll-1)*18+3;
            gthr = nh1-1+(ll-1)*18+6;
            tshr = nh1-1+(ll-1)*18+9;
            dbonhr = nh1-1+(ll-1)*18+12;
            dbond = hr(dbonhr);
            chat1nhr = nh2-1+(ll-1)*18+13;
            chat2nhr = nh2-1+(ll-1)*18+14;
            bonddbond1n = hr(chat1nhr);
            chattern = hr(chat2nhr);
            gtn_1 = [hr(gthr+1); hr(gthr+2); hr(gthr+3)];
            tT_n = [hr(tthr+1); hr(tthr+2); hr(tthr+3)];
                
            bnAdN1 = nvect*DmatL*BmatL;
            bnAdN2 = nvect*DmatR*BmatR;
            rhspulL = reshape(ulL,ndf*nelL,1);
            rhspulR = reshape(ulR,ndf*nelR,1);
            
            tvtr = 1/2*(bnAdN1*rhspulL + bnAdN2*rhspulR); %average stress
            sn = nvec'*tvtr;                              %normal stress
            tss = tvtr - sn*nvec;                         %shear stress
            jumpu = NmatR*rhspulR - NmatL*rhspulL;        %displacement jump
            gn = jumpu'*nvec;                             %normal gap
            tN_k = sn + eN*gn;                            %normal traction
            dddt = eye(3) - nmat;
            gt = jumpu - gn*nvec;                         %tangential gap
            
            if (dbond == 0 && (tN_k <= 0)) || (dbond == 1 && gn <= 0) || (iter > itchatn && chattern == 1) % compresion only
           
                
            ElemI(:,ll) = [gn
                           gt
                           sqrt(gt'*gt)
                           sn
                           tN_k
                           tT_n
                           sqrt(tT_n'*tT_n)
                           xint
                           yint
                           zint];
            
            else
                
                
            end
            
            end %lint
        
%         end %intt
%%
    case 26 % traction
        
        ElemS = zeros(nestr,1);
        
        t1 = zeros(3,1);
        t2 = t1;
        t3 = t1;
        
        dtol = 1e-11;
        itchatn = 8;100;12; %flag for iteration number to suppress chatter
        
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
        
%         eN = 100*muL;%/(0.001^2);
%         eN = 30*muR/((xlR(:,1)-xlR(:,7))'*(xlR(:,1)-xlR(:,7)));
        if length(mateprop)>4
            eT = mateprop(5);
        else
            eT = 1;
        end
        
        NmatL = zeros(3,nstL);
        BmatL = zeros(6,nstL);
        bnAdN1 = zeros(6,nstL);
        N1 = zeros(3,nstL);
        NmatR = zeros(3,nstR);
        BmatR = zeros(6,nstR);
        bnAdN2 = zeros(6,nstR);
        N2 = zeros(3,nstR);
        
        
        if nelL == 4
            lint = 3;
        elseif nelL == 10
            lint = 13;7;
        else
            lint = 4;
        end

        % Compute the tau integrals
        [tauL,intbL,volL] = TauE1_3d(xlL,DmatL,nelL);
        [tauR,intbR,volR] = TauE1_3d(xlR,DmatR,nelR);


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
        ib = 5;
        
        for ie = 1:lint            
% For separate bubble types on T and Q
           if nelL == 4 || nelL == 10
                [Wgt,ss] =  int3d_t(ie,lint,ib);
                ebeL = facebubbleT(ss,nelL);
                [shlL,shldL,shls,be] = shltt(ss,nelL,nelL,0,0);
           else
                [Wgt,ss] = intpntb(ie,lint,ib);
                [shlL,shldL,shls,be] = shlb(ss,nelL,nelL,0,0);
                ebeL = facebubbleQ(ss,nelL);
           end
                
           if nelR == 4 || nelR == 10
                [Wgt,ss] =  int3d_t(ie,lint,ib);
                ebeR = facebubbleT(ss,nelR);
           else
                [Wgt,ss] = intpntb(ie,lint,ib);
                ebeR = facebubbleQ(ss,nelR);
           end           
            
           if nelL == 4 || nelL == 10
                [Wgt,ss] =  int3d_t(ie,lint,ib);
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
        
        if nitvms == 1
        % VMS
        edgeK = tauL*ebL^2 + tauR*ebR^2;
        gamL = ebL^2*(edgeK\tauL);
        gamR = ebR^2*(edgeK\tauR);
        ep = pencoeff*intedge*inv(edgeK); %#ok<MINV>
        elseif nitvms == 2
        % Nitsche
        h = 2/(intedge/volR + intedge/volL);
        gamL = 1/2*eye(3);
        gamR = 1/2*eye(3);
        ep = 1000;pencoeff*eye(3)*max(muL,muR)/h;
        end
        eN = ep(1,1);
        
            lint = 1;
        ll = 0; % Counter for history variables

        
            for l = 1:lint

                ll = ll + 1;


                % Evaluate  basis functions at integration points
                if nelL == 4 || nelL == 10
                  [Wgt,ss] = int3d_t(l,lint,1);
                  [shlL,shldL,shls,be] = shltt(ss,nelL,nelL,0,0);
                  [PxyL,shgs,Jdet,bubble,xs] = shgtt(xlL(:,1:nelL),nelL,shldL,shls,nen,0,0,be);
                else
                  [Wgt,ss] = intpntb(l,lint,5);
                  [shlL,shldL,shls,be] = shlb(ss,nelL,nelL,0,0);
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
                  [shlR,shldR,shls,be] = shltt([rR sR tR],nelR,nelR,0,0);
                  PxyR = shgtt(xlR(:,1:nelR),nelR,shldR,shls,nen,0,0,be);
                else
                  [shlR,shldR,shls,be] = shlb([rR sR tR],nelR,nelR,0,0);
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
                nmat = nvec*nvec';

                c1 = Wgt*tm3;
                
            for i = 1:nelL
%                 NmatL(:,3*i-2:3*i) = shlL(i)*eye(3);
%                 BmatL(:,3*i-2:3*i) = [QxyL(i,1) 0 0 
%                                       0 QxyL(i,2) 0 
%                                       0 0 QxyL(i,3) 
%                                       QxyL(i,2) QxyL(i,1) 0 
%                                       0 QxyL(i,3) QxyL(i,2) 
%                                       QxyL(i,3) 0 QxyL(i,1) ];
                NmatL(1,(i-1)*ndf+1) = shlL(i);
                NmatL(2,(i-1)*ndf+2) = shlL(i);
                NmatL(3,(i-1)*ndf+3) = shlL(i);
                BmatL(Bcol1,(i-1)*ndf+1) = QxyL(i,col1);
                BmatL(Bcol2,(i-1)*ndf+2) = QxyL(i,col2);
                BmatL(Bcol3,(i-1)*ndf+3) = QxyL(i,col3);
            end
            
            for i = 1:nelR
%                 NmatR(:,3*i-2:3*i) = shlR(i)*eye(3);
%                 BmatR(:,3*i-2:3*i) = [QxyR(i,1) 0 0 
%                                       0 QxyR(i,2) 0 
%                                       0 0 QxyR(i,3) 
%                                       QxyR(i,2) QxyR(i,1) 0 
%                                       0 QxyR(i,3) QxyR(i,2) 
%                                       QxyR(i,3) 0 QxyR(i,1)];
                NmatR(1,(i-1)*ndf+1) = shlR(i);
                NmatR(2,(i-1)*ndf+2) = shlR(i);
                NmatR(3,(i-1)*ndf+3) = shlR(i);
                BmatR(Bcol1,(i-1)*ndf+1) = QxyR(i,col1);
                BmatR(Bcol2,(i-1)*ndf+2) = QxyR(i,col2);
                BmatR(Bcol3,(i-1)*ndf+3) = QxyR(i,col3);
            end
                
            bnAdN1 = nvect*DmatL*BmatL;
            bnAdN2 = nvect*DmatR*BmatR;
            rhspulL = reshape(ulL,ndf*nelL,1);
            rhspulR = reshape(ulR,ndf*nelR,1);
            
            tvtr = 1/2*(bnAdN1*rhspulL + bnAdN2*rhspulR); %average stress
            sn = nvec'*tvtr;                              %normal stress
            tss = tvtr - sn*nvec;                         %shear stress
            jumpu = NmatR*rhspulR - NmatL*rhspulL;        %displacement jump
            gn = jumpu'*nvec;                             %normal gap
            tN_k = sn + eN*gn;                            %normal traction
            gt = jumpu - gn*nvec;                         %tangential gap
            
            tT_tr = eT*gt;      %Elastic trial force
            
                tT_k = tT_tr+tN_k*nvec;
            
            
            % Force terms
            ElemS(13:16) = [tvtr; sn];
            ElemS(9:12) = [tT_k; tN_k];
            ElemS(5:8) = [jumpu; gn];
            
            end %lint
        
%         end %intt
        
end %Task Switch
