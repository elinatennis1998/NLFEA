% 02/24/2014

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

% Implementation copied from Interface5 folder, with modifications to make
% it work in NLFEA.
% Verified against the results in Interface5. Also, the maximum
% displacement level obtained was fairly large, like step 176 for the
% transverse displacement unit cell

% Modified 3/3/2014 (special day...) - according to Alfano and Crisfield
% IJNME 2001, the damage tangent rather than the unloading tangent should
% be used to initialize the next load step. This provided absolutely
% convergent results for the quasistatic tension problem (step=360), at
% least for rp = 1000*Hc; for 200*Hc, it got to step 138.

% New version 3/5/14: penalty parameter computed from bubble functions
% according to IVMS, but the matrices are converted into scalars by
% diagonalizing and taking the maximum. The flux weights are still computed
% and not set to 1/2. For transverse tension, the response still went all
% the way to step 360.

% This element also works for implicit dynamics now, at least
% CohAxialDynDGDcomp.m at tstep = 1e-7!!!!!!!!!!!
% Since the material properties are equal, the differences would only come
% from the sizes of elements making gamL~=gamR and that the penalty
% parameter would vary from interface segment to interface segment.
% Mat: CohAxialDynDGcompBubble500step.mat

Bcol1 = [1; 4; 6];
Bcol2 = [2; 4; 5];
Bcol3 = [3; 5; 6];
col1 = [1; 2; 3];
col2 = [2; 1; 3];
col3 = [3; 2; 1];

diagt = 1; % make taus diagonal
scalt = 1; % make taus scalar
minmaxt = 1;
pencoeff = 1;  %VMS

switch isw %Task Switch
    
    case 1
        
        if ndf > 3
            
            for i = 4:ndf
                lie(i,1) = 0;
            end
            
            
        end
        
        nh1 = 7*4*1; % 8 variables, 4 int points, 1
        nh3 = 1*4*1; % dstate variable
        iste = 21; % number of stresses per element

    case {3,6} % Interface stiffness
        
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
        itchat = 7;100;6;12;8; %flag for iteration number to suppress chatter
        
% %         beta = 1;
% %         beta2 = beta^2;
% %         sigmax = 100;
% %         dc = 0.2;
        sigmax = mateprop(5);
        dc = mateprop(6);
        beta = mateprop(7);
        if length(mateprop) < 8
            phitype = 2;
        else
            phitype = mateprop(8);
        end
        
        Hc = sigmax/dc;
        rp = 1000*Hc;
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
        
        if inter == 1
            numdbond = 0;
        end
        
        % Determin bounds of integration segment
        xlintL = xlL;
        xlintR = xlR;
        
        if nelL == 4 || nelL == 10
        lint = 3;4; % face
        lintTau = 4; % volume
        elseif nelL == 8 || nelL == 27
        lint = 4;3; % face
        lintTau = 8; % volume
        end
        ll = 0; % Counter for history variables

        % Compute tau from FS for linear elasticity
         [tauL,intb] = Tau3d41(xlL,DmatL,lintTau,nelL); %[Y^(-1)]       
         [tauR,intb] = Tau3d41(xlR,DmatR,lintTau,nelR);

            % Modifications to tau
            if exist('diagt','var') && diagt == 1 % make tauL and tauR diagonal
                if exist('scalt','var') && scalt == 1 % use a scalar of tau times the identity
                    if exist('minmaxt','var') && minmaxt == -1 % minimum entry
                        tau = diag(inv(tauL));
                        tau = min(tau);%*ones(ndm,1);
                        tauL = inv(diag(tau));
                        tau = diag(inv(tauR));
                        tau = min(tau);%*ones(ndm,1);
                        tauR = inv(diag(tau));
                    else %if minmaxt == 1 % maximum entry
                        tau = diag(inv(tauL));
                        tau = max(tau);%*ones(ndm,1);
                        tauL = inv(diag(tau));
                        tau = diag(inv(tauR));
                        tau = max(tau);%*ones(ndm,1);
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
                        tau = min(tau);%*ones(ndm,1);
                        tauL = inv(diag(tau));
                        tau = diag(inv(tauR));
                        tau = min(tau);%*ones(ndm,1);
                        tauR = inv(diag(tau));
                    else %if minmaxt == 1 % maximum entry
                        tau = diag(inv(tauL));
                        tau = max(tau);%*ones(ndm,1);
                        tauL = inv(diag(tau));
                        tau = diag(inv(tauR));
                        tau = max(tau);%*ones(ndm,1);
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
        
        % Integrate bubble function over face of element
        
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
        
        % Compute fine-scale stabilization tensors
        edgeK = tauL*ebL^2 + tauR*ebR^2;
        gamL = ebL^2*(edgeK\tauL);
        gamR = ebR^2*(edgeK\tauR);
        ep = pencoeff*intedge*inv(edgeK);
        rp = ep;
        
            for l = 1:lint

                ll = ll + 1;

                % Evaluate  basis functions at integration points
                if nelL == 4 || nelL == 10
                  [Wgt,ss] = int3d_t(l,lint,5);
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
            chat1hr = nh2-1+(ll-1)*7+6;
            chat2hr = nh2-1+(ll-1)*7+7;
            dvec = [hr(damhr+1); hr(damhr+2); hr(damhr+3)];
            dmax = hr(dmaxhr);
            bonddbond1 = hr(chat1hr);
            chatter = hr(chat2hr);
            dstate = hr(nh3-1+(ll-1)*1+1);
                
            bnAdN1 = gamL*nvect*DmatL*BmatL;
            bnAdN2 = gamR*nvect*DmatR*BmatR;
            rhspulL = reshape(ulL,ndf*nelL,1);
            rhspulR = reshape(ulR,ndf*nelR,1);
            
%             sigmax = sigmaxc;%*(1-.125+.125*cos(2*pi*zint/20000));
            
            tvtr = (bnAdN1*rhspulL + bnAdN2*rhspulR);
            jumpu = NmatR*rhspulR - NmatL*rhspulL;
            
            if step == 83
                step;
            end
            
            
            % Evaluate debonding interface model
            elast_init = transient == 0 && initia;
%             [dvec,dddt,dmax,dstate,bonddbond2,dbond,chatter,tvec] = RadialCZME(tvtr,jumpu,dmax, ...
%                 dvec,nvec,sigmax,Hc,dc,rp,dstate,iter,elast_init,itchat,chatter,bonddbond1,dtol,ndm);
            [dvec,dddt,dmax,dstate,bonddbond2,dbond,chatter] = GeneralCZM(tvtr,jumpu,dmax,dvec, ...
    nvec,[sigmax*dc/2 sigmax dc beta],phitype,rp,dstate,iter,elast_init,itchat,chatter,bonddbond1,dtol,ndm);

            
            ElemFL = ElemFL - c1*( - NmatL'*(tvtr + rp*(jumpu - dvec)) + bnAdN1'*(jumpu - dvec));
            ElemFR = ElemFR - c1*( + NmatR'*(tvtr + rp*(jumpu - dvec)) + bnAdN2'*(jumpu - dvec));

            ElemKLL = ElemKLL - c1*NmatL'*bnAdN1 - c1*bnAdN1'*NmatL;
            ElemKLR = ElemKLR - c1*NmatL'*bnAdN2 + c1*bnAdN1'*NmatR;
            ElemKRL = ElemKRL + c1*NmatR'*bnAdN1 - c1*bnAdN2'*NmatL;
            ElemKRR = ElemKRR + c1*NmatR'*bnAdN2 + c1*bnAdN2'*NmatR;

            ElemKLL = ElemKLL + c1*rp*(NmatL'*NmatL);
            ElemKLR = ElemKLR - c1*rp*(NmatL'*NmatR);
            ElemKRL = ElemKRL - c1*rp*(NmatR'*NmatL);
            ElemKRR = ElemKRR + c1*rp*(NmatR'*NmatR);

            % Damage terms
            
            ElemKLL = ElemKLL - c1*(bnAdN1'-rp*NmatL')*dddt*(bnAdN1-rp*NmatL);
            ElemKLR = ElemKLR - c1*(bnAdN1'-rp*NmatL')*dddt*(bnAdN2+rp*NmatR);
            ElemKRL = ElemKRL - c1*(bnAdN2'+rp*NmatR')*dddt*(bnAdN1-rp*NmatL);
            ElemKRR = ElemKRR - c1*(bnAdN2'+rp*NmatR')*dddt*(bnAdN2+rp*NmatR);
                    
            % Store history
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
            hr(nh3-1+(ll-1)*1+1) = dstate;
            numdbond = numdbond + abs(dbond);
            
            end %lint
        if inter == numSI
            numdbond
        end
        
%         end %intt
ElemKLL;
            ElemK = [ElemKLL ElemKLR; ElemKRL ElemKRR];
            ElemF = [ElemFL; ElemFR];

            
    case 51 % Surface strain homogenization
        
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
                  [shlL,shldL,shls,be] = shltt([rL sL tL],nelL,nelL,0,0);
                  PxyL = shgtt(xlL(:,1:nelL),nelL,shldL,shls,nen,0,0,be);
                else
                  [shlL,shldL,shls,be] = shlb([rL sL tL],nelL,nelL,0,0);
                  PxyL = shgb(xlL(:,1:nelL),nelL,shldL,shls,nen,0,0,be);
                end
                QxyL = PxyL;

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
                jumpu = -NmatL*reshape(ulL,ndf*nelL,1) + NmatR*reshape(ulR,ndf*nelR,1);
                epsili = -1/2*(nvec*jumpu' + jumpu*nvec');
                epsil = [epsili(1,1); epsili(2,2); epsili(3,3); epsili(1,2); epsili(2,3); epsili(3,1)];
                
                ElemSS(1:6) = ElemSS(1:6) + c1*epsil;
            
            end %lint
        
        end %intt
%%
    case 60
        
        numhr = 4;
        ElemI = zeros(10,numhr);
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
            ElemI(:,ll) = [theta
                           dn
                           jumpu'*mvec
                           normsig
                           tvtr'*mvec
                           (tvtr + rp*(jumpu-dvec))'*nvec
                           (tvtr + rp*(jumpu-dvec))'*mvec
                           normtvec
                           dvec'*nvec
                           dvec'*mvec];
            
            end %lint
        
%         end %intt

    case 26 % element stresses
        
        ElemS = zeros(nestr,1);
        t1 = zeros(3,1);
        t2 = t1;
        t3 = t1;
        
        dtol = 1e-11;
        
        sigmax = mateprop(5);
        dc = mateprop(6);
        beta = mateprop(7);
        if length(mateprop) < 8
            phitype = 2;
        else
            phitype = mateprop(8);
        end
        
        Hc = sigmax/dc;
        rp = 1000*Hc;
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
        
        % Determin bounds of integration segment
        xlintL = xlL;
        xlintR = xlR;
        
        if nelL == 4 || nelL == 10
        lint = 3;4; % face
        lintTau = 4; % volume
        elseif nelL == 8 || nelL == 27
        lint = 4;3; % face
        lintTau = 8; % volume
        end
        ll = 0; % Counter for history variables

        % Compute tau from FS for linear elasticity
         [tauL,intb] = Tau3d41(xlL,DmatL,lintTau,nelL); %[Y^(-1)]       
         [tauR,intb] = Tau3d41(xlR,DmatR,lintTau,nelR);

            % Modifications to tau
            if exist('diagt','var') && diagt == 1 % make tauL and tauR diagonal
                if exist('scalt','var') && scalt == 1 % use a scalar of tau times the identity
                    if exist('minmaxt','var') && minmaxt == -1 % minimum entry
                        tau = diag(inv(tauL));
                        tau = min(tau);%*ones(ndm,1);
                        tauL = inv(diag(tau));
                        tau = diag(inv(tauR));
                        tau = min(tau);%*ones(ndm,1);
                        tauR = inv(diag(tau));
                    else %if minmaxt == 1 % maximum entry
                        tau = diag(inv(tauL));
                        tau = max(tau);%*ones(ndm,1);
                        tauL = inv(diag(tau));
                        tau = diag(inv(tauR));
                        tau = max(tau);%*ones(ndm,1);
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
                        tau = min(tau);%*ones(ndm,1);
                        tauL = inv(diag(tau));
                        tau = diag(inv(tauR));
                        tau = min(tau);%*ones(ndm,1);
                        tauR = inv(diag(tau));
                    else %if minmaxt == 1 % maximum entry
                        tau = diag(inv(tauL));
                        tau = max(tau);%*ones(ndm,1);
                        tauL = inv(diag(tau));
                        tau = diag(inv(tauR));
                        tau = max(tau);%*ones(ndm,1);
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
        
        % Integrate bubble function over face of element
        
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
        
        % Compute fine-scale stabilization tensors
        edgeK = tauL*ebL^2 + tauR*ebR^2;
        gamL = ebL^2*(edgeK\tauL);
        gamR = ebR^2*(edgeK\tauR);
        ep = pencoeff*intedge*inv(edgeK);
        rp = ep;
        
            for l = 1:1

                ll = ll + 1;

                % Evaluate  basis functions at integration points
                if nelL == 4 || nelL == 10
                  [Wgt,ss] = int3d_t(l,lint,5);
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
            chat1hr = nh2-1+(ll-1)*7+6;
            chat2hr = nh2-1+(ll-1)*7+7;
            dvec = [hr(damhr+1); hr(damhr+2); hr(damhr+3)];
            dmax = hr(dmaxhr);
            bonddbond1 = hr(chat1hr);
            chatter = hr(chat2hr);
            dstate = hr(nh3-1+(ll-1)*1+1);
                
            bnAdN1 = gamL*nvect*DmatL*BmatL;
            bnAdN2 = gamR*nvect*DmatR*BmatR;
            rhspulL = reshape(ulL,ndf*nelL,1);
            rhspulR = reshape(ulR,ndf*nelR,1);
            
            tvtr = (bnAdN1*rhspulL + bnAdN2*rhspulR);
            jumpu = NmatR*rhspulR - NmatL*rhspulL;
            
            % Evaluate debonding interface model
            elast_init = transient == 0 && initia;
%             [dvec,dddt,dmax,dstate,bonddbond2,dbond,chatter,tvec] = RadialCZME(tvtr,jumpu,dmax, ...
%                 dvec,nvec,sigmax,Hc,dc,rp,dstate,iter,elast_init,itchat,chatter,bonddbond1,dtol,ndm);
            [dvec,dddt,dmax,dstate,bonddbond2,dbond,chatter] = GeneralCZM(tvtr,jumpu,dmax,dvec, ...
    nvec,[sigmax*dc/2 sigmax dc beta],phitype,rp,dstate,iter,elast_init,itchat,chatter,bonddbond1,dtol,ndm);

%       output stuff
       normsig = nvec(1)*tvec(1) + nvec(2)*tvec(2);
       dn = nvec(1)*jumpu(1) + nvec(2)*jumpu(2);
       dv = nvec(1)*dvec(1) + nvec(2)*dvec(2);
       tn = nvec(1)*tvtr(1) + nvec(2)*tvtr(2);
       if (ndm == 3) %then
          normsig = normsig + nvec(3)*tvec(3);
          dn = dn + nvec(3)*jumpu(3);
       dv = dv + nvec(3)*dvec(3);
       tn = tn + nvec(3)*tvtr(3);
       end %if
            ElemS(1 ) = nvec(1);
            ElemS(2 ) = nvec(2);
            ElemS(3 ) = nvec(3);
            ElemS(4 ) = ep(1,1);
            ElemS(5 ) = jumpu(1);
            ElemS(6 ) = jumpu(2);
            ElemS(7 ) = jumpu(3);
            ElemS(8 ) = dn;
            ElemS(9 ) = dvec(1);
            ElemS(10) = dvec(2);
            ElemS(11) = dvec(3);
            ElemS(12) = dv;
            ElemS(13) = tvec(1);
            ElemS(14) = tvec(2);
            ElemS(15) = tvec(3);
            ElemS(16) = normsig;
            ElemS(17) = tvtr(1);
            ElemS(18) = tvtr(2);
            ElemS(19) = tvtr(3);
            ElemS(20) = tn;
            ElemS(21) = dmax;
            
            end %lint
        
end %Task Switch
