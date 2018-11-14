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

% 06/24/12
% New version to use unequal tension/shear

Bcol1 = [1; 4; 6];
Bcol2 = [2; 4; 5];
Bcol3 = [3; 5; 6];
col1 = [1; 2; 3];
col2 = [2; 1; 3];
col3 = [3; 2; 1];

switch isw %Task Switch
    
    case 1
        
        if ndf > 3
            
            for i = 4:ndf
                lie(i,1) = 0;
            end
            
            
        end
        
        nh1 = 7*4*1; % 5 variables, 3 int points, 2 triangles
        iteri = 0; % set flag for initial elastic step at start of simulation
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
        
%         beta = 0.707;sqrt(1/2);
%         beta2 = beta^2;
%         sigmax = 100;
%         dc = 0.2;
        sigmax = mateprop(5);
        dc = mateprop(6);
        beta = mateprop(7);
        beta2 = beta^2;
        
        Hc = sigmax/dc;
        rp = 2000*Hc;
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
        
        if iter == 0 && iteri == 1 % two passes happen for iter=0, so only do prediction for the first one
            iteri = 0;
        elseif iter == 0 && iteri == 0
            iteri = 1;
        end
        
        if inter == 1
            numdbond = 0;
        end
        
        lint = 4;3;
        ll = 0; % Counter for history variables
        
            for l = 1:lint

                ll = ll + 1;

                % Evaluate  basis functions at integration points
                if nelL == 4 || nelL == 10
                  [shlL,shldL,shls,be] = shltt([rL sL tL],nelL,nelL,0,0);
                  PxyL = shgtt(xlL(:,1:nelL),nelL,shldL,shls,nen,0,0,be);
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
                         0 0 nLz 0 nLy nLx];
                nvec = [nLx; nLy; nLz];

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
            
            % Load history
            damhr = nh1-1+(ll-1)*7;
            dmaxhr = nh1-1+(ll-1)*7+4;
            chat1hr = nh2-1+(ll-1)*7+6;
            chat2hr = nh2-1+(ll-1)*7+7;
            dvec = [hr(damhr+1); hr(damhr+2); hr(damhr+3)];
            dmax = hr(dmaxhr);
            bonddbond1 = hr(chat1hr);
            chatter = hr(chat2hr);
                
            bnAdN1 = nvect*DmatL*BmatL;
            bnAdN2 = nvect*DmatR*BmatR;
            rhspulL = reshape(ulL,ndf*nelL,1);
            rhspulR = reshape(ulR,ndf*nelR,1);
            
%             sigmax = sigmaxc;%*(1-.125+.125*cos(2*pi*zint/20000));
            
            tvtr = 1/2*(bnAdN1*rhspulL + bnAdN2*rhspulR); %average stress
            sn = nvec'*tvtr;                              %normal stress
            jumpu = NmatR*rhspulR - NmatL*rhspulL;        %displacement jump
            un = jumpu'*nvec;                             %normal gap
            tn = sn+rp*un;                                %normal traction
            tss = tvtr + rp*jumpu - tn*nvec;              %shear traction v
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
                
                [dvec,dtdd,deq] = damageNR(tvec,nvec,dvec,normtvec,sigmax,dc,rp,beta2);
                if deq >= dc % opening/sliding
                    dvec = tvec/rp;
                    dddt = 1/rp*eye(3);
                    dn = dvec'*nvec;
                    dss = dvec - dn*nvec;
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
                dn = dvec'*nvec;
                dss = dvec - dn*nvec;
                dt = sqrt(dss'*dss);
                deq = sqrt(dn^2+beta2*dt^2);
                dmax = max(dmax,deq);
                bonddbond2 = -1;
                
              elseif dmax > 0 && dmax <= dc && iteri == 0 %elastic prediction
                
                [dvec,dddt,deq] = unloadingNR(tn,tt,tss,nvec,dmax,sigmax,dc,rp,beta2);
                bonddbond2 = -1;
                
              elseif dmax > 0 && rp*dmax <= tb && (tn^2/(rp*dmax + psik)^2 + beta2*tt^2/(rp*dmax + beta2*psik)^2 - 1) <= dtol %unloading
                
                [dvec,dddt,deq] = unloadingNR(tn,tt,tss,nvec,dmax,sigmax,dc,rp,beta2);
                bonddbond2 = -1;
                
              elseif (tn^2/(rp*dmax + psik)^2 + beta2*tt^2/(rp*dmax + beta2*psik) - 1)^2 > dtol % damage
                
                [dvec,dtdd,deq] = damageNR(tvec,nvec,dvec,normtvec,sigmax,dc,rp,beta2);
                if deq >= dc % opening/sliding
                    dvec = tvec/rp;
                    dddt = 1/rp*eye(3);
                    dn = dvec'*nvec;
                    dss = dvec - dn*nvec;
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
                dddt = 1/rp*(eye(3) - nvec*nvec');
                deq = beta*sqrt(dvec'*dvec);
                dmax = max(dmax,deq);
                bonddbond2 = -1;
                
              elseif dmax > 0 && dmax <= dc && iteri == 0 %elastic prediction
                
                ntilda = tss/tt;
                dvec = (dmax/beta)*ntilda;
                dddt = (dmax/beta)/tt*(eye(3) - ntilda*ntilda' - nvec*nvec');
%                 deq = sqrt(dvec'*dvec);
%                 dmax = max(dmax,deq);
                bonddbond2 = -1;
                
              elseif dmax > 0 && rp*dmax <= tb && (tb - (rp*dmax + beta2*psik)) <= dtol %unloading
                
                ntilda = tss/tt;
                dvec = (dmax/beta)*ntilda;
                dddt = (dmax/beta)/tt*(eye(3) - ntilda*ntilda' - nvec*nvec');
%                 deq = sqrt(dvec'*dvec);
%                 dmax = max(dmax,deq);
                bonddbond2 = -1;
                
              elseif dtol < (tb - (rp*dmax + beta2*psik)) % damage
                
                deq = beta*(tt - beta*sigmax)/(rp - beta2*Hc); % assume only damage
                ntilda = tss/tt;
                if deq >= dc % opening/sliding
                    dvec = tvec/rp;
                    deq = beta*sqrt(dvec'*dvec);
                    dddt = 1/rp*(eye(3) - nvec*nvec');
                else % only damage
                    dvec = (deq/beta)*ntilda;
                    dddt = 1/(1*(rp - beta2*Hc))*(ntilda*ntilda') + (deq/beta)/tt*(eye(3) - ntilda*ntilda' - nvec*nvec');
                end
                dmax = deq;
                bonddbond2 = -1;
                
              end
                
            end
            
            ElemFL = ElemFL - c1*( - NmatL'*(tvtr + rp*(jumpu - dvec)) + 1/2*bnAdN1'*(jumpu - dvec));
            ElemFR = ElemFR - c1*( + NmatR'*(tvtr + rp*(jumpu - dvec)) + 1/2*bnAdN2'*(jumpu - dvec));

            ElemKLL = ElemKLL - c1/2*NmatL'*bnAdN1 - c1/2*bnAdN1'*NmatL;
            ElemKLR = ElemKLR - c1/2*NmatL'*bnAdN2 + c1/2*bnAdN1'*NmatR;
            ElemKRL = ElemKRL + c1/2*NmatR'*bnAdN1 - c1/2*bnAdN2'*NmatL;
            ElemKRR = ElemKRR + c1/2*NmatR'*bnAdN2 + c1/2*bnAdN2'*NmatR;

            ElemKLL = ElemKLL + c1*rp*(NmatL'*NmatL);
            ElemKLR = ElemKLR - c1*rp*(NmatL'*NmatR);
            ElemKRL = ElemKRL - c1*rp*(NmatR'*NmatL);
            ElemKRR = ElemKRR + c1*rp*(NmatR'*NmatR);

            % Damage terms
            
%             dddu = rp*dddt;
%             ElemKLL = ElemKLL - c1*rp*(NmatL'*dddu*NmatL);
%             ElemKLR = ElemKLR + c1*rp*(NmatL'*dddu*NmatR);
%             ElemKRL = ElemKRL + c1*rp*(NmatR'*dddu*NmatL);
%             ElemKRR = ElemKRR - c1*rp*(NmatR'*dddu*NmatR);
%             
%             ElemKLL = ElemKLL + c1/2*rp*NmatL'*dddt*bnAdN1 + c1/2*bnAdN1'*dddu*NmatL;
%             ElemKLR = ElemKLR + c1/2*rp*NmatL'*dddt*bnAdN2 - c1/2*bnAdN1'*dddu*NmatR;
%             ElemKRL = ElemKRL - c1/2*rp*NmatR'*dddt*bnAdN1 + c1/2*bnAdN2'*dddu*NmatL;
%             ElemKRR = ElemKRR - c1/2*rp*NmatR'*dddt*bnAdN2 - c1/2*bnAdN2'*dddu*NmatR;
%             
%             ElemKLL = ElemKLL - c1/4*bnAdN1'*dddt*bnAdN1;
%             ElemKLR = ElemKLR - c1/4*bnAdN1'*dddt*bnAdN2;
%             ElemKRL = ElemKRL - c1/4*bnAdN2'*dddt*bnAdN1;
%             ElemKRR = ElemKRR - c1/4*bnAdN2'*dddt*bnAdN2;
            
            ElemKLL = ElemKLL - c1*(1/2*bnAdN1'-rp*NmatL')*dddt*(1/2*bnAdN1-rp*NmatL);
            ElemKLR = ElemKLR - c1*(1/2*bnAdN1'-rp*NmatL')*dddt*(1/2*bnAdN2+rp*NmatR);
            ElemKRL = ElemKRL - c1*(1/2*bnAdN2'+rp*NmatR')*dddt*(1/2*bnAdN1-rp*NmatL);
            ElemKRR = ElemKRR - c1*(1/2*bnAdN2'+rp*NmatR')*dddt*(1/2*bnAdN2+rp*NmatR);
                    
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
%             if iter <= itchat
%               if sign(bonddbond1) ~= sign(bonddbond2)
%                 chatter = 1;
%               else
%                 chatter = 0;
%               end
%             end
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
%             hr(chat1hr) = bonddbond2;
%             hr(chat2hr) = chatter;
            numdbond = numdbond + abs(dbond);
            
            end %lint
        
        if inter == numSI
            numdbond
        end
%         end %intt
ElemFL;
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
            sn = nvec'*tvtr;
            jumpu = NmatR*rhspulR - NmatL*rhspulL;
            un = jumpu'*nvec;
            tn = sn+rp*un;
            
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
%                            un
%                            jumpu'*mvec
%                            sn
%                            tvtr'*mvec
%                            tn
%                            (tvtr + rp*jumpu)'*mvec
%                            normtvec
%                            dvec'*nvec
%                            dvec'*mvec];
            ElemI(:,ll) = [theta
                           un
                           jumpu'*mvec
                           sn
                           tvtr'*mvec
                           (tvtr + rp*(jumpu-dvec))'*nvec
                           (tvtr + rp*(jumpu-dvec))'*mvec
                           normtvec
                           dvec'*nvec
                           dvec'*mvec];
            
            end %lint
        
%         end %intt
        
end %Task Switch
