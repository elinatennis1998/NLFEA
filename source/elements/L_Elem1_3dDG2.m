% 10/14/13

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


% Tim Truster
% 10/15/13
% DG for elasticity problem in 3D with nonconforming interface
% Hard-coded version which was benchmarked against the L_Elem3_3d.m in
% Interface3D folder. This element assumes the interface has elements
% oriented along a plane with constant z coordinate and having faces 5 & 6
% adjoining the interface. Penalty parameter and flux weights are
% constants.
% The DG element input should simply be ix = [ixL ixR maDG] with ixL and
% ixR the connectivity for the L and R elements without reordering (i.e.
% the bottom face doesn't have to point to the interface)

nitvms = 1;
if nitvms == 1 %VMS
pencoeff = 1;
elseif nitvms == 2 %Nitsche
pencoeff = 1;4;2;
else %RFB
pencoeff = 1;4;2;
end

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
        
        bf1 = 0;
        bf2 = 0;
        
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
        
        ulresL = reshape(ulL,ndf*nelL,1);
        ulresR = reshape(ulR,ndf*nelR,1);
        
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
            
            % Perform integration of various matrices
    % For separate bubble types on T and B
            if nelL == 8 || nelL == 27
            [tauL,intbL,volL] = TauE1_3dW(xlintL,DmatL,nelL,lint,2);
            else
            error('no tetrahedral tau yet')
            end
            if nelR == 8 || nelR == 27
            [tauR,intbR,volR] = TauE1_3dW(xlintR,DmatR,nelR,lint,2);
            else
            error('no tetrahedral tau yet')
            end

            der = 0;
            bf = 0;
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
            else
            % RFB
            end
            
            for ie = 1:lint

                %Integration point, weight, jacobian
                [Wgt,litr,lits] =  intpntt(ie,lint,0);
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
                [shlL,shldL,shls,be] = shlb([rL sL tL],nelL,nelL,0,0);
                QxyL = shgb(xlL,nelL,shldL,shls,nelL,0,0,be);

                xi = POU_Coord3(xint,yint,zint,xlR,1,nelR);
                rR = xi(1);
                sR = xi(2);
                tR = xi(3);

                % Evaluate  basis functions at integration points
                [shlR,shldR,shls,be] = shlb([rR sR tR],nelR,nelR,0,0);
                QxyR = shgb(xlR,nelR,shldR,shls,nelR,0,0,be);


                %Evaluate tangent and normal vectors
                t1 = xit*shld(:,1);%[xs(:,1); 0];
                [tm1, tu1] = VecNormalize(t1);
                t2 = xit*shld(:,2);%[xs(:,2); 0];
                [tm2, tu2] = VecNormalize(t2);
                t3 = VecCrossProd(t1,t2);
                [tm3, tu3] = VecNormalize(t3);
                % verify that unit vector points outward from L
                % elements are oriented so that nodes 1-4 are on bottom face of
                % the template element
                xcheck = xlL(:,1) + tu3';
                xi = POU_Coord3(xcheck(1),xcheck(2),xcheck(3),xlL,1,nelL);
                if ((nelL == 4 || nelL == 10) && xi(3) < 0) ...
                || ((nelL == 8 || nelL == 27) && xi(3) < -1)
                    nLx = tu3(1);
                    nLy = tu3(2);
                    nLz = tu3(3);
                else
                    nLx = -tu3(1);
                    nLy = -tu3(2);
                    nLz = -tu3(3);
                end
                nvect = [nLx 0 0 nLy 0   nLz 
                         0 nLy 0 nLx nLz 0 
                         0 0 nLz 0   nLy nLx];
                nvec = [nLx; nLy; nLz];

                c1 = Wgt*tm3;

                for i = 1:nelL
                    NmatL(1,(i-1)*ndf+1) = shlL(i);
                    NmatL(2,(i-1)*ndf+2) = shlL(i);
                    NmatL(3,(i-1)*ndf+3) = shlL(i);
                    BmatL(Bcol1,(i-1)*ndf+1) = QxyL(i,col1);
                    BmatL(Bcol2,(i-1)*ndf+2) = QxyL(i,col2);
                    BmatL(Bcol3,(i-1)*ndf+3) = QxyL(i,col3);
                end

                for i = 1:nelR
                    NmatR(1,(i-1)*ndf+1) = shlR(i);
                    NmatR(2,(i-1)*ndf+2) = shlR(i);
                    NmatR(3,(i-1)*ndf+3) = shlR(i);
                    BmatR(Bcol1,(i-1)*ndf+1) = QxyR(i,col1);
                    BmatR(Bcol2,(i-1)*ndf+2) = QxyR(i,col2);
                    BmatR(Bcol3,(i-1)*ndf+3) = QxyR(i,col3);
                end

                bnAdN1 = gamL*nvect*DmatL*BmatL;
                bnAdN2 = gamR*nvect*DmatR*BmatR;

                xint = xlL*shlL;

                tvtr = (bnAdN1*ulresL + bnAdN2*ulresR); %average stress
                jumpu = NmatR*ulresR - NmatL*ulresL;        %displacement jump

                % Not needed for linear problems
                ElemFL = ElemFL - c1*( - NmatL'*(tvtr + ep*jumpu) + bnAdN1'*jumpu);
                ElemFR = ElemFR - c1*( + NmatR'*(tvtr + ep*jumpu) + bnAdN2'*jumpu);

                ElemKLL = ElemKLL - c1*NmatL'*bnAdN1;
                ElemKLR = ElemKLR - c1*NmatL'*bnAdN2;
                ElemKRL = ElemKRL + c1*NmatR'*bnAdN1;
                ElemKRR = ElemKRR + c1*NmatR'*bnAdN2;

                ElemKLL = ElemKLL - c1*bnAdN1'*NmatL;
                ElemKLR = ElemKLR + c1*bnAdN1'*NmatR;
                ElemKRL = ElemKRL - c1*bnAdN2'*NmatL;
                ElemKRR = ElemKRR + c1*bnAdN2'*NmatR;

                ElemKLL = ElemKLL + c1*(NmatL'*ep*NmatL);
                ElemKLR = ElemKLR - c1*(NmatL'*ep*NmatR);
                ElemKRL = ElemKRL - c1*(NmatR'*ep*NmatL);
                ElemKRR = ElemKRR + c1*(NmatR'*ep*NmatR);

            end
        
        end
% ElemKLL
            ElemK = [ElemKLL ElemKLR; ElemKRL ElemKRR];
            ElemF = [ElemFL; ElemFR];
 %%
    case -1
        % Integrate over interface triangles, create surface traction to
        % subtract from the total surface traction, so that the remaining
        % traction in the global force vector only corresponds to the
        % exposed surfaces.
        ElemF = zeros(nst,1);
        
        lint = 7;13;

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
                [shlL,shldL,shls,be] = shlb([rL sL tL],nelL,nelL,0,0);
                shgL = shgb(xlL,nelL,shldL,shls,nelL,0,0,be);

                %Evaluate tangent and normal vectors
                t1 = xit*shld(:,1);%[xs(:,1); 0];
                [tm1, tu1] = VecNormalize(t1);
                t2 = xit*shld(:,2);%[xs(:,2); 0];
                [tm2, tu2] = VecNormalize(t2);
                t3 = VecCrossProd(t1,t2);
                [tm3, tu3] = VecNormalize(t3);
    %             t = [tu1' tu2' tu3'];
            
                if iprob == 5

                else
                    Traction = traction;
                end

                c1 = Wgt*tm3;

                for o=1:nelL

                    don = shlL(o);
                    F = don*Traction';

                    ElemF(ndf*(o-1)+1) = ElemF(ndf*(o-1)+1) + F(1)*c1;

                    ElemF(ndf*(o-1)+2) = ElemF(ndf*(o-1)+2) + F(2)*c1;

                    ElemF(ndf*(o-1)+3) = ElemF(ndf*(o-1)+3) + F(3)*c1;

                end %o

            end %ie
            
        end %intt
        ElemF;
%%
    case 11
        
        ElemE = zeros(numEn,1);

%%
    case 21 %interface stiffness
        
        ElemKRR = zeros(nstR,nstR);
        ElemKRL = zeros(nstR,nstL);
        ElemKLR = zeros(nstL,nstR);
        ElemKLL = zeros(nstL,nstL);
        t1 = zeros(3,1);
        t2 = t1;
        t3 = t1;
        
        bf1 = 0;
        bf2 = 0;
        
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
        
        ulresL = reshape(ulL,ndf*nelL,1);
        ulresR = reshape(ulR,ndf*nelR,1);
        
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
            
            % Perform integration of various matrices
    % For separate bubble types on T and B
            if nelL == 8 || nelL == 27
            [tauL,intbL,volL] = TauE1_3dW(xlintL,DmatL,nelL,lint,2);
            else
            error('no tetrahedral tau yet')
            end
            if nelR == 8 || nelR == 27
            [tauR,intbR,volR] = TauE1_3dW(xlintR,DmatR,nelR,lint,2);
            else
            error('no tetrahedral tau yet')
            end

            der = 0;
            bf = 0;
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
            else
            % RFB
            end
            
            for ie = 1:lint

                %Integration point, weight, jacobian
                [Wgt,litr,lits] =  intpntt(ie,lint,0);
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
                [shlL,shldL,shls,be] = shlb([rL sL tL],nelL,nelL,0,0);
                QxyL = shgb(xlL,nelL,shldL,shls,nelL,0,0,be);

                xi = POU_Coord3(xint,yint,zint,xlR,1,nelR);
                rR = xi(1);
                sR = xi(2);
                tR = xi(3);

                % Evaluate  basis functions at integration points
                [shlR,shldR,shls,be] = shlb([rR sR tR],nelR,nelR,0,0);
                QxyR = shgb(xlR,nelR,shldR,shls,nelR,0,0,be);


                %Evaluate tangent and normal vectors
                t1 = xit*shld(:,1);%[xs(:,1); 0];
                [tm1, tu1] = VecNormalize(t1);
                t2 = xit*shld(:,2);%[xs(:,2); 0];
                [tm2, tu2] = VecNormalize(t2);
                t3 = VecCrossProd(t1,t2);
                [tm3, tu3] = VecNormalize(t3);
                % verify that unit vector points outward from L
                % elements are oriented so that nodes 1-4 are on bottom face of
                % the template element
                xcheck = xlL(:,1) + tu3';
                xi = POU_Coord3(xcheck(1),xcheck(2),xcheck(3),xlL,1,nelL);
                if ((nelL == 4 || nelL == 10) && xi(3) < 0) ...
                || ((nelL == 8 || nelL == 27) && xi(3) < -1)
                    nLx = tu3(1);
                    nLy = tu3(2);
                    nLz = tu3(3);
                else
                    nLx = -tu3(1);
                    nLy = -tu3(2);
                    nLz = -tu3(3);
                end
                nvect = [nLx 0 0 nLy 0   nLz 
                         0 nLy 0 nLx nLz 0 
                         0 0 nLz 0   nLy nLx];
                nvec = [nLx; nLy; nLz];

                c1 = Wgt*tm3;

                for i = 1:nelL
                    NmatL(1,(i-1)*ndf+1) = shlL(i);
                    NmatL(2,(i-1)*ndf+2) = shlL(i);
                    NmatL(3,(i-1)*ndf+3) = shlL(i);
                    BmatL(Bcol1,(i-1)*ndf+1) = QxyL(i,col1);
                    BmatL(Bcol2,(i-1)*ndf+2) = QxyL(i,col2);
                    BmatL(Bcol3,(i-1)*ndf+3) = QxyL(i,col3);
                end

                for i = 1:nelR
                    NmatR(1,(i-1)*ndf+1) = shlR(i);
                    NmatR(2,(i-1)*ndf+2) = shlR(i);
                    NmatR(3,(i-1)*ndf+3) = shlR(i);
                    BmatR(Bcol1,(i-1)*ndf+1) = QxyR(i,col1);
                    BmatR(Bcol2,(i-1)*ndf+2) = QxyR(i,col2);
                    BmatR(Bcol3,(i-1)*ndf+3) = QxyR(i,col3);
                end

                bnAdN1 = gamL*nvect*DmatL*BmatL;
                bnAdN2 = gamR*nvect*DmatR*BmatR;

                ElemKLL = ElemKLL - c1*NmatL'*bnAdN1;
                ElemKLR = ElemKLR - c1*NmatL'*bnAdN2;
                ElemKRL = ElemKRL + c1*NmatR'*bnAdN1;
                ElemKRR = ElemKRR + c1*NmatR'*bnAdN2;

                ElemKLL = ElemKLL - c1*bnAdN1'*NmatL;
                ElemKLR = ElemKLR + c1*bnAdN1'*NmatR;
                ElemKRL = ElemKRL - c1*bnAdN2'*NmatL;
                ElemKRR = ElemKRR + c1*bnAdN2'*NmatR;

                ElemKLL = ElemKLL + c1*(NmatL'*ep*NmatL);
                ElemKLR = ElemKLR - c1*(NmatL'*ep*NmatR);
                ElemKRL = ElemKRL - c1*(NmatR'*ep*NmatL);
                ElemKRR = ElemKRR + c1*(NmatR'*ep*NmatR);

            end
        
        end
% ElemKLL
            ElemK = [ElemKLL ElemKLR; ElemKRL ElemKRR];

%%        
    case 25 %Stress Projection2

        
    case 60
        
        numhr = 3;
        ElemI = zeros(14,numhr);
        t1 = zeros(3,1);
        t2 = t1;
        t3 = t1;
        thick = 1;
        
        bf1 = 0;
        bf2 = 0;
        
        ElemEL = matepropL(1);
        ElemvL = matepropL(2);
        ElemER = matepropR(1);
        ElemvR = matepropR(2);
        lamdaR = ElemvR*ElemER/((1+ElemvR)*(1-2*ElemvR));
        muR = ElemER/(2*(1+ElemvR));
        lamdaL = ElemvL*ElemEL/((1+ElemvL)*(1-2*ElemvL));
        muL = ElemEL/(2*(1+ElemvL));
        if PSPS == 'n'
        DmatR = muR*diag([2 2 1]) + lamdaR*[1; 1; 0]*[1 1 0];
        DmatL = muL*diag([2 2 1]) + lamdaL*[1; 1; 0]*[1 1 0];
        else
        DmatL = ElemEL/(1-ElemvL^2)*[1      ElemvL  0
                                  ElemvL  1      0
                                  0      0      (1-ElemvL)/2];
        DmatR = ElemER/(1-ElemvR^2)*[1      ElemvR  0
                                  ElemvR  1      0
                                  0      0      (1-ElemvR)/2];
        end
        
        NmatL = zeros(2,nstL);
        BmatL = zeros(3,nstL);
        bnAdN1 = zeros(3,nstL);
        N1 = zeros(2,nstL);
        NmatR = zeros(2,nstR);
        BmatR = zeros(3,nstR);
        bnAdN2 = zeros(3,nstR);
        N2 = zeros(2,nstR);
        
        ulresL = reshape(ulL,ndf*nelL,1);
        ulresR = reshape(ulR,ndf*nelR,1);
        
        % Determin bounds of integration segment
        InterConn2D2 % InterConn2DT % 
        
%         h = 2/(hR + hL);

        % Set jacobian for integration space
        drdr = (eL2 - eL1)/drL;
        
        m = (eR2-eR1)/(eL1-eL2);
        
        % Perform integration of various matrices
        
%         [tauL,intbL] = TauE1_2dT(xlintL,DmatL,lintt6,nelL);
%         [tauR,intbR] = TauE1_2dT(xlintR,DmatR,lintt6,nelR);
% For separate bubble types on T and Q
        [tauL,intbL] = TauE1_2d(xlintL,DmatL,nelL,lintt6,lintq9);
        [tauR,intbR] = TauE1_2d(xlintR,DmatR,nelR,lintt6,lintq9);
%         % lumping
%         tauL(1,1) = sum(tauL(:,1));
%         tauL(2,2) = sum(tauL(:,2));
%         tauL(1,2) = 0;
%         tauL(2,1) = 0;
%         tauR(1,1) = sum(tauR(:,1));
%         tauR(2,2) = sum(tauR(:,2));
%         tauR(1,2) = 0;
%         tauR(2,1) = 0;
        tau = tauL + tauR;
        
        lint = 3;
        ib = 0;
        der = 0;
        bf = 0;
        ebL = 0;
        ebR = 0;
        intedge = 0;
        
        % THIS LOOP COULD BE EVALUATED IN CLOSED-FORM
        for ie = 1:lint

% % For triangle bubbles on T and Q
%             if nelL == 3 || nelL == 6
%                 [Wgt,r,s] = intpntt(ie,lint,1);
%                 rT = r;
%                 sT = 0;
%             elseif nelL == 4 || nelL == 9
%                 [Wgt,r,s] = intpntq(ie,lint,1);
%                 rT = (r + 1)/2;
%                 sT = 0;
%             end
%                     
%             ebeL = edgebubble(rT,sT);
%             ebeR = ebeL;
            
% For separate bubble types on T and Q
            if nelL == 3 || nelL == 6
                [Wgt,r,s] = intpntt(ie,lint,1);
                ebeL = edgebubble(r,s,nelL);
            elseif nelL == 4 || nelL == 9
                [Wgt,r,s] = intpntq(ie,lint,1);
                ebeL = edgebubbleQ(r,s,nelL);
            end
                    
            if nelL == 3 || nelL == 6
                [Wgt,rR,sR] = intpntt(ie,lint,1);
                ebeR = edgebubble(rR,sR,nelR);
            elseif nelL == 4 || nelL == 9
                [Wgt,rR,sR] = intpntq(ie,lint,1);
                ebeR = edgebubbleQ(rR,sR,nelR);
            end
            
            r = drdr*(r-roL)+eL1;
            
            if nelL == 3 || nelL == 6
                [shlL,shld,shls,be] = shlt(r,s,nelL,nelL,der,bf);
                [QxyL, shgs, Jdet, be, xsL] = shgt(xlL(:,1:nelL),nelL,shld,shls,nen,bf,der,be);
            elseif nelL == 4 || nelL == 9
                [shlL,shld,shls,be] = shlq(r,s,nelL,nelL,der,bf);
                [QxyL, shgs, Jdet, be, xsL] = shgq(xlL(:,1:nelL),nelL,shld,shls,nen,bf,der,be);
            end
                    
            %Evaluate tangent and normal vectors
            t1 = [xsL(:,1); 0];
            [tm1, tu1] = VecNormalize(t1);
            t2 = [0; 0; 1];
            tm2 = 1;
            tu2 = t2';
            t3 = VecCrossProd(t1,t2);
            [tm3, tu3] = VecNormalize(t3);
            nLx = tu3(1);
            nLy = tu3(2);
            nvect = [nLx nLy];
                    
            c1 = Wgt*tm3*drdr*thick;
            
            ebL = ebL + c1*ebeL;
            ebR = ebR + c1*ebeR;
            intedge = intedge + c1;
            
        end
        
        if nitvms == 1
        % VMS
%         volL = getvol(xlL,nelL);
%         volR = getvol(xlR,nelR);
%         volbL = getvol(xlintL,nelL);
%         volbR = getvol(xlintR,nelR);
%         tauL = tauL*(volL/volbL)^-1;
%         tauR = tauR*(volR/volbR)^-1;
        edgeK = tauL*ebL^2 + tauR*ebR^2;
        gamL = ebL^2*(edgeK\tauL);
        gamR = ebR^2*(edgeK\tauR);
        ep = pencoeff*intedge*inv(edgeK); %#ok<MINV>
        elseif nitvms == 2
        % Nitsche
        volL = getvol(xlL,nelL);
        volR = getvol(xlR,nelR);
        h = 2/(intedge/volR + intedge/volL);
        gamL = 1/2*eye(2);
        gamR = 1/2*eye(2);
        ep = pencoeff*eye(2)*max(muL,muR)/h;
        else
        % RFB
        tauL = RFBEfunc(xlintL,muL,lamdaL,8);
        tauR = RFBEfunc(xlintR,muR,lamdaR,8);
%         if inter == 65
%             tauL = 1.0e-008 *[
%    0.315216684586325   0.060542382391646
%    0.060542382391646   0.389745645413441];
%             tauR = 1.0e-009 *[
%    0.498227303287427   0.043693780768305
%    0.043693780768297   0.972566080948220];
%         end
%         tauL = RFBtaus(2*D/m1,1/((n1/m1)/(L/D)));
%         tauR = RFBtaus(2*D/m2,1/((n2/m2)/(L/D)));
        edgeK = tauL + tauR;
        gamL = (edgeK\tauL);
        gamR = (edgeK\tauR);
        ep = pencoeff*intedge*inv(edgeK); %#ok<MINV>
        end
        
        for ie = 1:lint
            
            if nelL == 3 || nelL == 6
                [Wgt,r,s] = intpntt(ie,lint,1);
            elseif nelL == 4 || nelL == 9
                [Wgt,r,s] = intpntq(ie,lint,1);
            end
            
            r = drdr*(r-roL)+eL1;
            
            if nelL == 3 || nelL == 6
                [shlL,shld,shls,be] = shlt(r,s,nelL,nelL,der,bf);
                [QxyL, shgs, Jdet, be, xsL] = shgt(xlL(:,1:nelL),nelL,shld,shls,nen,bf,der,be);
            elseif nelL == 4 || nelL == 9
                [shlL,shld,shls,be] = shlq(r,s,nelL,nelL,der,bf);
                [QxyL, shgs, Jdet, be, xsL] = shgq(xlL(:,1:nelL),nelL,shld,shls,nen,bf,der,be);
            end
            
            rR = m*(r-eL2) + eR1;
            xint = xlL(1,1:nelL)*shlL;
            yint = xlL(2,1:nelL)*shlL;
            
            if nelR == 3 || nelR == 6
                s = 0;
            else %if nelR == 4
                s = -1;
            end
            
            if nelR == 3 || nelR == 6
                [shlR,shld,shls,be] = shlt(rR,s,nelR,nelR,der,bf);
                [QxyR, shgs, Jdet, be, xsR] = shgt(xlR(:,1:nelR),nelR,shld,shls,nen,bf,der,be);
            elseif nelR == 4 || nelR == 9
                [shlR,shld,shls,be] = shlq(rR,s,nelR,nelR,der,bf);
                [QxyR, shgs, Jdet, be, xsR] = shgq(xlR(:,1:nelR),nelR,shld,shls,nen,bf,der,be);
            end
            
            %Evaluate tangent and normal vectors
            t1 = [xsL(:,1); 0];
            [tm1, tu1] = VecNormalize(t1);
            t2 = [0; 0; 1];
            tm2 = 1;
            tu2 = t2';
            t3 = VecCrossProd(t1,t2);
            [tm3, tu3] = VecNormalize(t3);
            nLx = tu3(1);
            nLy = tu3(2);
            nvect = [nLx 0 nLy
                     0 nLy nLx];
            nvec = [nLx; nLy];
                    
            c1 = Wgt*tm3*drdr*thick;
            
            for i = 1:nelL
                NmatL(1,(i-1)*2+1) = shlL(i);
                NmatL(2,(i-1)*2+2) = shlL(i);
                BmatL(Bcol1,(i-1)*2+1) = QxyL(i,col1);
                BmatL(Bcol2,(i-1)*2+2) = QxyL(i,col2);
            end
            
            for i = 1:nelR
                NmatR(1,(i-1)*2+1) = shlR(i);
                NmatR(2,(i-1)*2+2) = shlR(i);
                BmatR(Bcol1,(i-1)*2+1) = QxyR(i,col1);
                BmatR(Bcol2,(i-1)*2+2) = QxyR(i,col2);
            end
            
            bnAdN1 = gamL*nvect*DmatL*BmatL;
            bnAdN2 = gamR*nvect*DmatR*BmatR;
        
            xint = xlL*shlL;
            
            if iprob == 3
                T1 = h*q/12/kT*(6*xint(2)/h - 1) + 1/(2*h*rhoT*cT)*qtau;
                T2 = q/(4*kT*h)*xint(2)^2;
                T = T1 + T2;
                cvec = alph*T*[1; 1; 0];
            elseif iprob == 4
                cvec = alph*1*[1; 1; 0];
            end
            
            tvtr = (bnAdN1*ulresL + bnAdN2*ulresR) - (gamL*nvect*DmatL + gamR*nvect*DmatR)*cvec; %average stress
            jumpu = NmatR*ulresR - NmatL*ulresL;        %displacement jump

        %Compute value of exact fields at int. point
        if iprob == 1
        elseif iprob == 2
            [ue,duex,duey] = uexactbb(xint(1),xint(2),ElemEL,ElemvL,D);
            strvec = [duex(1); duey(2); duex(2)+duey(1)];
        elseif iprob == 3
            [ue,duex,duey] = uexacttb(xint(1),xint(2),ElemEL,ElemvL,alph,kT,rhoT,cT,q,qtau,h);
            T1 = h*q/12/kT*(6*xint(2)/h - 1) + 1/(2*h*rhoT*cT)*qtau;
            T2 = q/(4*kT*h)*xint(2)^2;
            T = T1 + T2;
            cvec = alph*T*[1; 1; 0];
            strvec = [duex(1); duey(2); duex(2)+duey(1)] - cvec;
            sxx2 = alph*ElemEL/(1-ElemvL)*(q/4/h/kT)*(h^2/3-xint(2)^2);
        end
        t_exact = nvect*DmatL*strvec;

        ElemI(:,ie) = [xint
                       jumpu
                       tvtr
                       (tvtr + ep*jumpu)
                       t_exact
                       (tvtr + ep*jumpu)-t_exact
                       (tvtr)-t_exact];
        end %lint
        
end %Task Switch
