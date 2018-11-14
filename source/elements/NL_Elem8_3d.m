% 04/30/2013
%
% -Example 3D linear elastic solid element
% -Demonstrates the behavior of shifting dof ordering in FEAP
% -In isw=1 case, the element turns off all dofs that exceed ndm=3 because
%  they are not used for pure-displacement.
% -The same solution will be obtained, regardless of DG elements in the
%  mesh, if the ordering of dofs is changed in MatTypeTable (see e.g.
%  CohShearAxialTestM3.m):
%      MatTypeTable = [1; 8; 0]; 
%      MatTypeTable = [1; 8; 0; 3; 1; 2; 0]; 
%      MatTypeTable = [1; 8; 0; 3; 2; 1; 0]; 
%      MatTypeTable = [1; 8; 0; 1; 3; 2; 0];
% -This is because the rows in ul and in the stiffness matrix are shuffled,
%  and this is recorded in iedof when the material is first loaded.
% -While there is no point to do this here, it shows how Matlab and FEAP
%  would act for complicated elements, e.g. thermomechanical.
% -Only isw=3 and -1 have been debugged as of 04/30

PatchE = mateprop(1);
Patchv = mateprop(2);

Bcol1 = [1; 4; 6];
Bcol2 = [2; 4; 5];
Bcol3 = [3; 5; 6];
col1 = [1; 2; 3];
col2 = [2; 1; 3];
col3 = [3; 2; 1];
iemat = ieFEAP(1:ndf,ma);

switch isw %Task Switch
    
    case 1
        
        if ndf > 3
            
            for i = 4:ndf
                lie(i,1) = 0;
            end
            
        end
        
        istv = 11;
        iste = 11;

    case 3
         
        ElemK = zeros(ndf*nel);
        ElemF = zeros(ndf*nel,1);

        %Set integration number
        if nel ~= 6
%         lint = IntPoint3(nel);
            if nel == 4
                lint = 1;4;11;5;16;
            elseif nel == 8
                lint = 8;1000; %1000 for body force problem 
            elseif nel == 10
                lint = 5;11;14;
    %             lint = 27;
            elseif nel == 18
                linttw = 13;
                lintlw = 3;
                lint = linttw*lintlw;
            else
                lint = 27;
            end
        else
            linttw = 7;
            lintlw = 2;
            lint = linttw*lintlw;
        end
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
        Nmat = zeros(3,ndf*nel);
        Bmat = zeros(6,ndf*nel);
        rhspul = reshape(ul,ndf*nen,1);

        for l = 1:lint                    

            % Evaluate 1-D basis functions at integration points
            if nel == 4 || nel == 10
              [w,ss] =  int3d_t(l,lint,ib);
              [shl,shld,shls,be] = shltt(ss,nel,nel,der,bf);
              [shg,shgs,Jdet,be] = shgtt(xl,nel,shld,shls,nel,bf,der,be);
            elseif nel == 6 || nel == 18 % wedge
              [w,ss] =  intpntw(l,linttw,lintlw,ib);
              [shl,shld,shls,be] = shlw(ss,nel,nel,der,bf);
              [shg,shgs,Jdet,be,xs] = shgw(xl,nel,shld,shls,nel,bf,der,be);
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

            Fvec = [fbx; fby; fbz];
            
            % Form B matrix
            for i = 1:nel
                
%             Nmat(:,(ie-1)*3+1:3*ie) = [shl(ie) 0       0      
%                                        0       shl(ie) 0      
%                                        0       0       shl(ie)];
              Nmat(1,(i-1)*ndf+iemat(1)) = shl(i);
              Nmat(2,(i-1)*ndf+iemat(2)) = shl(i);
              Nmat(3,(i-1)*ndf+iemat(3)) = shl(i);
                
%             Bmat(:,(ie-1)*3+1:3*ie) = [shg(ie,1) 0         0         
%                                        0         shg(ie,2) 0         
%                                        0         0         shg(ie,3) 
%                                        shg(ie,2) shg(ie,1) 0         
%                                        0         shg(ie,3) shg(ie,2) 
%                                        shg(ie,3) 0         shg(ie,1) ];
%             Bmat(1,(ie-1)*3+1) = shg(ie,1);
%             Bmat(4,(ie-1)*3+1) = shg(ie,2);
%             Bmat(6,(ie-1)*3+1) = shg(ie,3);
%             Bmat(2,(ie-1)*3+2) = shg(ie,2);
%             Bmat(4,(ie-1)*3+2) = shg(ie,1);
%             Bmat(5,(ie-1)*3+2) = shg(ie,3);
%             Bmat(3,3*ie      ) = shg(ie,3);
%             Bmat(5,3*ie      ) = shg(ie,2);
%             Bmat(6,3*ie      ) = shg(ie,1);
            Bmat(Bcol1,(i-1)*ndf+iemat(1)) = shg(i,col1);
            Bmat(Bcol2,(i-1)*ndf+iemat(2)) = shg(i,col2);
            Bmat(Bcol3,(i-1)*ndf+iemat(3)) = shg(i,col3);
                 
            end
            
%             ElemF = ElemF + c1*(Nmat'*Fve4 + BBmat'*Tmat*Fvec);
            ElemF = ElemF + c1*(Nmat'*Fvec - Bmat'*Dmat*(Bmat*rhspul(1:ndf*nel)));
            
% if unifelem == 0
    
%             ElemK = ElemK + c1*(Bmat'*Dmat*Bmat - BBmat'*Tmat*BBmat);
            ElemK = ElemK + c1*(Bmat'*Dmat*Bmat);

% end % unifelem

        end %je
ElemK;
% ElemF
% if unifelem == 1
%     
%     if elem <= 16
%         ElemK = ElemKB;
%     else
%         ElemK = ElemKT;
%     end
%     
% end

    case 11
        
        ElemE = zeros(numEn,1);
        %Set integration number
        lint = IntPoint3(nel);
        ib = 0;
        bf = 0;
        der = 0;

        lam = Patchv*PatchE/((1+Patchv)*(1-2*Patchv));
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

            c1 = Jdet*w*thick;

            xint = xl(1,1:nel)*shl;
            yint = xl(2,1:nel)*shl;
            zint = xl(3,1:nel)*shl;
            dux = ul(:,1:nel)*shg(:,1);
            duy = ul(:,1:nel)*shg(:,2);
            duz = ul(:,1:nel)*shg(:,3);
            u = ul(:,1:nel)*shl;

            %Compute value of exact fields at int. point
%             if iprob == 1
%                 [ue,duex,duey,duez] = uexact_bar3(xint,yint,zint,PatchE,Patchv,rho);
%             elseif iprob == 2
%                 [ue,duex,duey,duez] = uexact_fluid(xint,yint,zint);
%             else
                ue = zeros(4,1);duex=ue;duey=ue;duez=ue;
%             end

            %Add standard int. point error to element standard error
            for in = 1:3
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

    case 21
 
        ElemK = zeros(nst);

        %Set integration number
%         lint = IntPoint3(nel);
        if nel ~= 6
%         lint = IntPoint3(nel);
            if nel == 4
                lint = 1;4;11;5;16;
            elseif nel == 8
                lint = 8;1000; %1000 for body force problem 
            elseif nel == 10
                lint = 5;11;14;
    %             lint = 27;
            elseif nel == 18
                linttw = 13;
                lintlw = 3;
                lint = linttw*lintlw;
            else
                lint = 27;
            end
        else
            linttw = 7;
            lintlw = 2;
            lint = linttw*lintlw;
        end
        ib = 0;
        bf = 0;
        der = 0;

        lam = Patchv*PatchE/((1+Patchv)*(1-2*Patchv));
        mu = PatchE/(2*(1+Patchv));
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
            elseif nel == 6 || nel == 18 % wedge
              [w,ss] =  intpntw(l,linttw,lintlw,ib);
              [shl,shld,shls,be] = shlw(ss,nel,nel,der,bf);
              [shg,shgs,Jdet,be,xs] = shgw(xl,nel,shld,shls,nel,bf,der,be);
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
        
        % Reorder nodes for corresponding face of integration
        SurfOrientEtoS
        
        t1 = zeros(3,1);
        t2 = t1;
        t3 = t1;
        
        if nel == 4 || nel == 10
            lint = 13;
        else
            lint = 16;
        end
        der = 0;
        bf = 0;
        ib = 0;
        
        % Integration Loop
        for je = 1:lint

            if nel == 4 || nel == 10
              [Wgt,ss] =  int3d_t(je,lint,edge);
              [shl,shld,shls,be] = shltt(ss,nel,nel,der,bf);
              [QXY,shgs,Jdet,be,sx] = shgtt(xl,nel,shld,shls,nel,bf,der,be);
            else
              [Wgt,ss] =  intpntb(je,lint,5);
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

            %Force components are positive in positive coord. direction
            c1 = Wgt*tm3;
            for o=1:nel
                don = shl(o);
                F = don*traction';

                ElemF(ndf*(o-1)+iemat(1)) = ElemF(ndf*(o-1)+iemat(1)) + F(1)*c1;

                ElemF(ndf*(o-1)+iemat(2)) = ElemF(ndf*(o-1)+iemat(2)) + F(2)*c1;

                ElemF(ndf*(o-1)+iemat(3)) = ElemF(ndf*(o-1)+iemat(3)) + F(3)*c1;

            end %o
            
        end %je
        ElemF;
        
        % Reorder nodes back to the orientation of the element
        SurfOrientStoE
        ElemF(1:ndf*nel) = ElemF(ilist2);
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
    case 8 %interface stiffness
        
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
        eN = 30*muR/((xlR(:,1)-xlR(:,7))'*(xlR(:,1)-xlR(:,7)));
        
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
            
            tT_tr = tT_n + eT*(gt - gtn_1);      %Elastic trial force
            normtT_tr = sqrt(tT_tr'*tT_tr);
            etatr = tT_tr/normtT_tr;
            Psi_tr = normtT_tr + fmu*tN_k;          %Slip Condition
            
            if (dbond == 0 && (tN_k <= 0)) || (dbond == 1 && gn <= 0) || (iter > itchatn && chattern == 1) % compresion only
            
                dbond = 0;
                bonddbond2n = 0;
                
                if Psi_tr > dtol
                dele = max(0,Psi_tr/eT);                 %Plastic multiplier
                tT_k = tT_tr - eT*dele*etatr;    %Slip correction
                fdbond = 1;
                else
                dele = 0;
                tT_k = tT_tr;
                fdbond = 0;
                end

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

                if Psi_tr > dtol

                    dddt = eye(3) - nmat;
                    dddt2 = dddt - etatr*etatr';
                    dddt3 = eT*(-(fmu*tN_k/normtT_tr))*dddt2;

                    ElemKLL = ElemKLL + c1*NmatL'*dddt3*NmatL;
                    ElemKLR = ElemKLR - c1*NmatL'*dddt3*NmatR;
                    ElemKRL = ElemKRL - c1*NmatR'*dddt3*NmatL;
                    ElemKRR = ElemKRR + c1*NmatR'*dddt3*NmatR;

                    ElemKLL = ElemKLL + c1*fmu*NmatL'*etatr*nvec'*(1/2*bnAdN1-eN*NmatL);
                    ElemKLR = ElemKLR + c1*fmu*NmatL'*etatr*nvec'*(1/2*bnAdN2+eN*NmatR);
                    ElemKRL = ElemKRL - c1*fmu*NmatR'*etatr*nvec'*(1/2*bnAdN1-eN*NmatL);
                    ElemKRR = ElemKRR - c1*fmu*NmatR'*etatr*nvec'*(1/2*bnAdN2+eN*NmatR);
            numdbond = numdbond + 1;
                    
                else

                    dddt = eye(3) - nmat;

                    ElemKLL = ElemKLL + c1*eT*NmatL'*dddt*NmatL;
                    ElemKLR = ElemKLR - c1*eT*NmatL'*dddt*NmatR;
                    ElemKRL = ElemKRL - c1*eT*NmatR'*dddt*NmatL;
                    ElemKRR = ElemKRR + c1*eT*NmatR'*dddt*NmatR;
            
                end
            
            else
                dbond = 1;
                fdbond = -1;
                bonddbond2n = -1;
    %             bonddbond2f = bonddbond1f;
                tT_k = [0; 0; 0];
            end
            
            if iter <= itchatn
              if sign(bonddbond1n) ~= sign(bonddbond2n)
                chattern = 1;
              else
                chattern = 0;
              end
            end
                    
            % Store history
            gshr = nh2-1+(ll-1)*18;
            tthr = nh2-1+(ll-1)*18+3;
            gthr = nh2-1+(ll-1)*18+6;
            tshr = nh2-1+(ll-1)*18+9;
            chat1nhr = nh2-1+(ll-1)*18+13;
            chat2nhr = nh2-1+(ll-1)*18+14;
            fricthf = nh2-1+(ll-1)*18+15;
            hr(gthr+1) = gt(1);
            hr(gthr+2) = gt(2);
            hr(gthr+3) = gt(3);
            hr(tthr+1) = tT_k(1);
            hr(tthr+2) = tT_k(2);
            hr(tthr+3) = tT_k(3);
            hr(nh1-1+(ll-1)*18+12) = dbond;
            hr(nh2-1+(ll-1)*18+12) = dbond;
            hr(chat1nhr) = bonddbond2n;
            hr(chat2nhr) = chattern;
            hr(fricthf) = fdbond;
%             hr(chat1fhr) = bonddbond2f;
%             hr(chat2fhr) = chatterf;
%             numdbond = numdbond + dbond;
            
            end %lint
        
%         end %intt
ElemFL;
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
            
            epsil = Bmat*reshape(ul(1:ndf,1:nel),ndf*nel,1);
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
            
            epsil = Bmat*reshape(ul(:,1:nel),ndf*nel,1);
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
        
    case 26
        
        ElemS = zeros(nel,nestr);

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

            nint = 1;
        
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
            
            epsil = Bmat*reshape(ul(:,1:nel),ndf*nel,1);
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
            
            ElemS(stres) = sigmas;
            
            end

        end %je
        
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
    case 26 % Element Stress

        ElemS = zeros(npstr,1);

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
        ll = 1;
        
        der = 0;
        bf = 0;
        ib = 0;

        %Stress Loop
        for ll = 1:lint

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
            
            epsil = Bmat*reshape(ul(:,1:nel),ndf*nel,1);
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
            
            ElemS(stres) = sigmas;
            
            end

        end %je
        
end %Task Switch
