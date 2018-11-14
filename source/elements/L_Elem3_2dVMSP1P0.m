
% 07/22-23/2013
% P1P0 mixed elasticity, uses pressure jumps for stability, penalty is
% computed using IVMS
% NOTE: hard-coded for T3 only, for which nel comes in as 4 and you have to
% convert it to 3 in the way things are handled; I got around that in the
% traction load integration by adding 'badnum' to NL_FEA_Program

% Set Material Properties

PatchE = mateprop(1);
Patchv = mateprop(2);
thick = mateprop(3);
w = 6;

Bcol1 = [1; 3];
Bcol2 = [2; 3];
col1 = [1; 2];
col2 = [2; 1];

nitvms = 2;
if nitvms == 1 %VMS
pencoeff = 1;
elseif nitvms == 2 %Nitsche
pencoeff = 10;1;2;
else %RFB
pencoeff = 1;4;2;
end

switch isw %Task Switch
    
    case 1
        
            for j = 1:3
            for i = 3:ndf
                lie(i,1+j) = 0;
            end
            end
            for i = 1:2
                lie(i,5) = 0;
            end
            lie(:,6:end) = 0;
            
%%
    case 3 %Interior of PatchEnts: stiffness and force
        
        % for some reason, this block of code did not work for higher order
        % elements; this needs to be checked later.
%         ElemK = zeros(nst);
%         ElemF = zeros(nst,1);
%         Nmat = zeros(3,nst);
%         Bmat = zeros(4,nst);
% %         BBmat = zeros(2,nst);
%         BBmat = zeros(8,nst);
%         
%         % Load Gauss Points for quadrature
%         if enrich == 1
%             if nel == 3
%                 lint = lintt3;%13;
%             elseif nel == 4
%                 lint = lintq4;
%             elseif nel == 6
%                 lint = lintt6;%13;
%             elseif nel == 9
%                 lint = lintq9;
%             end
% %             [rlist, rWgts, rnum] = GaussPoints(pr+2);
% %             [slist, sWgts, snum] = GaussPoints(pr+1);
%         else
%             if nel == 3
%                 lint = lintt3;%13;
%             elseif nel == 4
%                 lint = lintq4;
%             elseif nel == 6
%                 lint = lintt6;%13;
%             elseif nel == 9
%                 lint = lintq9;
%             end
% %             [rlist, rWgts, rnum] = GaussPoints(pr+1);
% %             [slist, sWgts, snum] = GaussPoints(pr+1);
%         end
% 
%         lam = Patchv*PatchE/((1+Patchv)*(1-2*Patchv));
%         mu = PatchE/(2*(1+Patchv));
%         Dmat = mu*diag([2 2 1]);
%         cmat = Dmat;
%         Dmat = [Dmat  [1; 1; 0]
%                 [1 1 0] -1/lam];
%         I4 = [1 0 0 0 0 1
%               0 0 1 0 1 0];
%         P1 = [1 0 0 0 0 0
%               0 0 0 0 0 1
%               0 0 1 1 0 0
%               0 0 1 0 0 0
%               0 0 0 0 1 0
%               0 1 0 0 0 1];
%         Z2 = zeros(2);
%         
% %         thick = 1;
%         ib = 0;
%         der = 1;
%         bf = 1;
%         
%         fbx = 0;
%         if iprob == 5
%             grav = 9.81;
%             fby = -rho*grav;
%         else
%             fby = 0;
%         end
% 
%         if regelem == 1
%             [t11,t12,t21,t22] = Tausym(xl,mu,nel);
%         else
%             [t11,t12,t21,t22] = Tau(xl,mu,nel,nen,lint);
%             if nel == 3
%                 tcon = tcont3;%13;
%             elseif nel == 4
%                 tcon = tconq4;
%             elseif nel == 6
%                 tcon = tcont6;%13;
%             elseif nel == 9
%                 tcon = tconq9;
%             end
%             t11 = tcon*t11;
%             t12 = tcon*t12;
%             t21 = tcon*t21;
%             t22 = tcon*t22;
%         end
%         
%         for je = 1:lint
% 
%             if nel == 3 || nel == 6
%                 [Wgt,litr,lits] = intpntt(je,lint,ib);
%                 [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
%                 [Qxy, shgs, Jdet, be, xs] = shgt(xl,nel,shld,shls,nen,bf,der,be);
%             elseif nel == 4 || nel == 9
%                 [Wgt,litr,lits] = intpntq(je,lint,ib);
%                 [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
%                 [Qxy, shgs, Jdet, be, xs] = shgq(xl,nel,shld,shls,nen,bf,der,be);
%             end
%             c1 = Wgt*Jdet*thick;
% 
%             fb = [fbx; fby; 0];
%             
%             Tmat = be(3)*[t11 t12; t21 t22];
%             
%             % Form B matrix
%             for ie = 1:nel
%                 
%             Nmat(:,(ie-1)*3+1:3*ie) = shl(ie)*eye(3);
%                 
%             Bmat(:,(ie-1)*3+1:3*ie) = [Qxy(ie,1) 0         0
%                                        0         Qxy(ie,2) 0
%                                        Qxy(ie,2) Qxy(ie,1) 0
%                                        0         0         shl(ie)];
%                                    
% %             BBmat(:,(ie-1)*3+1:3*ie) = [mu*(2*cartd2(ie,1)+cartd2(ie,2)) mu*cartd2(ie,3) Qxy(ie,1)
% %                                         mu*cartd2(ie,3) mu*(cartd2(ie,1)+2*cartd2(ie,2)) Qxy(ie,2)];
%             BBmat(:,(ie-1)*3+1:3*ie) = [shgs(ie,1) 0         0       
%                      shgs(ie,2) 0         0       
%                      shgs(ie,3) 0         0       
%                      0         shgs(ie,1) 0       
%                      0         shgs(ie,2) 0       
%                      0         shgs(ie,3) 0      
%                      0         0         Qxy(ie,1)
%                      0         0         Qxy(ie,2)];
%                                     
%             end
%             term22 = eye(2);
%             diagc = [cmat(1:3,1:3) zeros(3); zeros(3) cmat(1:3,1:3)];
%             term8b = I4*diagc*P1; %page 8b
%             BBT = [term8b Z2; zeros(2,6) term22];
%             Tmat1 = [Tmat Z2; Z2 Tmat];
%             Tmat2 = [Z2 Tmat; Tmat Z2];
%             D22 = - BBT'*Tmat1*BBT - BBT'*Tmat2*BBT;
%             ElemTR = Tmat*(fb(1:2));
%             Tvec = -[ElemTR; ElemTR];
%             
% %             ElemF = ElemF + c1*(Nmat'*fb + BBmat'*Tmat*fb(1:2));
%             ElemF = ElemF + c1*(Nmat'*fb - BBmat'*BBT'*Tvec);
% %             ElemK = ElemK + c1*(Bmat'*Dmat*Bmat - BBmat'*Tmat*BBmat);
%             ElemK = ElemK + c1*(Bmat'*Dmat*Bmat + BBmat'*D22*BBmat);
%                 
%         end %je
%         
%         ElemK = 1/2*(ElemK + ElemK');
        
        ElemK = zeros(nst);
        ElemF = zeros(nst,1);
        Nmat = zeros(3,nst);
        Bmat = zeros(4,nst);
        BBmat = zeros(2,nst);
        
        % Load Gauss Points for quadrature
        lint = 1;3;lintt3;

        lam = Patchv*PatchE/((1+Patchv)*(1-2*Patchv));
        mu = PatchE/(2*(1+Patchv));
        Dmat = mu*diag([2 2 1]);
        Dmat = [Dmat  [1; 1; 0]
                [1 1 0] -1/lam];
        
%         thick = 1;
        ideriv = 1;
        
        fbx = 0;
        if iprob == 5
            grav = 9.81;
            fby = -rho*grav;
        else
            fby = 0;
        end

        der = 0;
        bf = 0;
        
        
        for je = 1:lint

%             if nel == 3 || nel == 6
                [Wgt,r,s] = intpntt(je,lint,0);
%             else %if nel == 4 || nel == 9
%                 [Wgt,r,s] = intpntq(je,lint,0);
%             end
            
%             if nel == 3 || nel == 6
                [shl,shld,shls,be] = shlt(r,s,3,3,der,bf);
                [shg, shgs, Jdet, be, xs] = shgt(xl,3,shld,shls,nen,bf,der,be);
%             elseif nel == 4 || nel == 9
%                 [shl,shld,shls,be] = shlq(r,s,nel,nel,der,bf);
%                 [shg, shgs, Jdet, be, xs] = shgq(xl,nel,shld,shls,nen,bf,der,be);
%             end
%             % Evaluate 1-D basis functions at integration points
%             [shp,shp2,be] = shpl_2d(r,s,nel,ideriv,ep,enrich);
%             %Evaluate first derivatives of basis functions at int. point
%             [Qxy, cartd2, Jdet] = shpg_2d(shp,shp2,xl,nel2,ideriv,be);
            c1 = Wgt*Jdet*thick;

%             if iprob == 4
%                 alp = 0.5;
%                 xint = xl(1,1:nel)*shl;
%                 yint = xl(2,1:nel)*shl;
%                 fbx = - 2*mu*(2*xint - 4) - (w*yint*cos(w*xint));
%                 fby = 2*mu*(w^2*sin(w*xint) - 2*xint + 4) - sin(w*xint) - (((5*((mu)^(1/2)/alp - 1))/(8*((2*(mu)^(1/2))/alp + 2)) - 17/48)*(2*yint + 1))/(2);
%             end
            fb = [fbx; fby; 0];
            
            
            % Form B matrix
            for i = 1:3
                
              Bmat(Bcol1,(i-1)*ndf+1) = shg(i,col1);
              Bmat(Bcol2,(i-1)*ndf+2) = shg(i,col2);
                                    
            end
              Bmat(4    ,(4-1)*ndf+3) = 1;
            
%             ElemF = ElemF + c1*(Nmat'*fb + BBmat'*Tmat*fb(1:2));
            ElemK = ElemK + c1*(Bmat'*Dmat*Bmat);
                
        end %je
        
        ElemK = 1/2*(ElemK + ElemK');
        
%%
    case 15 % body force for linear element
        
%         ElemF = zeros(nst,1); In FormFE
        
        if iprob == 5 || iprob == 4
            
        Nmat = zeros(3,nst);
        BBmat = zeros(2,nst);
        
        % Load Gauss Points for quadrature
%         if nel == 3
            lint = lintt3;%13;
%         elseif nel == 4
%             lint = lintq4;
%         elseif nel == 6
%             lint = lintt6;%13;
%         elseif nel == 9
%             lint = lintq9;
%         end

        lam = Patchv*PatchE/((1+Patchv)*(1-2*Patchv));
        mu = PatchE/(2*(1+Patchv));
        Dmat = mu*diag([2 2 1]);
        Dmat = [Dmat  [1; 1; 0]
                [1 1 0] -1/lam];
        
%         thick = 1;

        der = 0;
        bf = 0;
        
        if iprob == 5
            fbx = 0;
            grav = 9.81;
            fby = -rho*grav;
            fb = [fbx; fby; 0];
        else
            fb = bodyf;
        end
        
        for je = 1:lint

%             if nel == 3 || nel == 6
                [Wgt,r,s] = intpntt(je,lint,0);
%             else %if nel == 4 || nel == 9
%                 [Wgt,r,s] = intpntq(je,lint,0);
%             end
            
%             if nel == 3 || nel == 6
                [shl,shld,shls,be] = shlt(r,s,nel,nel,der,bf);
                [shg, shgs, Jdet, be, xs] = shgt(xl,nel,shld,shls,nen,bf,der,be);
%             elseif nel == 4 || nel == 9
%                 [shl,shld,shls,be] = shlq(r,s,nel,nel,der,bf);
%                 [shg, shgs, Jdet, be, xs] = shgq(xl,nel,shld,shls,nen,bf,der,be);
%             end
%             % Evaluate 1-D basis functions at integration points
%             [shp,shp2,be] = shpl_2d(r,s,nel,ideriv,ep,enrich);
%             %Evaluate first derivatives of basis functions at int. point
%             [Qxy, cartd2, Jdet] = shpg_2d(shp,shp2,xl,nel2,ideriv,be);
            c1 = Wgt*Jdet*thick;

%             if iprob == 4
%                 alp = 0.5;
%                 xint = xl(1,1:nel)*shl;
%                 yint = xl(2,1:nel)*shl;
%                 fbx = - 2*mu*(2*xint - 4) - (w*yint*cos(w*xint));
%                 fby = 2*mu*(w^2*sin(w*xint) - 2*xint + 4) - sin(w*xint) - (((5*((mu)^(1/2)/alp - 1))/(8*((2*(mu)^(1/2))/alp + 2)) - 17/48)*(2*yint + 1))/(2);
%             end
            
%             Tmat = be(3)*[t11 t12; t21 t22];
            
            % Form B matrix
            for i = 1:nel
                
              Nmat(1,(i-1)*ndf+1) = shl(i);
              Nmat(2,(i-1)*ndf+2) = shl(i);
              Nmat(3,(i-1)*ndf+3) = shl(i);
                                    
            end
            
            ElemF = ElemF + c1*(Nmat'*fb );
                
        end %je
        ElemF;
        end
        
%%        
    case -1 % boundary tractions
        
        ElemF = zeros(nst,1);

        % Reorder nodes for corresponding face of integration
        if edge == 1
            ilist = [1 2 3 4];
        elseif edge == 2
            ilist = [2 3 1 4];
        else %if edge == 3
            ilist = [3 1 2 4];
        end
        ElemFlag(1:nel) = ElemFlag(ilist);
        xl(1:ndm,1:nel) = xl(1:ndm,ilist);
        % Determine bounds of integration
        
%         if nel == 4 || nel == 9
%             
%         dr = 2;
%         ro = -1;
%         
%         % Upper Limit
%         if nodeA == ElemFlag(2)
%             eR2 = 1;
%         else %
%             eR2 = 0;
%         end
%         % Lower Limit
%         if nodeB == ElemFlag(1)
%             eR1 = -1;
%         else %nodeA == ElemFlag(5)
%             eR1 = 0;
%         end
        
%         elseif nel == 3 || nel == 6
            
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
        
%         end
        
        % Set jacobian for integration space
        drdr = (eR2 - eR1)/dr;
        
        lint = 4;
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

%         lamda = Patchv*PatchE/((1+Patchv)*(1-2*Patchv));
%         mu = PatchE/(2*(1+Patchv));
%         thick = 1;
        ib  = 1;
        der = 0;
        bf = 0;
        
        for je = 1:lint

%             if nel == 3 || nel == 6
                [Wgt,r,s] = intpntt(je,lint,1);
%             else %if nel == 4
%                 [Wgt,r,s] = intpntq(je,lint,1);
%             end
            r = drdr*(r-ro)+eR1;

%             if nel == 3 || nel == 6
                [shl,shld,shls,be] = shlt(r,s,nel,nel,der,bf);
                [Qxy, cartd2, Jdet, be, xs] = shgt(xl,nel,shld,shls,nen,bf,der,be);
%             elseif nel == 4 || nel == 9
%                 [shl,shld,shls,be] = shlq(r,s,nel,nel,der,bf);
%                 [Qxy, cartd2, Jdet, be, xs] = shgq(xl,nel,shld,shls,nen,bf,der,be);
%             end
                    
            %Evaluate tangent and normal vectors
            t1 = [xs(:,1); 0];
            [tm1, tu1] = VecNormalize(t1);
            t2 = [0; 0; 1];
            tm2 = 1;
            tu2 = t2';
            t3 = VecCrossProd(t1,t2);
            [tm3, tu3] = VecNormalize(t3);
            
            if iprob == 5
                x = 0;
                y = 0;

                for j = 1:nel
                   x = x + xl(1,j)*shl(j);
                   y = y + xl(2,j)*shl(j);
                end
                grav = 9.81d0;
                Cwid = D;
                Len = L/2;
                x2 = x-Len;
                s_xx = rho*grav*y - 1.d0/10.d0*rho*grav*y*(15.d0*Len^2-15.d0*x2^2 +4.d0*Cwid^2 +10.d0*y^2)/Cwid^2;
                s_xy = 3.d0/2.d0*rho*grav*x2*(Cwid^2-y^2)/Cwid^2;
                s_yy = rho*grav*y - 1.d0/2.d0*rho*grav*y*(3.d0*Cwid^2-y^2)/Cwid^2;
                Traction(1) = s_xx*tu3(1) + s_xy*tu3(2);
                Traction(2) = s_xy*tu3(1) + s_yy*tu3(2);
            elseif iprob == 4
                x = 0;
                y = 0;

                for j = 1:nel
                   x = x + xl(1,j)*shl(j);
                   y = y + xl(2,j)*shl(j);
                end
                Cwid = D;
                Len = L;
%                 Load = 2560;
                s_xx = 3.d0*(x-Len)*y*Load/(2*Cwid^3);
                s_xy = 3.d0/4.d0*(Cwid^2-y^2)*Load/Cwid^3;
                s_yy = 0;
                Traction(1) = s_xx*tu3(1) + s_xy*tu3(2);
                Traction(2) = s_xy*tu3(1) + s_yy*tu3(2);
            elseif iprob == 2
                x = 0;
                y = 0;

                for j = 1:nel
                   x = x + xl(1,j)*shl(j);
                   y = y + xl(2,j)*shl(j);
                end

                s_xx = 72*y;
                s_xy = 0;
                s_yy = 0;
                Traction(1) = s_xx*tu3(1) + s_xy*tu3(2);
                Traction(2) = s_xy*tu3(1) + s_yy*tu3(2);
            elseif iprob == 7
                x = 0;
                y = 0;

                for j = 1:nel
                   x = x + xl(1,j)*shl(j);
                   y = y + xl(2,j)*shl(j);
                end
                lam2 = 0.54448373678246398;
                Q = 0.54307557883673652;
                r = sqrt(x^2+y^2);
                if(y>=0.d0)
                  theta = acos(x/r);
                else
                  theta = -acos(x/r);
                end
                r = lam2*r^(lam2-1.d0);
                s_xx = r*((2.d0 - Q*(lam2+1.d0))*cos((lam2-1.d0)*theta) ...
                     - (lam2-1.d0)*cos((lam2-3.d0)*theta));
                s_yy = r*((2.d0 + Q*(lam2+1.d0))*cos((lam2-1.d0)*theta) ...
                     + (lam2-1.d0)*cos((lam2-3.d0)*theta));
                s_xy = r*((lam2-1.d0)*sin((lam2-3.d0)*theta) ...
                     + Q*(lam2+1.d0)*sin((lam2-1.d0)*theta));
                Traction(1) = s_xx*tu3(1) + s_xy*tu3(2);
                Traction(2) = s_xy*tu3(1) + s_yy*tu3(2);
            elseif iprob == 8
                x = 0;
                y = 0;

                for j = 1:nel
                   x = x + xl(1,j)*shl(j);
                   y = y + xl(2,j)*shl(j);
                end

                s_xx = q0*y^2;
                s_xy = 0;
                s_yy = -q0*x^2;
                Traction(1) = s_xx*tu3(1) + s_xy*tu3(2);
                Traction(2) = s_xy*tu3(1) + s_yy*tu3(2);
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

                ElemF(ndf*(o-1)+1) = ElemF(ndf*(o-1)+1) + F(1)*c1;

                ElemF(ndf*(o-1)+2) = ElemF(ndf*(o-1)+2) + F(2)*c1;

            end %o

        end %ie
        ElemF;
        
        % Reorder nodes back to the orientation of the element
        if edge == 1
            ilist = [1 2 3 4];
        elseif edge == 2
            ilist = [3 1 2 4];
        else %if edge == 3
            ilist = [2 3 1 4];
        end
        ilist2 = zeros(1,ndf*nel);
        for i = 1:ndf
            ilist2(i:ndf:ndf*nel-ndf+i) = ndf*ilist - ndf + i;
        end
        ElemF(1:ndf*nel) = ElemF(ilist2);
%%       
    case 11
        
        ElemE = zeros(numEn,1);
        
        % Load Gauss Points for quadrature
%         if nel == 3
            lint = lintt3;%13;
%         elseif nel == 4
%             lint = lintq4;
%         elseif nel == 6
%             lint = lintt6;%13;
%         elseif nel == 9
%             lint = lintq9;
%         end
        pfact = 1;

%         if PSPS == 's'
%             PatchE = PatchE*(1 + 2*Patchv)/(1 + Patchv)^2;
%             Patchv = Patchv/(1 + Patchv);
%         end
        lam = Patchv*PatchE/((1+Patchv)*(1-2*Patchv));
        mu = PatchE/(2*(1+Patchv));
        vis = mu;
%         thick = 1;
        ideriv = 1;
        fbx = 0;
        fby = 0;
        hx = 0;
        hy = 0;
        ib = 0;
        bf = 0;
        der = 0;

        el2el = zeros(3,1);
        eprixel = zeros(3,1);
        epriyel = zeros(3,1);
        el2fine = zeros(7,1);
        ue = zeros(3,1);
        duex = ue;
        duey = ue;

        if iprob == 3

            fbx = 2;

        elseif iprob == 5

            grav = 9.81;
            fby = -rho*grav;

        end

%         [t11,t12,t21,t22] = Tau3_2d(xl,mu,nel,nen,lint);
%         if nel == 3
%             tcon = tcont3;%13;
%         elseif nel == 4
%             tcon = tconq4;
%         elseif nel == 6
%             tcon = tcont6;%13;
%         elseif nel == 9
%             tcon = tconq9;
%         end
%         t11 = tcon*t11;
%         t12 = tcon*t12;
%         t21 = tcon*t21;
%         t22 = tcon*t22;

            for je = 1:lint

%                 if nel == 3 || nel == 6
                    [Wgt,r,s] = intpntt(je,lint,0);
%                 elseif nel == 4 || nel == 9
%                     [Wgt,r,s] = intpntq(je,lint,0);
%                 end

%                 if nel == 3 || nel == 6
                    [shl,shld,shls,be] = shlt(r,s,nel,nel,der,bf);
                    [shg, shgs, Jdet, be, xs] = shgt(xl,nel,shld,shls,nen,bf,der,be);
%                 elseif nel == 4 || nel == 9
%                     [shl,shld,shls,be] = shlq(r,s,nel,nel,der,bf);
%                     [shg, shgs, Jdet, be, xs] = shgq(xl,nel,shld,shls,nen,bf,der,be);
%                 end
                c1 = Wgt*Jdet*thick;
                b = be(3);

%                 if nelP == 3 || nelP == 6
%                     [shlp,shld,shls,bub] = shlt(r,s,nelP,nel,0,0);
%                     [shp, shp2, Jdet, bub, xs] = shgt(xl,nelP,shld,shls,nen,0,0,bub);
%                 elseif nelP == 4 || nelP == 9
%                     [shlp,shld,shls,bub] = shlq(r,s,nelP,nel,0,0);
%                     [shp, shp2, Jdet, bub, xs] = shgq(xl,nelP,shld,shls,nen,0,0,bub);
%                 end

                % Displacement Terms of Stiffness Matrix
%                 for jn = 1:nel2

                xint = xl(1,1:nel)*shl;
                yint = xl(2,1:nel)*shl;
                dux = ul(1:2,1:nel)*shg(:,1);
                duy = ul(1:2,1:nel)*shg(:,2);
                u = ul(1:2,1:nel)*shl;
                dux(3) = 0;
                duy(3) = 0;
                u(3) = ul(3,4);

                %Compute value of exact fields at int. point
                if iprob == 1
%                     [ue,duex,duey] = uexact_bar(xint,yint,PatchE,Patchv);
                elseif iprob == 2
                    [ue,duex,duey] = uexact_bend(xint,yint,PatchE,Patchv);
%                 elseif iprob == 3
%                     [ue,duex,duey] = uexact_faxe(xint,yint,PatchE,Patchv);
                elseif iprob == 4
%                     [ue,duex,duey] = uexact_beam(xint,yint,PatchE,Patchv,Load);
                    [ue,duex,duey] = uexact_vasil(xint,yint,mu,1,0.5,w);
                elseif iprob == 5
                    [ue,duex,duey] = uexact_selfw(xint,yint,PatchE,Patchv,rho);
                elseif iprob == 7
                    [ue,duex,duey] = uexact_lshp(xint,yint,PatchE,Patchv);
                elseif iprob == 8
                    [ue,duex,duey] = uexact_gh(xint,yint,PatchE,Patchv,q0);
                end

                %Add standard int. point error to element standard error
                for in = 1:3
                    un   = c1 * ( (u(in)-ue(in))^2 );
                    upnx   = c1 * ( (dux(in)-duex(in))^2 );
                    upny   = c1 * ( (duy(in)-duey(in))^2 );
                    el2el(in)   = el2el(in)   + un;
                    eprixel(in) = eprixel(in) + upnx;
                    epriyel(in) = epriyel(in) + upny;
                end

        end %je

        for in= 1:3
            ElemE(in) = el2el(in);
            ElemE(in+3) = eprixel(in);
            ElemE(in+6) = epriyel(in);
        end
%         ElemE(16) = ElemE(16) + el2fine(7);
        
        H1u = eprixel(1)+eprixel(2)+epriyel(1)+epriyel(2);
        Ieffvals(elem,1:3) = [sqrt(0/H1u) 0 H1u];
    
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
        DmatR = muR*diag([2 2 1]);
        DmatL = muL*diag([2 2 1]);
        
        ulresL = reshape(ulL,ndf*nelL,1);
        ulresR = reshape(ulR,ndf*nelR,1);
        
%         NmatL = zeros(2,nstL);
        NmatL = zeros(3,nstL);
        BmatL = zeros(4,nstL);
        bnAdN1 = zeros(4,nstL);
%         N1 = zeros(2,nstL);
        N1 = zeros(3,nstL);
%         NmatR = zeros(2,nstR);
        NmatR = zeros(3,nstR);
        BmatR = zeros(4,nstR);
        bnAdN2 = zeros(4,nstR);
%         N2 = zeros(2,nstR);
        N2 = zeros(3,nstR);
        
        % Determin bounds of integration segment
        InterConn2D2 % InterConn2DT % 
        
        h = 2/(hR + hL);
%         h = 0.125;
        
        eN = dNfac*efac/h;1e9;8;
        eT = dTfac*efac/h;
        ePP = pfac/efac*h;
        
        % Set jacobian for integration space
        drdr = (eL2 - eL1)/drL;
        
        m = (eR2-eR1)/(eL1-eL2);

        s = -1;
        
%         etauL = TauEE2d(xlintL,DmatL,lintt6);
%         etauR = TauEE2d(xlintR,DmatR,lintt6);
% For separate bubble types on T and Q
        [tauL,intbL] = TauE1_2d(xlintL,DmatL,nelL,lintt6,lintq9);
        [tauR,intbR] = TauE1_2d(xlintR,DmatR,nelR,lintt6,lintq9);
%         tau = tauL;
%         etau = etauL + etauR;
%         ep = tauR/tauL;
%         ep = 1;
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
        Kinv = [ep zeros(2,1); zeros(1,2) 0];
        elseif nitvms == 2
        % Nitsche
        volL = getvol(xlL,nelL);
        volR = getvol(xlR,nelR);
%         h = 2/(intedge/volR + intedge/volL);
        gamL = 1/2*eye(2);
        gamR = 1/2*eye(2);
        ep = pencoeff*eye(2)*max(muL,muR)/h;
        ePP = -pencoeff*2/max(muL,muR)*h;
        Kinv = [ep zeros(2,1); zeros(1,2) ePP];
        end

        for ie = 1:lint

            if nelL == 3 || nelL == 6
                [Wgt,r,s] = intpntt(ie,lint,1);
            elseif nelL == 4 || nelL == 9
                [Wgt,r,s] = intpntq(ie,lint,1);
            end
            
            r = drdr*(r-roL)+eL1;

            if nelLP == 3 || nelLP == 6
                [shpL,shld,shls,be] = shlt(r,s,nelLP,nelLP,0,0);
                [PxyL, shgs, Jdet] = shgt(xlL(:,1:nelLP),nelLP,shld,shls,nen,0,0,be);
            elseif nelLP == 4 || nelLP == 9
                [shpL,shld,shls,be] = shlq(r,s,nelLP,nelLP,0,0);
                [PxyL, shgs, Jdet] = shgq(xlL(:,1:nelLP),nelLP,shld,shls,nen,0,0,be);
            end
            
            if nelL == 3 || nelL == 6
                [shlL,shld,shls,be] = shlt(r,s,nelL,nelL,der,bf);
                [QxyL, shgs, Jdet, be, xsL] = shgt(xlL(:,1:nelL),nelL,shld,shls,nen,bf,der,be);
            elseif nelL == 4 || nelL == 9
                [shlL,shld,shls,be] = shlq(r,s,nelL,nelL,der,bf);
                [QxyL, shgs, Jdet, be, xsL] = shgq(xlL(:,1:nelL),nelL,shld,shls,nen,bf,der,be);
            end
            
            rR = m*(r-eL2) + eR1;
            
            if nelR == 3 || nelR == 6
                s = 0;
            else %if nelR == 4
                s = -1;
            end
            
            if nelRP == 3 || nelRP == 6
                [shpR,shld,shls,be] = shlt(rR,s,nelRP,nelRP,0,0);
                [PxyR, shgs, Jdet, be, xsR] = shgt(xlR(:,1:nelRP),nelRP,shld,shls,nen,0,0,be);
            elseif nelRP == 4 || nelRP == 9
                [shpR,shld,shls,be] = shlq(rR,s,nelRP,nelRP,0,0);
                [PxyR, shgs, Jdet, be, xsR] = shgq(xlR(:,1:nelRP),nelRP,shld,shls,nen,0,0,be);
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
                      0 nLy nLx]; %- ?
            
            c1 = Wgt*tm3*drdr*thick;
                
            for i = 1:nelL
%                 NmatL(:,3*i-2:3*i) = [shlL(i)*eye(2) zeros(2,1)];
                NmatL(:,3*i-2:3*i) = shlL(i)*eye(3);
                BmatL(:,3*i-2:3*i) = [QxyL(i,1) 0 0
                                      0 QxyL(i,2) 0
                                      QxyL(i,2) QxyL(i,1) 0
                                      0 0 shlL(i)];
            end
            
            for i = 1:nelR
%                 NmatR(:,3*i-2:3*i) = [shlR(i)*eye(2) zeros(2,1)];
                NmatR(:,3*i-2:3*i) = shlR(i)*eye(3);
                BmatR(:,3*i-2:3*i) = [QxyR(i,1) 0 0
                                      0 QxyR(i,2) 0
                                      QxyR(i,2) QxyR(i,1) 0
                                      0 0 shpR(i)];
            end
               
%%       
            if numSI > 0
            
            N1 = NmatL;
            N2 = NmatR;
            
            bnAdN1 = gamL*nvect*[DmatL [1; 1; 0]]*BmatL;
            bnAdN2 = gamR*nvect*[DmatR [1; 1; 0]]*BmatR;
        
            % Signs set to agree with L_Elem3_2d, original implementation
            ElemKLL = ElemKLL - c1*NmatL(1:2,:)'*bnAdN1 - c1*bnAdN1'*NmatL(1:2,:);
            ElemKLR = ElemKLR - c1*NmatL(1:2,:)'*bnAdN2 + c1*bnAdN1'*NmatR(1:2,:);
            ElemKRL = ElemKRL + c1*NmatR(1:2,:)'*bnAdN1 - c1*bnAdN2'*NmatL(1:2,:);
            ElemKRR = ElemKRR + c1*NmatR(1:2,:)'*bnAdN2 + c1*bnAdN2'*NmatR(1:2,:);
        
            ElemKLL = ElemKLL + c1*N1'*Kinv*N1 ;
            ElemKLR = ElemKLR - c1*N1'*Kinv*N2 ;
            ElemKRL = ElemKRL - c1*N2'*Kinv*N1 ;
            ElemKRR = ElemKRR + c1*N2'*Kinv*N2 ;
        
            end %numSI

        end %ie
        
        ElemKLL;

%%
    case 22 %Stress Projection
        
        ElemM = zeros(nst);
        ElemF = zeros(nst,1);
        
        % Load Gauss Points for quadrature
        if enrich == 1
            if nel == 3
                lint = lintt3;%13;
            elseif nel == 4
                lint = lintq4;
            elseif nel == 6
                lint = lintt6;%13;
            elseif nel == 9
                lint = lintq9;
            end
%             [rlist, rWgts, rnum] = GaussPoints(pr+2);
%             [slist, sWgts, snum] = GaussPoints(pr+1);
        else
            if nel == 3
                lint = lintt3;%13;
            elseif nel == 4
                lint = lintq4;
            elseif nel == 6
                lint = lintt6;%13;
            elseif nel == 9
                lint = lintq9;
            end
%             [rlist, rWgts, rnum] = GaussPoints(pr+1);
%             [slist, sWgts, snum] = GaussPoints(pr+1);
        end

%         if PSPS == 's'
%             PatchE = PatchE*(1 + 2*Patchv)/(1 + Patchv)^2;
%             Patchv = Patchv/(1 + Patchv);
%         end
        lam = Patchv*PatchE/((1+Patchv)*(1-2*Patchv));
        mu = PatchE/(2*(1+Patchv));
        vis = mu;
%         thick = 1;
        ib = 0;
        der = 0;
        bf = 0;
        fbx = 0;
        fby = 0;
        hx = 0;
        hy = 0;
        II = eye(2,2);
        inds = [1 2 1
                1 2 2];

        for je = 1:lint

            if nel == 3 || nel == 6
                [Wgt,litr,lits] = intpntt(je,lint,ib);
                [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
                [shg, cartd2, Jdet, be, xs] = shgt(xl,nel,shld,shls,nen,bf,der,be);
            elseif nel == 4 || nel == 9
                [Wgt,litr,lits] = intpntq(je,lint,ib);
                [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
                [shg, cartd2, Jdet, be, xs] = shgq(xl,nel,shld,shls,nen,bf,der,be);
            end
            c1 = Wgt*Jdet*thick;
            
            xint = xl(1,:)*shl;
            yint = xl(2,:)*shl;
            dux = ul(1:2,:)*shg(:,1);
            duy = ul(1:2,:)*shg(:,2);
            u = ul(1:2,:)*shl;
            strains = [2*dux(1) 2*duy(2) dux(2)+duy(1)];
            p = ul(3,:)*shl;
            
            sigmas = p*II(inds(1,stres),inds(2,stres)) + mu*strains(stres);

            for jn = 1:nel

                djn=shl(jn);

                ElemF(jn) = ElemF(jn) + c1*djn*sigmas;

                for in=1:nel2

                    din=shl(in);

                    ElemM(in,jn)=ElemM(in,jn) + c1*din*djn;

                end %in

            end %jn

        end %je
%%
    case 23
        
        ElemKRR = zeros(nstR,nstR);
        ElemKRL = zeros(nstR,nstL);
        ElemKLR = zeros(nstL,nstR);
        ElemKLL = zeros(nstL,nstL);
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
        
        if neqc > 0 %Lagrange Multipliers
            ElemKcL = zeros(nstC,nstL);
            ElemKcR = zeros(nstC,nstR);
            ElemKLc = zeros(nstL,nstC);
            ElemKRc = zeros(nstR,nstC);
        end

        % Determine bounds of integration, right
        
        if nelR == 4 || nelR == 9
            
            drR = 2;
            roR = -1;

            % Upper Limit
            if nodeAR == ElemFlagR(2)
                eR2 = 1;
            elseif nodeAR == ElemFlagR(nel2R)
                eR2 = epR;
            elseif nodeAR == -1 %no enrichment but DG instead
                xy = xlL(:,1);
                POUxi = POU_Coord(xy(1),xy(2),xlR,1,nelR);
                eR2 = POUxi(1);
            elseif nelR == 9 && nodeAR == ElemFlagR(5)
                eR2 = 0;
            end
            % Lower Limit
            if nodeBR == ElemFlagR(1)
                eR1 = -1;
            elseif nodeBR == ElemFlagR(nel2R)
                eR1 = epR;
            elseif nodeBR == -1 %no enrichment but DG instead
                xy = xlL(:,2);
                POUxi = POU_Coord(xy(1),xy(2),xlR,1,nelR);
                eR1 = POUxi(1);
            elseif nelR == 9 && nodeBR == ElemFlagR(5)
                eR1 = 0;
            end
        
        elseif nelR == 3 || nelR == 6
            
            drR = 1;
            roR = 0;

            % Upper Limit
            if nodeAR == ElemFlagR(2)
                eR2 = 1;
            elseif nodeAR == ElemFlagR(nel2R)
                eR2 = epR;
            elseif nodeAR == -1 %no enrichment but DG instead
                xy = xlL(:,1);
                POUxi = POU_Coord(xy(1),xy(2),xlR,1,nelR);
                eR2 = POUxi(1);
            elseif nelR == 6 && nodeAR == ElemFlagR(4)
                eR2 = 1/2;
            end
            % Lower Limit
            if nodeBR == ElemFlagR(1)
                eR1 = 0;
            elseif nodeBR == ElemFlagR(nel2R)
                eR1 = epR;
            elseif nodeBR == -1 %no enrichment but DG instead
                xy = xlL(:,2);
                POUxi = POU_Coord(xy(1),xy(2),xlR,1,nelR);
                eR1 = POUxi(1);
            elseif nelR == 6 && nodeBR == ElemFlagR(4)
                eR1 = 1/2;
            end
        
        end
        
        % Determine bounds of integration, left
        
        if nelL == 4 || nelL == 9
            
            drL = 2;
            roL = -1;

            % Upper Limit
            if nodeAL == ElemFlagL(1)
                eL1 = -1;
            elseif nodeAL == ElemFlagL(nel2L)
                eL1 = epL;
            elseif nodeAL == -1 %no enrichment but DG instead
                xy = xlR(:,2);
                POUxi = POU_Coord(xy(1),xy(2),xlL,1,nelL);
                eL1 = POUxi(1);
            elseif nelL == 9 && nodeAL == ElemFlagL(5)
                eL1 = 0;
            end
            % Lower Limit
            if nodeBL == ElemFlagL(2)
                eL2 = 1;
            elseif nodeBL == ElemFlagL(nel2L)
                eL2 = epL;
            elseif nodeBL == -1 %no enrichment but DG instead
                xy = xlR(:,1);
                POUxi = POU_Coord(xy(1),xy(2),xlL,1,nelL);
                eL2 = POUxi(1);
            elseif nelL == 9 && nodeB == ElemFlagL(5)
                eL2 = 0;
            end
        
        elseif nelL == 3 || nelL == 6
            
            drL = 1;
            roL = 0;

            % Upper Limit
            if nodeAL == ElemFlagL(1)
                eL1 = 0;
            elseif nodeAL == ElemFlagL(nel2L)
                eL1 = epL;
            elseif nodeAL == -1 %no enrichment but DG instead
                xy = xlR(:,2);
                POUxi = POU_Coord(xy(1),xy(2),xlL,1,nelL);
                eL1 = POUxi(1);
            elseif nelL == 6 && nodeAL == ElemFlagL(4)
                eL1 = 1/2;
            end
            % Lower Limit
            if nodeBL == ElemFlagL(2)
                eL2 = 1;
            elseif nodeBL == ElemFlagL(nel2L)
                eL2 = epL;
            elseif nodeBL == -1 %no enrichment but DG instead
                xy = xlR(:,1);
                POUxi = POU_Coord(xy(1),xy(2),xlL,1,nelL);
                eL2 = POUxi(1);
            elseif nelL == 6 && nodeBL == ElemFlagL(4)
                eL2 = 1/2;
            end
        
        end
        
        % Set jacobian for integration space
        drdr = (eL2 - eL1)/drL;
        
        m = (eR2-eR1)/(eL1-eL2);
        
        lint = 4;
        ib = 1;
        der = 0;
        bf = 0;

        s = 0;
        
        for ie = 1:lint

            if nelL == 3 || nelL == 6
                [Wgt,r,s] = intpntt(ie,lint,1);
            elseif nelL == 4 || nelL == 9
                [Wgt,r,s] = intpntq(ie,lint,1);
            end
            
            r = drdr*(r-roL)+eL1;
            
            if nelL == 3 || nelL == 6
                [shpL,shld,shls,be] = shlt(r,s,nelL,nelL,der,bf);
                [QxyL, shgs, Jdet, be, xsL] = shgt(xlL(:,1:nelL),nelL,shld,shls,nen,bf,der,be);
            elseif nelL == 4 || nelL == 9
                [shpL,shld,shls,be] = shlq(r,s,nelL,nelL,der,bf);
                [QxyL, shgs, Jdet, be, xsL] = shgq(xlL(:,1:nelL),nelL,shld,shls,nen,bf,der,be);
            end
            
            rR = m*(r-eL2) + eR1;
            
            if nelR == 3 || nelR == 6
                s = 0;
            else %if nelR == 4
                s = -1;
            end
            
            if nelR == 3 || nelR == 6
                [shpR,shld,shls,be] = shlt(rR,s,nelR,nelR,der,bf);
                [QxyR, shgs, Jdet, be, xsR] = shgt(xlR(:,1:nelR),nelR,shld,shls,nen,bf,der,be);
            elseif nelR == 4 || nelR == 9
                [shpR,shld,shls,be] = shlq(rR,s,nelR,nelR,der,bf);
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
            nRx = -tu3(1);
            nRy = -tu3(2);
            tLx = tu1(1);
            tLy = tu1(2);
            
            dux = ulL(1:2,:)*QxyL(:,1);
            duy = ulL(1:2,:)*QxyL(:,2);
            strainsL = [2*dux(1) 2*duy(2) dux(2)+duy(1)];
            pL = ulL(3,:)*shpL;
            sigmasL = pL*II(inds(1,stres),inds(2,stres)) + muL*strainsL(stres);
            dux = ulR(1:2,:)*QxyR(:,1);
            duy = ulR(1:2,:)*QxyR(:,2);
            strainsR = [2*dux(1) 2*duy(2) dux(2)+duy(1)];
            pR = ulR(3,:)*shpR;
            sigmasR = pR*II(inds(1,stres),inds(2,stres)) + muR*strainsR(stres);
            
            sigmas = sigmasR-sigmasL;
                    
            c1 = Wgt*tm3*drdr*thick;
                
            for a = 1:nel2L
            
                NLa = shpL(a)*c1;

                for c = 1:nel2L

                    NLc = shpL(c);

                    ElemKLL(c,a) = ElemKLL(c,a) - ...
                                             10*(NLc*NLa);
                                         
                end
                                         
                for d = 1:nel2R

                    NRd = shpR(d);

                    ElemKRL(d,a) = ElemKRL(d,a) + ...
                                             10*(NRd*NLa);
                                         
                end

            end
            
            for b = 1:nel2R
            
                NRb = shpR(b)*c1;

                for c = 1:nel2L

                    NLc = shpL(c);

                    ElemKLR(c,b) = ElemKLR(c,b) + ...
                                             10*(NLc*NRb);
                                         
                end
                                         
                for d = 1:nel2R

                    NRd = shpR(d);

                    ElemKRR(d,b) = ElemKRR(d,b) - ...
                                             10*(NRd*NRb);
                                         
                end

            end
            
        end %lint
        
        ElemKLL = 1/10*ElemKLL;
        ElemKRL = 1/10*ElemKRL;
        ElemKLR = 1/10*ElemKLR;
        ElemKRR = 1/10*ElemKRR;
%%        

    case 60
        
        numhr = 3;
        ElemI = zeros(14,numhr);
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
        DmatR = muR*diag([2 2 1]);
        DmatL = muL*diag([2 2 1]);
        
        ulresL = reshape(ulL,ndf*nelL,1);
        ulresR = reshape(ulR,ndf*nelR,1);
        
%         NmatL = zeros(2,nstL);
        NmatL = zeros(3,nstL);
        BmatL = zeros(4,nstL);
        bnAdN1 = zeros(4,nstL);
%         N1 = zeros(2,nstL);
        N1 = zeros(3,nstL);
%         NmatR = zeros(2,nstR);
        NmatR = zeros(3,nstR);
        BmatR = zeros(4,nstR);
        bnAdN2 = zeros(4,nstR);
%         N2 = zeros(2,nstR);
        N2 = zeros(3,nstR);
        
        % Determin bounds of integration segment
        InterConn2D2 % InterConn2DT % 
        
        h = 2/(hR + hL);
%         h = 0.125;
        
        eN = dNfac*efac/h;1e9;8;
        eT = dTfac*efac/h;
        ePP = pfac/efac*h;
        
        % Set jacobian for integration space
        drdr = (eL2 - eL1)/drL;
        
        m = (eR2-eR1)/(eL1-eL2);

        s = -1;
        
%         etauL = TauEE2d(xlintL,DmatL,lintt6);
%         etauR = TauEE2d(xlintR,DmatR,lintt6);
% For separate bubble types on T and Q
        [tauL,intbL] = TauE1_2d(xlintL,DmatL,nelL,lintt6,lintq9);
        [tauR,intbR] = TauE1_2d(xlintR,DmatR,nelR,lintt6,lintq9);
%         tau = tauL;
%         etau = etauL + etauR;
%         ep = tauR/tauL;
%         ep = 1;
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
        Kinv = [ep zeros(2,1); zeros(1,2) 0];
        elseif nitvms == 2
        % Nitsche
        volL = getvol(xlL,nelL);
        volR = getvol(xlR,nelR);
%         h = 2/(intedge/volR + intedge/volL);
        gamL = 1/2*eye(2);
        gamR = 1/2*eye(2);
        ep = pencoeff*eye(2)*max(muL,muR)/h;
        ePP = -pencoeff*2/max(muL,muR)*h;
        Kinv = [ep zeros(2,1); zeros(1,2) ePP];
        end

        for ie = 1:lint

            if nelL == 3 || nelL == 6
                [Wgt,r,s] = intpntt(ie,lint,1);
            elseif nelL == 4 || nelL == 9
                [Wgt,r,s] = intpntq(ie,lint,1);
            end
            
            r = drdr*(r-roL)+eL1;
            xint = xlL(1,1:nelL)*shlL;
            yint = xlL(2,1:nelL)*shlL;

            if nelLP == 3 || nelLP == 6
                [shpL,shld,shls,be] = shlt(r,s,nelLP,nelLP,0,0);
                [PxyL, shgs, Jdet] = shgt(xlL(:,1:nelLP),nelLP,shld,shls,nen,0,0,be);
            elseif nelLP == 4 || nelLP == 9
                [shpL,shld,shls,be] = shlq(r,s,nelLP,nelLP,0,0);
                [PxyL, shgs, Jdet] = shgq(xlL(:,1:nelLP),nelLP,shld,shls,nen,0,0,be);
            end
            
            if nelL == 3 || nelL == 6
                [shlL,shld,shls,be] = shlt(r,s,nelL,nelL,der,bf);
                [QxyL, shgs, Jdet, be, xsL] = shgt(xlL(:,1:nelL),nelL,shld,shls,nen,bf,der,be);
            elseif nelL == 4 || nelL == 9
                [shlL,shld,shls,be] = shlq(r,s,nelL,nelL,der,bf);
                [QxyL, shgs, Jdet, be, xsL] = shgq(xlL(:,1:nelL),nelL,shld,shls,nen,bf,der,be);
            end
            
            rR = m*(r-eL2) + eR1;
            
            if nelR == 3 || nelR == 6
                s = 0;
            else %if nelR == 4
                s = -1;
            end
            
            if nelRP == 3 || nelRP == 6
                [shpR,shld,shls,be] = shlt(rR,s,nelRP,nelRP,0,0);
                [PxyR, shgs, Jdet, be, xsR] = shgt(xlR(:,1:nelRP),nelRP,shld,shls,nen,0,0,be);
            elseif nelRP == 4 || nelRP == 9
                [shpR,shld,shls,be] = shlq(rR,s,nelRP,nelRP,0,0);
                [PxyR, shgs, Jdet, be, xsR] = shgq(xlR(:,1:nelRP),nelRP,shld,shls,nen,0,0,be);
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
                      0 nLy nLx]; %- ?
            
            c1 = Wgt*tm3*drdr*thick;
                
            for i = 1:nelL
%                 NmatL(:,3*i-2:3*i) = [shlL(i)*eye(2) zeros(2,1)];
                NmatL(:,3*i-2:3*i) = shlL(i)*eye(3);
                BmatL(:,3*i-2:3*i) = [QxyL(i,1) 0 0
                                      0 QxyL(i,2) 0
                                      QxyL(i,2) QxyL(i,1) 0
                                      0 0 shlL(i)];
            end
            
            for i = 1:nelR
%                 NmatR(:,3*i-2:3*i) = [shlR(i)*eye(2) zeros(2,1)];
                NmatR(:,3*i-2:3*i) = shlR(i)*eye(3);
                BmatR(:,3*i-2:3*i) = [QxyR(i,1) 0 0
                                      0 QxyR(i,2) 0
                                      QxyR(i,2) QxyR(i,1) 0
                                      0 0 shpR(i)];
            end
               
%             xint = xlL*shlL;
            
            N1 = NmatL;
            N2 = NmatR;
            
            bnAdN1 = gamL*nvect*[DmatL [1; 1; 0]]*BmatL;
            bnAdN2 = gamR*nvect*[DmatR [1; 1; 0]]*BmatR;
            
            tvtr = (bnAdN1*ulresL + bnAdN2*ulresR); %average stress
            jumpu = NmatR*ulresR - NmatL*ulresL;        %displacement jump
            jumpu = jumpu(1:2);

        %Compute value of exact fields at int. point
        if iprob == 1
        elseif iprob == 5
            [ue,duex,duey] = uexact_selfw(xint,yint,PatchE,Patchv,rho);
            strvec = [duex(1); duey(2); duex(2)+duey(1)];
        end
        t_exact = nvect*(DmatL*strvec + [1; 1; 0]*ue(3));

        ElemI(:,ie) = [xint
                       yint
                       jumpu
                       tvtr
                       (tvtr + ep*jumpu)
                       t_exact
                       (tvtr + ep*jumpu)-t_exact
                       (tvtr)-t_exact];

        end %ie
        
end %Task Switch
