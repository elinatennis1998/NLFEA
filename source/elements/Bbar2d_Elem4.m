% Tim Truster
% CEE 577
% 04/21/2017
% Bbar for triangle patch

% Set Material Properties

ElemYM = mateprop(2);
Elemv = mateprop(3);
thick = mateprop(1);
Khard = mateprop(4);
Hhard = 0;%mateprop(5); % no kinematic hardening allowed
sigy = mateprop(6);

switch isw %Task Switch
%%    
    case 1
        
        if ndf > 2
            
            for i = 3:ndf
                lie(i,1) = 0;
            end
            
        end
        
        nh1 = (4+1)*4;
        
%%
    case 3
        
        bulk = ElemYM/(3*(1-2*Elemv));%Elemv*ElemYM/((1+Elemv)*(1-2*Elemv));
        mu = ElemYM/(2*(1+Elemv));
        One = [1; 1; 0];
        ind3to2 = [1 2 4];
        I4 = diag([1.0d0 1.0d0 1.0d0 0.5d0 0.5d0 0.5d0]);
        OneOne = ([1.0d0; 1.0d0; 1.0d0; 0.0d0; 0.0d0; 0.0d0]*[1.0d0; 1.0d0; 1.0d0; 0.0d0; 0.0d0; 0.0d0]');

        % Initialize Matrix and Vector

        nst = nel*ndf;
        ElemK = zeros(nst);
        ElemF = zeros(nst,1);
        % Load Guass Integration Points

%         if nel == 3 || nel == 6
%             lint = 13;
%         elseif nel == 4
            lint = 4;
% lint = 16;
%         else
%             lint = 9;
%         end
        der = 0;
        bf = 0;
        ib = 0;
%         if iprob == 5
%             fy = -1000*9.81;
%         else
%             fy = 0;
%         end

        triagNodes = [1 4 6
                      2 5 4
                      3 6 5
                      4 5 6];
%         Bmat = zeros(3,2*nel);
%         bvecT = zeros(1,2*nel);
        bbarvecT = zeros(1,2*nel);
        % [rlist, rWgts, rnum] = GaussPoints(2);
        % [slist, sWgts, snum] = GaussPoints(2);
        Gamvec = 1;
        H = 0;
%         Pdev = diag([1 1 1 1/2 1/2 1/2]) - 1/3*([1; 1; 1; 0; 0; 0]*[1 1 1 0 0 0]);
        
        % Compute bbarvec
        for l = 1:lint

            nodeT = triagNodes(l,1:3);
            xlT = xl(1:2,nodeT);
%                 if nel == 3 || nel == 6
              [Wgt,litr,lits] =  intpntt(1,1,ib);
              [shl,shld,shls,be] = shlt(litr,lits,3,3,der,bf);
              [Qxy, shgs, Jdet] = shgt(xlT,3,shld,shls,3,bf,der,be);
%                 else
%                   [Wgt,litr,lits] =  intpntq(l,lint,ib);
%                   [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
%                   [Qxy, shgs, Jdet] = shgq(xl,nel,shld,shls,nel,bf,der,be);
%                 end

            % Form b matrix
            bvecT = zeros(1,2*nel);
            for i = 1:3
                ie = nodeT(i);
                bvecT(1,(ie-1)*2+1:2*ie) = [Qxy(i,1) Qxy(i,2)];
            end

            % Update integration weighting factor
            W = Wgt*Jdet*thick;

            bbarvecT = bbarvecT + W*Gamvec*bvecT;
            H = H + W*(Gamvec*Gamvec');

        %     end %ie
        end %j
        
        % Loop over integration points
        for l = 1:lint

            nodeT = triagNodes(l,1:3);
            xlT = xl(1:2,nodeT);

              [Wgt,litr,lits] =  intpntt(1,1,ib);
              [shl,shld,shls,be] = shlt(litr,lits,3,3,der,bf);
              [Qxy, shgs, Jdet] = shgt(xlT,3,shld,shls,3,bf,der,be);

            % Form B matrix
            Bmat = zeros(3,2*nel);
            bvecT = zeros(1,2*nel);
            for i = 1:3
                ie = nodeT(i);
                Bmat(:,(ie-1)*2+1:2*ie) = [Qxy(i,1) 0       
                                           0        Qxy(i,2)
                                           Qxy(i,2) Qxy(i,1)];
                bvecT(1,(ie-1)*2+1:2*ie) = [Qxy(i,1) Qxy(i,2)];
            end
                
                Bbar = Bmat - 1/3*One*bvecT + 1/3*One*Gamvec'/H*bbarvecT;
                pn_1 = bulk*Gamvec'/H*bbarvecT*reshape(ul,nst,1);
                
                % Compute input for Radial Return
                deps2d = Bbar*reshape(ul-ul_n,nst,1); %3D enhanced strain
                deps_n1 = [deps2d(1); deps2d(2); deps2d(3); 0];
                deps3d = [deps2d(1); deps2d(2); 0; deps2d(3); 0; 0];
                ephr = nh1-1+(l-1)*5; %pointer for plastic strain at pt l
%                 betahr = nh1-1+(l-1)*7+3;
                ahr = nh1-1+l*5;
                ee_n = [hr(ephr+1); hr(ephr+2); hr(ephr+3); hr(ephr+4)];
                ee_n3d = [hr(ephr+1); hr(ephr+2); hr(ephr+4); hr(ephr+3); 0; 0];
%                 beta_n = [hr(betahr+1); hr(betahr+2); 0; hr(betahr+3); 0; 0];
%                 beta_n(3) = -(beta_n(1) + beta_n(2));
                a_n = hr(ahr);
                eps_n1 = ee_n + deps_n1;
                
%                 [DGAMA,LALGVA,RSTAVA,STRES] = SUVM(2,2,[ElemYM,Elemv],[0,100;sigy,Khard*100+sigy],[ee_n; a_n],eps_n1);
%                 C_n1 = CTVM(DGAMA,LALGVA(1),2,2,[ElemYM,Elemv],[0,100;sigy,Khard*100+sigy],RSTAVA,STRES);
% %                 [DGAM,LALGVA,RSTAVA,STRES] = SUTR(2,2,[ElemYM,Elemv],[0,100;sigy,Khard*100+sigy],[ee_n; a_n],eps_n1);
% %                 C_n1 = CTTR(LALGVA(1),2,LALGVA,2,[ElemYM,Elemv],[0,100;sigy,Khard*100+sigy],RSTAVA,eps_n1,STRES);
%                 
%                 % Convert output from 3D to 2D
%                 sig_n1 = [STRES(1); STRES(2); STRES(4); STRES(3); 0; 0];
%                 Cmat = C_n1(1:3,1:3);
                [sig_n1,C_n1,ee_n1,beta_n1,a_n1] = J2RadialReturnE(deps3d,ee_n3d,0*ee_n3d,a_n,mu,bulk,Khard,Hhard,sigy);
                Cmat = C_n1(ind3to2,ind3to2);
                sdev3 = sig_n1 - 1/3*OneOne*sig_n1;
                sdev2 = sdev3(ind3to2);
                
                % Store history variables
                ephr = nh2-1+(l-1)*5; %pointer for plastic strain at pt l
%                 betahr = nh2-1+(l-1)*7+3;
                ahr = nh2-1+l*5;
                hr(ephr+1) = ee_n1(1);
                hr(ephr+2) = ee_n1(2);
                hr(ephr+3) = ee_n1(4);
                hr(ephr+4) = ee_n1(3);
%                 hr(ephr+1) = RSTAVA(1);
%                 hr(ephr+2) = RSTAVA(2);
%                 hr(ephr+3) = RSTAVA(3);
%                 hr(ephr+4) = RSTAVA(4);
%                 hr(betahr+1) = beta_n1(1);
%                 hr(betahr+2) = beta_n1(2);
%                 hr(betahr+3) = beta_n1(4);
                hr(ahr) = a_n1;%RSTAVA(5);

                % Update integration weighting factor
                W = Wgt*Jdet*thick;

                ElemF = ElemF - W*Bbar'*(pn_1*One + sdev2);
                ElemK = ElemK + W*Bbar'*(Cmat - bulk*(One*One'))*Bbar;

        %     end %ie
        end %je
        
        Gambbar = One*bbarvecT;
        ElemK = ElemK + bulk/ndm*Gambbar'/H*Gambbar; %%% NOTE: this factor of one half does not work for general loading
        % You need the expanded arrays discussed by Hughes FEA Static and
        % Dynamic
% ElemK
    case 24
        
        ElemP = zeros(12,nel);
                
        bulk = ElemYM/(3*(1-2*Elemv));%Elemv*ElemYM/((1+Elemv)*(1-2*Elemv));
        mu = ElemYM/(2*(1+Elemv));
        One = [1; 1; 0];
        ind3to2 = [1 2 4];
        I4 = diag([1.0d0 1.0d0 1.0d0 0.5d0 0.5d0 0.5d0]);
        OneOne = ([1.0d0; 1.0d0; 1.0d0; 0.0d0; 0.0d0; 0.0d0]*[1.0d0; 1.0d0; 1.0d0; 0.0d0; 0.0d0; 0.0d0]');

        % Initialize Matrix and Vector

        nst = nel*ndf;
        % Load Guass Integration Points

%         if nel == 3 || nel == 6
%             lint = 13;
%         elseif nel == 4
            lint = 4;
% lint = 16;
%         else
%             lint = 9;
%         end
        der = 0;
        bf = 0;
        ib = 0;
%         if iprob == 5
%             fy = -1000*9.81;
%         else
%             fy = 0;
%         end

        Bmat = zeros(3,2*nel);
        bvecT = zeros(1,2*nel);
        bbarvecT = zeros(1,2*nel);
        % [rlist, rWgts, rnum] = GaussPoints(2);
        % [slist, sWgts, snum] = GaussPoints(2);
        Gamvec = 1;
        H = 0;
        Pdev = diag([1 1 1 1/2 1/2 1/2]) - 1/3*([1; 1; 1; 0; 0; 0]*[1 1 1 0 0 0]);
        
        % Compute bbarvec
        for l = 1:lint

                if nel == 3 || nel == 6
                  [Wgt,litr,lits] =  intpntt(l,lint,ib);
                  [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
                  [Qxy, shgs, Jdet] = shgt(xl,nel,shld,shls,nel,bf,der,be);
                else
                  [Wgt,litr,lits] =  intpntq(l,lint,ib);
                  [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
                  [Qxy, shgs, Jdet] = shgq(xl,nel,shld,shls,nel,bf,der,be);
                end

                % Form b matrix
                for ie = 1:nel
                bvecT(1,(ie-1)*2+1:2*ie) = [Qxy(ie,1) Qxy(ie,2)];
                end

                % Update integration weighting factor
                W = Wgt*Jdet*thick;

                bbarvecT = bbarvecT + W*Gamvec*bvecT;
                H = H + W*(Gamvec*Gamvec');

        %     end %ie
        end %j
        
        % Loop over integration points
        for l = 1:lint

                if nel == 3 || nel == 6
                  [Wgt,litr,lits] =  intpntt(l,lint,ib);
                  [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
                  [Qxy, shgs, Jdet] = shgt(xl,nel,shld,shls,nel,bf,der,be);
                else
                  [Wgt,litr,lits] =  intpntq(l,lint,ib);
                  [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
                  [Qxy, shgs, Jdet] = shgq(xl,nel,shld,shls,nel,bf,der,be);
                end

                % Form B matrix
                for ie = 1:nel
                Bmat(:,(ie-1)*2+1:2*ie) = [Qxy(ie,1) 0       
                                           0        Qxy(ie,2)
                                           Qxy(ie,2) Qxy(ie,1)];
                bvecT(1,(ie-1)*2+1:2*ie) = [Qxy(ie,1) Qxy(ie,2)];
                end
                
                Bbar = Bmat - 1/3*One*bvecT + 1/3*One*Gamvec'/H*bbarvecT;
                pn_1 = bulk*Gamvec'/H*bbarvecT*reshape(ul,nst,1);
                
                % Compute input for Radial Return
                deps2d = Bbar*reshape(ul-ul_n,nst,1); %3D enhanced strain
                deps_n1 = [deps2d(1); deps2d(2); deps2d(3); 0];
                deps3d = [deps2d(1); deps2d(2); 0; deps2d(3); 0; 0];
                ephr = nh1-1+(l-1)*5; %pointer for plastic strain at pt l
%                 betahr = nh1-1+(l-1)*7+3;
                ahr = nh1-1+l*5;
                ee_n = [hr(ephr+1); hr(ephr+2); hr(ephr+3); hr(ephr+4)];
                ee_n3d = [hr(ephr+1); hr(ephr+2); hr(ephr+4); hr(ephr+3); 0; 0];
%                 beta_n = [hr(betahr+1); hr(betahr+2); 0; hr(betahr+3); 0; 0];
%                 beta_n(3) = -(beta_n(1) + beta_n(2));
                a_n = hr(ahr);
                eps_n1 = ee_n + deps_n1;
                
%                 [DGAMA,LALGVA,RSTAVA,STRES] = SUVM(2,2,[ElemYM,Elemv],[0,100;sigy,Khard*100+sigy],[ee_n; a_n],eps_n1);
%                 C_n1 = CTVM(DGAMA,LALGVA(1),2,2,[ElemYM,Elemv],[0,100;sigy,Khard*100+sigy],RSTAVA,STRES);
% %                 [DGAM,LALGVA,RSTAVA,STRES] = SUTR(2,2,[ElemYM,Elemv],[0,100;sigy,Khard*100+sigy],[ee_n; a_n],eps_n1);
% %                 C_n1 = CTTR(LALGVA(1),2,LALGVA,2,[ElemYM,Elemv],[0,100;sigy,Khard*100+sigy],RSTAVA,eps_n1,STRES);
%                 
%                 % Convert output from 3D to 2D
%                 sig_n1 = [STRES(1); STRES(2); STRES(4); STRES(3); 0; 0];
%                 Cmat = C_n1(1:3,1:3);
                [sig_n1,C_n1,ee_n1,beta_n1,a_n1] = J2RadialReturnE(deps3d,ee_n3d,0*ee_n3d,a_n,mu,bulk,Khard,Hhard,sigy);
                Cmat = C_n1(ind3to2,ind3to2);
                sdev3 = sig_n1 - 1/3*OneOne*sig_n1;
                sdev2 = sdev3(ind3to2);
                
                ElemP(1,l) = pn_1 + sdev2(1);
                ElemP(2,l) = pn_1 + sdev2(2);
                ElemP(3,l) = sdev2(3);
                ElemP(4,l) = pn_1 + sdev3(3);
                ElemP(5,l) = a_n1;%RSTAVA(5);
%                 ElemP(6,l) = beta_n1(1);
%                 ElemP(7,l) = beta_n1(2);
%                 ElemP(8,l) = beta_n1(4);
%                 ElemP(9,l) = beta_n1(3);
                ElemP(10,l) = Cmat(1,1);
                ElemP(11,l) = Cmat(2,2);
                ElemP(12,l) = Cmat(3,3);

        end %je
%%
    case -1 % boundary tractions
        
        ElemF = zeros(nst,1);
        
        % Reorder nodes for corresponding face of integration
        SurfOrientEtoS

        triagNodes = [1 4 6
                      4 2 5];

        % Determine bounds of integration
        
        if nel == 4 || nel == 9
            
        dr = 2;
        ro = -1;
        
        % Upper Limit
%         if nodeA == ElemFlag(2)
            eR2 = 1;
%         else %
%             eR2 = 0;
%         end
        % Lower Limit
%         if nodeB == ElemFlag(1)
            eR1 = -1;
%         else %nodeA == ElemFlag(5)
%             eR1 = 0;
%         end
        
        elseif nel == 3 || nel == 6
            
        dr = 1;
        ro = 0;
            
        % Upper Limit
%         if nodeA == ElemFlag(2)
            eR2 = 1;
%         else %nodeA == ElemFlagR(5)
%             eR2 = 1/2;
%         end
        % Lower Limit
%         if nodeB == ElemFlag(1)
            eR1 = 0;
%         else %nodeA == ElemFlag(5)
%             eR1 = 1/2;
%         end
        
        end
        
        % Set jacobian for integration space
        drdr = (eR2 - eR1)/dr;
        
        lint = 2;
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
        ideriv = 0;
        der = 0;
        bf = 0;
        
        for je = 1:lint

            nodeT = triagNodes(1,1:3);
            xlT = xl(1:2,nodeT);

            if nel == 3 || nel == 6
                [Wgt,r,s] = intpntt(je,lint,1);
            else %if nel == 4
                [Wgt,r,s] = intpntq(je,lint,1);
            end
            r = drdr*(r-ro)+eR1;

              [shl,shld,shls,be] = shlt(r,s,3,3,der,bf);
              [shg, shgs, Jdet, be, xs] = shgt(xlT,3,shld,shls,3,bf,der,be);
                    
            %Evaluate tangent and normal vectors
            t1 = [xs(:,1); 0];
            [tm1, tu1] = VecNormalize(t1);
            t2 = [0; 0; 1];
            tm2 = 1;
            tu2 = t2';
            t3 = VecCrossProd(t1,t2);
            [tm3, tu3] = VecNormalize(t3);
            
            if edge < 0 % negative edge denotes pressure
                pressur = -traction(1); % pressure is defined positive in compression
                Traction = pressur*tu3(1:2);
            else
                Traction = traction(1:2);
            end
            
            c1 = Wgt*tm3*drdr*thick;

            Nmat = zeros(2,ndf*nel);
            
            for i = 1:3
                ie = nodeT(i);
                
              Nmat(1,(ie-1)*ndf+1) = shl(i);
              Nmat(2,(ie-1)*ndf+2) = shl(i);
                 
            end
            
            ElemF = ElemF + c1*Nmat'*Traction';

            nodeT = triagNodes(2,1:3);
            xlT = xl(1:2,nodeT);

            if nel == 3 || nel == 6
                [Wgt,r,s] = intpntt(je,lint,1);
            else %if nel == 4
                [Wgt,r,s] = intpntq(je,lint,1);
            end
            r = drdr*(r-ro)+eR1;

              [shl,shld,shls,be] = shlt(r,s,3,3,der,bf);
              [shg, shgs, Jdet, be, xs] = shgt(xlT,3,shld,shls,3,bf,der,be);
                    
            %Evaluate tangent and normal vectors
            t1 = [xs(:,1); 0];
            [tm1, tu1] = VecNormalize(t1);
            t2 = [0; 0; 1];
            tm2 = 1;
            tu2 = t2';
            t3 = VecCrossProd(t1,t2);
            [tm3, tu3] = VecNormalize(t3);
            
            if edge < 0 % negative edge denotes pressure
                pressur = -traction(1); % pressure is defined positive in compression
                Traction = pressur*tu3(1:2);
            else
                Traction = traction(1:2);
            end
            
            c1 = Wgt*tm3*drdr*thick;

            Nmat = zeros(2,ndf*nel);
            
            for i = 1:3
                ie = nodeT(i);
                
              Nmat(1,(ie-1)*ndf+1) = shl(i);
              Nmat(2,(ie-1)*ndf+2) = shl(i);
                 
            end
            
            ElemF = ElemF + c1*Nmat'*Traction';

        end %ie
        
        % Reorder nodes back to the orientation of the element
        SurfOrientStoE
        ElemF(1:ndf*nel) = ElemF(ilist2);
%%
    case 26
        
        bulk = ElemYM/(3*(1-2*Elemv));%Elemv*ElemYM/((1+Elemv)*(1-2*Elemv));
        mu = ElemYM/(2*(1+Elemv));
        One = [1; 1; 0];
        ind3to2 = [1 2 4];
        I4 = diag([1.0d0 1.0d0 1.0d0 0.5d0 0.5d0 0.5d0]);
        OneOne = ([1.0d0; 1.0d0; 1.0d0; 0.0d0; 0.0d0; 0.0d0]*[1.0d0; 1.0d0; 1.0d0; 0.0d0; 0.0d0; 0.0d0]');
        I1 = [1; 1; 0];

        % Initialize Matrix and Vector

        ElemS = zeros(nestr,1);
        % Load Guass Integration Points

%         if nel == 3 || nel == 6
%             lint = 13;
%         elseif nel == 4
            lint = 4;
% lint = 16;
%         else
%             lint = 9;
%         end
        der = 0;
        bf = 0;
        ib = 0;
        stresID = [1 2 3 0 1 3];
%         if iprob == 5
%             fy = -1000*9.81;
%         else
%             fy = 0;
%         end

        triagNodes = [1 4 6
                      2 5 4
                      3 6 5
                      4 5 6];
%         Bmat = zeros(3,2*nel);
%         bvecT = zeros(1,2*nel);
        bbarvecT = zeros(1,2*nel);
        % [rlist, rWgts, rnum] = GaussPoints(2);
        % [slist, sWgts, snum] = GaussPoints(2);
        Gamvec = 1;
        H = 0;
%         Pdev = diag([1 1 1 1/2 1/2 1/2]) - 1/3*([1; 1; 1; 0; 0; 0]*[1 1 1 0 0 0]);
        
        % Compute bbarvec
        for l = 1:lint

            nodeT = triagNodes(l,1:3);
            xlT = xl(1:2,nodeT);
%                 if nel == 3 || nel == 6
              [Wgt,litr,lits] =  intpntt(1,1,ib);
              [shl,shld,shls,be] = shlt(litr,lits,3,3,der,bf);
              [Qxy, shgs, Jdet] = shgt(xlT,3,shld,shls,3,bf,der,be);
%                 else
%                   [Wgt,litr,lits] =  intpntq(l,lint,ib);
%                   [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
%                   [Qxy, shgs, Jdet] = shgq(xl,nel,shld,shls,nel,bf,der,be);
%                 end

            % Form b matrix
            bvecT = zeros(1,2*nel);
            for i = 1:3
                ie = nodeT(i);
                bvecT(1,(ie-1)*2+1:2*ie) = [Qxy(i,1) Qxy(i,2)];
            end

            % Update integration weighting factor
            W = Wgt*Jdet*thick;

            bbarvecT = bbarvecT + W*Gamvec*bvecT;
            H = H + W*(Gamvec*Gamvec');

        %     end %ie
        end %j
        
        % Loop over integration points
        for l = 1:lint

            nodeT = triagNodes(l,1:3);
            xlT = xl(1:2,nodeT);

              [Wgt,litr,lits] =  intpntt(1,1,ib);
              [shl,shld,shls,be] = shlt(litr,lits,3,3,der,bf);
              [Qxy, shgs, Jdet] = shgt(xlT,3,shld,shls,3,bf,der,be);

            % Form B matrix
            Bmat = zeros(3,2*nel);
            bvecT = zeros(1,2*nel);
            for i = 1:3
                ie = nodeT(i);
                Bmat(:,(ie-1)*2+1:2*ie) = [Qxy(i,1) 0       
                                           0        Qxy(i,2)
                                           Qxy(i,2) Qxy(i,1)];
                bvecT(1,(ie-1)*2+1:2*ie) = [Qxy(i,1) Qxy(i,2)];
            end
                
                Bbar = Bmat - 1/3*One*bvecT + 1/3*One*Gamvec'/H*bbarvecT;
                pn_1 = bulk*Gamvec'/H*bbarvecT*reshape(ul,nst,1);
                
                % Compute input for Radial Return
                deps2d = Bbar*reshape(ul-ul_n,nst,1); %3D enhanced strain
                deps_n1 = [deps2d(1); deps2d(2); deps2d(3); 0];
                deps3d = [deps2d(1); deps2d(2); 0; deps2d(3); 0; 0];
                ephr = nh1-1+(l-1)*5; %pointer for plastic strain at pt l
%                 betahr = nh1-1+(l-1)*7+3;
                ahr = nh1-1+l*5;
                ee_n = [hr(ephr+1); hr(ephr+2); hr(ephr+3); hr(ephr+4)];
                ee_n3d = [hr(ephr+1); hr(ephr+2); hr(ephr+4); hr(ephr+3); 0; 0];
%                 beta_n = [hr(betahr+1); hr(betahr+2); 0; hr(betahr+3); 0; 0];
%                 beta_n(3) = -(beta_n(1) + beta_n(2));
                a_n = hr(ahr);
                eps_n1 = ee_n + deps_n1;
                
%                 [DGAMA,LALGVA,RSTAVA,STRES] = SUVM(2,2,[ElemYM,Elemv],[0,100;sigy,Khard*100+sigy],[ee_n; a_n],eps_n1);
%                 C_n1 = CTVM(DGAMA,LALGVA(1),2,2,[ElemYM,Elemv],[0,100;sigy,Khard*100+sigy],RSTAVA,STRES);
% %                 [DGAM,LALGVA,RSTAVA,STRES] = SUTR(2,2,[ElemYM,Elemv],[0,100;sigy,Khard*100+sigy],[ee_n; a_n],eps_n1);
% %                 C_n1 = CTTR(LALGVA(1),2,LALGVA,2,[ElemYM,Elemv],[0,100;sigy,Khard*100+sigy],RSTAVA,eps_n1,STRES);
%                 
%                 % Convert output from 3D to 2D
%                 sig_n1 = [STRES(1); STRES(2); STRES(4); STRES(3); 0; 0];
%                 Cmat = C_n1(1:3,1:3);
                [sig_n1,C_n1,ee_n1,beta_n1,a_n1] = J2RadialReturnE(deps3d,ee_n3d,0*ee_n3d,a_n,mu,bulk,Khard,Hhard,sigy);
                Cmat = C_n1(ind3to2,ind3to2);
                sdev3 = sig_n1 - 1/3*OneOne*sig_n1;
                sdev2 = sdev3(ind3to2);
                sigma = (pn_1*One + sdev2);
            
            for stres = 1:npstr
            
            if stres <= 3 % stress components
                sigmas = sigma(stresID(stres));
            elseif stres >= 5
                if stres <= 6 % principal stresses
%                     if PSPS == 'n'
                        sigz = bulk*(epsil(1)+epsil(2));
%                     else
%                         sigz = 0;
%                     end
                    sigma2 = [sigma(1) sigma(3) 0; sigma(3) sigma(2) 0; 0 0 sigz];
                    psig = eig(sigma2);
                    sigmas = psig(stresID(stres));
                else % hydrostatic stress
%                     if PSPS == 'n'
                        sigz = bulk*(epsil(1)+epsil(2));
%                     else
%                         sigz = 0;
%                     end
                    sigmas = 1/3*(sigma'*I1 + sigz);
                end
            else % von Mises stress
%                 if PSPS == 'n'
                    sigz = bulk*(epsil(1)+epsil(2));
%                 else
%                     sigz = 0;
%                 end
                trs = sigma'*I1 + sigz;
                dsig = sigma - 1/3*trs*I1;
                sigmas = sqrt(3/2*(dsig'*dsig));
            end
            
            ElemS(stres) = ElemS(stres) + sigmas;
            
            end

        %     end %ie
        end %je
        
        ElemS = ElemS/4;
        
end %Task Switch