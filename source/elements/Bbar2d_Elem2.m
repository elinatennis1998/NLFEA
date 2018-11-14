% Tim Truster
% CEE 577
% 04/18/2012
% Midterm Exam
% Plane strain B-bar element for small strain J2-Plasticity with linear
% isotrop hardening
% Bilinear quadrilateral element with consistent tangent
% Built using deSouza functions; agrees with Hughes version exactly
% No modifications in SUVM were necessary; only the input had to be
% specified properly

% Does not work for general loading; see comments below.
% Shear seems to work

% Verified Tresca from deSouza works for single or double return; using
% large steps for the shear problem, the return to corner is encountered.

% Verified against NL Mixed Inelasticity on 6/30/2013 using IME_Q2_2d.m

% Set Material Properties

ElemYM = mateprop(2);
Elemv = mateprop(3);
thick = mateprop(1);
Khard = mateprop(4);
Hhard = mateprop(5);
sigy = mateprop(6);

plasmodel = 1;

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

        Bmat = zeros(3,2*nel);
        bvecT = zeros(1,2*nel);
        bbarvecT = zeros(1,2*nel);
        % [rlist, rWgts, rnum] = GaussPoints(2);
        % [slist, sWgts, snum] = GaussPoints(2);
        Gamvec = 1;
        H = 0;
%         Pdev = diag([1 1 1 1/2 1/2 1/2]) - 1/3*([1; 1; 1; 0; 0; 0]*[1 1 1 0 0 0]);
        
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
                ephr = nh1-1+(l-1)*5; %pointer for plastic strain at pt l
%                 betahr = nh1-1+(l-1)*7+3;
                ahr = nh1-1+l*5;
                ee_n = [hr(ephr+1); hr(ephr+2); hr(ephr+3); hr(ephr+4)];
%                 beta_n = [hr(betahr+1); hr(betahr+2); 0; hr(betahr+3); 0; 0];
%                 beta_n(3) = -(beta_n(1) + beta_n(2));
                a_n = hr(ahr);
                eps_n1 = ee_n + deps_n1;
                
                if plasmodel == 1
                [DGAMA,LALGVA,RSTAVA,STRES] = SUVM(2,2,[ElemYM,Elemv],[0,100;sigy,Khard*100+sigy],[ee_n; a_n],eps_n1);
                C_n1 = CTVM(DGAMA,LALGVA(1),2,2,[ElemYM,Elemv],[0,100;sigy,Khard*100+sigy],RSTAVA,STRES);
                elseif plasmodel == 2
                [DGAM,LALGVA,RSTAVA,STRES] = SUTR(2,2,[ElemYM,Elemv],[0,100;sigy,Khard*100+sigy],[ee_n; a_n],eps_n1);
                C_n1 = CTTR(LALGVA(1),2,LALGVA,2,[ElemYM,Elemv],[0,100;sigy,Khard*100+sigy],RSTAVA,eps_n1,STRES);
                end
                
                % Convert output from 3D to 2D
                sig_n1 = [STRES(1); STRES(2); STRES(4); STRES(3); 0; 0];
                sdev3 = sig_n1 - 1/3*OneOne*sig_n1;
                sdev2 = sdev3(ind3to2);
                Cmat = C_n1(1:3,1:3);
                
                % Store history variables
                ephr = nh2-1+(l-1)*5; %pointer for plastic strain at pt l
%                 betahr = nh2-1+(l-1)*7+3;
                ahr = nh2-1+l*5;
                hr(ephr+1) = RSTAVA(1);
                hr(ephr+2) = RSTAVA(2);
                hr(ephr+3) = RSTAVA(3);
                hr(ephr+4) = RSTAVA(4);
%                 hr(betahr+1) = beta_n1(1);
%                 hr(betahr+2) = beta_n1(2);
%                 hr(betahr+3) = beta_n1(4);
                hr(ahr) = RSTAVA(5);

                % Update integration weighting factor
                W = Wgt*Jdet*thick;

                ElemF = ElemF - W*Bmat'*(pn_1*One + sdev2);
                ElemK = ElemK + W*Bmat'*(Cmat - bulk*(One*One'))*Bmat;

        %     end %ie
        end %je
        
        Gambbar = One*bbarvecT;
        ElemK = ElemK + bulk/ndm*Gambbar'/H*Gambbar; %%% NOTE: this factor of one half does not work for general loading
        % You need the expanded arrays discussed by Hughes FEA Static and
        % Dynamic
% ElemK
    case 22
        
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
                ephr = nh1-1+(l-1)*5; %pointer for plastic strain at pt l
%                 betahr = nh1-1+(l-1)*7+3;
                ahr = nh1-1+l*5;
                ee_n = [hr(ephr+1); hr(ephr+2); hr(ephr+3); hr(ephr+4)];
%                 beta_n = [hr(betahr+1); hr(betahr+2); 0; hr(betahr+3); 0; 0];
%                 beta_n(3) = -(beta_n(1) + beta_n(2));
                a_n = hr(ahr);
                eps_n1 = ee_n + deps_n1;
                
                if plasmodel == 1
                [DGAMA,LALGVA,RSTAVA,STRES] = SUVM(2,2,[ElemYM,Elemv],[0,100;sigy,Khard*100+sigy],[ee_n; a_n],eps_n1);
                C_n1 = CTVM(DGAMA,LALGVA(1),2,2,[ElemYM,Elemv],[0,100;sigy,Khard*100+sigy],RSTAVA,STRES);
                elseif plasmodel == 2
                [DGAM,LALGVA,RSTAVA,STRES] = SUTR(2,2,[ElemYM,Elemv],[0,100;sigy,Khard*100+sigy],[ee_n; a_n],eps_n1);
                C_n1 = CTTR(LALGVA(1),2,LALGVA,2,[ElemYM,Elemv],[0,100;sigy,Khard*100+sigy],RSTAVA,eps_n1,STRES);
                end
                
                % Convert output from 3D to 2D
                sig_n1 = [STRES(1); STRES(2); STRES(4); STRES(3); 0; 0];
                sdev3 = sig_n1 - 1/3*OneOne*sig_n1;
                sdev2 = sdev3(ind3to2);
                Cmat = C_n1(1:3,1:3);
                
                ElemP(1,l) = pn_1 + sdev2(1);
                ElemP(2,l) = pn_1 + sdev2(2);
                ElemP(3,l) = sdev2(3);
                ElemP(4,l) = pn_1 + sdev3(3);
                ElemP(5,l) = RSTAVA(5);
%                 ElemP(6,l) = beta_n1(1);
%                 ElemP(7,l) = beta_n1(2);
%                 ElemP(8,l) = beta_n1(4);
%                 ElemP(9,l) = beta_n1(3);
                ElemP(10,l) = Cmat(1,1);
                ElemP(11,l) = Cmat(2,2);
                ElemP(12,l) = Cmat(3,3);

        end %je
        
end %Task Switch