% Tim Truster
% 05/11/2014
% Plane strain Q2Q1 for small strain J2-Plasticity with linear
% isotropic/kinematic hardening
% Biquadratic-bilinear mixed quadrilateral element with consistent tangent

% Tested using MT52U9M1, and results agree with MT51U4M1

% Element is hard-coded as a Q2Q1 element, namely the number of pressure
% nodes is assumed to be 4.

% Set Material Properties

ElemYM = mateprop(2);
Elemv = mateprop(3);
thick = mateprop(1);
Khard = mateprop(4);
Hhard = mateprop(5);
sigy = mateprop(6);
if length(mateprop)>6
    plasmodel = mateprop(7);
    if length(mateprop)>7
        plasversion = mateprop(8);
    else
        plasversion = 1;2;
    end
else
    plasmodel = 1;
    plasversion = 1;2;
end

switch isw %Task Switch
    
    case 1
        
        if ndf > 3
            
            for i = 4:ndf
                lie(i,1) = 0;
            end
            
        end
        
        nh1 = 7*13;(nen+1);
        istv = 12;8;

    case 3
        
        bulk = ElemYM/(3*(1-2*Elemv));%Elemv*ElemYM/((1+Elemv)*(1-2*Elemv));
        mu = ElemYM/(2*(1+Elemv));
        One = [1; 1; 0];
        I2 = [1; 1; 1; 0; 0; 0];
        ind3to2 = [1 2 4];
        I4 = diag([1.0d0 1.0d0 1.0d0 0.5d0 0.5d0 0.5d0]);
        OneOne = ([1.0d0; 1.0d0; 1.0d0; 0.0d0; 0.0d0; 0.0d0]*[1.0d0; 1.0d0; 1.0d0; 0.0d0; 0.0d0; 0.0d0]');

        % Initialize Matrix and Vector
        ElemK = zeros(nst);
        ElemF = zeros(nst,1);
        
        % Load Guass Integration Points
        if nel == 9
            lint = 9;
        elseif nel == 6
            lint = 7;
        else
            error('bad element')
        end

        der = 0;
        bf = 0;
        ib = 0;

        Bmat = zeros(4,nst);
        Pdev = diag([1 1 1 1/2 1/2 1/2]) - 1/3*([1; 1; 1; 0; 0; 0]*[1 1 1 0 0 0]);
        
        % Loop over integration points
        for l = 1:lint

                if nel == 3 || nel == 6
                  [Wgt,litr,lits] =  intpntt(l,lint,ib);
                  [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
                  [Qxy, shgs, Jdet] = shgt(xl,nel,shld,shls,nel,bf,der,be);
                  [shlp,shld,shls,bub] = shlt(litr,lits,3,nel,0,0);
                else
                  [Wgt,litr,lits] =  intpntq(l,lint,ib);
                  [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
                  [Qxy, shgs, Jdet] = shgq(xl,nel,shld,shls,nel,bf,der,be);
                  [shlp,shld,shls,bub] = shlq(litr,lits,4,nel,0,0);
                end

                % Form B matrix
                for ie = 1:nel
                Bmat(:,(ie-1)*3+1:3*ie) = [Qxy(ie,1) 0        0
                                           0        Qxy(ie,2) 0
                                           Qxy(ie,2) Qxy(ie,1) 0
                                           0 0 shlp(ie)];
                end
                
                xint = xl(:,1:nel)*shl;
                
                
                if plasversion == 1
                % Compute input for Radial Return
                du = Bmat*reshape(ul,nst,1);
                pn_1 = du(4);
                
                eps2d = du(1:3); %2D enhanced strain
                theta_n1 = One'*eps2d - pn_1/bulk;
                eps3d = [eps2d(1); eps2d(2); 0; eps2d(3); 0; 0];
                ephr = nh1-1+(l-1)*7; %pointer for plastic strain at pt l
                betahr = nh1-1+(l-1)*7+3;
                ahr = nh1-1+l*7;
                ep_n = [hr(ephr+1); hr(ephr+2); 0; hr(ephr+3); 0; 0];
                ep_n(3) = -(ep_n(1) + ep_n(2));
                beta_n = [hr(betahr+1); hr(betahr+2); 0; hr(betahr+3); 0; 0];
                beta_n(3) = -(beta_n(1) + beta_n(2));
                a_n = hr(ahr);
                
%                 [sig_n1,C_n1,ep_n1,beta_n1,a_n1,sdev3,Cdev_n1] = J2RadialReturn(eps3d,ep_n,beta_n,a_n,mu,bulk,Khard,Hhard,sigy);
                [sig_n1,C_n1,ep_n1,beta_n1,a_n1,sdev3,Cdev_n1] = J2RadialReturn0(eps3d,ep_n,beta_n,a_n,mu,bulk,Khard,Hhard,sigy);
                
                % Convert output from 3D to 2D
                sdev2 = sdev3(ind3to2);
                Cmat = Cdev_n1(ind3to2,ind3to2);
                
                % Store history variables
                ephr = nh2-1+(l-1)*7; %pointer for plastic strain at pt l
                betahr = nh2-1+(l-1)*7+3;
                ahr = nh2-1+l*7;
                hr(ephr+1) = ep_n1(1);
                hr(ephr+2) = ep_n1(2);
                hr(ephr+3) = ep_n1(4);
                hr(betahr+1) = beta_n1(1);
                hr(betahr+2) = beta_n1(2);
                hr(betahr+3) = beta_n1(4);
                hr(ahr) = a_n1;
                
                else
                du = Bmat*reshape(ul,nst,1);
                pn_1 = du(4);
                
                eps2d = du(1:3); %2D enhanced strain
                theta_n1 = One'*eps2d - pn_1/bulk;
                deps2d = Bmat*reshape(ul-ul_n,nst,1); %3D enhanced strain
                deps_n1 = [deps2d(1); deps2d(2); deps2d(3); 0];
                ephr = nh1-1+(l-1)*7; %pointer for plastic strain at pt l
%                 betahr = nh1-1+(l-1)*7+3;
                ahr = nh1-1+l*7;
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
%                 if DGAM ~= 0 && elem < 8
%                     elem, xint
%                 end
                
                % Convert output from 3D to 2D
                sig_n1 = [STRES(1); STRES(2); STRES(4); STRES(3); 0; 0];
                sdev3 = sig_n1 - 1/3*OneOne*sig_n1;
                sdev2 = sdev3(ind3to2);
                Cmat = C_n1(1:3,1:3) - bulk*[1 1 0]'*[1 1 0];
                
                % Store history variables
                ephr = nh2-1+(l-1)*7; %pointer for plastic strain at pt l
%                 betahr = nh2-1+(l-1)*7+3;
                ahr = nh2-1+l*7;
                hr(ephr+1) = RSTAVA(1);
                hr(ephr+2) = RSTAVA(2);
                hr(ephr+3) = RSTAVA(3);
                hr(ephr+4) = RSTAVA(4);
%                 hr(betahr+1) = beta_n1(1);
%                 hr(betahr+2) = beta_n1(2);
%                 hr(betahr+3) = beta_n1(4);
                hr(ahr) = RSTAVA(5);
                end

                % Update integration weighting factor
                W = Wgt*Jdet*thick;

                ElemF = ElemF - W*Bmat'*[pn_1*One + sdev2; theta_n1];
                ElemK = ElemK + W*Bmat'*[[Cmat One]; [One' -1/bulk]]*Bmat;

        %     end %ie
        end %je
ElemK;

    case 6
        
        bulk = ElemYM/(3*(1-2*Elemv));%Elemv*ElemYM/((1+Elemv)*(1-2*Elemv));
        mu = ElemYM/(2*(1+Elemv));
        One = [1; 1; 0];
        ind3to2 = [1 2 4];
        I4 = diag([1.0d0 1.0d0 1.0d0 0.5d0 0.5d0 0.5d0]);
        OneOne = ([1.0d0; 1.0d0; 1.0d0; 0.0d0; 0.0d0; 0.0d0]*[1.0d0; 1.0d0; 1.0d0; 0.0d0; 0.0d0; 0.0d0]');

        % Initialize Matrix and Vector
        ElemF = zeros(nst,1);
        
        % Load Guass Integration Points
        if nel == 9
            lint = 9;
        elseif nel == 6
            lint = 7;
        else
            error('bad element')
        end

        der = 0;
        bf = 0;
        ib = 0;

        Bmat = zeros(4,nst);
        Pdev = diag([1 1 1 1/2 1/2 1/2]) - 1/3*([1; 1; 1; 0; 0; 0]*[1 1 1 0 0 0]);
        
        % Loop over integration points
        for l = 1:lint

                if nel == 3 || nel == 6
                  [Wgt,litr,lits] =  intpntt(l,lint,ib);
                  [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
                  [Qxy, shgs, Jdet] = shgt(xl,nel,shld,shls,nel,bf,der,be);
                  [shlp,shld,shls,bub] = shlt(litr,lits,3,nel,0,0);
                else
                  [Wgt,litr,lits] =  intpntq(l,lint,ib);
                  [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
                  [Qxy, shgs, Jdet] = shgq(xl,nel,shld,shls,nel,bf,der,be);
                  [shlp,shld,shls,bub] = shlq(litr,lits,4,nel,0,0);
                end

                % Form B matrix
                for ie = 1:nel
                Bmat(:,(ie-1)*3+1:3*ie) = [Qxy(ie,1) 0        0
                                           0        Qxy(ie,2) 0
                                           Qxy(ie,2) Qxy(ie,1) 0
                                           0 0 shlp(ie)];
                end
                
                
                if plasversion == 1
                % Compute input for Radial Return
                du = Bmat*reshape(ul,nst,1);
                pn_1 = du(4);
                
                eps2d = du(1:3); %2D enhanced strain
                theta_n1 = One'*eps2d - pn_1/bulk;
                eps3d = [eps2d(1); eps2d(2); 0; eps2d(3); 0; 0];
                ephr = nh1-1+(l-1)*7; %pointer for plastic strain at pt l
                betahr = nh1-1+(l-1)*7+3;
                ahr = nh1-1+l*7;
                ep_n = [hr(ephr+1); hr(ephr+2); 0; hr(ephr+3); 0; 0];
                ep_n(3) = -(ep_n(1) + ep_n(2));
                beta_n = [hr(betahr+1); hr(betahr+2); 0; hr(betahr+3); 0; 0];
                beta_n(3) = -(beta_n(1) + beta_n(2));
                a_n = hr(ahr);
                
%                 [sig_n1,C_n1,ep_n1,beta_n1,a_n1,sdev3,Cdev_n1] = J2RadialReturn(eps3d,ep_n,beta_n,a_n,mu,bulk,Khard,Hhard,sigy);
                [sig_n1,C_n1,ep_n1,beta_n1,a_n1,sdev3,Cdev_n1] = J2RadialReturn0(eps3d,ep_n,beta_n,a_n,mu,bulk,Khard,Hhard,sigy);
                
                % Convert output from 3D to 2D
                sdev2 = sdev3(ind3to2);
                
                
                % Store history variables
                ephr = nh2-1+(l-1)*7; %pointer for plastic strain at pt l
                betahr = nh2-1+(l-1)*7+3;
                ahr = nh2-1+l*7;
                hr(ephr+1) = ep_n1(1);
                hr(ephr+2) = ep_n1(2);
                hr(ephr+3) = ep_n1(4);
                hr(betahr+1) = beta_n1(1);
                hr(betahr+2) = beta_n1(2);
                hr(betahr+3) = beta_n1(4);
                hr(ahr) = a_n1;
                
                else
                du = Bmat*reshape(ul,nst,1);
                pn_1 = du(4);
                
                eps2d = du(1:3); %2D enhanced strain
                theta_n1 = One'*eps2d - pn_1/bulk;
                deps2d = Bmat*reshape(ul-ul_n,nst,1); %3D enhanced strain
                deps_n1 = [deps2d(1); deps2d(2); deps2d(3); 0];
                ephr = nh1-1+(l-1)*7; %pointer for plastic strain at pt l
%                 betahr = nh1-1+(l-1)*7+3;
                ahr = nh1-1+l*7;
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
                
                % Store history variables
                ephr = nh2-1+(l-1)*7; %pointer for plastic strain at pt l
%                 betahr = nh2-1+(l-1)*7+3;
                ahr = nh2-1+l*7;
                hr(ephr+1) = RSTAVA(1);
                hr(ephr+2) = RSTAVA(2);
                hr(ephr+3) = RSTAVA(3);
                hr(ephr+4) = RSTAVA(4);
%                 hr(betahr+1) = beta_n1(1);
%                 hr(betahr+2) = beta_n1(2);
%                 hr(betahr+3) = beta_n1(4);
                hr(ahr) = RSTAVA(5);
                end

                % Update integration weighting factor
                W = Wgt*Jdet*thick;

                ElemF = ElemF - W*Bmat'*[pn_1*One + sdev2; theta_n1];

        %     end %ie
        end %je
%%
    case -1 % boundary tractions
        
        ElemF = zeros(nst,1);
        
        % Reorder nodes for corresponding face of integration
        SurfOrientEtoS

        % Determine bounds of integration
        
        if nel == 4 || nel == 9
            
        dr = 2;
        ro = -1;
        
        % Upper Limit
        if nodeA == ElemFlag(2)
            eR2 = 1;
        else %
            eR2 = 0;
        end
        % Lower Limit
        if nodeB == ElemFlag(1)
            eR1 = -1;
        else %nodeA == ElemFlag(5)
            eR1 = 0;
        end
        
        elseif nel == 3 || nel == 6
            
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
        
        end
        
        % Set jacobian for integration space
        drdr = (eR2 - eR1)/dr;
        
%         if exist('iprob','var') == 1 && iprob == 6
%             lint = 10;
%         else
            lint = 4; % Use 10 for body force BF2U4M0.m
%         end
        
        ideriv = 0;
        der = 0;
        bf = 0;
        thick = 1;
        
        for je = 1:lint

            if nel == 3 || nel == 6
                [Wgt,r,s] = intpntt(je,lint,1);
            elseif nel == 4 || nel == 9
                [Wgt,r,s] = intpntq(je,lint,1);
            end
            r = drdr*(r-ro)+eR1;

            if nel == 3 || nel == 6
              [shl,shld,shls,be] = shlt(r,s,nel,nel,der,bf);
              [QXY, shgs, Jdet, be, xs] = shgt(xl,nel,shld,shls,nel,bf,der,be);
            else
              [shl,shld,shls,be] = shlq(r,s,nel,nel,der,bf);
              [QXY, shgs, Jdet, be, xs] = shgq(xl,nel,shld,shls,nel,bf,der,be);
            end
                    
            %Evaluate tangent and normal vectors
            t1 = [xs(:,1); 0];
            [tm1, tu1] = VecNormalize(t1);
            t2 = [0; 0; 1];
            tm2 = 1;
            tu2 = t2';
            t3 = VecCrossProd(t1,t2);
            [tm3, tu3] = VecNormalize(t3);
            
            if exist('iprob','var')
                if iprob == 6
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
                elseif edge < 0 % negative edge denotes pressure
                    pressur = -traction(1); % pressure is defined positive in compression
                    Traction = pressur*tu3;
                else
                    Traction = traction;
                end
            elseif edge < 0 % negative edge denotes pressure
                pressur = -traction(1); % pressure is defined positive in compression
                Traction = pressur*tu3;
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

                ElemF(ndf*o-2) = ElemF(ndf*o-2) + F(1)*c1;

                ElemF(ndf*o-1)   = ElemF(ndf*o-1)   + F(2)*c1;

            end %o

        end %ie
        
        ElemF;
        
        % Reorder nodes back to the orientation of the element
        SurfOrientStoE
        ElemF(1:ndf*nel) = ElemF(ilist2);
        
    case 13 % plastic dissipation
        
        ElemD = 0;
        
        bulk = ElemYM/(3*(1-2*Elemv));%Elemv*ElemYM/((1+Elemv)*(1-2*Elemv));
        mu = ElemYM/(2*(1+Elemv));
        One = [1; 1; 1; 0; 0; 0];
        I1 = [1; 1; 0];
        ind3to2 = [1 2 4];
        I4 = diag([1.0d0 1.0d0 1.0d0 0.5d0 0.5d0 0.5d0]);
        OneOne = ([1.0d0; 1.0d0; 1.0d0; 0.0d0; 0.0d0; 0.0d0]*[1.0d0; 1.0d0; 1.0d0; 0.0d0; 0.0d0; 0.0d0]');
        
        % Load Guass Integration Points
        if nel == 9
            lint = 9;
        elseif nel == 6
            lint = 7;
        else
            error('bad element')
        end
        der = 0;
        bf = 0;
        ib = 0;

        Bmat = zeros(4,nst);
        
        % Loop over integration points
        for l = 1:lint

                if nel == 3 || nel == 6
                  [Wgt,litr,lits] =  intpntt(l,lint,ib);
                  [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
                  [Qxy, shgs, Jdet] = shgt(xl,nel,shld,shls,nel,bf,der,be);
                  [shlp,shld,shls,bub] = shlt(litr,lits,3,nel,0,0);
                else
                  [Wgt,litr,lits] =  intpntq(l,lint,ib);
                  [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
                  [Qxy, shgs, Jdet] = shgq(xl,nel,shld,shls,nel,bf,der,be);
                  [shlp,shld,shls,bub] = shlq(litr,lits,4,nel,0,0);
                end
                
                % Update integration weighting factor
                W = Wgt*Jdet*thick;
                
                % Form B matrix
                for ie = 1:nel
                Bmat(:,(ie-1)*3+1:3*ie) = [Qxy(ie,1) 0        0
                                           0        Qxy(ie,2) 0
                                           Qxy(ie,2) Qxy(ie,1) 0
                                           0 0 shlp(ie)];
                end

                if plasversion == 1
                % Compute input for Radial Return
                du = Bmat*reshape(ul,nst,1);
                pn_1 = du(4);
                
                eps2d = du(1:3); %2D enhanced strain
                theta_n1 = I1'*eps2d - pn_1/bulk;
                eps3d = [eps2d(1); eps2d(2); 0; eps2d(3); 0; 0];
                ephr = nh1-1+(l-1)*7; %pointer for plastic strain at pt l
                betahr = nh1-1+(l-1)*7+3;
                ahr = nh1-1+l*7;
                ep_n = [hr(ephr+1); hr(ephr+2); 0; hr(ephr+3); 0; 0];
                ep_n(3) = -(ep_n(1) + ep_n(2));
                beta_n = [hr(betahr+1); hr(betahr+2); 0; hr(betahr+3); 0; 0];
                beta_n(3) = -(beta_n(1) + beta_n(2));
                a_n = hr(ahr);
                
                [sig_n1,C_n1,ep_n1,beta_n1,a_n1,sdev3,Cdev_n1] = J2RadialReturn0(eps3d,ep_n,beta_n,a_n,mu,bulk,Khard,Hhard,sigy);

                % mixed stress tensor
                sig_n1 = pn_1*One + sdev3;
                
                % compute incremental plastic dissipation
                eps_n1 = [ep_n1(1); ep_n1(2); -(ep_n1(1) + ep_n1(2)); ep_n1(4); 0; 0]; %step n+1
                eps_n = [ep_n(1); ep_n(2); -(ep_n(1) + ep_n(2)); ep_n(4); 0; 0]; % step n
                beta_n1(3) = -(beta_n1(1) + beta_n1(2));
                q_n1 = [-Khard*a_n1; beta_n1; beta_n1(4:6)];
                q_n = [-Khard*a_n; beta_n; beta_n1(4:6)];
                ElemD = ElemD + W*((eps_n1 - eps_n)'*sig_n1);
                if Khard > 0
                ElemD = ElemD + W*( - q_n1(1)*(q_n1(1) - q_n(1))/Khard);
                end
                if Hhard > 0
                ElemD = ElemD + W*(- sum(q_n1(2:10).*(q_n1(2:10) - q_n(2:10)))/(2/3*Hhard));
                end

                else
                du = Bmat*reshape(ul,nst,1);
                pn_1 = du(4);
                
                eps2d = du(1:3); %2D enhanced strain
                theta_n1 = I1'*eps2d - pn_1/bulk;
                deps2d = Bmat*reshape(ul-ul_n,nst,1); %3D enhanced strain
                deps_n1 = [deps2d(1); deps2d(2); deps2d(3); 0];
                ephr = nh1-1+(l-1)*7; %pointer for plastic strain at pt l
%                 betahr = nh1-1+(l-1)*7+3;
                ahr = nh1-1+l*7;
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
                sig_n1 = pn_1*One + sdev3;
                
                
                % compute incremental plastic dissipation
                epp = eps2d - RSTAVA(1:3);
                eps_n1 = [epp(1); epp(2); -RSTAVA(4); epp(3); 0; 0]; %step n+1
                es_n = Bmat*reshape(ul_n,nst,1);
                epp = es_n(1:3) - ee_n(1:3);
                eps_n = [epp(1); epp(2); -ee_n(4); epp(3); 0; 0]; % step n
                q_n1 = [-Khard*RSTAVA(5)];
                q_n = [-Khard*a_n];
                % VERIFY these factors for the equivalent plastic strain
                % for deSouza; I did not check it on 6/4/14 because the
                % problem is perfectly plastic.
                % ALSO: the PLFUN seems to give the wrong value for the
                % initial data point in SUVM; check that out.
                ElemD = ElemD + W*((eps_n1 - eps_n)'*sig_n1);
                if Khard > 0
                ElemD = ElemD + W*( - q_n1(1)*(q_n1(1) - q_n(1))/Khard);
                end
                
                end

        %     end %ie
        end %je
%%       
    case 11
        
        ElemE = zeros(numEn,1);
        
        % Load Gauss Points for quadrature
        if nel == 9
            lint = 9;
        elseif nel == 6
            lint = 7;
        else
            error('bad element')
        end
        pfact = 1;

        ideriv = 1;
        fbx = 0;
        fby = 0;
        hx = 0;
        hy = 0;
        ib = 0;
        bf = 1;
        der = 1;

        el2el = zeros(3,1);
        eprixel = zeros(3,1);
        epriyel = zeros(3,1);
        el2fine = zeros(7,1);
        ue = zeros(3,1);
        duex = ue;
        duey = ue;
        
        if elem == 1 && exist('iprob','var') && iprob == 1
            if  sigy < inf
            % Determine c
            options = optimoptions('fsolve','Display','off');
            cPPC = fsolve(@(c)(sigy/2*(1-c^2/Ro^2 + log(c^2/Ri^2))-pPPC/stepmax*step),Ri,options); % eq (4.3.27)
            else
                cPPC = Ro;
            end
        end

        for je = 1:lint

                if nel == 3 || nel == 6
                    [Wgt,r,s] = intpntt(je,lint,0);
                elseif nel == 4 || nel == 9
                    [Wgt,r,s] = intpntq(je,lint,0);
                end

                if nel == 3 || nel == 6
                    [shl,shld,shls,be] = shlt(r,s,nel,nel,der,bf);
                    [shg, shgs, Jdet, be, xs] = shgt(xl,nel,shld,shls,nen,bf,der,be);
                elseif nel == 4 || nel == 9
                    [shl,shld,shls,be] = shlq(r,s,nel,nel,der,bf);
                    [shg, shgs, Jdet, be, xs] = shgq(xl,nel,shld,shls,nen,bf,der,be);
                end
                c1 = Wgt*Jdet*thick;
                b = be(3);

                if nelP == 3 || nelP == 6
                    [shlp,shld,shls,bub] = shlt(r,s,nelP,nel,0,0);
                    [shp, shp2, Jdet, bub, xs] = shgt(xl,nelP,shld,shls,nen,0,0,bub);
                elseif nelP == 4 || nelP == 9
                    [shlp,shld,shls,bub] = shlq(r,s,nelP,nel,0,0);
                    [shp, shp2, Jdet, bub, xs] = shgq(xl,nelP,shld,shls,nen,0,0,bub);
                end

                % Displacement Terms of Stiffness Matrix
%                 for jn = 1:nel2

                xint = xl(1,1:nel)*shl;
                yint = xl(2,1:nel)*shl;
                dux = ul(1:2,1:nel)*shg(:,1);
                duy = ul(1:2,1:nel)*shg(:,2);
                u = ul(1:2,1:nel)*shl;
                dux(3) = ul(3,1:nelP)*shp(:,1)*pfact;
                duy(3) = ul(3,1:nelP)*shp(:,2)*pfact;
                u(3) = ul(3,1:nelP)*shlp*pfact;
                px = dux(3);
                py = duy(3);

                %Compute value of exact fields at int. point
                if exist('iprob','var') && iprob == 1
                    if mod(elem,nu/2) <= nv/4 && mod(elem,nu/2) > 0 %%|| 1%mod(elem,nu/2) <= nv/16 %  
    %                 [ue,duex,duey] = uexact_thickcylE(xint,yint,0,ElemYM,Elemv,5,15,.32*lamda,0);
                    [ue,duex,duey] = uexact_thickcylEP(xint,yint,0,ElemYM,Elemv,sigy/2,Ri,Ro,cPPC,pPPC/stepmax*step,0);
                    duey(1) = 0;
                    dux(2) = 1/2*(dux(2) + duy(1)); % strain component
                    duy(1) = 0;
                    else
                        u = 0*u;
                        dux = 0*dux;
                        duy = 0*duy;
                    end
                else
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
        Ieffvals(elem,1:3) = [sqrt(1/H1u) 0 H1u];

    case 24
        
        ElemP = zeros(12,nel);
                
        bulk = ElemYM/(3*(1-2*Elemv));%Elemv*ElemYM/((1+Elemv)*(1-2*Elemv));
        mu = ElemYM/(2*(1+Elemv));
        One = [1; 1; 0];
        ind3to2 = [1 2 4];
        I4 = diag([1.0d0 1.0d0 1.0d0 0.5d0 0.5d0 0.5d0]);
        OneOne = ([1.0d0; 1.0d0; 1.0d0; 0.0d0; 0.0d0; 0.0d0]*[1.0d0; 1.0d0; 1.0d0; 0.0d0; 0.0d0; 0.0d0]');
        
        % Load Guass Integration Points
        if nel == 9
            lint = 9;
        elseif nel == 6
            lint = 7;
        else
            error('bad element')
        end

        der = 0;
        bf = 0;
        ib = 0;

        Bmat = zeros(4,2*nel);
        Pdev = diag([1 1 1 1/2 1/2 1/2]) - 1/3*([1; 1; 1; 0; 0; 0]*[1 1 1 0 0 0]);
        
        % Loop over integration points
        for l = 1:lint

                if nel == 3 || nel == 6
                  [Wgt,litr,lits] =  intpntt(l,lint,ib);
                  [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
                  [Qxy, shgs, Jdet] = shgt(xl,nel,shld,shls,nel,bf,der,be);
                  [shlp,shld,shls,bub] = shlt(litr,lits,3,nel,0,0);
                else
                  [Wgt,litr,lits] =  intpntq(l,lint,ib);
                  [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
                  [Qxy, shgs, Jdet] = shgq(xl,nel,shld,shls,nel,bf,der,be);
                  [shlp,shld,shls,bub] = shlq(litr,lits,4,nel,0,0);
                end

                % Form B matrix
                for ie = 1:nel
                Bmat(:,(ie-1)*3+1:3*ie) = [Qxy(ie,1) 0        0
                                           0        Qxy(ie,2) 0
                                           Qxy(ie,2) Qxy(ie,1) 0
                                           0 0 shlp(ie)];
                end
                
                
%                 if plasversion == 1
                % Compute input for Radial Return
                du = Bmat*reshape(ul,nst,1);
                pn_1 = du(4);
                
                eps2d = du(1:3); %2D enhanced strain
                theta_n1 = One'*eps2d - pn_1/bulk;
                eps3d = [eps2d(1); eps2d(2); 0; eps2d(3); 0; 0];
                ephr = nh1-1+(l-1)*7; %pointer for plastic strain at pt l
                betahr = nh1-1+(l-1)*7+3;
                ahr = nh1-1+l*7;
                ep_n = [hr(ephr+1); hr(ephr+2); 0; hr(ephr+3); 0; 0];
                ep_n(3) = -(ep_n(1) + ep_n(2));
                beta_n = [hr(betahr+1); hr(betahr+2); 0; hr(betahr+3); 0; 0];
                beta_n(3) = -(beta_n(1) + beta_n(2));
                a_n = hr(ahr);
                
%                 [sig_n1,C_n1,ep_n1,beta_n1,a_n1,sdev3,Cdev_n1] = J2RadialReturn(eps3d,ep_n,beta_n,a_n,mu,bulk,Khard,Hhard,sigy);
                [sig_n1,C_n1,ep_n1,beta_n1,a_n1,sdev3,Cdev_n1] = J2RadialReturn0(eps3d,ep_n,beta_n,a_n,mu,bulk,Khard,Hhard,sigy);
                
                % Convert output from 3D to 2D
                sdev2 = sdev3(ind3to2);
                Cmat = Cdev_n1(ind3to2,ind3to2);
                
                ElemP(1,l) = pn_1 + sdev2(1);
                ElemP(2,l) = pn_1 + sdev2(2);
                ElemP(3,l) = sdev2(3);
                ElemP(4,l) = pn_1 + sdev3(3);
                ElemP(5,l) = a_n1;
                ElemP(6,l) = beta_n1(1);
                ElemP(7,l) = beta_n1(2);
                ElemP(8,l) = beta_n1(4);
                ElemP(9,l) = beta_n1(3);
                ElemP(10,l) = Cmat(1,1);
                ElemP(11,l) = Cmat(2,2);
                ElemP(12,l) = Cmat(3,3);

        end %je
        
%%        
    case 25 %Stress Projection2

% stress calculation updated 06/03/14 to use stresses at integration points
% for extrapolating
% Stresses for T6 will be zero because of lint=7; don't know how to handle
% shape functions from integration points.
        
        ElemS = zeros(nel,npstr+1);
        ElemS2 = zeros(nel,npstr);
        stresID = [1 2 3 0 1 3];

        bulk = ElemYM/(3*(1-2*Elemv));%Elemv*ElemYM/((1+Elemv)*(1-2*Elemv));
        mu = ElemYM/(2*(1+Elemv));
        thick = 1;
        fbx = 0;
        fby = 0;
        fbz = 0;
        One = [1; 1; 1; 0; 0; 0];
        I1 = [1; 1; 0];
        ind3to2 = [1 2 4];
        
        % Load Guass Integration Points

        if nel == 3
            lint = 1;
            nint = 1;
        elseif nel == 4
%             lint = 4;
            lint = 4;
            nint = 4;
        elseif nel == 6
            lint = 7;
            nint = 6;
        else
            lint = 9;
            nint = 9;
        end
        
        der = 0;
        bf = 0;
        ib = 0;
        
        ulres = reshape(ul,nst,1);

        Bmat = zeros(4,nst);
        
        % Loop over integration points
        for l = 1:nint

                if nel == 3 || nel == 6
                  [Wgt,litr,lits] =  intpntt(l,lint,ib);
                  [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
                  [Qxy, shgs, Jdet] = shgt(xl,nel,shld,shls,nel,bf,der,be);
                  [shlp,shld,shls,bub] = shlt(litr,lits,3,nel,0,0);
                else
                  [Wgt,litr,lits] =  intpntq(l,lint,ib);
                  [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
                  [Qxy, shgs, Jdet] = shgq(xl,nel,shld,shls,nel,bf,der,be);
                  [shlp,shld,shls,bub] = shlq(litr,lits,4,nel,0,0);
                end

                % Form B matrix
                for ie = 1:nel
                Bmat(:,(ie-1)*3+1:3*ie) = [Qxy(ie,1) 0        0
                                           0        Qxy(ie,2) 0
                                           Qxy(ie,2) Qxy(ie,1) 0
                                           0 0 shlp(ie)];
                end
                
                if plasversion == 1
                    
                % Compute input for Radial Return
                du = Bmat*reshape(ul,nst,1);
                pn_1 = du(4);
                
                eps2d = du(1:3); %2D enhanced strain
                theta_n1 = I1'*eps2d - pn_1/bulk;
                eps3d = [eps2d(1); eps2d(2); 0; eps2d(3); 0; 0];
                ephr = nh1-1+(l-1)*7; %pointer for plastic strain at pt l
                betahr = nh1-1+(l-1)*7+3;
                ahr = nh1-1+l*7;
                ep_n = [hr(ephr+1); hr(ephr+2); 0; hr(ephr+3); 0; 0];
                ep_n(3) = -(ep_n(1) + ep_n(2));
                beta_n = [hr(betahr+1); hr(betahr+2); 0; hr(betahr+3); 0; 0];
                beta_n(3) = -(beta_n(1) + beta_n(2));
                a_n = hr(ahr);
                
%                 [sig_n1,C_n1,ep_n1,beta_n1,a_n1,sdev3,Cdev_n1] = J2RadialReturn(eps3d,ep_n,beta_n,a_n,mu,bulk,Khard,Hhard,sigy);
                [sig_n1,C_n1,ep_n1,beta_n1,a_n1,sdev3,Cdev_n1] = J2RadialReturn0(eps3d,ep_n,beta_n,a_n,mu,bulk,Khard,Hhard,sigy);
                
                % Convert output from 3D to 2D
                sig_n1 = pn_1*One + sdev3;
                sigma = sig_n1(ind3to2);
                
                else
                % Compute input for Radial Return
                du = Bmat*reshape(ul,nst,1);
                pn_1 = du(4);
                
                eps2d = du(1:3); %2D enhanced strain
                deps2d = Bmat*reshape(ul-ul_n,nst,1); %3D enhanced strain
                deps_n1 = [deps2d(1); deps2d(2); deps2d(3); 0];
                ephr = nh1-1+(l-1)*7; %pointer for plastic strain at pt l
%                 betahr = nh1-1+(l-1)*7+3;
                ahr = nh1-1+l*7;
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
                a_n1 = RSTAVA(5);
                
                % Convert output from 3D to 2D
                sig_n1 = [STRES(1); STRES(2); STRES(4); STRES(3); 0; 0];
                sdev3 = sig_n1 - 1/3*OneOne*sig_n1;
                sig_n1 = pn_1*One + sdev3;
                sigma = sig_n1(ind3to2);
                end

        %Stress Loop
            
            for stres = 1:npstr
            
            if stres <= 3 % stress components
                sigmas = sigma(stresID(stres));
            elseif stres >= 5
                if stres <= 6 % principal stresses
                    sigz = sig_n1(3);
                    sigma2 = [sigma(1) sigma(3) 0; sigma(3) sigma(2) 0; 0 0 sigz];
                    psig = eig(sigma2);
                    sigmas = psig(stresID(stres));
                elseif stres == 8 % equivalent plastic strain
                    sigmas = a_n1;
                elseif stres == 7 % hydrostatic stress
                    sigz = sig_n1(3);
                    sigmas = 1/3*(sigma'*I1 + sigz);
                else
                    if plasversion == 1
                    if stres >= 9 && stres <=11
                        sigmas = beta_n1(stresID(stres-8));
                    else % von Mises stress
                        sigz = beta_n1(3);
                        trs = beta_n1(ind3to2)'*I1 + sigz;
                        dsig = beta_n1(ind3to2) - 1/3*trs*I1;
                        dsig(4) = sigz - 1/3*trs;
                        sigmas = sqrt(3/2*(dsig'*dsig + dsig(3)^2)); % Correct calculation of vonMises stress
                    end
                    end
                end
            else % von Mises stress
                sigz = sig_n1(3);
                trs = sigma'*I1 + sigz;
                dsig = sigma - 1/3*trs*I1;
                dsig(4) = sigz - 1/3*trs;
                sigmas = sqrt(3/2*(dsig'*dsig + dsig(3)^2)); % Correct calculation of vonMises stress
            end
            
            ElemS2(l,stres) = sigmas;
            
            end

        end %je
        
        % interpolate stress at nodes
        if nel == 3
            plist = [0 1 0
                     0 0 1];
        elseif nel == 4
            sqr3 = sqrt(3);
            plist = [-sqr3 sqr3 sqr3 -sqr3
                     -sqr3 -sqr3 sqr3 sqr3];
        elseif nel == 6
            plist = [-1/3 5/3 -1/3 2/3 2/3 -1/3
                     -1/3 -1/3 5/3 -1/3 2/3 2/3];
        else
            sqr3 = sqrt(5/3);
            plist = [-sqr3 sqr3 sqr3 -sqr3 0 sqr3 0 -sqr3 0
                     -sqr3 -sqr3 sqr3 sqr3 -sqr3 0 sqr3 0 0];
        end
        
        for ll = 1:nel
            
            r = plist(1,ll);
            s = plist(2,ll);
            shpS = sshp2d(r,s,nint);
            
%             for stres = 1:npstr
%                 
%                 sigmas = ElemS2(1:nint,stres)'*shpS;
%                 ElemS(ll,stres) = sigmas;
%                 
%             end
            ElemS(ll,1:npstr) = (ElemS2(1:nint,:)'*shpS)';
            
        end
        
%         %Integration Loop
%         Vol = 0;
%         for ll = 1:lint
% 
%             %Evaluate first derivatives of basis functions at int. point
%             if nel == 3 || nel == 6
%               [Wgt,litr,lits] =  intpntt(ll,lint,ib);
%               [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
%               [Qxy, shgs, Jdet,be] = shgt(xl+ul(1:2,:),nel,shld,shls,nel,bf,der,be);
%             else
%               [Wgt,litr,lits] =  intpntq(ll,lint,ib);
%               [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
%               [Qxy, shgs, Jdet,be] = shgq(xl+ul(1:2,:),nel,shld,shls,nel,bf,der,be);
%             end
% 
%             w = Wgt*Jdet*thick;
%             
%             Vol = Vol + w;
% 
%         end %je

        for i = 1:nel
        ElemS(i,npstr+1) = 1; % use this for simple average Vol; % use this for weighted average 
        end

end %Task Switch