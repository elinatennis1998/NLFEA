% Tim Truster
% 05/10/2015
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

switch isw %Task Switch
    
    case 1
        
        if ndf > 4
            
            for i = 5:ndf
                lie(i,1) = 0;
            end
            
        end
        
        nh1 = 13*27;(nen+1);
        istv = 12;8;

    case 3
        
        bulk = ElemYM/(3*(1-2*Elemv));%Elemv*ElemYM/((1+Elemv)*(1-2*Elemv));
        mu = ElemYM/(2*(1+Elemv));
        One = [1; 1; 1; 0; 0; 0];
        I2 = [1; 1; 1; 0; 0; 0];
        I4 = diag([1.0d0 1.0d0 1.0d0 0.5d0 0.5d0 0.5d0]);
        OneOne = ([1.0d0; 1.0d0; 1.0d0; 0.0d0; 0.0d0; 0.0d0]*[1.0d0; 1.0d0; 1.0d0; 0.0d0; 0.0d0; 0.0d0]');

        % Initialize Matrix and Vector
        ElemK = zeros(nel*ndf);
        ElemF = zeros(nel*ndf,1);
        
        % Load Guass Integration Points
        if nel == 27
            lint = 27;
        elseif nel == 10
            lint = 14;
        else
            error('bad element')
        end

        der = 0;
        bf = 0;
        ib = 0;

        Bmat = zeros(7,nel*ndf);
        Pdev = diag([1 1 1 1/2 1/2 1/2]) - 1/3*([1; 1; 1; 0; 0; 0]*[1 1 1 0 0 0]);
        ulrshp = reshape(ul,nst,1);
        
        % Loop over integration points
        for l = 1:lint

                if nel == 4 || nel == 10
                  [Wgt,ss] =  int3d_t(l,lint,ib);
                  [shl,shld,shls,be] = shltt(ss,nel,nel,der,bf);
                  [Qxy, shgs, Jdet] = shgtt(xl,nel,shld,shls,nel,bf,der,be);
                  [shlp,shld,shls,bub] = shltt(ss,4,nel,0,0);
                else
                  [Wgt,ss] =  intpntb(l,lint,ib);
                  [shl,shld,shls,be] = shlb(ss,nel,nel,der,bf);
                  [Qxy, shgs, Jdet] = shgb(xl,nel,shld,shls,nel,bf,der,be);
                  [shlp,shld,shls,bub] = shlb(ss,8,nel,0,0);
                end

                % Form B matrix
                for mm = 1:nel
                Bmat(:,4*mm-3:4*mm) = [Qxy(mm,1) 0         0        0 
                        0         Qxy(mm,2) 0         0
                        0         0         Qxy(mm,3) 0
                        Qxy(mm,2) Qxy(mm,1) 0         0
                        0         Qxy(mm,3) Qxy(mm,2) 0
                        Qxy(mm,3) 0         Qxy(mm,1) 0
                                           0 0 0 shlp(mm)];
                end
                
                xint = xl(:,1:nel)*shl;
                
                % Compute input for Radial Return
                du = Bmat*ulrshp(1:nel*ndf);
                pn_1 = du(7);
                
                eps3d = du(1:6); %2D enhanced strain
                theta_n1 = One'*eps3d - pn_1/bulk;
                ephr = nh1-1+(l-1)*13; %pointer for plastic strain at pt l
                betahr = nh1-1+(l-1)*13+6;
                ahr = nh1-1+l*13;
                ep_n = [hr(ephr+1); hr(ephr+2); hr(ephr+3); hr(ephr+4); hr(ephr+5); hr(ephr+6)];
                beta_n = [hr(betahr+1); hr(betahr+2); hr(betahr+3); hr(betahr+4); hr(betahr+5); hr(betahr+6)];
                a_n = hr(ahr);
                
%                 [sig_n1,C_n1,ep_n1,beta_n1,a_n1,sdev3,Cdev_n1] = J2RadialReturn(eps3d,ep_n,beta_n,a_n,mu,bulk,Khard,Hhard,sigy);
                [sig_n1,C_n1,ep_n1,beta_n1,a_n1,sdev3,Cmat] = J2RadialReturn0(eps3d,ep_n,beta_n,a_n,mu,bulk,Khard,Hhard,sigy);
                
                % Store history variables
                ephr = nh2-1+(l-1)*13; %pointer for plastic strain at pt l
                betahr = nh2-1+(l-1)*13+6;
                ahr = nh2-1+l*13;
                hr(ephr+1:ephr+6) = ep_n1(1:6);
                hr(betahr+1:betahr+6) = beta_n1(1:6);
                hr(ahr) = a_n1;

                % Update integration weighting factor
                W = Wgt*Jdet*thick;

                ElemF = ElemF - W*Bmat'*[pn_1*One + sdev3; theta_n1];
                ElemK = ElemK + W*Bmat'*[[Cmat One]; [One' -1/bulk]]*Bmat;

        %     end %ie
        end %je
ElemK;

    case 6
        
        bulk = ElemYM/(3*(1-2*Elemv));%Elemv*ElemYM/((1+Elemv)*(1-2*Elemv));
        mu = ElemYM/(2*(1+Elemv));
        One = [1; 1; 1; 0; 0; 0];
        I2 = [1; 1; 1; 0; 0; 0];
        I4 = diag([1.0d0 1.0d0 1.0d0 0.5d0 0.5d0 0.5d0]);
        OneOne = ([1.0d0; 1.0d0; 1.0d0; 0.0d0; 0.0d0; 0.0d0]*[1.0d0; 1.0d0; 1.0d0; 0.0d0; 0.0d0; 0.0d0]');

        % Initialize Matrix and Vector
        ElemK = zeros(nst);
        ElemF = zeros(nst,1);
        
        % Load Guass Integration Points
        if nel == 27
            lint = 27;
        elseif nel == 10
            lint = 14;
        else
            error('bad element')
        end

        der = 0;
        bf = 0;
        ib = 0;

        Bmat = zeros(7,nst);
        Pdev = diag([1 1 1 1/2 1/2 1/2]) - 1/3*([1; 1; 1; 0; 0; 0]*[1 1 1 0 0 0]);
        
        % Loop over integration points
        for l = 1:lint

                if nel == 4 || nel == 10
                  [Wgt,ss] =  int3d_t(l,lint,ib);
                  [shl,shld,shls,be] = shltt(ss,nel,nel,der,bf);
                  [Qxy, shgs, Jdet] = shgtt(xl,nel,shld,shls,nel,bf,der,be);
                  [shlp,shld,shls,bub] = shltt(ss,4,nel,0,0);
                else
                  [Wgt,ss] =  intpntb(l,lint,ib);
                  [shl,shld,shls,be] = shlb(ss,nel,nel,der,bf);
                  [Qxy, shgs, Jdet] = shgb(xl,nel,shld,shls,nel,bf,der,be);
                  [shlp,shld,shls,bub] = shlb(ss,8,nel,0,0);
                end

                % Form B matrix
                for mm = 1:nel
                Bmat(:,4*mm-3:4*mm) = [Qxy(mm,1) 0         0        0 
                        0         Qxy(mm,2) 0         0
                        0         0         Qxy(mm,3) 0
                        Qxy(mm,2) Qxy(mm,1) 0         0
                        0         Qxy(mm,3) Qxy(mm,2) 0
                        Qxy(mm,3) 0         Qxy(mm,1) 0
                                           0 0 0 shlp(mm)];
                end
                
                xint = xl(:,1:nel)*shl;
                
                % Compute input for Radial Return
                du = Bmat*reshape(ul,nst,1);
                pn_1 = du(7);
                
                eps3d = du(1:6); %2D enhanced strain
                theta_n1 = One'*eps3d - pn_1/bulk;
                ephr = nh1-1+(l-1)*13; %pointer for plastic strain at pt l
                betahr = nh1-1+(l-1)*13+6;
                ahr = nh1-1+l*13;
                ep_n = [hr(ephr+1); hr(ephr+2); hr(ephr+3); hr(ephr+4); hr(ephr+5); hr(ephr+6)];
                beta_n = [hr(betahr+1); hr(betahr+2); hr(betahr+3); hr(betahr+4); hr(betahr+5); hr(betahr+6)];
                a_n = hr(ahr);
                
%                 [sig_n1,C_n1,ep_n1,beta_n1,a_n1,sdev3,Cdev_n1] = J2RadialReturn(eps3d,ep_n,beta_n,a_n,mu,bulk,Khard,Hhard,sigy);
                [sig_n1,C_n1,ep_n1,beta_n1,a_n1,sdev3,Cmat] = J2RadialReturn0(eps3d,ep_n,beta_n,a_n,mu,bulk,Khard,Hhard,sigy);
                
                % Store history variables
                ephr = nh2-1+(l-1)*13; %pointer for plastic strain at pt l
                betahr = nh2-1+(l-1)*13+6;
                ahr = nh2-1+l*13;
                hr(ephr+1:ephr+6) = ep_n1(1:6);
                hr(betahr+1:betahr+6) = beta_n1(1:6);
                hr(ahr) = a_n1;

                % Update integration weighting factor
                W = Wgt*Jdet*thick;

                ElemF = ElemF - W*Bmat'*[pn_1*One + sdev3; theta_n1];
                ElemK = ElemK + W*Bmat'*[[Cmat One]; [One' -1/bulk]]*Bmat;

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
        I2 = [1; 1; 1; 0; 0; 0];
        I4 = diag([1.0d0 1.0d0 1.0d0 0.5d0 0.5d0 0.5d0]);
        OneOne = ([1.0d0; 1.0d0; 1.0d0; 0.0d0; 0.0d0; 0.0d0]*[1.0d0; 1.0d0; 1.0d0; 0.0d0; 0.0d0; 0.0d0]');
        
        % Load Guass Integration Points
        if nel == 27
            lint = 27;
        elseif nel == 10
            lint = 14;
        else
            error('bad element')
        end

        der = 0;
        bf = 0;
        ib = 0;

        Bmat = zeros(7,nst);
        Pdev = diag([1 1 1 1/2 1/2 1/2]) - 1/3*([1; 1; 1; 0; 0; 0]*[1 1 1 0 0 0]);
        
        % Loop over integration points
        for l = 1:lint

                if nel == 4 || nel == 10
                  [Wgt,ss] =  int3d_t(l,lint,ib);
                  [shl,shld,shls,be] = shltt(ss,nel,nel,der,bf);
                  [Qxy, shgs, Jdet] = shgtt(xl,nel,shld,shls,nel,bf,der,be);
                  [shlp,shld,shls,bub] = shltt(ss,4,nel,0,0);
                else
                  [Wgt,ss] =  intpntb(l,lint,ib);
                  [shl,shld,shls,be] = shlb(ss,nel,nel,der,bf);
                  [Qxy, shgs, Jdet] = shgb(xl,nel,shld,shls,nel,bf,der,be);
                  [shlp,shld,shls,bub] = shlb(ss,8,nel,0,0);
                end

                % Form B matrix
                for mm = 1:nel
                Bmat(:,4*mm-3:4*mm) = [Qxy(mm,1) 0         0        0 
                        0         Qxy(mm,2) 0         0
                        0         0         Qxy(mm,3) 0
                        Qxy(mm,2) Qxy(mm,1) 0         0
                        0         Qxy(mm,3) Qxy(mm,2) 0
                        Qxy(mm,3) 0         Qxy(mm,1) 0
                                           0 0 0 shlp(mm)];
                end
                
                xint = xl(:,1:nel)*shl;
                
                % Compute input for Radial Return
                du = Bmat*reshape(ul,nst,1);
                pn_1 = du(7);
                
                eps3d = du(1:6); %2D enhanced strain
                theta_n1 = One'*eps3d - pn_1/bulk;
                ephr = nh1-1+(l-1)*13; %pointer for plastic strain at pt l
                betahr = nh1-1+(l-1)*13+6;
                ahr = nh1-1+l*13;
                ep_n = [hr(ephr+1); hr(ephr+2); hr(ephr+3); hr(ephr+4); hr(ephr+5); hr(ephr+6)];
                beta_n = [hr(betahr+1); hr(betahr+2); hr(betahr+3); hr(betahr+4); hr(betahr+5); hr(betahr+6)];
                a_n = hr(ahr);
                
%                 [sig_n1,C_n1,ep_n1,beta_n1,a_n1,sdev3,Cdev_n1] = J2RadialReturn(eps3d,ep_n,beta_n,a_n,mu,bulk,Khard,Hhard,sigy);
                [sig_n1,C_n1,ep_n1,beta_n1,a_n1,sdev3,Cmat] = J2RadialReturn0(eps3d,ep_n,beta_n,a_n,mu,bulk,Khard,Hhard,sigy);
                
                
                % compute incremental plastic dissipation
                eps_n1 = ep_n1; %step n+1
                eps_n1(3) = -(eps_n1(1) + eps_n1(2));
                eps_n = ep_n; % step n
                beta_n1(3) = -(beta_n1(1) + beta_n1(2));
                q_n1 = [-Khard*a_n1; beta_n1];
                q_n = [-Khard*a_n; beta_n];
                ElemD = ElemD + W*((eps_n1 - eps_n)'*sig_n1);
                if Khard > 0
                ElemD = ElemD + W*( - q_n1(1)*(q_n1(1) - q_n(1))/Khard);
                end
                if Hhard > 0
                ElemD = ElemD + W*(- sum(q_n1(2:10).*(q_n1(2:10) - q_n(2:10)))/(2/3*Hhard));
                end
                

        %     end %ie
        end %je
        
%%       
    case 11
        
        ElemE = zeros(numEn,1);
        
        % Load Gauss Points for quadrature
        if nel == 27
            lint = 27;
        elseif nel == 10
            lint = 14;
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

        el2el = zeros(4,1);
        eprixel = zeros(4,1);
        epriyel = zeros(4,1);
        eprizel = zeros(4,1);
        el2fine = zeros(13,1);
        ue = zeros(4,1);
        duex = ue;
        duey = ue;
        
        if elem == 1 && exist('iprob','var') && iprob == 1
            pPPC = presses(step);
            if  sigy < inf
            % Determine c
            options = optimoptions('fsolve','Display','off','TolFun',1e-16);
            cPPC = fsolve(@(c)(sigy*2/3*(1-c^3/Ro^3 + log(c^3/Ri^3))-pPPC),Ri,options); % eq (4.3.6), with typo fixed
            else
                cPPC = Ro;
            end
        end
        
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

            if nelP == 4 || nelP == 10
                [shlp,shld,shls,be] = shltt(ss,nelP,nel,der,bf);
                [shp, shp2, Jdet2, bub, xs] = shgtt(xl,nelP,shld,shls,nen,0,0,be);
            elseif nelP == 8 || nelP == 27
                [shlp,shld,shls,be] = shlb(ss,nelP,nel,0,0);
                [shp, shp2, Jdet2, bub, xs] = shgq(xl,nelP,shld,shls,nen,0,0,be);
            end

            c1 = Jdet*w*thick;

            xint = xl(1,1:nel)*shl;
            yint = xl(2,1:nel)*shl;
            zint = xl(3,1:nel)*shl;
            dux = ul(:,1:nel)*shg(:,1);
            duy = ul(:,1:nel)*shg(:,2);
            duz = ul(:,1:nel)*shg(:,3);
            u = ul(:,1:nel)*shl;
            dux(4) = ul(4,1:nelP)*shp(:,1)*pfact;
            duy(4) = ul(4,1:nelP)*shp(:,2)*pfact;
            duz(4) = ul(4,1:nelP)*shp(:,3)*pfact;
            u(4) = ul(4,1:nelP)*shlp*pfact;
            px = dux(4);
            py = duy(4);
            pz = duz(4);

            %Compute value of exact fields at int. point
            if exist('iprob','var') && iprob == 1
%                 [ue,duex,duey] = uexact_thickcylE(xint,yint,0,ElemYM,Elemv,5,15,.32*lamda,0);
                [ue,duex,duey,duez] = uexact_thicksphEP(xint,yint,zint,ElemYM,Elemv,sigy,Ri,Ro,cPPC,pPPC,0);
                duey(1) = 0;
                dux(2) = 1/2*(dux(2) + duy(1)); % strain component
                duy(1) = 0;
                duez(2) = 0;
                duy(3) = 1/2*(duy(3) + duz(2)); % strain component
                duz(2) = 0;
                duex(3) = 0;
                duz(1) = 1/2*(duz(1) + dux(3)); % strain component
                dux(3) = 0;
            else
            ue = zeros(4,1);duex=ue;duey=ue;duez=ue;
            end

            %Add standard int. point error to element standard error
            for in = 1:4
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

        for in= 1:4
            ElemE(in) = el2el(in);
            ElemE(in+4) = eprixel(in);
            ElemE(in+8) = epriyel(in);
            ElemE(in+12) = eprizel(in);
        end
        
        H1u = eprixel(1)+eprixel(2)+eprixel(3)+epriyel(1)+epriyel(2)+epriyel(3)+eprizel(1)+eprizel(2)+eprizel(3);
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
        I1 = One;
        
        % Load Guass Integration Points

        if nel == 3
            lint = 1;
            nint = 1;
        elseif nel == 4
%             lint = 4;
            lint = 4;
            nint = 4;
        elseif nel == 10
            lint = 14;
            nint = 10;
        else
            lint = 27;
            nint = 27;
        end
        
        der = 0;
        bf = 0;
        ib = 0;
        
        ulres = reshape(ul,nst,1);

        Bmat = zeros(7,nst);
        
        % Loop over integration points
        for l = 1:nint

                if nel == 4 || nel == 10
                  [Wgt,ss] =  int3d_t(l,lint,ib);
                  [shl,shld,shls,be] = shltt(ss,nel,nel,der,bf);
                  [Qxy, shgs, Jdet] = shgtt(xl,nel,shld,shls,nel,bf,der,be);
                  [shlp,shld,shls,bub] = shltt(ss,4,nel,0,0);
                else
                  [Wgt,ss] =  intpntb(l,lint,ib);
                  [shl,shld,shls,be] = shlb(ss,nel,nel,der,bf);
                  [Qxy, shgs, Jdet] = shgb(xl,nel,shld,shls,nel,bf,der,be);
                  [shlp,shld,shls,bub] = shlb(ss,8,nel,0,0);
                end

                % Form B matrix
                for mm = 1:nel
                Bmat(:,4*mm-3:4*mm) = [Qxy(mm,1) 0         0        0 
                        0         Qxy(mm,2) 0         0
                        0         0         Qxy(mm,3) 0
                        Qxy(mm,2) Qxy(mm,1) 0         0
                        0         Qxy(mm,3) Qxy(mm,2) 0
                        Qxy(mm,3) 0         Qxy(mm,1) 0
                                           0 0 0 shlp(mm)];
                end
                
                xint = xl(:,1:nel)*shl;
                
                % Compute input for Radial Return
                du = Bmat*reshape(ul,nst,1);
                pn_1 = du(7);
                
                eps3d = du(1:6); %2D enhanced strain
                theta_n1 = One'*eps3d - pn_1/bulk;
                ephr = nh1-1+(l-1)*13; %pointer for plastic strain at pt l
                betahr = nh1-1+(l-1)*13+6;
                ahr = nh1-1+l*13;
                ep_n = [hr(ephr+1); hr(ephr+2); hr(ephr+3); hr(ephr+4); hr(ephr+5); hr(ephr+6)];
                beta_n = [hr(betahr+1); hr(betahr+2); hr(betahr+3); hr(betahr+4); hr(betahr+5); hr(betahr+6)];
                a_n = hr(ahr);
                
%                 [sig_n1,C_n1,ep_n1,beta_n1,a_n1,sdev3,Cdev_n1] = J2RadialReturn(eps3d,ep_n,beta_n,a_n,mu,bulk,Khard,Hhard,sigy);
                [sig_n1,C_n1,ep_n1,beta_n1,a_n1,sdev3,Cmat] = J2RadialReturn0(eps3d,ep_n,beta_n,a_n,mu,bulk,Khard,Hhard,sigy);
                
                sigma = pn_1*One + sdev3;

        %Stress Loop
            
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
            
            ElemS2(l,stres) = sigmas;
            
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
%         Vol = 0;
%         for ll = 1:lint
% 
%             %Evaluate first derivatives of basis functions at int. point
%             if nel == 4 || nel == 10
%               [Wgt,ss] =  int3d_t(ll,lint,ib);
%               [shl,shld,shls,be] = shltt(ss,nel,nel,der,bf);
%               [Qxy,shgs,Jdet,be] = shgtt(xl+ul(1:3,:),nel,shld,shls,nel,bf,der,be);
%             else
%               [Wgt,ss] =  intpntb(ll,lint,ib);
%               [shl,shld,shls,be] = shlb(ss,nel,nel,der,bf);
%               [Qxy,shgs,Jdet,be] = shgb(xl+ul(1:3,:),nel,shld,shls,nel,bf,der,be);
%             end
%             
%             [fi,JxX,F] = kine3d(Qxy,-ul,nel,0); %this is equivalent to ikine2d
%             JxX = 1/JxX; %this is equivalent to ikine2d
%     %         [F,JxX,fi] = ikine2d(Qxy,ul,nel); %this is equivalent to above two lines
%             Jdet = Jdet/JxX;
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