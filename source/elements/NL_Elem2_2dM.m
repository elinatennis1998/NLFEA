% function ElemK = NL_Elem2_2d(ElemXYZ, nel, ndf, R, S, T, pr, ps, PatchE, Patchv)
%
% Tim Truster
% 10/2011
% UIUC

% Master pure-displacement element
% body force problem re-run for Pinlei on 5/28/13, no other changes made

switch isw %Task Switch
%%    
    case 1
        
        if ndf > 2
            
            for i = 3:ndf
                lie(i,1) = 0;
            end
            
        end
        
%%
    case 3
        
%         ElemK = zeros(nst);
%         ElemF = zeros(nst,1);
        ElemK = zeros(ndf*nel);
        ElemF = zeros(ndf*nel,1);

        lam = getlam(mateprop);
        
        thick = 1;
        
        % Load Guass Integration Points

        if nel == 3 || nel == 6
            if exist('iprob','var') == 1 && iprob == 6
                lint = 91;
            else
                lint = 13;
            end
        elseif nel == 4
            if exist('iprob','var') == 1 && iprob == 6
                lint = 100;4;
            else
                lint = 4;
            end
        else
            if exist('iprob','var') == 1 && iprob == 6
                lint = 100;
            else
                lint = 9;
            end
        end
        der = 0;
        bf = 0;
        ib = 0;

        %Integration Loop
        for ll = 1:lint

            %Evaluate first derivatives of basis functions at int. point
            if nel == 3 || nel == 6
              [Wgt,litr,lits] =  intpntt(ll,lint,ib);
              [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
              [Qxy, shgs, Jdet] = shgt(xl+ul(1:2,:),nel,shld,shls,nel,bf,der,be);
            else
              [Wgt,litr,lits] =  intpntq(ll,lint,ib);
              [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
              [Qxy, shgs, Jdet] = shgq(xl+ul(1:2,:),nel,shld,shls,nel,bf,der,be);
            end
                
            if nel == 3
            % Form B matrix
            Bmat = [Qxy(1,1)  0        Qxy(2,1)  0        Qxy(3,1)  0       
                    0         Qxy(1,2) 0         Qxy(2,2) 0         Qxy(3,2)
                    Qxy(1,2)  Qxy(1,1) Qxy(2,2)  Qxy(2,1) Qxy(3,2)  Qxy(3,1)
                    Qxy(1,2) -Qxy(1,1) Qxy(2,2) -Qxy(2,1) Qxy(3,2) -Qxy(3,1)];
            elseif nel == 4
            % Form B matrix
            Bmat = [Qxy(1,1)  0        Qxy(2,1)  0        Qxy(3,1)  0        Qxy(4,1)  0       
                    0         Qxy(1,2) 0         Qxy(2,2) 0         Qxy(3,2) 0         Qxy(4,2) 
                    Qxy(1,2)  Qxy(1,1) Qxy(2,2)  Qxy(2,1) Qxy(3,2)  Qxy(3,1) Qxy(4,2)  Qxy(4,1)
                    Qxy(1,2) -Qxy(1,1) Qxy(2,2) -Qxy(2,1) Qxy(3,2) -Qxy(3,1) Qxy(4,2) -Qxy(4,1)];
            elseif nel == 6
            % Form B matrix
            Bmat = [Qxy(1,1)  0        Qxy(2,1)  0        Qxy(3,1)  0        Qxy(4,1)  0        Qxy(5,1)  0        Qxy(6,1)  0       
                    0         Qxy(1,2) 0         Qxy(2,2) 0         Qxy(3,2) 0         Qxy(4,2) 0         Qxy(5,2) 0         Qxy(6,2)
                    Qxy(1,2)  Qxy(1,1) Qxy(2,2)  Qxy(2,1) Qxy(3,2)  Qxy(3,1) Qxy(4,2)  Qxy(4,1) Qxy(5,2)  Qxy(5,1) Qxy(6,2)  Qxy(6,1)
                    Qxy(1,2) -Qxy(1,1) Qxy(2,2) -Qxy(2,1) Qxy(3,2) -Qxy(3,1) Qxy(4,2) -Qxy(4,1) Qxy(5,2) -Qxy(5,1) Qxy(6,2) -Qxy(6,1)];
            elseif nel == 9
            % Form B matrix
            Bmat = [Qxy(1,1)  0        Qxy(2,1)  0        Qxy(3,1)  0        Qxy(4,1)  0        Qxy(5,1)  0        Qxy(6,1)  0        Qxy(7,1)  0        Qxy(8,1)  0        Qxy(9,1)  0       
                    0         Qxy(1,2) 0         Qxy(2,2) 0         Qxy(3,2) 0         Qxy(4,2) 0         Qxy(5,2) 0         Qxy(6,2) 0         Qxy(7,2) 0         Qxy(8,2) 0         Qxy(9,2)
                    Qxy(1,2)  Qxy(1,1) Qxy(2,2)  Qxy(2,1) Qxy(3,2)  Qxy(3,1) Qxy(4,2)  Qxy(4,1) Qxy(5,2)  Qxy(5,1) Qxy(6,2)  Qxy(6,1) Qxy(7,2)  Qxy(7,1) Qxy(8,2)  Qxy(8,1) Qxy(9,2)  Qxy(9,1)
                    Qxy(1,2) -Qxy(1,1) Qxy(2,2) -Qxy(2,1) Qxy(3,2) -Qxy(3,1) Qxy(4,2) -Qxy(4,1) Qxy(5,2) -Qxy(5,1) Qxy(6,2) -Qxy(6,1) Qxy(7,2) -Qxy(7,1) Qxy(8,2) -Qxy(8,1) Qxy(9,2) -Qxy(9,1)];
            end
                
            [fi,JxX,F] = kine2d(Qxy,-ul,nel,0); %this is equivalent to ikine2d
            JxX = 1/JxX; %this is equivalent to ikine2d
    %         [F,JxX,fi] = ikine2d(Qxy,ul,nel); %this is equivalent to above two lines
            Jdet = Jdet/JxX;
                
            c1 = Wgt*Jdet*thick;

            [sigma3, cmat] = SigmaCmat2(F,JxX,mateprop,lam);
            
            Smat = [sigma3(1) 0  sigma3(3)/2 sigma3(3)/2
                    0 sigma3(2)  sigma3(3)/2 -sigma3(3)/2
                    sigma3(3)/2  sigma3(3)/2 (sigma3(2)+sigma3(1))/4 (sigma3(2)-sigma3(1))/4
                    sigma3(3)/2 -sigma3(3)/2 (sigma3(2)-sigma3(1))/4 (sigma3(2)+sigma3(1))/4]+cmat;
            
            ElemF = ElemF - c1*(Bmat'*(sigma3));
            
            ElemK = ElemK + c1*(Bmat'*Smat*Bmat);

        end %je
    ElemK;
%%
    case 6
        
%         ElemF = zeros(nst,1);
        ElemF = zeros(ndf*nel,1);

        lam = getlam(mateprop);
        
        thick = 1;
        
        % Load Guass Integration Points

        if nel == 3 || nel == 6
            if exist('iprob','var') == 1 && iprob == 6
                lint = 91;
            else
                lint = 13;
            end
        elseif nel == 4
            if exist('iprob','var') == 1 && iprob == 6
                lint = 100;4;
            else
                lint = 4;
            end
        else
            if exist('iprob','var') == 1 && iprob == 6
                lint = 100;
            else
                lint = 9;
            end
        end
        der = 0;
        bf = 0;
        ib = 0;

        %Integration Loop
        for ll = 1:lint

            %Evaluate first derivatives of basis functions at int. point
            if nel == 3 || nel == 6
              [Wgt,litr,lits] =  intpntt(ll,lint,ib);
              [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
              [Qxy, shgs, Jdet] = shgt(xl+ul(1:2,:),nel,shld,shls,nel,bf,der,be);
            else
              [Wgt,litr,lits] =  intpntq(ll,lint,ib);
              [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
              [Qxy, shgs, Jdet] = shgq(xl+ul(1:2,:),nel,shld,shls,nel,bf,der,be);
            end
                
            if nel == 3            
            % Form B matrix    
            Bmat = [Qxy(1,1)  0        Qxy(2,1)  0        Qxy(3,1)  0       
                    0         Qxy(1,2) 0         Qxy(2,2) 0         Qxy(3,2)
                    Qxy(1,2)  Qxy(1,1) Qxy(2,2)  Qxy(2,1) Qxy(3,2)  Qxy(3,1)
                    Qxy(1,2) -Qxy(1,1) Qxy(2,2) -Qxy(2,1) Qxy(3,2) -Qxy(3,1)];
            elseif nel == 4            
            % Form B matrix
            Bmat = [Qxy(1,1)  0        Qxy(2,1)  0        Qxy(3,1)  0        Qxy(4,1)  0       
                    0         Qxy(1,2) 0         Qxy(2,2) 0         Qxy(3,2) 0         Qxy(4,2) 
                    Qxy(1,2)  Qxy(1,1) Qxy(2,2)  Qxy(2,1) Qxy(3,2)  Qxy(3,1) Qxy(4,2)  Qxy(4,1)
                    Qxy(1,2) -Qxy(1,1) Qxy(2,2) -Qxy(2,1) Qxy(3,2) -Qxy(3,1) Qxy(4,2) -Qxy(4,1)];
            elseif nel == 6
            % Form B matrix
            Bmat = [Qxy(1,1)  0        Qxy(2,1)  0        Qxy(3,1)  0        Qxy(4,1)  0        Qxy(5,1)  0        Qxy(6,1)  0       
                    0         Qxy(1,2) 0         Qxy(2,2) 0         Qxy(3,2) 0         Qxy(4,2) 0         Qxy(5,2) 0         Qxy(6,2)
                    Qxy(1,2)  Qxy(1,1) Qxy(2,2)  Qxy(2,1) Qxy(3,2)  Qxy(3,1) Qxy(4,2)  Qxy(4,1) Qxy(5,2)  Qxy(5,1) Qxy(6,2)  Qxy(6,1)
                    Qxy(1,2) -Qxy(1,1) Qxy(2,2) -Qxy(2,1) Qxy(3,2) -Qxy(3,1) Qxy(4,2) -Qxy(4,1) Qxy(5,2) -Qxy(5,1) Qxy(6,2) -Qxy(6,1)];
            elseif nel == 9
            % Form B matrix
            Bmat = [Qxy(1,1)  0        Qxy(2,1)  0        Qxy(3,1)  0        Qxy(4,1)  0        Qxy(5,1)  0        Qxy(6,1)  0        Qxy(7,1)  0        Qxy(8,1)  0        Qxy(9,1)  0       
                    0         Qxy(1,2) 0         Qxy(2,2) 0         Qxy(3,2) 0         Qxy(4,2) 0         Qxy(5,2) 0         Qxy(6,2) 0         Qxy(7,2) 0         Qxy(8,2) 0         Qxy(9,2)
                    Qxy(1,2)  Qxy(1,1) Qxy(2,2)  Qxy(2,1) Qxy(3,2)  Qxy(3,1) Qxy(4,2)  Qxy(4,1) Qxy(5,2)  Qxy(5,1) Qxy(6,2)  Qxy(6,1) Qxy(7,2)  Qxy(7,1) Qxy(8,2)  Qxy(8,1) Qxy(9,2)  Qxy(9,1)
                    Qxy(1,2) -Qxy(1,1) Qxy(2,2) -Qxy(2,1) Qxy(3,2) -Qxy(3,1) Qxy(4,2) -Qxy(4,1) Qxy(5,2) -Qxy(5,1) Qxy(6,2) -Qxy(6,1) Qxy(7,2) -Qxy(7,1) Qxy(8,2) -Qxy(8,1) Qxy(9,2) -Qxy(9,1)];
            end
                
            [fi,JxX,F] = kine2d(Qxy,-ul,nel,0); %this is equivalent to ikine2d
            JxX = 1/JxX; %this is equivalent to ikine2d
    %         [F,JxX,fi] = ikine2d(Qxy,ul,nel); %this is equivalent to above two lines
            Jdet = Jdet/JxX;
                
            c1 = Wgt*Jdet*thick;

            sigma3 = SigmaCmat2(F,JxX,mateprop,lam);
            
            Smat = [sigma3(1) 0  sigma3(3)/2 sigma3(3)/2
                    0 sigma3(2)  sigma3(3)/2 -sigma3(3)/2
                    sigma3(3)/2  sigma3(3)/2 (sigma3(2)+sigma3(1))/4 (sigma3(2)-sigma3(1))/4
                    sigma3(3)/2 -sigma3(3)/2 (sigma3(2)-sigma3(1))/4 (sigma3(2)+sigma3(1))/4]+cmat;
            
            ElemF = ElemF - c1*(Bmat'*(sigma3));

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
        
        if exist('iprob','var') == 1 && iprob == 6
            lint = 10;
        else
            lint = 4; % Use 10 for body force BF2U4M0.m
        end
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

%         lamda = Elemv*ElemE/((1+Elemv)*(1-2*Elemv));
%         mu = ElemE/(2*(1+Elemv));
%         thick = 1;
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
                    mu = 40;
                    lam = 40;
                    X = xl(1,1:nel)*shl;
                    Y = xl(2,1:nel)*shl;
                    P = [ (101*X*lam*(101*X + 100))/10000, (101*Y*mu)/(101*X + 100) - (10201*X*Y*lam)/10000, 0
                         (101*Y*mu)/100, mu - (100*mu)/(101*X + 100) + X*((101*lam)/100 + (101*mu)/100),     0
                              0,                                            0, (101*X*lam*(101*X + 100))/10000];
                    Traction = P*tu3';
                else
                    Traction = traction;
                end
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

                ElemF(ndf*o-1) = ElemF(ndf*o-1) + F(1)*c1;

                ElemF(ndf*o-0)   = ElemF(ndf*o-0)   + F(2)*c1;

            end %o

        end %ie
        
        ElemF;
        
        % Reorder nodes back to the orientation of the element
        SurfOrientStoE
        ElemF(1:ndf*nel) = ElemF(ilist2);
        
        
%%
    case 15 % Body Force
        
        ElemF = zeros(nst,1);

        lam = getlam(mateprop);
        
        thick = 1;
        
        % Load Guass Integration Points

        if nel == 3 || nel == 6
            if exist('iprob','var') == 1 && iprob == 6
                lint = 91;
            else
                lint = 13;
            end
        elseif nel == 4
            if exist('iprob','var') == 1 && iprob == 6
                lint = 100;4;
            else
                lint = 4;
            end
        else
            if exist('iprob','var') == 1 && iprob == 6
                lint = 100;
            else
                lint = 9;
            end
        end
        der = 0;
        bf = 0;
        ib = 0;

        %Integration Loop
        for ll = 1:lint

            %Evaluate first derivatives of basis functions at int. point
            if nel == 3 || nel == 6
              [Wgt,litr,lits] =  intpntt(ll,lint,ib);
              [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
              [Qxy, shgs, Jdet] = shgt(xl,nel,shld,shls,nel,bf,der,be);
            else
              [Wgt,litr,lits] =  intpntq(ll,lint,ib);
              [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
              [Qxy, shgs, Jdet] = shgq(xl,nel,shld,shls,nel,bf,der,be);
            end
                
            if nel == 3
%             % Form N matrix
            Nmat = [shl(1)  0        shl(2)  0        shl(3)  0       
                    0         shl(1) 0         shl(2) 0         shl(3)];
            elseif nel == 4
%             % Form N matrix
            Nmat = [shl(1)  0        shl(2)  0        shl(3)  0        shl(4)  0       
                    0         shl(1) 0         shl(2) 0         shl(3) 0         shl(4)];
            elseif nel == 6
%             % Form N matrix
            Nmat = [shl(1)  0        shl(2)  0        shl(3)  0        shl(4)  0        shl(5)  0        shl(6)  0       
                    0         shl(1) 0         shl(2) 0         shl(3) 0         shl(4) 0         shl(5) 0         shl(6)];
            elseif nel == 9
%             % Form N matrix
            Nmat = [shl(1)  0        shl(2)  0        shl(3)  0        shl(4)  0        shl(5)  0        shl(6)  0        shl(7)  0        shl(8)  0        shl(9)  0       
                    0         shl(1) 0         shl(2) 0         shl(3) 0         shl(4) 0         shl(5) 0         shl(6) 0         shl(7) 0         shl(8) 0         shl(9)];
            end
                
            c1 = Wgt*Jdet*thick;
            
            if exist('iprob','var')
                if iprob == 6
                    mu = 40;
                    X = xl(1,1:nel)*shl;
                    Y = xl(2,1:nel)*shl;
                    fb = [ - (101*mu)/(101*X + 100) - (101*lam*(101*X + 100))/10000
                                                      0];
                else
                    fb = bodyf(1:2)';
                end
            else
                fb = bodyf(1:2)';
            end
            
            ElemF = ElemF + c1*Nmat'*fb;

        end %je
    ElemF;
        
%%
    case 7 % user boundary tractions
        
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
        
        if exist('iprob','var') == 1 && iprob == 6
            lint = 10;
        else
            lint = 4; % Use 10 for body force BF2U4M0.m
        end
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

%         lamda = Elemv*ElemE/((1+Elemv)*(1-2*Elemv));
%         mu = ElemE/(2*(1+Elemv));
%         thick = 1;
        ideriv = 0;
        der = 0;
        bf = 0;
        thick = 1;
        
        if load == 1
        if exist('iprob','var') && iprob ==8
            if lamda > 0
                psi = BendAngle*lamda;
                syms aaa bbb positive
                eq1 = (bbb^2-aaa^2)/(H)*psi/2 - L; % psi only subtends half of beam, so (Ro^2-Ri^2)*psi/2 = L*H
                eq2 = psi^2*aaa*bbb - L^2;
                solv = solve(eq1, eq2);
                Ro = double(solv.bbb);
                Ri = double(solv.aaa);
            else
                Ro = 0;
                Ri = 0;
            end
        end
        end
        
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
            if lamda > 0 && tu3(1) ~= 0
               lamda; 
            end
            
            if exist('iprob','var')
                if iprob == 6
                    mu = 40;
                    lam = 40;
                    X = xl(1,1:nel)*shl;
                    Y = xl(2,1:nel)*shl;
                    Traction = PureBendBF(X,Y,mu,lam,BFdelta*lamda,tu3(1:2)');
                elseif iprob == 8
                    PatchE = mateprop(4);
                    Patchv = mateprop(5);
                    mu = PatchE/(2*(1+Patchv));%80.19;
                    X = xl(1,1:nel)*shl;
                    Y = xl(2,1:nel)*shl;
                    Traction = PureBendTB(X,Y,mu,BendAngle*lamda,Ro,Ri,8,1,tu3(1:2)');
                elseif iprob == 9
                    PatchE = mateprop(4);
                    Patchv = mateprop(5);
                    X = xl(1,1:nel)*shl;
                    Y = xl(2,1:nel)*shl;
                    Traction = PureBend3TB(X,Y,PatchE,Patchv,BendAngle*lamda,8,tu3(1:2)');
                else
                    Traction = traction;
                end
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

                ElemF(ndf*o-1) = ElemF(ndf*o-1) + F(1)*c1;

                ElemF(ndf*o-0)   = ElemF(ndf*o-0)   + F(2)*c1;

            end %o

        end %ie
        
        ElemF;
        
        % Reorder nodes back to the orientation of the element
        SurfOrientStoE
        ElemF(1:ndf*nel) = ElemF(ilist2);
        
%%
    case -15 % User Body Force
        
        ElemF = zeros(nst,1);

        lam = getlam(mateprop);
        
        thick = 1;
        
        % Load Guass Integration Points

        if nel == 3 || nel == 6
            if exist('iprob','var') == 1 && iprob == 8
                lint = 91;
            else
                lint = 13;
            end
        elseif nel == 4
            if exist('iprob','var') == 1 && iprob == 8
                lint = 100;4;
            else
                lint = 4;
            end
        else
            if exist('iprob','var') == 1 && iprob == 8
                lint = 100;
            else
                lint = 9;
            end
        end
        der = 0;
        bf = 0;
        ib = 0;

        %Integration Loop
        for ll = 1:lint

            %Evaluate first derivatives of basis functions at int. point
            if nel == 3 || nel == 6
              [Wgt,litr,lits] =  intpntt(ll,lint,ib);
              [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
              [Qxy, shgs, Jdet] = shgt(xl,nel,shld,shls,nel,bf,der,be);
            else
              [Wgt,litr,lits] =  intpntq(ll,lint,ib);
              [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
              [Qxy, shgs, Jdet] = shgq(xl,nel,shld,shls,nel,bf,der,be);
            end
                
            if nel == 3
%             % Form N matrix
            Nmat = [shl(1)  0        shl(2)  0        shl(3)  0       
                    0         shl(1) 0         shl(2) 0         shl(3)];
            elseif nel == 4
%             % Form N matrix
            Nmat = [shl(1)  0        shl(2)  0        shl(3)  0        shl(4)  0       
                    0         shl(1) 0         shl(2) 0         shl(3) 0         shl(4)];
            elseif nel == 6
%             % Form N matrix
            Nmat = [shl(1)  0        shl(2)  0        shl(3)  0        shl(4)  0        shl(5)  0        shl(6)  0       
                    0         shl(1) 0         shl(2) 0         shl(3) 0         shl(4) 0         shl(5) 0         shl(6)];
            elseif nel == 9
%             % Form N matrix
            Nmat = [shl(1)  0        shl(2)  0        shl(3)  0        shl(4)  0        shl(5)  0        shl(6)  0        shl(7)  0        shl(8)  0        shl(9)  0       
                    0         shl(1) 0         shl(2) 0         shl(3) 0         shl(4) 0         shl(5) 0         shl(6) 0         shl(7) 0         shl(8) 0         shl(9)];
            end
                
            c1 = Wgt*Jdet*thick;
            
            if exist('iprob','var')
                if iprob == 6
                    mu = 40;
                    lam = 40;
                    X = xl(1,1:nel)*shl;
                    Y = xl(2,1:nel)*shl;
                    [Traction,fb] = PureBendBF(X,Y,mu,lam,BFdelta*lamda,[1;0]);
                elseif iprob == 8
                    PatchE = mateprop(4);
                    Patchv = mateprop(5);
                    mu = PatchE/(2*(1+Patchv));%80.19;
                    X = xl(1,1:nel)*shl;
                    Y = xl(2,1:nel)*shl;
                    [Traction,fb] = PureBendTB(X,Y,mu,BendAngle*lamda,Ro,Ri,8,1,[1;0]);
                elseif iprob == 9
                    PatchE = mateprop(4);
                    Patchv = mateprop(5);
                    X = xl(1,1:nel)*shl;
                    Y = xl(2,1:nel)*shl;
                    [Traction,fb] = PureBend3TB(X,Y,PatchE,Patchv,BendAngle*lamda,8,[1;0]);
                else
                    fb = bodyf(1:2)';
                end
            else
                fb = bodyf(1:2)';
            end
            
            ElemF = ElemF + c1*Nmat'*fb;

        end %je
    ElemF;
    
%%        
    case 11 % Error estimation
        
        ElemE = zeros(numEn,1);
        %Set integration number
        if nel == 3 || nel == 6
            if exist('iprob','var') == 1 && iprob == 6
                lint = 91;
            else
                lint = 13;
            end
        elseif nel == 4
            if exist('iprob','var') == 1 && iprob == 6
                lint = 100;4;
            else
                lint = 4;
            end
        else
            if exist('iprob','var') == 1 && iprob == 6
                lint = 100;
            else
                lint = 9;
            end
        end
        ib = 0;
        bf = 1;
        der = 1;

        lam = getlam(mateprop);
        
        thick = 1;
        el2el = zeros(3,1);
        eprixel = zeros(3,1);
        epriyel = zeros(3,1);
        
        for ll = 1:lint
                    
            %Evaluate first derivatives of basis functions at int. point
            if nel == 3 || nel == 6
                [Wgt,litr,lits] =  intpntt(ll,lint,ib);
              [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
              [Qxy, shgs, Jdet, be, sx] = shgt(xl+ul(1:2,:),nel,shld,shls,nen,bf,der,be);
            else
                [Wgt,litr,lits] =  intpntq(ll,lint,ib);
              [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
              [Qxy, shgs, Jdet, be, sx] = shgq(xl+ul(1:2,:),nel,shld,shls,nen,bf,der,be);
            end

            [fi,JxX,F,QXY] = kine2d(Qxy,-ul,nel,1); %this is equivalent to ikine2d
            JxX = 1/JxX; %this is equivalent to ikine2d
    %         [F,JxX,fi] = ikine2d(Qxy,ul,nel); %this is equivalent to above two lines
    %         xs = sx;
            Jdet = Jdet/JxX;

            c1 = Wgt*Jdet*thick;

            xint = xl(1,1:nel)*shl;
            yint = xl(2,1:nel)*shl;
            dux = ul(:,1:nel)*QXY(1:nel,1);
            duy = ul(:,1:nel)*QXY(:,2);
            u = ul(:,1:nel)*shl;

%             Compute value of exact fields at int. point
            if iprob == 6
                [ue,duex,duey] = uexact_BF(xint,yint,lam);
            elseif iprob == 8
                [ue,due] = PureBendU(xint,yint,BendAngle*lamda,Ro,Ri,8,1);
                duex = due(:,1);
                duey = due(:,2);
            elseif iprob == 9
                [ue,due] = PureBend3U(xint,yint,BendAngle*lamda,8);
                duex = due(:,1);
                duey = due(:,2);
            else
                ue = zeros(3,1);
                duex = zeros(3,1);
                duey = zeros(3,1);
            end

            %Add standard int. point error to element standard error
            for in = 1:ndf
                un   = c1 * ( (u(in)-ue(in))^2 );
                upnx   = c1 * ( (dux(in)-duex(in))^2 );
                upny   = c1 * ( (duy(in)-duey(in))^2 );
                el2el(in)   = el2el(in)   + un;
                eprixel(in) = eprixel(in) + upnx;
                epriyel(in) = epriyel(in) + upny;
            end

        end %je

        for in= 1:ndf
            ElemE(in) = el2el(in);
            ElemE(in+2) = eprixel(in);
            ElemE(in+4) = epriyel(in);
        end
        
%         H1up = el2fine(3)+el2fine(4)+el2fine(5)+el2fine(6);
        H1u = eprixel(1)+eprixel(2)+epriyel(1)+epriyel(2);
%         Ieffvals(elem,:) = [sqrt(H1up/H1u) H1up H1u];
        Ieffvals(elem,1:3) = [el2el(1) el2el(2) H1u];% JxXm divu];
        
%%        
    case 25 %Stress Projection2

        ElemS = zeros(nel,npstr+1);
        ElemS2 = zeros(nel,npstr);

        lam = getlam(mateprop);
        
        I1 = [1; 1; 1; 0; 0; 0; 0; 0; 0];
        stresID = [1 2 4 0 1 3];
        spvec0 = I1;
        
        thick = 1;
        
        % Load Guass Integration Points

        if nel == 3
            lint = 1;
            nint = 1;
        elseif nel == 4
%             lint = 4;
            lint = 4;
            nint = 1;
        elseif nel == 6
            lint = 7;
            nint = 3;
        else
            lint = 9;
            nint = 4;
        end
        
        der = 0;
        bf = 1;
        ib = 0;

        %Stress Loop
        for ll = 1:nint

            %Evaluate first derivatives of basis functions at int. point
            if nel == 3 || nel == 6
              [Wgt,litr,lits] =  intpntt(ll,nint,ib);
              [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
              [Qxy, shgs, Jdet,be] = shgt(xl+ul(1:2,:),nel,shld,shls,nel,bf,der,be);
            else
              [Wgt,litr,lits] =  intpntq(ll,nint,ib);
              [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
              [Qxy, shgs, Jdet,be] = shgq(xl+ul(1:2,:),nel,shld,shls,nel,bf,der,be);
            end
            
            [fi,JxX,F] = kine2d(Qxy,-ul,nel,0); %this is equivalent to ikine2d
            JxX = 1/JxX; %this is equivalent to ikine2d
    %         [F,JxX,fi] = ikine2d(Qxy,ul,nel); %this is equivalent to above two lines
            Jdet = Jdet/JxX;
            
            sigma = SigmaCmat3([F zeros(2,1); zeros(1,2) 1],JxX,mateprop,lam);
            sigma = sigma/JxX;
            
            for stres = 1:npstr
            
            if stres <= 3 % stress components
                sigmas = sigma(stresID(stres));
            elseif stres >= 5
                if stres <= 6 % principal stresses
                sigma2 = [sigma(1) sigma(4) sigma(6); sigma(4) sigma(2) sigma(5); sigma(6) sigma(5) sigma(3)];
                psig = eig(sigma2);
                sigmas = psig(stresID(stres));
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
        if nel == 3
            plist = [0 1 0
                     0 0 1];
        elseif nel == 4
            plist = [-1 1 1 -1
                     -1 -1 1 1];
        elseif nel == 6
            plist = [-1/3 5/3 -1/3 2/3 2/3 -1/3
                     -1/3 -1/3 5/3 -1/3 2/3 2/3];
        else
            sqr3 = sqrt(3);
            plist = [-sqr3 sqr3 sqr3 -sqr3 0 sqr3 0 -sqr3 0
                     -sqr3 -sqr3 sqr3 sqr3 -sqr3 0 sqr3 0 0];
        end
        
        for ll = 1:nelS
            
            r = plist(1,ll);
            s = plist(2,ll);
            shpS = sshp2d(r,s,nint);
            
            for stres = 1:npstr
                
                sigmas = ElemS2(1:nint,stres)'*shpS;
                ElemS(ll,stres) = sigmas;
                
            end
            
        end
        
        %Integration Loop
        Vol = 0;
        for ll = 1:lint

            %Evaluate first derivatives of basis functions at int. point
            if nel == 3 || nel == 6
              [Wgt,litr,lits] =  intpntt(ll,lint,ib);
              [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
              [Qxy, shgs, Jdet,be] = shgt(xl+ul(1:2,:),nel,shld,shls,nel,bf,der,be);
            else
              [Wgt,litr,lits] =  intpntq(ll,lint,ib);
              [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
              [Qxy, shgs, Jdet,be] = shgq(xl+ul(1:2,:),nel,shld,shls,nel,bf,der,be);
            end
            
            [fi,JxX,F] = kine2d(Qxy,-ul,nel,0); %this is equivalent to ikine2d
            JxX = 1/JxX; %this is equivalent to ikine2d
    %         [F,JxX,fi] = ikine2d(Qxy,ul,nel); %this is equivalent to above two lines
            Jdet = Jdet/JxX;

            w = Wgt*Jdet*thick;
            
            Vol = Vol + w;

        end %je

        for i = 1:nel
        ElemS(i,npstr+1) = Vol; % use this for weighted average 1; % use this for simple average 
        end
%%        
    case 26 % Element Stress

        ElemS = zeros(1,nestr);

        lam = getlam(mateprop);
        
        I1 = [1; 1; 1; 0; 0; 0; 0; 0; 0];
        stresID = [1 2 4 0 1 3];
        spvec0 = I1;
        
        % Load Guass Integration Points

        
        der = 0;
        bf = 0;
        ib = 0;


            %Evaluate first derivatives of basis functions at int. point
            if nel == 3 || nel == 6
              ss = 1/3*[1 1];
              [shl,shld,shls,be] = shlt(ss(1),ss(2),nel,nel,der,bf);
              [Qxy,shgs,Jdet,be] = shgt(xl+ul(1:2,:),nel,shld,shls,nel,bf,der,be);
            else
              ss = [0 0];
              [shl,shld,shls,be] = shlq(ss(1),ss(2),nel,nel,der,bf);
              [Qxy,shgs,Jdet,be] = shgq(xl+ul(1:2,:),nel,shld,shls,nel,bf,der,be);
            end
            
            if nelP == 3 || nelP == 6
              shlp = shlt(ss(1),ss(2),nelP,nel,0,0);
            else
              shlp = shlq(ss(1),ss(2),nelP,nel,0,0);
            end
            
            [fi,JxX,F] = kine2d(Qxy,-ul,nel,0); %this is equivalent to ikine2d
            JxX = 1/JxX; %this is equivalent to ikine2d
    %         [F,JxX,fi] = ikine2d(Qxy,ul,nel); %this is equivalent to above two lines
            Jdet = Jdet/JxX;
            
            sigma = SigmaCmat3([F zeros(2,1); zeros(1,2) 1],JxX,mateprop,lam);
            sigma = sigma/JxX;
            
            for stres = 1:nestr
            
            if stres <= 3 % stress components
                sigmas = sigma(stresID(stres));
            elseif stres >= 5
                if stres <= 6 % principal stresses
                sigma2 = [sigma(1) sigma(4) sigma(6); sigma(4) sigma(2) sigma(5); sigma(6) sigma(5) sigma(3)];
                psig = eig(sigma2);
                sigmas = psig(stresID(stres));
                else % hydrostatic stress
                sigmas = 1/3*sigma'*I1;
                end
            else % von Mises stress
                trs = sigma'*I1;
                dsig = sigma - 1/3*trs*I1;
                sigmas = sqrt(3/2*(dsig'*dsig));
            end
            
            ElemS(1,stres) = sigmas;
            
            end
end