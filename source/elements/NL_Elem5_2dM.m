% function ElemK = NL_Elem2_2d(ElemXYZ, nel, ndf, R, S, T, pr, ps, PatchE, Patchv)
%
% Tim Truster
% 10/27/2011
% UIUC
%

% Master non-symmetric mixed formulation element, without stabilization

switch isw %Task Switch
%%
    case 3
        
        ElemK = zeros(nst);
        ElemF = zeros(nst,1);
%         ElemP = zeros(2,1);
%         ElemR = zeros(2,1);
%         Nmat = zeros(3,3*nel);
        Bmat = zeros(4,3*nel);
%         BBmat = zeros(6,3*nel);

        lam = getlam(mateprop);
        
        I1 = [1; 1; 0; 0];
        spvec0 = I1;
        cpmat1 = [1 1 0 0
                  1 1 0 0
                  0 0 0 0
                  0 0 0 0];
        cpmat2 = [-2 0 0 0
                  0 -2 0 0
                  0 0 -1 0
                  0 0 0  0];
        cpmat0 = cpmat1 + cpmat2;
        
        thick = 1;
        
        % Load Guass Integration Points

        if nel == 3
            if exist('iprob','var') == 1 && iprob == 6
                lint = 25;
            else
                lint = 7; %minimum of 7 for all integrals in deformed state
            end
        elseif nel == 6
            if exist('iprob','var') == 1 && iprob == 6
                lint = 25;
            else
                lint = 13; %minimum of 13 for all integrals in deformed state
            end
        elseif nel == 4
            if exist('iprob','var') == 1 && iprob == 6
                lint = 100;4;
            else
                lint = 16;
            end
        else
            if exist('iprob','var') == 1 && iprob == 6
                lint = 100;
            else
                lint = 16;
            end
        end
        bf = 0;
        ib = 0;
        
        %Integration Loop
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

            if nelP == 3 || nelP == 6
              [shlp,shld,shls] = shlt(litr,lits,nelP,nel,0,0);
              [Pxy, shgp] = shgt(xl+ul(1:2,:),nelP,shld,shls,nelP,0,0,be);
            else
              [shlp,shld,shls] = shlq(litr,lits,nelP,nel,0,0);
              [Pxy, shgp] = shgq(xl+ul(1:2,:),nelP,shld,shls,nelP,0,0,be);
            end

            if Jdet < 0
                Jdet
            end
            
            if nel == 3
%             % Form N matrix
            Nmat = [shl(1)  0        0 shl(2)  0        0 shl(3)  0  0     
                    0         shl(1) 0  0         shl(2) 0 0         shl(3) 0
                    0        0         shlp(1,1) 0        0         shlp(2,1) 0        0         shlp(3,1)];
            % Form B matrix
            Bmat = [Qxy(1,1)  0        0         Qxy(2,1)  0        0         Qxy(3,1)  0        0
                    0         Qxy(1,2) 0         0         Qxy(2,2) 0         0         Qxy(3,2) 0         
                    Qxy(1,2)  Qxy(1,1) 0         Qxy(2,2)  Qxy(2,1) 0         Qxy(3,2)  Qxy(3,1) 0         
                    Qxy(1,2) -Qxy(1,1) 0         Qxy(2,2) -Qxy(2,1) 0         Qxy(3,2) -Qxy(3,1) 0         
                    0        0         shlp(1,1) 0        0         shlp(2,1) 0        0         shlp(3,1)];
            elseif nel == 4
%             % Form N matrix
            Nmat = [shl(1)  0        0 shl(2)  0        0 shl(3)  0        0 shl(4)  0        0
                    0         shl(1) 0 0         shl(2) 0 0         shl(3) 0 0         shl(4) 0
                    0        0         shlp(1,1) 0        0         shlp(2,1) 0        0         shlp(3,1) 0        0         shlp(4,1)]; % Add back for BF2U4M0.m
            % Form B matrix
            Bmat = [Qxy(1,1)  0        0         Qxy(2,1)  0        0         Qxy(3,1)  0        0         Qxy(4,1)  0        0
                    0         Qxy(1,2) 0         0         Qxy(2,2) 0         0         Qxy(3,2) 0         0         Qxy(4,2) 0         
                    Qxy(1,2)  Qxy(1,1) 0         Qxy(2,2)  Qxy(2,1) 0         Qxy(3,2)  Qxy(3,1) 0         Qxy(4,2)  Qxy(4,1) 0         
                    Qxy(1,2) -Qxy(1,1) 0         Qxy(2,2) -Qxy(2,1) 0         Qxy(3,2) -Qxy(3,1) 0         Qxy(4,2) -Qxy(4,1) 0         
                    0        0         shlp(1,1) 0        0         shlp(2,1) 0        0         shlp(3,1) 0        0         shlp(4,1)];
            elseif nel == 6
%             % Form N matrix
            Nmat = [shl(1)  0        0 shl(2)  0        0 shl(3)  0        0 shl(4)  0        0 shl(5)  0        0 shl(6)  0        0
                    0         shl(1) 0 0         shl(2) 0 0         shl(3) 0 0         shl(4) 0 0         shl(5) 0 0         shl(6) 0
                    0        0         shlp(1,1) 0        0         shlp(2,1) 0        0         shlp(3,1) 0        0         shlp(4,1) 0        0         shlp(5,1) 0        0         shlp(6,1)];
            % Form B matrix
            Bmat = [Qxy(1,1)  0        0         Qxy(2,1)  0        0         Qxy(3,1)  0        0         Qxy(4,1)  0        0         Qxy(5,1)  0        0         Qxy(6,1)  0        0
                    0         Qxy(1,2) 0         0         Qxy(2,2) 0         0         Qxy(3,2) 0         0         Qxy(4,2) 0         0         Qxy(5,2) 0         0         Qxy(6,2) 0
                    Qxy(1,2)  Qxy(1,1) 0         Qxy(2,2)  Qxy(2,1) 0         Qxy(3,2)  Qxy(3,1) 0         Qxy(4,2)  Qxy(4,1) 0         Qxy(5,2)  Qxy(5,1) 0         Qxy(6,2)  Qxy(6,1) 0
                    Qxy(1,2) -Qxy(1,1) 0         Qxy(2,2) -Qxy(2,1) 0         Qxy(3,2) -Qxy(3,1) 0         Qxy(4,2) -Qxy(4,1) 0         Qxy(5,2) -Qxy(5,1) 0         Qxy(6,2) -Qxy(6,1) 0
                    0        0         shlp(1,1) 0        0         shlp(2,1) 0        0         shlp(3,1) 0        0         shlp(4,1) 0        0         shlp(5,1) 0        0         shlp(6,1)];
            elseif nel == 9
%             % Form N matrix
            Nmat = [shl(1)  0      0 shl(2)  0      0 shl(3)  0      0 shl(4)  0      0 shl(5)  0      0 shl(6)  0      0 shl(7)  0      0 shl(8)  0      0 shl(9)  0      0     
                    0       shl(1) 0 0       shl(2) 0 0       shl(3) 0 0       shl(4) 0 0       shl(5) 0 0       shl(6) 0 0       shl(7) 0 0       shl(8) 0 0       shl(9) 0
                    0       0      shlp(1,1) 0    0 shlp(2,1) 0    0 shlp(3,1) 0    0 shlp(4,1) 0    0 shlp(5,1) 0    0 shlp(6,1) 0    0 shlp(7,1) 0    0 shlp(8,1) 0    0 shlp(9,1)];
            % Form B matrix
            Bmat = [Qxy(1,1)  0        0         Qxy(2,1)  0        0         Qxy(3,1)  0        0         Qxy(4,1)  0        0         Qxy(5,1)  0        0         Qxy(6,1)  0        0         Qxy(7,1)  0        0         Qxy(8,1)  0        0         Qxy(9,1)  0        0
                    0         Qxy(1,2) 0         0         Qxy(2,2) 0         0         Qxy(3,2) 0         0         Qxy(4,2) 0         0         Qxy(5,2) 0         0         Qxy(6,2) 0         0         Qxy(7,2) 0         0         Qxy(8,2) 0         0         Qxy(9,2) 0
                    Qxy(1,2)  Qxy(1,1) 0         Qxy(2,2)  Qxy(2,1) 0         Qxy(3,2)  Qxy(3,1) 0         Qxy(4,2)  Qxy(4,1) 0         Qxy(5,2)  Qxy(5,1) 0         Qxy(6,2)  Qxy(6,1) 0         Qxy(7,2)  Qxy(7,1) 0         Qxy(8,2)  Qxy(8,1) 0         Qxy(9,2)  Qxy(9,1) 0
                    Qxy(1,2) -Qxy(1,1) 0         Qxy(2,2) -Qxy(2,1) 0         Qxy(3,2) -Qxy(3,1) 0         Qxy(4,2) -Qxy(4,1) 0         Qxy(5,2) -Qxy(5,1) 0         Qxy(6,2) -Qxy(6,1) 0         Qxy(7,2) -Qxy(7,1) 0         Qxy(8,2) -Qxy(8,1) 0         Qxy(9,2) -Qxy(9,1) 0
                    0        0         shlp(1,1) 0        0         shlp(2,1) 0        0         shlp(3,1) 0        0         shlp(4,1) 0        0         shlp(5,1) 0        0         shlp(6,1) 0        0         shlp(7,1) 0        0         shlp(8,1) 0        0         shlp(9,1)];
            end

            [fi,JxX,F] = kine2d(Qxy,-ul,nel,0); %this is equivalent to ikine2d
            JxX = 1/JxX; %this is equivalent to ikine2d
    %         [F,JxX,fi] = ikine2d(Qxy,ul,nel); %this is equivalent to above two lines
            Jdet = Jdet/JxX;
            [theta1,theta2,theta3] = ThetaNS(JxX,mateprop);
            spvec = JxX*spvec0;
            cpmat = JxX*cpmat0;
                
            c1 = Wgt*Jdet*thick;

            [sigmai, cmati] = SigmaCmatNSCST2i(F,JxX,mateprop);
            press = ul(3,:)*shlp;
            sigmap = press*spvec;
            cmatp = press*cpmat;
            sigma = sigmai + sigmap;
            cmat = cmati + cmatp;
            
%             if exist('iprob','var')
%                 if iprob == 6
%                     mu = 40;
%                     X = xl(1,:)*shl;
%                     Y = xl(2,:)*shl;
%                     fb = [ - (101*mu)/(101*X + 100) - (101*lam*(101*X + 100))/10000
%                                                       0; 0];
%                 else
%                     fb = zeros(3,1);
%                 end
%             else
%                 fb = zeros(3,1);
%             end
            
            Svec = [sigma; (theta1-press/lam)];
            sigma3 = sigma;
            Smat = [sigma3(1) 0  sigma3(3)/2 sigma3(3)/2
                    0 sigma3(2)  sigma3(3)/2 -sigma3(3)/2
                    sigma3(3)/2  sigma3(3)/2 (sigma3(2)+sigma3(1))/4 (sigma3(2)-sigma3(1))/4
                    sigma3(3)/2 -sigma3(3)/2 (sigma3(2)-sigma3(1))/4 (sigma3(2)+sigma3(1))/4];
            Cmat = cmat;
            Cvec12 = JxX*I1;
            Cvec21 = theta2*JxX*I1';
            
            D11 = [Smat+Cmat Cvec12
                   Cvec21 -1/(lam)];
            
            ElemF = ElemF - c1*(Bmat'*(Svec));% - Nmat'*fb  Add back for BF2U4M0.m
            
            ElemK = ElemK + c1*(Bmat'*D11*Bmat);

        end %je
ElemK;
%         if elem == 48
%             ElemK(1,1)
%         end
%%
    case 6
        
        ElemF = zeros(nst,1);
%         ElemP = zeros(2,1);
%         ElemR = zeros(2,1);
%         Nmat = zeros(3,3*nel);
        Bmat = zeros(4,3*nel);
%         BBmat = zeros(6,3*nel);

        lam = getlam(mateprop);
        
        I1 = [1; 1; 0; 0];
        spvec0 = I1;
        
        thick = 1;
        
        % Load Guass Integration Points

        if nel == 3
            if exist('iprob','var') == 1 && iprob == 6
                lint = 25;
            else
                lint = 7; %minimum of 7 for all integrals in deformed state
            end
        elseif nel == 6
            if exist('iprob','var') == 1 && iprob == 6
                lint = 25;
            else
                lint = 13; %minimum of 13 for all integrals in deformed state
            end
        elseif nel == 4
            if exist('iprob','var') == 1 && iprob == 6
                lint = 100;4;
            else
                lint = 16;
            end
        else
            if exist('iprob','var') == 1 && iprob == 6
                lint = 100;
            else
                lint = 16;
            end
        end
        bf = 0;
        ib = 0;

        [Qxy, shgs, Jdet0,be] = shgt(xl,nel,shld,shls,nel,bf,der,be);
        
        %Integration Loop
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

            if nelP == 3 || nelP == 6
              [shlp,shld,shls] = shlt(litr,lits,nelP,nel,0,0);
              [Pxy, shgp] = shgt(xl+ul(1:2,:),nelP,shld,shls,nelP,0,0,be);
            else
              [shlp,shld,shls] = shlq(litr,lits,nelP,nel,0,0);
              [Pxy, shgp] = shgq(xl+ul(1:2,:),nelP,shld,shls,nelP,0,0,be);
            end

            if Jdet < 0.01*Jdet0
                Jdet
            end
%             Nmat(:,(ie-1)*3+1:3*ie) = [shlp(ie,1) 0         0
%                                        0         shlp(ie,1) 0
%                                        0         0         shlp(ie,1)];
            if nel == 3            
            % Form B matrix    
            Bmat = [Qxy(1,1)  0        0         Qxy(2,1)  0        0         Qxy(3,1)  0        0
                    0         Qxy(1,2) 0         0         Qxy(2,2) 0         0         Qxy(3,2) 0         
                    Qxy(1,2)  Qxy(1,1) 0         Qxy(2,2)  Qxy(2,1) 0         Qxy(3,2)  Qxy(3,1) 0         
                    Qxy(1,2) -Qxy(1,1) 0         Qxy(2,2) -Qxy(2,1) 0         Qxy(3,2) -Qxy(3,1) 0         
                    0        0         shlp(1,1) 0        0         shlp(2,1) 0        0         shlp(3,1)];
            elseif nel == 4            
            % Form B matrix
            Bmat = [Qxy(1,1)  0        0         Qxy(2,1)  0        0         Qxy(3,1)  0        0         Qxy(4,1)  0        0
                    0         Qxy(1,2) 0         0         Qxy(2,2) 0         0         Qxy(3,2) 0         0         Qxy(4,2) 0         
                    Qxy(1,2)  Qxy(1,1) 0         Qxy(2,2)  Qxy(2,1) 0         Qxy(3,2)  Qxy(3,1) 0         Qxy(4,2)  Qxy(4,1) 0         
                    Qxy(1,2) -Qxy(1,1) 0         Qxy(2,2) -Qxy(2,1) 0         Qxy(3,2) -Qxy(3,1) 0         Qxy(4,2) -Qxy(4,1) 0         
                    0        0         shlp(1,1) 0        0         shlp(2,1) 0        0         shlp(3,1) 0        0         shlp(4,1)];
            elseif nel == 6
            % Form B matrix
            Bmat = [Qxy(1,1)  0        0         Qxy(2,1)  0        0         Qxy(3,1)  0        0         Qxy(4,1)  0        0         Qxy(5,1)  0        0         Qxy(6,1)  0        0
                    0         Qxy(1,2) 0         0         Qxy(2,2) 0         0         Qxy(3,2) 0         0         Qxy(4,2) 0         0         Qxy(5,2) 0         0         Qxy(6,2) 0
                    Qxy(1,2)  Qxy(1,1) 0         Qxy(2,2)  Qxy(2,1) 0         Qxy(3,2)  Qxy(3,1) 0         Qxy(4,2)  Qxy(4,1) 0         Qxy(5,2)  Qxy(5,1) 0         Qxy(6,2)  Qxy(6,1) 0
                    Qxy(1,2) -Qxy(1,1) 0         Qxy(2,2) -Qxy(2,1) 0         Qxy(3,2) -Qxy(3,1) 0         Qxy(4,2) -Qxy(4,1) 0         Qxy(5,2) -Qxy(5,1) 0         Qxy(6,2) -Qxy(6,1) 0
                    0        0         shlp(1,1) 0        0         shlp(2,1) 0        0         shlp(3,1) 0        0         shlp(4,1) 0        0         shlp(5,1) 0        0         shlp(6,1)];
            elseif nel == 9
            % Form B matrix
            Bmat = [Qxy(1,1)  0        0         Qxy(2,1)  0        0         Qxy(3,1)  0        0         Qxy(4,1)  0        0         Qxy(5,1)  0        0         Qxy(6,1)  0        0         Qxy(7,1)  0        0         Qxy(8,1)  0        0         Qxy(9,1)  0        0
                    0         Qxy(1,2) 0         0         Qxy(2,2) 0         0         Qxy(3,2) 0         0         Qxy(4,2) 0         0         Qxy(5,2) 0         0         Qxy(6,2) 0         0         Qxy(7,2) 0         0         Qxy(8,2) 0         0         Qxy(9,2) 0
                    Qxy(1,2)  Qxy(1,1) 0         Qxy(2,2)  Qxy(2,1) 0         Qxy(3,2)  Qxy(3,1) 0         Qxy(4,2)  Qxy(4,1) 0         Qxy(5,2)  Qxy(5,1) 0         Qxy(6,2)  Qxy(6,1) 0         Qxy(7,2)  Qxy(7,1) 0         Qxy(8,2)  Qxy(8,1) 0         Qxy(9,2)  Qxy(9,1) 0
                    Qxy(1,2) -Qxy(1,1) 0         Qxy(2,2) -Qxy(2,1) 0         Qxy(3,2) -Qxy(3,1) 0         Qxy(4,2) -Qxy(4,1) 0         Qxy(5,2) -Qxy(5,1) 0         Qxy(6,2) -Qxy(6,1) 0         Qxy(7,2) -Qxy(7,1) 0         Qxy(8,2) -Qxy(8,1) 0         Qxy(9,2) -Qxy(9,1) 0
                    0        0         shlp(1,1) 0        0         shlp(2,1) 0        0         shlp(3,1) 0        0         shlp(4,1) 0        0         shlp(5,1) 0        0         shlp(6,1) 0        0         shlp(7,1) 0        0         shlp(8,1) 0        0         shlp(9,1)];
            end

            [fi,JxX,F] = kine2d(Qxy,-ul,nel,0); %this is equivalent to ikine2d
            JxX = 1/JxX; %this is equivalent to ikine2d
    %         [F,JxX,fi] = ikine2d(Qxy,ul,nel); %this is equivalent to above two lines
            Jdet = Jdet/JxX;
            [theta1,theta2,theta3] = ThetaNS(JxX,mateprop);
            spvec = JxX*spvec0;
                
            c1 = Wgt*Jdet*thick;

            sigmai = SigmaCmatNSCST2i(F,JxX,mateprop);
            press = ul(3,:)*shlp;
            sigmap = press*spvec;
            sigma = sigmai + sigmap;
            
            Svec = [sigma; (theta1-press/lam)];
            
            ElemF = ElemF - c1*(Bmat'*(Svec));

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
            lint = 10; % Use 10 for body force BF2U4M0.m
        else
            lint = 4;
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
                    X = xl(1,:)*shl;
                    Y = xl(2,:)*shl;
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

                ElemF(ndf*o-2) = ElemF(ndf*o-2) + F(1)*c1;

                ElemF(ndf*o-1)   = ElemF(ndf*o-1)   + F(2)*c1;

            end %o

        end %ie
        
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
                lint = 25;
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
            shlp = shl;
                
            if nel == 3
            % Form N matrix 
            Nmat = [shl(1)  0        0 shl(2)  0        0 shl(3)  0  0     
                    0         shl(1) 0  0         shl(2) 0 0         shl(3) 0
                    0        0         shlp(1,1) 0        0         shlp(2,1) 0        0         shlp(3,1)];
            elseif nel == 4
            % Form N matrix
            Nmat = [shl(1)  0        0 shl(2)  0        0 shl(3)  0        0 shl(4)  0        0
                    0         shl(1) 0 0         shl(2) 0 0         shl(3) 0 0         shl(4) 0
                    0        0         shlp(1,1) 0        0         shlp(2,1) 0        0         shlp(3,1) 0        0         shlp(4,1)]; % Added back for body force BF2U4M0.m
            elseif nel == 6
            % Form N matrix
            Nmat = [shl(1)  0        0 shl(2)  0        0 shl(3)  0        0 shl(4)  0        0 shl(5)  0        0 shl(6)  0        0
                    0         shl(1) 0 0         shl(2) 0 0         shl(3) 0 0         shl(4) 0 0         shl(5) 0 0         shl(6) 0
                    0        0         shlp(1,1) 0        0         shlp(2,1) 0        0         shlp(3,1) 0        0         shlp(4,1) 0        0         shlp(5,1) 0        0         shlp(6,1)];
            elseif nel == 9
            % Form N matrix
            Nmat = [shl(1)  0      0 shl(2)  0      0 shl(3)  0      0 shl(4)  0      0 shl(5)  0      0 shl(6)  0      0 shl(7)  0      0 shl(8)  0      0 shl(9)  0      0     
                    0       shl(1) 0 0       shl(2) 0 0       shl(3) 0 0       shl(4) 0 0       shl(5) 0 0       shl(6) 0 0       shl(7) 0 0       shl(8) 0 0       shl(9) 0
                    0       0      shlp(1,1) 0    0 shlp(2,1) 0    0 shlp(3,1) 0    0 shlp(4,1) 0    0 shlp(5,1) 0    0 shlp(6,1) 0    0 shlp(7,1) 0    0 shlp(8,1) 0    0 shlp(9,1)];
            end
                
            c1 = Wgt*Jdet*thick;
            
            if exist('iprob','var')
                if iprob == 6
                    mu = 40;
                    X = xl(1,1:nel)*shl;
                    Y = xl(2,1:nel)*shl;
                    fb = [ - (101*mu)/(101*X + 100) - (101*lam*(101*X + 100))/10000
                                                      0
                                                      0];
                else
                    fb = [bodyf(1:2)'; 0];
                end
            else
                fb = [bodyf(1:2)'; 0];
            end
            
            ElemF = ElemF + c1*Nmat'*fb;

        end %je
    ElemF;
%%
    case 22
        
        ElemM = zeros(nst);
        ElemF = zeros(nst,1);

        lam = getlam(mateprop);
        
        I1 = [1; 1; 1; 0; 0; 0; 0; 0; 0];
        stresID = [1 2 4 0 1 3];
        spvec0 = I1;
        
        thick = 1;
        
        % Load Guass Integration Points

        if nel == 3
            lint = 7; %minimum of 7 for all integrals in deformed state
            der = 0;
        elseif nel == 4
%             lint = 4;
            lint = 16;
            der = 1;
        elseif nel == 6
            lint = 13; %minimum of 13 for all integrals in deformed state
            der = 1;
        else
%             lint = 9;
            lint = 25;
            der = 0;
        end

        %Integration Loop
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

            if nelP == 3 || nelP == 6
              [shlp,shld,shls] = shlt(litr,lits,nelP,nel,0,0);
              [Pxy, shgp] = shgt(xl+ul(1:2,:),nelP,shld,shls,nelP,0,0,be);
            else
              [shlp,shld,shls] = shlq(litr,lits,nelP,nel,0,0);
              [Pxy, shgp] = shgq(xl+ul(1:2,:),nelP,shld,shls,nelP,0,0,be);
            end

            if nelS == 3 || nelS == 6
              [shlS,shld,shls] = shlt(litr,lits,nelS,nel,0,0);
            else
              [shlS,shld,shls] = shlq(litr,lits,nelS,nel,0,0);
            end
            
            [fi,JxX,F] = kine2d(Qxy,-ul,nel,0); %this is equivalent to ikine2d
            JxX = 1/JxX; %this is equivalent to ikine2d
    %         [F,JxX,fi] = ikine2d(Qxy,ul,nel); %this is equivalent to above two lines
            Jdet = Jdet/JxX;
            [theta1,theta2,theta3] = ThetaNS(JxX,mateprop);
            spvec = JxX*spvec0;

            sigmai = SigmaCmatNSCST3i([F zeros(2,1); zeros(1,2) 1],JxX,mateprop);
            press = ul(3,:)*shlp;
            sigmap = press*spvec;
            sigma = (sigmai + sigmap)/JxX;
            
            % Form B matrix
            Nmat = shlS';

            w = Wgt*Jdet*thick;
            
            if stres <= 3 % stress components
                sigmas = sigma(stresID(stres));
            elseif stres >= 5
                if stres <= 6 % principal stresses
                sigma = [sigma(1) sigma(4) sigma(6); sigma(4) sigma(2) sigma(5); sigma(6) sigma(5) sigma(3)];
                psig = eig(sigma);
                sigmas = psig(stresID(stres));
                else % hydrostatic stress
                sigmas = 1/3*sigma'*I1;
                end
            else % von Mises stress
                trs = sigma'*I1;
                dsig = sigma - 1/3*trs*I1;
                sigmas = sqrt(3/2*(dsig'*dsig));
            end
            
            ElemF = ElemF + w*Nmat'*sigmas;
            
            ElemM = ElemM + w*(Nmat'*Nmat);

        end %je
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

            if nelP == 3 || nelP == 6
              shlp = shlt(litr,lits,nelP,nel,0,0);
            else
              shlp = shlq(litr,lits,nelP,nel,0,0);
            end
            
            [fi,JxX,F] = kine2d(Qxy,-ul,nel,0); %this is equivalent to ikine2d
            JxX = 1/JxX; %this is equivalent to ikine2d
    %         [F,JxX,fi] = ikine2d(Qxy,ul,nel); %this is equivalent to above two lines
            Jdet = Jdet/JxX;
            [theta1,theta2,theta3] = ThetaNS(JxX,mateprop);
            spvec = JxX*spvec0;

            sigmai = SigmaCmatNSCST3i([F zeros(2,1); zeros(1,2) 1],JxX,mateprop);
            press = ul(3,:)*shlp;
            sigmap = press*spvec;
            sigma = (sigmai + sigmap)/JxX;
            
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
        ElemS(i,npstr+1) = 1; % use this for simple average Vol; % use this for weighted average 
        end
end