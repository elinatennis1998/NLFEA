
% Tim Truster
% 04/07/2014
% 2-D element for crystal plasticity 
% Bilinear quadrilateral element with consistent tangent
% Uses elastic strain as history variable

% Elastic tangent verified 4/7/14; error was in incorrect arrangement of
% Bmat

% Uses deSouza implementation

% Tested using CrysPlasTest, 2, and 3.m; seems to converge quadratically or
% at least pretty well. Results for applied displacement match with the
% book results for shear.

switch isw %Task Switch
%%    
    case 1
        
        if ndf > 2
            
            for i = 3:ndf
                lie(i,1) = 0;
            end
            
        end
        
        nh1 = 5*4;
        % Sets number of Lie variables for interpolation, later
        numLie = 18; % tensors Re and Ue
        
%%
    case 3
        
        % Get material data from mateprop
        IPROPS = mateprop(1:3);
        NTYPE = 2;
        RPROPS = mateprop(4:end);
        
        ElemK = zeros(nst);
        ElemF = zeros(nst,1);
        Bmat = zeros(4,2*nel);
        
        % Load Guass Integration Points

            lint = 4;
        der = 0;
        bf = 0;
        ib = 0;


        %Integration Loop
        for l = 1:lint

                  [Wgt,litr,lits] =  intpntq(l,lint,ib);
                  [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
                  [QXY, shgs, Jdet] = shgq(xl,nel,shld,shls,nel,bf,der,be);

            % Push forward kinematics
            [Fn1,JxX,fi,df,Fn,Qxy] = kine2m(QXY,ul,uld,nel,1);
            df = inv(Fn*fi);
            JxXi = det(df);
            
%                 Fn1 = [Fn1 zeros(2,1); zeros(1,2) 1]; %#ok<AGROW>
%                 Fn = [Fn zeros(2,1); zeros(1,2) 1]; %#ok<AGROW>
%                 df = [df zeros(2,1); zeros(1,2) 1]; %#ok<AGROW>
        
            for mm = 1:nel    
 Bmat(:,2*mm-1:2*mm) = [Qxy(mm,1) 0        
                        0 Qxy(mm,1)
                        Qxy(mm,2) 0
                        0         Qxy(mm,2)];
            end

            % Compute input for Radial Return
            ephr = nh1-1+(l-1)*5; %pointer for plastic strain at pt l
            RSTAVN = hr(ephr+1:ephr+5);
                
            [DGAM,LALGVA,RSTAVA,STRES] = SUPDSC(df,IPROPS,NTYPE,RPROPS,RSTAVN);
            if max(DGAM) > 0
                EPFLAG = 1;
            else
                EPFLAG = 0;
            end
            AMATX = CSTPDS(DGAM,EPFLAG,df,IPROPS,LALGVA,NTYPE,RPROPS,RSTAVA,RSTAVN,STRES);
            
            % Store history variables
            ephr = nh2-1+(l-1)*5; %pointer for plastic strain at pt l
            hr(ephr+1:ephr+5) = RSTAVA;
            
            sigma = STRES([1 3 3 2])'*JxX(1);
            cmat = AMATX*JxX(1);

            % Update integration weighting factor
            W = Wgt*Jdet;

            ElemF = ElemF - W*Bmat'*(sigma);
            ElemK = ElemK + W*Bmat'*(cmat)*Bmat;
            
        end %je
   ElemK;     
%%
    case 6
        
        ElemF = zeros(nst,1);
        Bmat = zeros(4,2*nel);
        
        % Get material data from mateprop
        IPROPS = mateprop(1:3);
        NTYPE = 2;
        RPROPS = mateprop(4:end);
        
        % Load Guass Integration Points

            lint = 4;
        der = 0;
        bf = 0;
        ib = 0;


        %Integration Loop
        for l = 1:lint

                  [Wgt,litr,lits] =  intpntq(l,lint,ib);
                  [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
                  [QXY, shgs, Jdet] = shgq(xl,nel,shld,shls,nel,bf,der,be);

            % Push forward kinematics
            [Fn1,JxX,fi,df,Fn,Qxy] = kine2m(QXY,ul,uld,nel,1);
            df = inv(Fn*fi);
            JxXi = det(df);
            
%                 Fn1 = [Fn1 zeros(2,1); zeros(1,2) 1]; %#ok<AGROW>
%                 Fn = [Fn zeros(2,1); zeros(1,2) 1]; %#ok<AGROW>
%                 df = [df zeros(2,1); zeros(1,2) 1]; %#ok<AGROW>
        
            for mm = 1:nel    
 Bmat(:,2*mm-1:2*mm) = [Qxy(mm,1) 0        
                        0         Qxy(mm,2)
                        Qxy(mm,1) 0
                        0         Qxy(mm,2)];
            end

            % Compute input for Radial Return
            ephr = nh1-1+(l-1)*5; %pointer for plastic strain at pt l
            RSTAVN = hr(ephr+1:ephr+5);
                
            [DGAM,LALGVA,RSTAVA,STRES] = SUPDSC(df,IPROPS,NTYPE,RPROPS,RSTAVN);
            if max(DGAM) > 0
                EPFLAG = 1;
            else
                EPFLAG = 0;
            end
            
            % Store history variables
            ephr = nh2-1+(l-1)*5; %pointer for plastic strain at pt l
            hr(ephr+1:ephr+5) = RSTAVA;
            
            sigma = STRES([1 3 3 2])';

            % Update integration weighting factor
            W = Wgt*Jdet;

            ElemF = ElemF - W*Bmat'*(sigma);
            
        end %je
   ElemF;
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
    case 24 % Output plasticity results, at integration points
        
        ElemP = zeros(12,4);
        
        % Initialize Matrix and Vector

        nst = nel*ndf;
        
        % Load Guass Integration Points
        lint = 4;

        der = 0;
        bf = 0;
        ib = 0;
        
        % Loop over integration points
        for l = 1:lint

                  [Wgt,litr,lits] =  intpntq(l,lint,ib);
                  [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
                  [QXY, shgs, Jdet] = shgq(xl,nel,shld,shls,nel,bf,der,be);

            % Push forward kinematics
            [Fn1,JxX,fi,df,Fn,Qxy] = kine2m(QXY,ul,uld,nel,1);
            df = inv(Fn*fi);
            JxXi = det(df);
            

            % Compute input for Radial Return
            ephr = nh1-1+(l-1)*5; %pointer for plastic strain at pt l
            RSTAVN = hr(ephr+1:ephr+5);
                
            [DGAM,LALGVA,RSTAVA,STRES] = SUPDSC(df,IPROPS,NTYPE,RPROPS,RSTAVN);
                
                % store outputs: stress and total hardening
                
                ElemP(1,l) = STRES(1);
                ElemP(2,l) = STRES(2);
                ElemP(3,l) = STRES(4);
                ElemP(4,l) = STRES(3);
                ElemP(5,l) = RSTAVA(5);

        end %je
                
        
    case 40 % Initialize Intermediate Configuration
        
        lint = 4;
        
        % Loop over integration points
        for l = 1:lint
                
                % Store history variables
                ephr = nh1-1+(l-1)*5; %pointer for plastic strain at pt l
                hr(ephr+1) = 1.d0;
                hr(ephr+4) = 1.d0;

        end %je
        
    case 98 % Nye tensor field post-processing

        ElemL = zeros(nel,numLie+1);
        ElemL2 = zeros(numLie,nel);

        % Loop over nodes, convert Lie group to Lie algebra
        for node = 1:nel
        
            ephr = nh2-1+(node-1)*5; %pointer for elastic deformation gradient n+1
            Fe2d = reshape(hr(ephr+1:ephr+4),2,2);
            Fe = [Fe2d zeros(2,1); zeros(1,2) 1.0];
            
            % Polar decomposition of Fe
            [Re_grp, Ue_grp] = poldec(Fe);
            
            % convert rotation
            Re_alg = logm(Re_grp);
            ElemL2(1:9,node) = reshape(Re_alg,1,9);

            % convert stretch
            Ue_alg = logm(Ue_grp);
            ElemL2(10:18,node) = reshape(Ue_alg,1,9);
            
        end
        
        % Extrapolation parametric locations for nodes
        if nel == 4
            lint = 4;
            nelS = 4;
            sqr3 = sqrt(3);
            plist = [-sqr3 sqr3 sqr3 -sqr3
                     -sqr3 -sqr3 sqr3 sqr3];
            %          n1    n2   n3   n4
            % Map history ordering to node ordering
            ElemL2 = ElemL2(:,[1 2 4 3]);
%         elseif nel == 6
%             lint = 7;
%             nint = 3;
%             plist = [-1/3 5/3 -1/3 2/3 2/3 -1/3
%                      -1/3 -1/3 5/3 -1/3 2/3 2/3];
%         else
%             lint = 9;
%             nint = 4;
%             plist = [-sqr3 sqr3 sqr3 -sqr3 0 sqr3 0 -sqr3 0
%                      -sqr3 -sqr3 sqr3 sqr3 -sqr3 0 sqr3 0 0];
        end
        
        der = 0;
        bf = 0;
        
        % Extrapolate Lie algebra values from integration points to nodes
        for ll = 1:nelS
            
            r = plist(1,ll);
            s = plist(2,ll);
            
            %Evaluate first derivatives of basis functions at int. point
            if nel == 3 || nel == 6
              shl = shlt(r,s,nel,nel,der,bf);
            else
              shl = shlq(r,s,nel,nel,der,bf);
            end
            
            % Extrapolation is a product of shape functions and integration
            % point values
            Extra_vals = ElemL2*shl;
            ElemL(ll,1:numLie) = Extra_vals';
            
        end
        
        %Integration Loop for area
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

            w = Wgt*Jdet;
            
            Vol = Vol + w;

        end %je

        for i = 1:nel
        ElemL(i,numLie+1) = 1; % use this for simple average Vol; % use this for weighted average 
        end
        
    case 99 % Convert Lie algebra to Lie group
        
        Re_alg = ElemL(1:9);
        Re_alg = reshape(Re_alg,3,3);
        Re_grp = expm(Re_alg);
        ElemL(1:9) = reshape(Re_grp,1,9);
        
        Ue_alg = ElemL(10:18);
        Ue_alg = reshape(Ue_alg,3,3);
        Ue_grp = expm(Ue_alg);
        ElemL(10:18) = reshape(Ue_grp,1,9);
        
end