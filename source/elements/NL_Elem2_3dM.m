% function ElemK = NL_Elem2_2d(ElemXYZ, nel, ndf, R, S, T, pr, ps, PatchE, Patchv)
%
% Tim Truster
% 01/04/2012
% UIUC

% Master pure-displacement element, 3D
% Verified against NL_Elem2_3d.m on 01/04/2012 using PT2A10L1b.m,
% BC2U8L1.m, and BC2U4L1.m.
% Stress post-processing functions and gives meaningful results, not yet
% verified.


switch isw %Task Switch
%%

    case 1
        
        if ndf > 3
            
            for i = 4:ndf
                lie(i,1) = 0;
            end
            
        end
        
    case 3
        
        ElemK = zeros(nst);
        ElemF = zeros(nst,1);
        Bmat = zeros(9,3*nel);

        lam = getlam(mateprop);
        
        % Load Guass Integration Points

        if nel == 4
            lint = 1;4;11;5;16;
        elseif nel == 8
            lint = 8;1000; %1000 for body force problem 
        elseif nel == 10
            lint = 14;
%             lint = 27;
        else
            lint = 27;
        end
        der = 0;
        bf = 0;
        ib = 0;

        %Integration Loop
        for ll = 1:lint

                %Evaluate first derivatives of basis functions at int. point
                if nel == 4 || nel == 10
                  [Wgt,ss] =  int3d_t(ll,lint,ib);
                  [shl,shld,shls,be] = shltt(ss,nel,nel,der,bf);
                  [Qxy, shgs, Jdet] = shgtt(xl+ul,nel,shld,shls,nel,bf,der,be);
                else
                  [Wgt,ss] =  intpntb(ll,lint,ib);
                  [shl,shld,shls,be] = shlb(ss,nel,nel,der,bf);
                  [Qxy, shgs, Jdet] = shgb(xl+ul,nel,shld,shls,nel,bf,der,be);
                end
                
            for mm = 1:nel  
 Nmat(:,3*mm-2:3*mm) = [shl(mm,1)     0          0
                           0        shl(mm,1)     0
                           0            0       shl(mm,1) ];  % for body force part
 Bmat(:,3*mm-2:3*mm) = [Qxy(mm,1) 0         0         
                        0         Qxy(mm,2) 0         
                        0         0         Qxy(mm,3) 
                        Qxy(mm,2) Qxy(mm,1) 0         
                        0         Qxy(mm,3) Qxy(mm,2) 
                        Qxy(mm,3) 0         Qxy(mm,1) 
                        Qxy(mm,2) -Qxy(mm,1) 0         
                        0         Qxy(mm,3) -Qxy(mm,2) 
                        -Qxy(mm,3) 0         Qxy(mm,1) ];
            end
                
            [fi,JxX,F] = kine3d(Qxy,-ul(:,1:nel),nel,0); %this is equivalent to ikine2d
            JxX = 1/JxX; %this is equivalent to ikine2d
    %         [F,JxX,fi] = ikine2d(Qxy,ul,nel); %this is equivalent to above two lines
            Jdet = Jdet/JxX;
                
            c1 = Wgt*Jdet;
            
            [sigma2, cmat] = SigmaCmat3(F,JxX,mateprop,lam);
            
            Smat = ...
[    sigma2(1),        0,        0,           sigma2(4)/2,                 0,           sigma2(6)/2,           sigma2(4)/2,                 0,          -sigma2(6)/2
         0,    sigma2(2),        0,           sigma2(4)/2,           sigma2(5)/2,                 0,          -sigma2(4)/2,           sigma2(5)/2,                 0
         0,        0,    sigma2(3),                 0,           sigma2(5)/2,           sigma2(6)/2,                 0,          -sigma2(5)/2,           sigma2(6)/2
   sigma2(4)/2,  sigma2(4)/2,        0, sigma2(1)/4 + sigma2(2)/4,           sigma2(6)/4,           sigma2(5)/4, sigma2(2)/4 - sigma2(1)/4,           sigma2(6)/4,          -sigma2(5)/4
         0,  sigma2(5)/2,  sigma2(5)/2,           sigma2(6)/4, sigma2(2)/4 + sigma2(3)/4,           sigma2(4)/4,          -sigma2(6)/4, sigma2(3)/4 - sigma2(2)/4,           sigma2(4)/4
   sigma2(6)/2,        0,  sigma2(6)/2,           sigma2(5)/4,           sigma2(4)/4, sigma2(1)/4 + sigma2(3)/4,           sigma2(5)/4,          -sigma2(4)/4, sigma2(1)/4 - sigma2(3)/4
   sigma2(4)/2, -sigma2(4)/2,        0, sigma2(2)/4 - sigma2(1)/4,          -sigma2(6)/4,           sigma2(5)/4, sigma2(1)/4 + sigma2(2)/4,          -sigma2(6)/4,          -sigma2(5)/4
         0,  sigma2(5)/2, -sigma2(5)/2,           sigma2(6)/4, sigma2(3)/4 - sigma2(2)/4,          -sigma2(4)/4,          -sigma2(6)/4, sigma2(2)/4 + sigma2(3)/4,          -sigma2(4)/4
  -sigma2(6)/2,        0,  sigma2(6)/2,          -sigma2(5)/4,           sigma2(4)/4, sigma2(1)/4 - sigma2(3)/4,          -sigma2(5)/4,          -sigma2(4)/4, sigma2(1)/4 + sigma2(3)/4];
            
            Smat = Smat + cmat;          

            ElemF(1:ndf*nel) = ElemF(1:ndf*nel) - c1*Bmat'*(sigma2);
            
            ElemK(1:ndf*nel,1:ndf*nel) = ElemK(1:ndf*nel,1:ndf*nel) + c1*(Bmat'*Smat*Bmat);
            
        end %je
   ElemK;  
%%
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
            
            if exist('iprob','var')
                if iprob == 6
                    mu = 40;
                    lam = getlam(mateprop);
                    X = xl(1,1:nel)*shl;
                    Y = xl(2,1:nel)*shl;
                    P = [ (101*X*lam*(101*X + 100))/10000, (101*Y*mu)/(101*X + 100) - (10201*X*Y*lam)/10000, 0
                         (101*Y*mu)/100, mu - (100*mu)/(101*X + 100) + X*((101*lam)/100 + (101*mu)/100),     0
                              0,                                            0, (101*X*lam*(101*X + 100))/10000];
                    Traction = P*tu3';
                elseif iprob == 7
                    X = xl(1,1:nel)*shl;
                    Y = xl(2,1:nel)*shl;
                    Z = xl(3,1:nel)*shl;
                    if exist('bimat','var') && bimat == 1 
                        if ma == 1
                            Traction = TensTorsTB(X,Y,Z,mateprop(4),mateprop(5),Twist/L1,0,tu3',eye(3));
                        else
                            Z = Z - L1;
                            Q = [cos(pi*Twist) -sin(pi*Twist)
                                 sin(pi*Twist)  cos(pi*Twist)];
                            XY = Q*[X; Y];
                            Twist2 = MateT(1,4)/MateT(2,4)*Twist;
                            NQ = [Q zeros(2,1); 0 0 1]*tu3';
                            Traction = TensTorsTB(XY(1),XY(2),Z,mateprop(4),mateprop(5),Twist2/L2,0,NQ,eye(3));
                        end
                    else
                        Traction = TensTorsTB(X,Y,Z,mateprop(4),mateprop(5),Twist/((L1+L2)*(1+Axialdelta)),Axialdelta,tu3',eye(3));
                    end
%                     [Traction,fb] = TensTorsTB2(X,Y,Z,mateprop(4),mateprop(5),Twist,tu3');
                elseif iprob == 9
                    X = xl(1,1:nel)*shl;
                    Y = xl(2,1:nel)*shl;
                    Z = xl(3,1:nel)*shl;
                    Traction = BlockJ1TB(X,Y,Z,mateprop(4),mateprop(5),tu3',eye(3));
                else
                    Traction = traction;
                end
            else
                Traction = traction;
            end

            %Force components are positive in positive coord. direction
            c1 = Wgt*tm3;
            for o=1:nel
                don = shl(o);
                F = don*Traction';  %traction for the one without body force

                ElemF(ndf*o-2) = ElemF(ndf*o-2) + F(1)*c1;

                ElemF(ndf*o-1) = ElemF(ndf*o-1) + F(2)*c1;

                ElemF(ndf*o-0) = ElemF(ndf*o-0) + F(3)*c1;

            end %o
            
        end %je
        
        % Reorder nodes back to the orientation of the element
        SurfOrientStoE
        ElemF(1:ndf*nel) = ElemF(ilist2);
%%
    case 6
        
        ElemF = zeros(nst,1);
        Bmat = zeros(9,3*nel);

        lam = getlam(mateprop);
        
        % Load Guass Integration Points

        if nel == 4
            lint = 1;11;5;16;
        elseif nel == 8
            lint = 8;
        elseif nel == 10
            lint = 14;
%             lint = 27;
        else
            lint = 27;
        end
        der = 0;
        bf = 0;
        ib = 0;

        %Integration Loop
        for ll = 1:lint

                %Evaluate first derivatives of basis functions at int. point
                if nel == 4 || nel == 10
                  [Wgt,ss] =  int3d_t(ll,lint,ib);
                  [shl,shld,shls,be] = shltt(ss,nel,nen,der,bf);
                  [Qxy, shgs, Jdet] = shgtt(xl+ul,nen,shld,shls,nen,bf,der,be);
                else
                  [Wgt,ss] =  intpntb(ll,lint,ib);
                  [shl,shld,shls,be] = shlb(ss,nel,nen,der,bf);
                  [Qxy, shgs, Jdet] = shgb(xl+ul,nen,shld,shls,nen,bf,der,be);
                end
                

            for mm = 1:nel    
 Bmat(:,3*mm-2:3*mm) = [Qxy(mm,1) 0         0         
                        0         Qxy(mm,2) 0         
                        0         0         Qxy(mm,3) 
                        Qxy(mm,2) Qxy(mm,1) 0         
                        0         Qxy(mm,3) Qxy(mm,2) 
                        Qxy(mm,3) 0         Qxy(mm,1) 
                        Qxy(mm,2) -Qxy(mm,1) 0         
                        0         Qxy(mm,3) -Qxy(mm,2) 
                        -Qxy(mm,3) 0         Qxy(mm,1) ];
            end
                
                
            [fi,JxX,F] = kine3d(Qxy(1:nel,:),-ul,nel,0); %this is equivalent to ikine2d
            JxX = 1/JxX; %this is equivalent to ikine2d
    %         [F,JxX,fi] = ikine2d(Qxy,ul,nel); %this is equivalent to above two lines
            Jdet = Jdet/JxX;
                
            c1 = Wgt*Jdet;
            
            sigma2 = SigmaCmat3(F,JxX,mateprop,lam);
            
            ElemF(1:ndf*nel) = ElemF(1:ndf*nel) - c1*(Bmat'*(sigma2));
            
        end %je
        
%%
    case 15 % Body Force
        
        ElemF = zeros(nst,1);
        Nmat = zeros(3,nst);

        lam = getlam(mateprop);
        
        % Load Guass Integration Points

        if nel == 4
            lint = 11;5;16;
        elseif nel == 8
            lint = 8;1000; %1000 for body force
        elseif nel == 10
            lint = 14;
%             lint = 27;
        else
            lint = 64;27;
        end
        der = 0;
        bf = 0;
        ib = 0;

        %Integration Loop
        for ll = 1:lint
            
            %Evaluate first derivatives of basis functions at int. point
            if nel == 4 || nel == 10
              [Wgt,ss] =  int3d_t(ll,lint,ib);
              [shl,shld,shls,be] = shltt(ss,nel,nel,der,bf);
              [Qxy, shgs, Jdet] = shgtt(xl,nel,shld,shls,nel,bf,der,be);
            else
              [Wgt,ss] =  intpntb(ll,lint,ib);
              [shl,shld,shls,be] = shlb(ss,nel,nel,der,bf);
              [Qxy, shgs, Jdet] = shgb(xl,nel,shld,shls,nel,bf,der,be);
            end
                
            for mm = 1:nel    
                Nmat(:,3*mm-2:3*mm) = shl(mm)*eye(3);
            end
                
            c1 = Wgt*Jdet;
            
            if exist('iprob','var')
                if iprob == 6
                    mu = 40;
                    X = xl(1,1:nel)*shl;
                    Y = xl(2,1:nel)*shl;
                    fb = [ - (101*mu)/(101*X + 100) - (101*lam*(101*X + 100))/10000
                                                      0; 0];
                elseif iprob == 7
                    X = xl(1,1:nel)*shl;
                    Y = xl(2,1:nel)*shl;
                    Z = xl(3,1:nel)*shl;
                    if exist('bimat','var') && bimat == 1 
                        if ma == 1
                            [Traction,fb] = TensTorsTB(X,Y,Z,mateprop(4),mateprop(5),Twist/L1,0,[0 0 1]',eye(3));
                        else
                            Z = Z - L1;
                            Q = [cos(pi*Twist) -sin(pi*Twist)
                                 sin(pi*Twist)  cos(pi*Twist)];
                            XY = Q*[X; Y];
                            Twist2 = MateT(1,4)/MateT(2,4)*Twist;
                            [Traction,fb] = TensTorsTB(XY(1),XY(2),Z,mateprop(4),mateprop(5),Twist2/L2,0,[0 0 1]',eye(3));
                        end
                    else
                        [Traction,fb] = TensTorsTB(X,Y,Z,mateprop(4),mateprop(5),Twist/((L1+L2)*(1+Axialdelta)),Axialdelta,[0 0 1]',eye(3));
                    end
%                     [Traction,fb] = TensTorsTB2(X,Y,Z,mateprop(4),mateprop(5),Twist,[0 0 1]');
                elseif iprob == 9
                    X = xl(1,1:nel)*shl;
                    Y = xl(2,1:nel)*shl;
                    Z = xl(3,1:nel)*shl;
                    [Traction,fb] = BlockJ1TB(X,Y,Z,mateprop(4),mateprop(5),[0 0 1]',eye(3));
                else
                    fb = bodyf(1:3)';
                end
            else
                fb = bodyf(1:3)';
            end
            
            ElemF = ElemF + c1*Nmat'*fb;

        end %je
    ElemF;
      
%%
    case 7 % User defined loading, updated every step/iteration

        ElemF = zeros(nst,1);
        
        % Reorder nodes for corresponding face of integration
        SurfOrientEtoS
        
        t1 = zeros(3,1);
        t2 = t1;
        t3 = t1;
        
        if nel == 4 || nel == 10
            lint = 13;
        else
            lint = 16;4;16;
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
            
            if exist('iprob','var')
                if iprob == 7
                    X = xl(1,1:nel)*shl;
                    Y = xl(2,1:nel)*shl;
                    Z = xl(3,1:nel)*shl;
                    if exist('bimat','var') && bimat == 1 
                        if ma == 1
                            Traction = TensTorsTB(X,Y,Z,mateprop(4),mateprop(5),lamda*Twist/L1,0,tu3',eye(3));
                        else
%                             Z = Z - L1;
%                             Q = [cos(pi*lamda*Twist) -sin(pi*lamda*Twist)
%                                  sin(pi*lamda*Twist)  cos(pi*lamda*Twist)];
%                             XY = Q*[X; Y];
%                             Twist2 = MateT(1,4)/MateT(2,4)*Twist;
%                             NQ = [Q zeros(2,1); 0 0 1]*tu3';
%                             Traction = TensTorsTB(XY(1),XY(2),Z,mateprop(4),mateprop(5),lamda*Twist2/L2,0,NQ,eye(3));
                            Z = Z - L1;
                            Q = [cos(pi*lamda*Twist) -sin(pi*lamda*Twist)
                                 sin(pi*lamda*Twist)  cos(pi*lamda*Twist)];
                            Twist2 = MateT(1,4)/MateT(2,4)*Twist;
                            Traction = TensTorsTB(X,Y,Z,mateprop(4),mateprop(5),lamda*Twist2/L2,0,tu3',[Q zeros(2,1); 0 0 1]);
                        end
                    else
                        Traction = TensTorsTB(X,Y,Z,mateprop(4),mateprop(5),lamda*Twist/((L1+L2)*(1+lamda*Axialdelta)),lamda*Axialdelta,tu3',eye(3));
                    end
%                     Traction = TensTorsTB2(X,Y,Z,mateprop(4),mateprop(5),lamda*Twist,tu3');
                else
                    Traction = zeros(3,1);
                end
            else
                Traction = zeros(3,1);
            end

            %Force components are positive in positive coord. direction
            c1 = Wgt*tm3;
            for o=1:nel
                don = shl(o);
                F = don*Traction';

                ElemF(ndf*o-2) = ElemF(ndf*o-2) + F(1)*c1;

                ElemF(ndf*o-1) = ElemF(ndf*o-1) + F(2)*c1;

                ElemF(ndf*o-0) = ElemF(ndf*o-0) + F(3)*c1;

            end %o
            
        end %je
%         if load == 33
%             load;
%         end
        % Reorder nodes back to the orientation of the element
        SurfOrientStoE
        ElemF(1:ndf*nel) = ElemF(ilist2);
%%
    case -15 % User Body Force
        
        ElemF = zeros(nst,1);
        Nmat = zeros(3,nst);

        lam = getlam(mateprop);
        
        % Load Guass Integration Points

        if nel == 4
            lint = 11;5;16;
        elseif nel == 8
            lint = 8;64;
        elseif nel == 10
            lint = 14;
%             lint = 27;
        else
            lint = 64;27;
        end
        der = 0;
        bf = 0;
        ib = 0;

        %Integration Loop
        for ll = 1:lint
            
            %Evaluate first derivatives of basis functions at int. point
            if nel == 4 || nel == 10
              [Wgt,ss] =  int3d_t(ll,lint,ib);
              [shl,shld,shls,be] = shltt(ss,nel,nel,der,bf);
              [Qxy, shgs, Jdet] = shgtt(xl,nel,shld,shls,nel,bf,der,be);
            else
              [Wgt,ss] =  intpntb(ll,lint,ib);
              [shl,shld,shls,be] = shlb(ss,nel,nel,der,bf);
              [Qxy, shgs, Jdet] = shgb(xl,nel,shld,shls,nel,bf,der,be);
            end
                
            for mm = 1:nel    
                Nmat(:,3*mm-2:3*mm) = shl(mm)*eye(3);
            end
                
            c1 = Wgt*Jdet;
            
            if exist('iprob','var')
                if iprob == 7
                    X = xl(1,1:nel)*shl;
                    Y = xl(2,1:nel)*shl;
                    Z = xl(3,1:nel)*shl;
                    if exist('bimat','var') && bimat == 1 
                        if ma == 1
                            [Traction,fb] = TensTorsTB(X,Y,Z,mateprop(4),mateprop(5),lamda*Twist/L1,0,[0 0 1]',eye(3));
                        else
%                             Z = Z - L1;
%                             Q = [cos(pi*lamda*Twist) -sin(pi*lamda*Twist)
%                                  sin(pi*lamda*Twist)  cos(pi*lamda*Twist)];
%                             XY = Q*[X; Y];
%                             Twist2 = MateT(1,4)/MateT(2,4)*Twist;
%                             [Traction,fb] = TensTorsTB(XY(1),XY(2),Z,mateprop(4),mateprop(5),lamda*Twist2/L2,0,[0 0 1]',eye(3));
                            Z = Z - L1;
                            Q = [cos(pi*lamda*Twist) -sin(pi*lamda*Twist)
                                 sin(pi*lamda*Twist)  cos(pi*lamda*Twist)];
                            Twist2 = MateT(1,4)/MateT(2,4)*Twist;
                            [Traction,fb] = TensTorsTB(X,Y,Z,mateprop(4),mateprop(5),lamda*Twist2/L2,0,[0 0 1]',[Q zeros(2,1); 0 0 1]);
                        end
                    else
                        [Traction,fb] = TensTorsTB(X,Y,Z,mateprop(4),mateprop(5),lamda*Twist/((L1+L2)*(1+lamda*Axialdelta)),lamda*Axialdelta,[0 0 1]',eye(3));
                    end
%                     [Traction,fb] = TensTorsTB2(X,Y,Z,mateprop(4),mateprop(5),lamda*Twist,[0 0 1]');
                else
                    fb = zeros(3,1);
                end
            else
                fb = zeros(3,1);
            end
            
            ElemF = ElemF + c1*Nmat'*fb;

        end %je
    ElemF;

    case 11 % Error evaluation
        
        ElemE = zeros(numEn,1);

        lam = getlam(mateprop);
        
        % Load Guass Integration Points

        if nel == 4
            lint = 16;11;5;16;
        elseif nel == 8
            lint = 27;8;1000; %1000 for body force problem 
        elseif nel == 10
            lint = 14;
%             lint = 27;
        else
            lint = 27;
        end
        der = 0;
        bf = 0;
        ib = 0;
        
        el2el = zeros(3,1);
        eprixel = zeros(3,1);
        epriyel = zeros(3,1);
        eprizel = zeros(3,1);
        ue = zeros(3,1);
        duex = ue;
        duey = ue;

        %Integration Loop
        for ll = 1:lint

                %Evaluate first derivatives of basis functions at int. point
                if nel == 4 || nel == 10
                  [Wgt,ss] =  int3d_t(ll,lint,ib);
                  [shl,shld,shls,be] = shltt(ss,nel,nel,der,bf);
                  [Qxy, shgs, Jdet] = shgtt(xl,nel,shld,shls,nel,bf,der,be);
                else
                  [Wgt,ss] =  intpntb(ll,lint,ib);
                  [shl,shld,shls,be] = shlb(ss,nel,nel,der,bf);
                  [Qxy, shgs, Jdet] = shgb(xl,nel,shld,shls,nel,bf,der,be);
                end
                
                xint = xl(1,1:nel)*shl;
                yint = xl(2,1:nel)*shl;
                zint = xl(3,1:nel)*shl;
                dux = ul(1:3,1:nel)*Qxy(:,1);
                duy = ul(1:3,1:nel)*Qxy(:,2);
                duz = ul(1:3,1:nel)*Qxy(:,3);
                u = ul(1:3,1:nel)*shl;
                
            [F,JxX,fi] = kine3d(Qxy,ul(:,1:nel),nel,0);
                
            c1 = Wgt*Jdet;
            
            if iprob == 7
                if exist('bimat','var') && bimat == 1 
                    if ma == 1
                        [ue,due] = TensTorsU(xint,yint,zint,mateprop(5),lamda*Twist/L1,0);
                    else
                        zint = zint - L1;
                        Q = [cos(pi*lamda*Twist) -sin(pi*lamda*Twist)
                             sin(pi*lamda*Twist)  cos(pi*lamda*Twist)];
                        Twist2 = MateT(1,4)/MateT(2,4)*Twist;
                        [ue,due] = TensTorsU(xint,yint,zint,mateprop(5),lamda*Twist2/L2,0);
                        Q2 = due(1:2,1:2)+eye(2); % rotation from interface to current cross-section
                        XY = [xint; yint];
                        xyQ = Q2*Q*XY; % get rotated current coordinates
                        ue(1:2) = xyQ - [xint; yint]; % compute displacement wrt undeformed coordinates
                        due(1:2,1:2) = Q2*Q-eye(2);
                        due(1:2,3) = Q*due(1:2,3);
                    end
                else
                    [ue,due] = TensTorsU(xint,yint,zint,mateprop(5),lamda*Twist/((L1+L2)*(1+lamda*Axialdelta)),lamda*Axialdelta);
                end
%                 [ue,due] = TensTorsU2(xint,yint,zint,lamda*Twist);
                duex = due(:,1);
                duey = due(:,2);
                duez = due(:,3);
            elseif iprob == 9
                [ue,due] = BlockJ1U(xint,yint,zint);
                duex = due(:,1);
                duey = due(:,2);
                duez = due(:,3);
            else
                ue = zeros(3,1);
                duex = ue;
                duey = ue;
                duez = ue;
            end
            
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
            
        end %je

        for in= 1:3
            ElemE(in) = el2el(in);
            ElemE(in+3) = eprixel(in);
            ElemE(in+6) = epriyel(in);
            ElemE(in+9) = eprizel(in);
        end
        
    case 22
        
        ElemM = zeros(nst);
        ElemF = zeros(nst,1);

        lam = getlam(mateprop);
        
        I1 = [1; 1; 1; 0; 0; 0; 0; 0; 0];
        spvec0 = I1;
        
        % Load Guass Integration Points

        if nel == 4
            lint = 11; %minimum of 7 for all integrals in deformed state
            der = 0;
        elseif nel == 8
%             lint = 4;
            lint = 8;
            der = 1;
        elseif nel == 10
            lint = 14; %minimum of 13 for all integrals in deformed state
            der = 1;
        else
%             lint = 9;
            lint = 27;
            der = 0;
        end

        %Integration Loop
        for ll = 1:lint

            %Evaluate first derivatives of basis functions at int. point
            if nel == 4 || nel == 10
              [Wgt,ss] =  int3d_t(ll,lint,ib);
              [shl,shld,shls,be] = shltt(ss,nel,nen,der,bf);
              [Qxy, shgs, Jdet] = shgtt(xl+ul,nel,shld,shls,nel,bf,der,be);
            else
              [Wgt,ss] =  intpntb(ll,lint,ib);
              [shl,shld,shls,be] = shlb(ss,nel,nen,der,bf);
              [Qxy, shgs, Jdet] = shgb(xl+ul,nel,shld,shls,nel,bf,der,be);
            end

            if nelP == 4 || nelP == 10
              [shlS,shld,shls] = shltt(ss,nelS,nel,0,0);
            else
              [shlS,shld,shls] = shlb(ss,nelS,nel,0,0);
            end
            
            [fi,JxX,F] = kine3d(Qxy,-ul,nel,0); %this is equivalent to ikine2d
            JxX = 1/JxX; %this is equivalent to ikine2d
    %         [F,JxX,fi] = ikine2d(Qxy,ul,nel); %this is equivalent to above two lines
            Jdet = Jdet/JxX;

            sigma = SigmaCmat3(F,JxX,mateprop,lam);
            sigma = sigma/JxX;
            
            % Form B matrix
            Nmat = shlS';

            w = Wgt*Jdet;
            
            if stres <= 6 % stress components
                sigmas = sigma(stres);
            elseif stres >= 8
                if stres <= 10 % principal stresses
                sigma = [sigma(1) sigma(4) sigma(6); sigma(4) sigma(2) sigma(5); sigma(6) sigma(5) sigma(3)];
                psig = eig(sigma);
                sigmas = psig(stres-7);
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
        spvec0 = I1;
        
        thick = 1;
        
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
              [Qxy,shgs,Jdet,be] = shgtt(xl+ul(1:3,:),nel,shld,shls,nel,bf,der,be);
            else
              [Wgt,ss] =  intpntb(ll,nint,ib);
              [shl,shld,shls,be] = shlb(ss,nel,nel,der,bf);
              [Qxy,shgs,Jdet,be] = shgb(xl+ul(1:3,:),nel,shld,shls,nel,bf,der,be);
            end
            
            if nelP == 4 || nelP == 10
              shlp = shltt(ss,nelP,nel,0,0);
            else
              shlp = shlb(ss,nelP,nel,0,0);
            end
            
            [fi,JxX,F] = kine3d(Qxy,-ul,nel,0); %this is equivalent to ikine2d
            JxX = 1/JxX; %this is equivalent to ikine2d
    %         [F,JxX,fi] = ikine2d(Qxy,ul,nel); %this is equivalent to above two lines
            Jdet = Jdet/JxX;
            
            sigma = SigmaCmat3(F,JxX,mateprop,lam);
            sigma = sigma/JxX;
            
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
              [Qxy,shgs,Jdet,be] = shgtt(xl+ul(1:3,:),nel,shld,shls,nel,bf,der,be);
            else
              [Wgt,ss] =  intpntb(ll,lint,ib);
              [shl,shld,shls,be] = shlb(ss,nel,nel,der,bf);
              [Qxy,shgs,Jdet,be] = shgb(xl+ul(1:3,:),nel,shld,shls,nel,bf,der,be);
            end
            
            [fi,JxX,F] = kine3d(Qxy,-ul,nel,0); %this is equivalent to ikine2d
            JxX = 1/JxX; %this is equivalent to ikine2d
    %         [F,JxX,fi] = ikine2d(Qxy,ul,nel); %this is equivalent to above two lines
            Jdet = Jdet/JxX;

            w = Wgt*Jdet*thick;
            
            Vol = Vol + w;

        end %je

        for i = 1:nel
        ElemS(i,npstr+1) = 1; % use this for simple average Vol; % use this for weighted average 
        end
%%        
    case 26 % Element Stress

        ElemS = zeros(1,nestr);

        lam = getlam(mateprop);
        
        I1 = [1; 1; 1; 0; 0; 0; 0; 0; 0];
        spvec0 = I1;
        
        % Load Guass Integration Points

        
        der = 0;
        bf = 0;
        ib = 0;


            %Evaluate first derivatives of basis functions at int. point
            if nel == 4 || nel == 10
              ss = .25*[1 1 1];
              [shl,shld,shls,be] = shltt(ss,nel,nel,der,bf);
              [Qxy,shgs,Jdet,be] = shgtt(xl+ul(1:3,:),nel,shld,shls,nel,bf,der,be);
            else
              ss = [0 0 0];
              [shl,shld,shls,be] = shlb(ss,nel,nel,der,bf);
              [Qxy,shgs,Jdet,be] = shgb(xl+ul(1:3,:),nel,shld,shls,nel,bf,der,be);
            end
            
            if nelP == 4 || nelP == 10
              shlp = shltt(ss,nelP,nel,0,0);
            else
              shlp = shlb(ss,nelP,nel,0,0);
            end
            
            [fi,JxX,F] = kine3d(Qxy,-ul,nel,0); %this is equivalent to ikine2d
            JxX = 1/JxX; %this is equivalent to ikine2d
    %         [F,JxX,fi] = ikine2d(Qxy,ul,nel); %this is equivalent to above two lines
            Jdet = Jdet/JxX;
            
            sigma = SigmaCmat3(F,JxX,mateprop,lam);
            sigma = sigma/JxX;
            
            for stres = 1:nestr
            
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
            
            ElemS(1,stres) = sigmas;
            
            end
        
end