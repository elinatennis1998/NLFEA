% function ElemK = NL_Elem2_2d(ElemXYZ, nel, ndf, R, S, T, pr, ps, PatchE, Patchv)
%
% Tim Truster
% 7/2009
% UIUC
%
%   nel:             = number of nodes on current element
%
%   ndf:             = max number of DOF per node
%
%   ndm:             = space dimension of mesh
%
%   nst:             = size of element arrays (normally ndf*nel)
%
%
%   wf:              = weight factor; adds additional dimension to
%                      coordinate arrays (0 = Lagrangian elements, 
%                      1 = NURBS elements)
%
%   ElemXYZ:          = Local array containing (x,y) coordinates of nodes
%                      forming the element passed to the ComElemStif
%                      function
%
%                       X       Y
%                       1.54    -2.367
%
%   pr,ps,pt=0:      = degree of polynomial for elements in
%                      uvw-directions
%
%   mr,ms,mt:        = number of konts on element in uvw-directions
%
%   nr,ns,nt:        = number of nodes on element in uvw-directions
%
%   R,S,T:           = Element knot vectors
%

% Set Material Properties

PatchE = mateprop(2);
Patchv = mateprop(3);
rho0 = mateprop(4);

switch isw %Task Switch

    case 3 %Compute Stiffness and Residual
        
        ElemK = zeros(nst);
        ElemF = zeros(nst,1);
        ElemFn = zeros(nst,1);

        lam = Patchv*PatchE/((1+Patchv)*(1-2*Patchv));
        mu = PatchE/(2*(1+Patchv));
        
        ul12 = 1/2*(ul + ul_n);
        one = [1; 1; 0; 0];
        mat1 = one*one';
        matE = diag([2,2,1,0]);
        Cijkl = zeros(2,2,2,2);
        Cijkl(1,1,1,1) = lam + 2*mu;
        Cijkl(2,2,2,2) = lam + 2*mu;
        Cijkl(1,1,2,2) = lam;
        Cijkl(2,2,1,1) = lam;
        Cijkl(1,2,1,2) = mu;
        Cijkl(1,2,2,1) = mu;
        Cijkl(2,1,1,2) = mu;
        Cijkl(2,1,2,1) = mu;
        
        thick = 1;
        
        % Load Guass Integration Points

        if nel == 3 || nel == 6
            lint = 13;
        elseif nel == 4
            lint = 4;
            % lint = 16;
        else
            lint = 9;
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
                  [QXY, shgs, Jdet0] = shgt(xl,nel,shld,shls,nel,bf,der,be);
                else
                  [Wgt,litr,lits] =  intpntq(ll,lint,ib);
                  [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
                  [QXY, shgs, Jdet0] = shgq(xl,nel,shld,shls,nel,bf,der,be);
                end
                
                %Compute rho, accel
                x = xl(1,:)*shl;
                y = xl(2,:)*shl;
                rho = rho0;%*(1-0.1*x); %rho = f(x,y);
                
                if transient == 4
                    
                accel = 2/tstep^2*(ul-ul_n-tstep*vl_n)*shl;
                
                Fn1 = ul(1:2,:)*QXY+eye(2);
                Fn = ul_n(1:2,:)*QXY+eye(2);
                [F,JxX,fi,Qxy] = kine2d(QXY,ul12,nel,1);
                F12 = 1/2*(Fn1);% + Fn);
                Jdet = Jdet0;%*JxX;
            
                Cn = Fn'*Fn;
                Cn1 = Fn1'*Fn1;
                Ch = 1/2*(Cn1 + Cn); %%%% VERY IMPORTANT: not C(F_n+1/2)
                Ev = 1/2*[Ch(1,1)-1; Ch(2,2)-1; 2*Ch(1,2); 0];
                trE = one'*Ev;
                Sh = lam*trE*one + mu*matE*Ev;
                Dh = lam*mat1 + mu*matE;
                
                P = [F(1,1)^2 F(1,2)^2 2*F(1,1)*F(1,2) 0
                     F(2,1)^2 F(2,2)^2 2*F(2,1)*F(2,2) 0
                     F(1,1)*F(2,1) F(1,2)*F(2,2) F(1,1)*F(2,2)+F(1,2)*F(2,1) 0
                    zeros(1,4)];
                sigma = P*Sh;
                Dmat = P*Dh*P'; % This doesn't match the true linearization
                % NOTE: Simo's tangent is incorrect; you need beta factors
                % on the geometric and material terms; this is where the
                % 1/2 factors shown below and on F12 come from. I found it
                % by trial and error, and one day I'll write it on paper,
                % but it makes sense, because it's S(1/2C_n+1 + 1/2*C_n).
                
                % If I teach this in a class, I should not expect them to
                % be able to figure out or implement the exact consistent
                % tangent; instead, get them to be happy with the rhs, and
                % if they desire to know one-on-one then I can show them.
                
                c1 = Wgt*Jdet*thick;
                c2 = Wgt*Jdet0*thick*rho;
                %Evaluate Stiffness terms
                for o=1:nel
                    dox = Qxy(o,1)*c1;
                    doy = Qxy(o,2)*c1;
                    doN = shl(o)*c2;
                    BTsigma = [dox*sigma(1) + doy*sigma(3); dox*sigma(3) + doy*sigma(2)];

                    ElemF(ndf*o-1) = ElemF(ndf*o-1) - doN*accel(1) - coeff1*BTsigma(1);
                    ElemF(ndf*o) =   ElemF(ndf*o)   - doN*accel(2) - coeff1*BTsigma(2);
                    
                    ElemFn(ndf*o-1) = ElemFn(ndf*o-1) - Nalpha*BTsigma(1);
                    ElemFn(ndf*o) =   ElemFn(ndf*o)   - Nalpha*BTsigma(2);

%                     for l=1:nel
                    for l=1:nel
                        dlx = Qxy(l,1);
                        dly = Qxy(l,2);
                        BTsigmaB = BTsigma(1)*dlx + BTsigma(2)*dly;
                        %(wx,ux)
                        ElemK(ndf*o-1,ndf*l-1) = ElemK(ndf*o-1,ndf*l-1) + 1/2*BTsigmaB;
                        %(wy,uy)
                        ElemK(ndf*o,ndf*l)     = ElemK(ndf*o,ndf*l) + 1/2*BTsigmaB;
                    end %l
                    do = c1*[QXY(o,1); QXY(o,2)];
                    for q=1:nel
                        dq = [QXY(q,1); QXY(q,2)];
                        for i = 1:2
                            for k = 1:2
                                term = 0;
                                for I = 1:2
                                    for J = 1:2
                                        for K = 1:2
                                            for L = 1:2
                                                term = term + F(i,I)*do(J)*Cijkl(I,J,K,L)*F12(k,K)*dq(L);
                                            end
                                        end
                                    end
                                end
                                ElemK(ndf*(o-1)+i,ndf*(q-1)+k) = ElemK(ndf*(o-1)+i,ndf*(q-1)+k) + term;
                            end
                        end
                    end %l
                end %o
                
                else

                accel = al*shl;
                
                [F,JxX,fi,Qxy] = kine2d(QXY,ul,nel,1);
                Jdet = Jdet0;%*JxX;
                [Sh,Dh] = SigmaCmat8(F,mu,lam);
                P = [F(1,1)^2 F(1,2)^2 2*F(1,1)*F(1,2) 0
                     F(2,1)^2 F(2,2)^2 2*F(2,1)*F(2,2) 0
                     F(1,1)*F(2,1) F(1,2)*F(2,2) F(1,1)*F(2,2)+F(1,2)*F(2,1) 0
                    zeros(1,4)];
                sigma = P*Sh;
                Dmat = P*Dh*P';
                
                c1 = Wgt*Jdet*thick;
                c2 = Wgt*Jdet0*thick*rho;
                %Evaluate Stiffness terms
                for o=1:nel
                    dox = Qxy(o,1)*c1;
                    doy = Qxy(o,2)*c1;
                    doN = shl(o)*c2;
                    BTsigma = [dox*sigma(1) + doy*sigma(3); dox*sigma(3) + doy*sigma(2)];

                    ElemF(ndf*o-1) = ElemF(ndf*o-1) - doN*accel(1) - coeff1*BTsigma(1);
                    ElemF(ndf*o) =   ElemF(ndf*o)   - doN*accel(2) - coeff1*BTsigma(2);
                    
                    ElemFn(ndf*o-1) = ElemFn(ndf*o-1) - Nalpha*BTsigma(1);
                    ElemFn(ndf*o) =   ElemFn(ndf*o)   - Nalpha*BTsigma(2);

%                     for l=1:nel
                    for l=1:nel
                        dlx = Qxy(l,1);
                        dly = Qxy(l,2);
                        BTsigmaB = BTsigma(1)*dlx + BTsigma(2)*dly;
                        %(wx,ux)
                        ElemK(ndf*o-1,ndf*l-1) = ElemK(ndf*o-1,ndf*l-1) + ...
                                                 (dox*Dmat(1,1)+doy*Dmat(3,1))*dlx + ...
                                                 (dox*Dmat(1,3)+doy*Dmat(3,3))*dly + ...
                                                 BTsigmaB;
                        %(wx,uy)
                        ElemK(ndf*o-1,ndf*l)   = ElemK(ndf*o-1,ndf*l) + ...
                                                 (dox*Dmat(1,2)+doy*Dmat(3,2))*dly + ...
                                                 (dox*Dmat(1,3)+doy*Dmat(3,3))*dlx;
%                         %(wy,ux)
                        ElemK(ndf*o,ndf*l-1)   = ElemK(ndf*o,ndf*l-1) + ...
                                                 (doy*Dmat(2,1)+dox*Dmat(3,1))*dlx + ...
                                                 (doy*Dmat(2,3)+dox*Dmat(3,3))*dly;
                        %(wy,uy)
                        ElemK(ndf*o,ndf*l)     = ElemK(ndf*o,ndf*l) + ...
                                                 (doy*Dmat(2,2)+dox*Dmat(3,2))*dly + ...
                                                 (doy*Dmat(2,3)+dox*Dmat(3,3))*dlx + ...
                                                 BTsigmaB;
                    end %l
                end %o
                    
                end

        end %je
        
    case 5 %Compute Mass
        
        ElemM = zeros(nst);
        
        thick = 1;
        
        % Load Guass Integration Points

        if nel == 3 || nel == 6
            lint = 13;
        elseif nel == 4
            lint = 4;
            % lint = 16;
        else
            lint = 9;
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
                
                x = xl(1,:)*shl;
                y = xl(2,:)*shl;
                rho = rho0;%*(1-0.1*x); %rho = f(x,y);
%                 rho = rho0;
                
                c1 = Wgt*Jdet*thick;
                for o=1:nel
                    don = shl(o)*c1;

%                     for l=1:nel
                    for l=o:nel
                        dln = shl(l);
                        Mmat = don*rho*dln;
                        %(wx,ux)
                        ElemM(ndf*o-1,ndf*l-1) = ElemM(ndf*o-1,ndf*l-1) + ...
                                                 Mmat;
                        %(wy,uy)
                        ElemM(ndf*o,ndf*l)     = ElemM(ndf*o,ndf*l) + ...
                                                 Mmat;
                    end %l
                end %o

        end %je
        
        %Copy symmetric values
        
        for o=1:nel

            for l=o+1:nel
                
                Mmat = ElemM(ndf*o-1,ndf*l-1);
                
                %(wx,ux)
                ElemM(ndf*l-1,ndf*o-1) = Mmat;
                %(wy,uy)
                ElemM(ndf*l,ndf*o)     = Mmat;
                
            end %l
        end %o
% ElemM
    case 6 %Compute Residual
        
        ElemF = zeros(nst,1);
        ElemFn = zeros(nst,1);

        lam = Patchv*PatchE/((1+Patchv)*(1-2*Patchv));
        mu = PatchE/(2*(1+Patchv));
        
        thick = 1;
        
        % Load Guass Integration Points

        if nel == 3 || nel == 6
            lint = 13;
        elseif nel == 4
            lint = 4;
            % lint = 16;
        else
            lint = 9;
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
                  [QXY, shgs, Jdet0] = shgt(xl,nel,shld,shls,nel,bf,der,be);
                else
                  [Wgt,litr,lits] =  intpntq(ll,lint,ib);
                  [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
                  [QXY, shgs, Jdet0] = shgq(xl,nel,shld,shls,nel,bf,der,be);
                end
                
                %Compute rho, accel
                x = xl(1,:)*shl;
                y = xl(2,:)*shl;
                rho = rho0;%*(1-0.1*x); %rho = f(x,y);
                accel = al*shl;
                
                %Compute D
%                 F = zeros(2,2);
%                 for ii = 1:2
%                     for II = 1:2
%                         for kk = 1:nel
%                             F(ii,II) = F(ii,II) + Elemxyz(ii,kk)*QXY(II,kk);
%                         end
%                     end
%                 end
% 
%                 JxX = det(F);
                [F,JxX,fi,Qxy] = kine2d(QXY,ul,nel,1);
                Jdet = Jdet0*JxX;
                sigma = ComSigma09(JxX,F,mu,lam);
                Dmat = ComCijkl09(JxX,F,mu,lam);
                
                c1 = Wgt*Jdet*thick;
                c2 = Wgt*Jdet0*thick*rho;
                %Evaluate Stiffness terms
                for o=1:nel
                    dox = Qxy(o,1)*c1;
                    doy = Qxy(o,2)*c1;
                    doN = shl(o)*c2;
                    BTsigma = [dox*sigma(1) + doy*sigma(3); dox*sigma(3) + doy*sigma(2)];

                    ElemF(ndf*o-1) = ElemF(ndf*o-1) - doN*accel(1) - coeff1*BTsigma(1);
                    ElemF(ndf*o) =   ElemF(ndf*o)   - doN*accel(2) - coeff1*BTsigma(2);
                    
                    ElemFn(ndf*o-1) = ElemFn(ndf*o-1) - Nalpha*BTsigma(1);
                    ElemFn(ndf*o) =   ElemFn(ndf*o)   - Nalpha*BTsigma(2);

                end %o

        end %je

    case 21 %Compute Stiffness
        
        ElemK = zeros(nst);

        lam = Patchv*PatchE/((1+Patchv)*(1-2*Patchv));
        mu = PatchE/(2*(1+Patchv));
        
        thick = 1;
        
        % Load Guass Integration Points

        if nel == 3 || nel == 6
            lint = 13;
        elseif nel == 4
            lint = 4;
            % lint = 16;
        else
            lint = 9;
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
                  [QXY, shgs, Jdet0] = shgt(xl,nel,shld,shls,nel,bf,der,be);
                else
                  [Wgt,litr,lits] =  intpntq(ll,lint,ib);
                  [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
                  [QXY, shgs, Jdet0] = shgq(xl,nel,shld,shls,nel,bf,der,be);
                end
                
                %Compute D
%                 F = zeros(2,2);
%                 for ii = 1:2
%                     for II = 1:2
%                         for kk = 1:nel
%                             F(ii,II) = F(ii,II) + Elemxyz(ii,kk)*QXY(II,kk);
%                         end
%                     end
%                 end
% 
%                 JxX = det(F);
                [F,JxX,fi,Qxy] = kine2d(QXY,ul,nel,1);
                Jdet = Jdet0*JxX;
                sigma = ComSigma09(JxX,F,mu,lam);
                Dmat = ComCijkl09(JxX,F,mu,lam);
                
                c1 = Wgt*Jdet*thick;
                %Evaluate Stiffness terms
                for o=1:nel
                    dox = Qxy(o,1)*c1;
                    doy = Qxy(o,2)*c1;

%                     for l=1:nel
                    for l=o:nel
                        dlx = Qxy(l,1);
                        dly = Qxy(l,2);
                        BTsigmaB = BTsigma(1)*dlx + BTsigma(2)*dly;
                        %(wx,ux)
                        ElemK(ndf*o-1,ndf*l-1) = ElemK(ndf*o-1,ndf*l-1) + ...
                                                 (dox*Dmat(1,1)+doy*Dmat(3,1))*dlx + ...
                                                 (dox*Dmat(1,3)+doy*Dmat(3,3))*dly + ...
                                                 BTsigmaB;
                        %(wx,uy)
                        ElemK(ndf*o-1,ndf*l)   = ElemK(ndf*o-1,ndf*l) + ...
                                                 (dox*Dmat(1,2)+doy*Dmat(3,2))*dly + ...
                                                 (dox*Dmat(1,3)+doy*Dmat(3,3))*dlx;
%                         %(wy,ux)
                        ElemK(ndf*o,ndf*l-1)   = ElemK(ndf*o,ndf*l-1) + ...
                                                 (doy*Dmat(2,1)+dox*Dmat(3,1))*dlx + ...
                                                 (doy*Dmat(2,3)+dox*Dmat(3,3))*dly;
                        %(wy,uy)
                        ElemK(ndf*o,ndf*l)     = ElemK(ndf*o,ndf*l) + ...
                                                 (doy*Dmat(2,2)+dox*Dmat(3,2))*dly + ...
                                                 (doy*Dmat(2,3)+dox*Dmat(3,3))*dlx + ...
                                                 BTsigmaB;
                    end %l
                end %o

        end %je
        
        %Copy symmetric values
        
        for o=1:nst

            for l=o+1:nst
                
                ElemK(l,o) = ElemK(o,l);
                
            end %l
        end %o

    case 22
        
        ElemP = zeros(2*(3*ndf-3),nel);

        lam = Patchv*PatchE/((1+Patchv)*(1-2*Patchv));
        mu = PatchE/(2*(1+Patchv));
        
        thick = 1;
        
        % Load Guass Integration Points

        if nel == 3 || nel == 6
            lint = 13;
        elseif nel == 4
            lint = 4;
            % lint = 16;
        else
            lint = 9;
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
                  [QXY, shgs, Jdet] = shgt(xl,nel,shld,shls,nel,bf,der,be);
                else
                  [Wgt,litr,lits] =  intpntq(ll,lint,ib);
                  [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
                  [QXY, shgs, Jdet] = shgq(xl,nel,shld,shls,nel,bf,der,be);
                end
                
                %Compute D
%                 F = zeros(2,2);
%                 for ii = 1:2
%                     for II = 1:2
%                         for kk = 1:nel
%                             F(ii,II) = F(ii,II) + Elemxyz(ii,kk)*QXY(II,kk);
%                         end
%                     end
%                 end
% 
%                 JxX = det(F);
                [F,JxX,fi,Qxy] = kine2d(QXY,ul,nel,1);
                Jdet = Jdet*JxX;
                Estrain = 1/2*(F'*F - eye(ndm));
                sigma = ComSigma09(JxX,F,mu,lam);
                
                ElemP(6*(ll-1)+1) = Estrain(1,1);
                ElemP(6*(ll-1)+2) = Estrain(2,2);
                ElemP(6*(ll-1)+3) = Estrain(1,2);
                ElemP(6*(ll-1)+4) = sigma(1);
                ElemP(6*(ll-1)+5) = sigma(2);
                ElemP(6*(ll-1)+6) = sigma(3);

        end %je
        
    case 12 % energy
        
        ElemE = 0;
        Nmat = zeros(2,nst);
        Bmat = zeros(3,nst);
        
        thick = 1;
        
        % Load Guass Integration Points

        if nel == 3 || nel == 6
            lint = 13;
        elseif nel == 4
            lint = 4;
            % lint = 16;
        else
            lint = 9;
        end
        der = 0;
        bf = 0;
        ib = 0;
        
        one = [1; 1; 0; 0];
        mat1 = one*one';
        matE = diag([2,2,1,0]);

        %Integration Loop
        for ll = 1:lint

                %Evaluate first derivatives of basis functions at int. point
                if nel == 3 || nel == 6
                  [Wgt,litr,lits] =  intpntt(ll,lint,ib);
                  [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
                  [QXY, shgs, Jdet0] = shgt(xl,nel,shld,shls,nel,bf,der,be);
                else
                  [Wgt,litr,lits] =  intpntq(ll,lint,ib);
                  [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
                  [QXY, shgs, Jdet0] = shgq(xl,nel,shld,shls,nel,bf,der,be);
                end
                
                F = ul(1:2,:)*QXY+eye(2);
            
                C = F'*F;
                Ev = 1/2*[C(1,1)-1; C(2,2)-1; 2*C(1,2); 0];
                Dmat = lam*mat1 + mu*matE;
                
                Ue = 1/2*(Ev'*Dmat*Ev);
                
                c1 = Wgt*Jdet0*thick;
                ElemE = ElemE + Ue*c1;

        end %je
        
end