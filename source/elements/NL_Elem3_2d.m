% Linear dynamics 2D element
% Updated and verified 1/30/2014 for a-form and d-form with mass and
% stiffness combined in one routine/call; alpha can also be nonzero

% Set Material Properties

PatchE = mateprop(2);
Patchv = mateprop(3);
rho0 = mateprop(4);

switch isw %Task Switch

    case 3 %Compute Stiffness and Residual
        
        ElemK = zeros(nst);
        ElemF = zeros(nst,1);
        ElemM = zeros(nst);
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
                  [Qxy, shgs, Jdet] = shgt(xl,nel,shld,shls,nel,bf,der,be);
                else
                  [Wgt,litr,lits] =  intpntq(ll,lint,ib);
                  [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
                  [Qxy, shgs, Jdet] = shgq(xl,nel,shld,shls,nel,bf,der,be);
                end
                
                %Compute rho, accel
                rho = rho0; %rho = f(x,y);
                accel = al*shl;
                
                %Compute F, sigma, D
                [F] = kine2d(Qxy,ul,nel,1);
                F = F - eye(2);
%                 sigma = ComSigma10(JxX,F,mu,Kmat);
%                 Dmat = ComCijkl10(JxX,F,mu,Kmat);
                sigma = ComSigma570(F,mu,lam);
                %Compute F, sigma, D
                [Fn] = kine2d(Qxy,ul_n,nel,1);
                Fn = Fn - eye(2);
%                 sigma = ComSigma10(JxX,F,mu,Kmat);
%                 Dmat = ComCijkl10(JxX,F,mu,Kmat);
                sigman = ComSigma570(Fn,mu,lam);
                
                c1 = Wgt*Jdet*thick;
                c2 = Wgt*Jdet*thick*rho;
                %Evaluate Stiffness terms
                for o=1:nel
                    dox = Qxy(o,1)*c1;
                    doy = Qxy(o,2)*c1;
                    doN = shl(o)*c2;
                    don = shl(o)*c1;
                    BTsigma = [dox*sigma(1) + doy*sigma(3); dox*sigma(3) + doy*sigma(2)];
                    BTsigman = [dox*sigman(1) + doy*sigman(3); dox*sigman(3) + doy*sigman(2)];

                    ElemF(ndf*o-1) = ElemF(ndf*o-1) - doN*accel(1) - coeffkl2*BTsigma(1) + Nalpha*BTsigman(1);
                    ElemF(ndf*o) =   ElemF(ndf*o)   - doN*accel(2) - coeffkl2*BTsigma(2) + Nalpha*BTsigman(2);

                    for l=1:nel
                        dlx = Qxy(l,1);
                        dly = Qxy(l,2);
                        %(wx,ux)
                        ElemK(ndf*o-1,ndf*l-1) = ElemK(ndf*o-1,ndf*l-1) + ...
                                                 dox*dlx*lam + 2*dox*mu*dlx + ...
                                                 doy*mu*dly;
                        %(wx,uy)
                        ElemK(ndf*o-1,ndf*l)   = ElemK(ndf*o-1,ndf*l) + ...
                                                 dox*lam*dly + doy*mu*dlx;
                        %(wy,ux)
                        ElemK(ndf*o,ndf*l-1)   = ElemK(ndf*o,ndf*l-1) + ...
                                                 doy*lam*dlx + dox*mu*dly;
                        %(wy,uy)
                        ElemK(ndf*o,ndf*l)     = ElemK(ndf*o,ndf*l) + ...
                                                 doy*dly*lam + 2*doy*mu*dly + ...
                                                 dox*mu*dlx;
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
        
        if lumping == 1
            ElemM = diag(sum(ElemM));
            ElemF = -ElemM*reshape(al,ndf*nel,1) -ElemK*reshape(ul,ndf*nel,1);
        end
        
        ElemK = coeffm*ElemM + coeffk*ElemK;
        
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
%                 rho = rho0*(1-1/30*x); %rho = f(x,y);
                rho = rho0;
                
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
        
        if lumping == 1
            ElemM = diag(sum(ElemM));
        end
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
                  [Qxy, shgs, Jdet] = shgt(xl,nel,shld,shls,nel,bf,der,be);
                else
                  [Wgt,litr,lits] =  intpntq(ll,lint,ib);
                  [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
                  [Qxy, shgs, Jdet] = shgq(xl,nel,shld,shls,nel,bf,der,be);
                end
                
                %Compute rho, accel
                rho = rho0; %rho = f(x,y);
                accel = al*shl;
                
                %Compute F, sigma, D
                [F] = kine2d(Qxy,ul,nel,1);
                F = F - eye(2);
%                 sigma = ComSigma10(JxX,F,mu,Kmat);
%                 Dmat = ComCijkl10(JxX,F,mu,Kmat);
                sigma = ComSigma570(F,mu,lam);
                
                c1 = Wgt*Jdet*thick;
                c2 = Wgt*Jdet*thick*rho;
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
        ElemM = zeros(nst);

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
                  [Qxy, shgs, Jdet] = shgt(xl,nel,shld,shls,nel,bf,der,be);
                else
                  [Wgt,litr,lits] =  intpntq(ll,lint,ib);
                  [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
                  [Qxy, shgs, Jdet] = shgq(xl,nel,shld,shls,nel,bf,der,be);
                end
                
                c1 = Wgt*Jdet*thick;
                
                %Evaluate Stiffness terms
                for o=1:nel
                    dox = Qxy(o,1)*c1;
                    doy = Qxy(o,2)*c1;
                    don = shl(o)*c1;

                    for l=1:nel
                        dlx = Qxy(l,1);
                        dly = Qxy(l,2);
                        %(wx,ux)
                        ElemK(ndf*o-1,ndf*l-1) = ElemK(ndf*o-1,ndf*l-1) + ...
                                                 dox*dlx*lam + 2*dox*mu*dlx + ...
                                                 doy*mu*dly;
                        %(wx,uy)
                        ElemK(ndf*o-1,ndf*l)   = ElemK(ndf*o-1,ndf*l) + ...
                                                 dox*lam*dly + doy*mu*dlx;
                        %(wy,ux)
                        ElemK(ndf*o,ndf*l-1)   = ElemK(ndf*o,ndf*l-1) + ...
                                                 doy*lam*dlx + dox*mu*dly;
                        %(wy,uy)
                        ElemK(ndf*o,ndf*l)     = ElemK(ndf*o,ndf*l) + ...
                                                 doy*dly*lam + 2*doy*mu*dly + ...
                                                 dox*mu*dlx;
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
                  [Qxy, shgs, Jdet] = shgt(xl,nel,shld,shls,nel,bf,der,be);
                else
                  [Wgt,litr,lits] =  intpntq(ll,lint,ib);
                  [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
                  [Qxy, shgs, Jdet] = shgq(xl,nel,shld,shls,nel,bf,der,be);
                end
                
                %Compute strain and stress
                E = zeros(3,1);
                for o = 1:nel

                    E(1) = E(1) + Qxy(o,1)*ul(1,o);
                    E(2) = E(2) + Qxy(o,2)*ul(2,o);
                    E(3) = E(3) + 1/2*(Qxy(o,1)*ul(2,o) + Qxy(o,2)*ul(1,o));

                end
                Ekk = E(1) + E(2);
                sigma = lam*Ekk*[1 1 0]' + 2*mu*E;
                
                ElemP(6*(ll-1)+1) = E(1);
                ElemP(6*(ll-1)+2) = E(2);
                ElemP(6*(ll-1)+3) = E(3);
                ElemP(6*(ll-1)+4) = sigma(1);
                ElemP(6*(ll-1)+5) = sigma(2);
                ElemP(6*(ll-1)+6) = sigma(3);

        end %je
        
    case 12 % energy
        
        ElemE = 0;
        Nmat = zeros(2,nst);
        Bmat = zeros(3,nst);

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
                D = [lam+2*mu lam 0 
                     lam lam+2*mu 0
                     0 0  mu];

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
                
                for i = 1:nel
                    Bmat(:,2*i-1:2*i) = [Qxy(i,1) 0
                                         0 Qxy(i,2)
                                         Qxy(i,2) Qxy(i,1)];
                end
                %Compute F, sigma, D
                E = Bmat*reshape(ul,nst,1);
                velo = vl*shl;
                
                c1 = Wgt*Jdet*thick;
                ElemE = ElemE + c1/2*E'*D*E + c1/2*velo'*rho*velo;

        end %je

end