
% Tim Truster
% 2/10/2014
% This version now gives the same results as in Harari's paper
% However, I do not agree philosophically with what they have done.
% They have integrated by parts the predictor term on the displacement
% field. This ruins the consistency of the method, because that term is
% supposed to vanish for linear elements. Namely, you cannot integrate by
% parts ( div(grad(w) , div(grad(u))
% In reality, they are only using the mass part of L* applied to f_tilde
% Thus, I think we can drop pursuing this further.

% Set Material Properties

cwave = mateprop(1);
t_on = 1;0;

Bcol1 = [1; 3];
Bcol2 = [2; 3];
col1 = [1; 2];
col2 = [2; 1];

switch isw %Task Switch
%%    
    case 1
        
        if ndf > 1
            
            for i = 2:ndf
                lie(i,1) = 0;
            end
            
        end
        
%         nh1 = nen*ndm*3;
        
%%
    case 3 %Compute Stiffness and Residual
        
        ElemK = zeros(ndf*nel);
        ElemF = zeros(ndf*nel,1);

        %Set integration number
        lint = nel;
        
        ib = 0;
        bf = 1;
        der = 1;
        
        fb = 0;
        Dmat = cwave^2;
        Nmat = zeros(1,ndf*nel);
        Bmat = zeros(1,ndf*nel);
        BBmat = zeros(1,nst);
        
        if t_on
            % tau from Harari
        tau = 1+6*Nbeta*CFL^2*(1-cosh(1/sqrt(Nbeta)/CFL))/(2+cosh(1/sqrt(Nbeta)/CFL));
        else
        tau = 0;
        end

        % Compute kinematic fields at intermediate time levels
        ul_n_af = (1 - Nalphaf)*ul + Nalphaf*ul_n;
        al_n_am = (1 - Nalpham)*al + Nalpham*al_n;
        ulres = reshape(ul_n_af,ndf*nel,1);
        alres = reshape(al_n_am,ndf*nel,1);
        
        ul_pred = ul_n + vl_n*tstep + tstep^2/2*(1-2*Nbeta)*al_n;
        ulpres = reshape(ul_pred,ndf*nel,1);
        
        sw = int1d(lint);
        [shl,shld,shls,be] = shl1d(sw(1,:),lint,nel-1);
        
            if(step==0.6*L/cwave) && elem == 1
                nind = 0;
            end
            
        for ll = 1:lint                    

            % Evaluate 1-D basis functions at integration points
            Wgt = sw(2,ll);
            [shg, shgs, Jdet, be] = shg1d(xl,ndm,nel,shld(ll,:),shls(ll,:),nen,bf,der,be);
            
            c1 = Wgt*Jdet;
            
            Nmat(1:ndf:ndf*(nel-1)+1) = shl(ll,:);
                
            Bmat(1:ndf:ndf*(nel-1)+1) = shg';
            
            BBmat(1:ndf:ndf*(nel-1)+1) = shgs';

            %Compute rho, accel
            accel = Nmat*alres;
            xint = Nmat*xl';
            if t_on
                sigma = Dmat*(Bmat*ulpres);
            else
                sigma = Dmat*(Bmat*ulres);
            end
            resid = -(Nmat*alres - Dmat*BBmat*ulres);
            
            Fvec = fb;
            
            if t_on
            % This version is eq (49) and (52) with integration by parts on
            % the predictor term. I do not agree that this stabilization is
            % consistent. Also, it leads to a different value of sigma that
            % is scaled by tau, so the external loads seen by the system
            % are completely different
            ElemF = ElemF + c1*(Nmat'*Fvec - (1-tau)*Bmat'*sigma);
            else
            % This version is the standard pure displacement plus consistent stabilization
            ElemF = ElemF + c1*(Nmat'*Fvec - Nmat'*accel - Bmat'*sigma ...
                                - (Nmat - (Nbeta*tstep^2)*Dmat*BBmat)'*tau*resid);
            end
            
            ElemK = ElemK + c1*(Nmat'*Nmat + Bmat'*(Nbeta*tstep^2)*Dmat*Bmat ...
                               - (Nmat - (Nbeta*tstep^2)*Dmat*BBmat)'*tau*(Nmat - (Nbeta*tstep^2)*Dmat*BBmat));
            
            if(step==stepmax)
                nind = nind + 1;
                Sig1(nind,:)=[xint/L (1-tau)*sigma/NodeLoadnp(3)];
                Sig1(nind,:)=[xint/L (1-tau)*Dmat*(Bmat*ulres)/NodeLoadnp(3)];
                
            end

        end %je
        
        % Put the acceleration part of the force vector on the RHS so that
        % NR converges (i.e. the corrector part)
        if t_on
        ElemF = ElemF - ElemK*alres;
        end
        
    case 5 %Compute Mass
        
        ElemM = zeros(nst);
        
        thick = 1;
        
        % Load Guass Integration Points
        lint = nel;
        der = 0;
        bf = 1;
        ib = 0;
        Nmat = zeros(1,nst);
        
        if t_on
        tau = 1+6*Nbeta*CFL^2*(1-cosh(1/sqrt(Nbeta)/CFL))/(2+cosh(1/sqrt(Nbeta)/CFL));
        else
        tau = 0;
        end
        
        sw = int1d(lint);
        [shl,shld,shls,be] = shl1d(sw(1,:),lint,nel-1);
        
        for ll = 1:lint                    

            % Evaluate 1-D basis functions at integration points
            Wgt = sw(2,ll);
            [shg, shgs, Jdet, be] = shg1d(xl,ndm,nel,shld(ll,:),shls(ll,:),nen,bf,der,be);
            
            c1 = Wgt*Jdet;
            
            Nmat(1:ndf:ndf*(nel-1)+1) = shl(ll,:);
            
            ElemM = ElemM + c1*(Nmat'*Nmat)*(1-tau);
%             ElemM = ElemM + c1*rho*(Nmat'*Nmat);

        end %je
        ElemM;

    case 6 %Compute Residual
        
        ElemF = zeros(ndf*nel,1);

        %Set integration number
        lint = nel;
        
        ib = 0;
        bf = 1;
        der = 1;
        
        fb = 0;
        Dmat = cwave^2;
        Nmat = zeros(1,ndf*nel);
        Bmat = zeros(1,ndf*nel);
        BBmat = zeros(1,nst);
        
        if t_on
        tau = 1+6*Nbeta*CFL^2*(1-cosh(1/sqrt(Nbeta)/CFL))/(2+cosh(1/sqrt(Nbeta)/CFL));
        else
        tau = 0;
        end

        % Compute kinematic fields at intermediate time levels
        ul_n_af = (1 - Nalphaf)*ul + Nalphaf*ul_n;
        al_n_am = (1 - Nalpham)*al + Nalpham*al_n;
        ulres = reshape(ul_n_af,ndf*nel,1);
        alres = reshape(al_n_am,ndf*nel,1);
        
        sw = int1d(lint);
        [shl,shld,shls,be] = shl1d(sw(1,:),lint,nel-1);
        
        for ll = 1:lint                    

            % Evaluate 1-D basis functions at integration points
            Wgt = sw(2,ll);
            [shg, shgs, Jdet, be] = shg1d(xl,ndm,nel,shld(ll,:),shls(ll,:),nen,bf,der,be);
            
            c1 = Wgt*Jdet;
            
            Nmat(1:ndf:ndf*(nel-1)+1) = shl(ll,:);
                
            Bmat(1:ndf:ndf*(nel-1)+1) = shg';
            
            BBmat(1:ndf:ndf*(nel-1)+1) = shgs';

            %Compute rho, accel
            accel = Nmat*alres;
            sigma = Dmat*(Bmat*ulres);
            resid = -(Nmat*alres - Dmat*BBmat*ulres);
%             resid = -(Nmat*alres + Dmat*Bmat*ulres);
            
            Fvec = fb;
            
            ElemF = ElemF + c1*(Nmat'*Fvec - Nmat'*accel - Bmat'*sigma ...
                                - (Nmat - (Nbeta*tstep^2)*Dmat*BBmat)'*tau*resid);

        end %je

    case 21 %Compute Stiffness
        
        ElemK = zeros(nst);
        
    case 40 % Initialize FS acceleration
        
        
    case 12 % energy
        
        ElemE = 0;
        
        
        % Load Guass Integration Points

        lint = nel;
        der = 1;
        bf = 1;
        ib = 0;
        
        Dmat = cwave^2;
        D = Dmat;
        Nmat = zeros(1,nst);
        Bmat = zeros(1,nst);
        BBmat = zeros(1,nst);
        
        sw = int1d(lint);
        [shl,shld,shls,be] = shl1d(sw(1,:),lint,nel-1);

        %Integration Loop
        for ll = 1:lint

                % Evaluate 1-D basis functions at integration points
                Wgt = sw(2,ll);
                [shg, shgs, Jdet, be] = shg1d(xl,ndm,nel,shld(ll,:),shls(ll,:),nen,bf,der,be);

                c1 = Wgt*Jdet;

                Nmat(1:ndf:ndf*(nel-1)+1) = shl(ll,:);

                Bmat(1:ndf:ndf*(nel-1)+1) = shg';

                BBmat(1:ndf:ndf*(nel-1)+1) = shgs';
                
                %Compute F, sigma, D
                E = Bmat*reshape(ul,nst,1);
                velo = Nmat*reshape(vl,nst,1);
                
                ElemE = ElemE + c1/2*(E'*D*E + velo'*velo);

        end %je
ElemE;
end