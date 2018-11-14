function [ElemE,Ieffvals] = NL_Elem2_2dF11(mateprop,ul,xl,ElemFlag,hr,nh1,nh2,nh3,ndf,ndm,nst,nel,nelP,nen,numEn,iprob,Ieffvals)
        
% Purpose: Compute error norms for element
% Called by: Explicit_Error_Estimation.m
% Notes: General routine for assembling scalar quantities across all
%        elements, typically the element error norms but other quantities
%        could be used instead. The usual ordering is by derivative and
%        then by dof; see L_Elem3_2dVMS.m for an example.
% Example: L_Elem3_2dVMS.m
        
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
        