% Tim Truster
% 10/16/2016
%
% Master pure-displacement element
% body force problem re-run for Pinlei on 5/28/13, no other changes made


%% Global declarations
% Place here any commands that should be executed
% whenever the element routine is called.
if exist('iprob','var') == 0
    iprob = 0;
end


%% Task Switch - Perform various element functions
switch isw 
%%    
    case 1 % Setup up elemental degrees of freedom (OPTIONAL)
        
% Purpose: Rearrange order of dofs on element and/or flag certain
%          dofs to be presecribed away
% Called by: pmatin.m
% Notes: Setting values of lie(i,1) zero for i = 1:ndf flags those
%        dofs on all nodes to be unused; global dofs will be
%        prescribed as zero if all elements connected to that node
%        have flagged that dof.
%        Individual dofs on each node are handled by slots lie(i,j)
%        for j=2:nen+1.
%        History variables are allocated using nh1 and nh3.
%        Dof reordering is handled in MatTypeTable in the input
%        file; see NL_Elem8_3d.m for an example.
%        See FEAP pmanual.pdf and pmatin.m for more details.
% Example: NL_Elem8_3d.m
        
        % Example
        if ndf > 2
            
            for i = 3:ndf
                lie(i,1) = 0;
            end
            
        end

%         nh1 = 0; % number of time dependent history variables
%         nh3 = 0; % number of time independent history variables
%         istv = 0; % number of stresses per node for post-processing
%         iste = 0; % number of stresses per element (total number, including all integration points)
%         istee = 0; % number of error norm quantities

%%
    case 3 % Stiffness and internal force vector (REQUIRED)
        
        [ElemK,ElemF,hr] = NL_Elem2_2dF03(mateprop,ul,xl,ElemFlag,hr,nh1,nh2,nh3,ndf,ndm,nst,nel,nelP,nen,iprob);
        

%%        
    case -1 % Boundary tractions (RECOMMENDED)
        
        ElemF = NL_Elem2_2dFm1(mateprop,nodeA,nodeB,elem,edge,traction,xl,ElemFlag,hr,nh1,nh2,nh3,ndf,ndm,nst,nel,nelP,nen,iprob);
        
        
%%
    case 6 % Internal force vector (RECOMMENDED)
        
        [ElemF,hr] = NL_Elem2_2dF06(mateprop,ul,xl,ElemFlag,hr,nh1,nh2,nh3,ndf,ndm,nst,nel,nelP,nen,iprob);


%%
    case 7 % user boundary tractions (OPTIONAL)

% Purpose: Compute user surface traction for an element; these forces are
%          recomputed at every step (iteration)
% Called by: ploadu.m
% Notes: The ordering of dofs in the external load vector are as follows:
%        ElemF = [node1_dof1 node1_dof2 ... node1_dofndf node2_dof1 ...
%                 node2_dofndf node3_dof1 ... nodenen_dofndf]';
%        Integration of the force can be handled in a standard fashion by
%        rotating the parameteric space of the current element to match a
%        template face over which integration is always performed.
% Example: NL_Elem2_2dM.m    
        
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
    case 11 % Error estimation (OPTIONAL)
        
        [ElemE,Ieffvals] = NL_Elem2_2dF11(mateprop,ul,xl,ElemFlag,hr,nh1,nh2,nh3,ndf,ndm,nst,nel,nelP,nen,numEn,iprob,Ieffvals);
        
        
%%
    case 15 % Body force calculation (OPTIONAL)
        
        ElemF = NL_Elem2_2dF15(mateprop,bodyf,xl,ElemFlag,hr,nh1,nh2,nh3,ndf,ndm,nst,nel,nelP,nen,iprob);


%%
    case -15 % User Body Force (OPTIONAL)
        
% Purpose: Compute body force for element; these forces are recomputed at 
%          every step (iteration)
% Called by: pbodyfu.m
% Notes: The ordering of dofs in the external load vector are as follows:
%        ElemF = [node1_dof1 node1_dof2 ... node1_dofndf node2_dof1 ...
%                 node2_dofndf node3_dof1 ... nodenen_dofndf]';
%        Body forces are treated as external loads, which can be scaled by
%        a proportional loading parameter.
%        ElemF is initialized as zero inside pbodyf.m and does not need to
%        be reinitialized here.
% Example: NL_Elem2_2dM.m
        
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
    case 25 % Stress projection to nodes (RECOMMENDED)
        
        ElemS = NL_Elem2_2dF25(mateprop,ul,xl,ElemFlag,hr,nh1,nh2,nh3,ndf,ndm,nst,nel,nelP,nelS,nen,npstr);
        

%%        
    case 26 % Element Stress (OPTIONAL)
        
        ElemS = NL_Elem2_2dF26(mateprop,ul,xl,ElemFlag,hr,nh1,nh2,nh3,ndf,ndm,nst,nel,nelP,nen,nestr);
        

end %Task Switch
