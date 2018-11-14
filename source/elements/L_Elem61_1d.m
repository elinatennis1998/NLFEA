
% 1D element for bar problem
PatchE = mateprop(1);% Elastic modulus
PatchA = mateprop(2);% Crossectional Area of the bar
%% Task Switch - Perform various element functions
switch isw 
    %%    
    case 1 % Setup up elemental degrees of freedom (OPTIONAL)
        
% Purpose: Rearrange order of dofs on element and/or flag certain
%          dofs to be presecribed away
        if ndf > 1
            
            for i = 2:ndf
                lie(i,1) = 0;
            end
            
        end
        istv = 1; % number of stresses per node for post-processing
        iste = nen; % number of stresses per element (total number, including all integration points)
 %%
    case 3 % Stiffness and force vector (REQUIRED)
        
% Purpose: Compute stiffness matrix and force vector for element
% Called by: FormFE.m
        ElemK = zeros(ndf*nel,1); %initiate Elemental Stiffness matrix
        ElemF= zeros(ndf*nel,1); %initiate Elemental Force vector
        Nmat = zeros(1,ndf*nel); %initiate Elemental Shape function
        Dmat = PatchE;
        Bmat = zeros(1,ndf*nel);
        ulres = reshape(ul,ndf*nen,1);

        % Load Gauss Points for quadrature
        lint = nel - 1;
        der = 0;
        bf = 0;
        ib = 0;
        sw = int1d(lint);
        [shl,shld,shls,be] = shl1d(sw(1,:),lint,nel-1);
        
        for je = 1:lint

            Wgt = sw(2,je);
            [shg, shgs, Jdet, be] = shg1d(xl,ndm,nel,shld(je,:),shls(je,:),nen,bf,der,be);
                    
            c1 = Wgt*Jdet*PatchA;
                
            Nmat(1:ndf:ndf*(nel-1)+1) = shl(je,:);
                
            Bmat(1:ndf:ndf*(nel-1)+1) = shg';
            
%            sigma = Dmat*Bmat*ulres(1:ndf*nel);
%             
%             ElemF = ElemF - c1*Bmat'*sigma;
            
            ElemK = c1*Bmat'*Dmat*Bmat;
            
        end %je           
%%
   case 25 %Stress Projection2

        ElemS = zeros(nel,npstr+1);
        ElemS2 = zeros(nel,npstr);

        Dmat = PatchE;
        Bmat = zeros(1,ndf*nel);
        
        % Load Guass Integration Points

        if nel == 2
            lint = 2;
            nint = 1;
        elseif nel == 3
            lint = 3;
            nint = 2;
        else
            lint = 3;
            nint = 2;
        end
        
        der = 0;
        bf = 0;
        ib = 0;
        
        ulres = reshape(ul,ndf*nen,1);

        %Stress Loop
        
        sw = int1d(nint);
        [shl,shld,shls,be] = shl1d(sw(1,:),lint,nel-1);
        
        
        for ll = 1:nint                    

            % Evaluate 1-D basis functions at integration points
            Wgt = sw(2,ll);
            [shg, shgs, Jdet, be] = shg1d(xl,ndm,nel,shld(ll,:),shls(ll,:),nen,bf,der,be);
            
            % Form B matrix
            Bmat(1:ndf:ndf*(nel-1)+1) = shg';
            
            epsil = Bmat*ulres(1:ndf*nel);
            sigma = Dmat*epsil;
            
            ElemS2(ll,1) = sigma;

        end %je
        
        % interpolate stress at nodes
        if nel == 2
            plist = [-1 1];
        elseif nel == 3
            sqr3 = sqrt(3);
            plist = [-sqr3 0 sqr3];
        else
            sqr3 = sqrt(3);
            plist = [-sqr3 0 sqr3];
        end
        
        for ll = 1:lint
            
            r = plist(1,ll);
            shpS = shl1d(r,1,nint-1); % polynomial order p = nint - 1
            
%             for stres = 1:npstr
%                 
%                 sigmas = ElemS2(1:nint,stres)'*shpS;
%                 ElemS(ll,stres) = sigmas;
%                 
%             end
            ElemS(ll,1:npstr) = (shpS*ElemS2(1:nint,:))';
            
        end
        
%         %Integration Loop
%         Vol = xl(2) - xl(1);

        for i = 1:nel
        ElemS(i,npstr+1) = 1; % use this for simple average Vol; % use this for weighted average 
        end
        
end %Task Switch