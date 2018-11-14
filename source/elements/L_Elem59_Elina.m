


%Elina Geut. Last modified 11/10/2018
%Subroutin for 1D solid element. The results for all cases are
%verified with Patch Test for p =1,2. Used in DG input file as well

 PatchE = mateprop(1);                                                     %Young's Modulus
 PatchA = mateprop(2);                                                     %Cross-Sectional Area
 
switch isw                                                                 %Task switch 
%%
    case 1
        if ndf > ndm 
            for i = ndm+1:ndf
                lie(i,1) = 0;
            end
            
        end
        nh1 = 0;                                                           %# of time dependent history vars
        nh3 = 0;                                                           %# of time indep history vars
        istv = 0;                                                          %# of stresses per node for post-proc
        iste = 0;                                                          %# of stresses per element
        istee = 0;                                                         %# of error norm quantities
        
%%
    case 3 %Element stiffness matrix
        %Case 3 was checked and runs properly. Last modified 11/07/2018
        
        ElemK = zeros(ndf*nel,ndf*nel);
        ElemF = zeros(ndf*nel,ndm);
        Nmat = zeros(ndm,nel);
        Bmat = zeros(ndm,ndf*nel);
        Dmat = PatchE*PatchA;                                              %Constitutive Relationship
         
        nel;p;nen;lint;ndm;
        
        [r,w] = gauss_Elina(lint);
        [Nmat,DerN] = shf1d_Elina_M(r,ndm,ndf,p,lint);
        
        for i = 1:lint 
            w1 = w(:,i);
            [Bmat,DetJ] = numint_Elina(xl,DerN(i,:),nel);
            c = w1*DetJ;
            ElemK = ElemK+c.*Bmat'*Dmat*Bmat;
        end
                                                                                        
%%  
    case 6 %Element Residual Vector
%         [ElemF] = BlankElem06(mateprop,ul,xl,ElemFlag,ndf,ndm,nst,nel,nen)
%Case 6 was checked for 1 case only and passed that test.
%Residual forces fro 1D DG element with p=1 works properly 
%Last Modified 11/11/2018
        ElemF = zeros(ndf*numel,1);
        Nmat = zeros(ndm,nel);
        Bmat = zeros(ndm,ndf*nel);
        Dmat = PatchE*PatchA;  
        
        ul = ul(ndf,ndf:nen);                                              %Local Displacements
        ulres = reshape(ul,ndf*nen,1);
        
        nel;p;nen;lint;ndm;
        
        [r,w] = gauss_Elina(lint);
        [Nmat,DerN] = shf1d_Elina_M(r,ndm,ndf,p,lint);
        
        for i = 1:lint 
            w1 = w(:,i);
            [Bmat,DetJ] = numint_Elina(xl,DerN(i,:),nel);
            c = w1*DetJ*PatchA; 
            
            sigma = Dmat*Bmat*ulres;
            ElemF = ElemF - c*Bmat'*sigma;
        end
                
        
%%
    case 15 %Include Body Forces
%         ElemF = BlankElem15(mateprop,ul,xl,ElemFlag,ndf,ndm,nst,nel,nen)

%The case 15 was checked and runs properly (CG & p=1 DG). Last modified 11/11/2018

        nel;p;nen;lint;ndm;
     
        [r,w] = gauss_Elina(lint);
        [Nmat,DerN] = shf1d_Elina_M(r,ndm,ndf,p,lint);
        
        bf = bodyf;
        dint = lint-1;
        
        for i = 1:lint 
            w1 = w(:,i);
            [Bmat,DetJ] = numint_Elina(xl,DerN(i,:),nel);
            c = w1*DetJ*PatchA; 
            %Error within the bf calculations. Nmat causing error
            ElemF = ElemF + c*(Nmat(i,:)*bf)';

        end
        
        
%%
    case 25 %Stress Projections 
     %The issues are fixed. Case 25 runs properly. Last modified 11/11/2018
     %The issue with DG stress projection was fixed
     %Talk to Dr. Truster about the stress projection polinomial order of
     %interpolation
     
        PatchE = mateprop(1);                                              %Young's Modulus
        PatchA = mateprop(2);                                              %Cross-Sectional Area
        Sigma1 = zeros(nel,npstr+1);
        Sigma2 = zeros(nel,npstr+1);
        Sigma = zeros(nel,npstr+1);
        Dmat = PatchE*PatchA;

        
        ul = ul(ndf,ndf:nen);                                              
        ulres = reshape(ul,ndf*nen,1);      
        nel;p;

        [r,w] = gauss_Elina(dint);
        [Nmat,DerN] = shf1d_Elina_M(r,ndm,ndf,p,dint);
        
         if nel == 2
            ri = [-1 1];
        elseif nel > 2
            ri = [-sqrt(3) 0 sqrt(3)];
            lint = 3;
         end 
         
        dint = lint-1;                                                     %# of int.pts. for displ. approx.
        
        for i = 1:dint 
            w1 = w(:,i);
            [Bmat,DetJ] = numint_Elina(xl,DerN(i,:),nel);
            Epsil = Bmat*ulres(1:ndf*nel);
            Sigma1(i,1) = Dmat*Epsil;
        end
        
        for j = 1:lint 
            [Nmat,DerN] = shf1d_Elina_M(ri(j),ndm,ndf,p-1,1);
            Sigma(j,1) = (Nmat*Sigma1(1:dint,1))';
            
        end
        for i = 1:nel
            Sigma2(i,npstr+1) = 1;
        end
        
         ElemS = Sigma+Sigma2;
        
end

 

