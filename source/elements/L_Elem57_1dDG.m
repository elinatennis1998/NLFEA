%Elina Geut. Last modified 11/11/2018
%Element subroutine for DG 1D element. 

%pL - polynomial order on the left
%pR - plynomial order on the right
%lintL, lintR - # of integrational points on the Left and Right
%respectively
%Patch(A)(E) - material properties such as cross section and Young's
%Modulus
%Dmat - constitutve relationship L-left, R-right
%r - integration point, w - integration weight
%xl - coordiates of the current element
%nel - nodes on element 
%Bmat - derivative of shape functions 
 
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
    case 3 %Element Stiffness matrix for DG (1D)
%Modification of 11/12/2018, The set of if statements is implemented for a
%case of connection of p=1 and p=2 shape functions
   
        pL = mateprop(3);                                                  
        pR = mateprop(4);
        lintL = 1; lintR = 1;
        nelL = pL+1;
        nelR = pR+1;
        
        PatchEL = MateT(mateprop(1),1);
        PatchAL = MateT(mateprop(1),2);
        PatchER = MateT(mateprop(2),1);
        PatchAR = MateT(mateprop(2),2);

        DmatL = PatchEL*PatchAL;
        DmatR = PatchER*PatchAR;
       
         % Use one point integration rule
        
        lintL;lintR;ndm;ndf;pL;pR;nel;
        
        %Calculations for left side. Created 11/11/2018
        
        rL = 1;                                                          
        rR = -1;
        wL = 1;
        wR = 1;
        
        [NmatL,DerNL,NmatR,DerNR] = shf1dDG_Elina_M(rL,rR,ndm,ndf,pL,pR,lintL,lintR);
        xlL = xl(1:nelL);
        Bmat = zeros(lintL,nelL);
        
        for i = 1:lintL 
            w1 = wL(:,i);
            [Bmat(i,:),DetJ] = numint_Elina(xlL,DerNL(i,:),nelL);
            cL = w1*DetJ;
        end
        
        BmatL = Bmat;
        DetJL = DetJ;
        
        %Calculations for right side. Created 11/11/2018
        
        xlR = xl(:,(nelL+1:nelL+nelR));
        Bmat = zeros(lintR,nelR);
        
        for i = 1:lintR 
            w1 = wR(:,i);
            [Bmat(i,:),DetJ] = numint_Elina(xlR,DerNR(i,:),nelR);
            cR = w1*DetJ;
        end
        
        DetJR = DetJ;
        BmatR = Bmat;
        
        KLL2 = (-0.5*wL).*NmatL'*DmatL*BmatL;
        KRL2 = (0.5*wR).*NmatR'*DmatL*BmatL;
        KLR2 = (-0.5*wL).*NmatL'*DmatR*BmatR;
        KRR2 = (0.5*wR).*NmatR'*DmatR*BmatR;

        
        K2 = [KLL2 KLR2                                                      
              KRL2 KRR2];
          
        %Third Term of Weak Form
        
        K3 = K2';

        %Fourth Term of Weak Form
        
        tau = 10^(1);  
        KLL4 = (0.5*wL*tau).*(NmatL'*NmatL);
        KLR4 = (-0.5*wL*tau).*(NmatL'*NmatR);
        KRL4 = (-0.5*wR*tau).*(NmatR'*NmatL);
        KRR4 = (0.5*wR*tau).*(NmatR'*NmatR);

        K4 = [KLL4 KLR4                                                      
              KRL4 KRR4];
       
        ElemK = K2+K3+K4;
        
end

        
        
        
         