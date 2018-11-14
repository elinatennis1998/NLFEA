



%Element Subroutine for a Quadrelataral Element with Linear Shape Function

PatchE = mateprop(1,1);                                                    %Young's Modulus 
Patchv = mateprop(1,2);                                                    %Poisson's Ratio
Patcht = mateprop(1,3);                                                    %Thickness

switch isw 
%%
    case 1 %Setup up elemental degrees of freedom
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
    case 3 %Stiffness Matrix %Works for 1 element, but not for 2. the problem is in coordinate changes
        ElemK = zeros(ndf*nen,ndf*nen);                                    %Space Preallocation
        Nmat = zeros(ndm,nen);
        Bmat = zeros(ndm,ndf*nen);
        DerN = zeros(ndm,ndf*nen);
        xl = Coordinates; 
        
%         if rem(elem/2,1) == 0
%             r(1) = 1/(sqrt(3));                                            %Isoparametric Coordinates
%             r(2) = 1/(sqrt(3));                                            % N1 = (r1,s1)  
%             r(3) = 1/(sqrt(3));                                            % N2 = (r2,s2)                          
%             r(4) = -1/(sqrt(3));                                           % N3 = (r3,s3)
%             s(1) = -1/(sqrt(3));                                           % N4 = (r4,s4)
%             s(2) = 1/(sqrt(3));
%             s(3) = 1/(sqrt(3));
%             s(4) = 1/(sqrt(3));
%         else
            r(1) = -1/(sqrt(3));                                           %Isoparametric Coordinates
            r(2) = 1/(sqrt(3));                                            % N1 = (r1,s1)  
            r(3) = 1/(sqrt(3));                                            % N2 = (r2,s2)                          
            r(4) = -1/(sqrt(3));                                           % N3 = (r3,s3)
            s(1) = -1/(sqrt(3));                                           % N4 = (r4,s4)
            s(2) = -1/(sqrt(3));
            s(3) = 1/(sqrt(3));
            s(4) = 1/(sqrt(3));
%         end
                                                                           %Constitutive Relation
        T1 = PatchE/((1+Patchv)*(1-2*Patchv));                             %Term 1 of Dmat
        T2 = [1-Patchv Patchv 0];                                          %Term 2 of Dmat
        T3 = [Patchv 1-Patchv 0];                                          %Term 3 of Dmat
        T4 = [0 0 0.5*(1-2*Patchv)];                                       %Term 4 of Dmat
        
        Dmat = T1.*[T2;T3;T4];
        
         for i=1:4
            
            w = 1;
        
            Nmat1(1,1) = 0.25*(1-r(i))*(1-s(i));
            Nmat1(1,2) = 0.25*(1+r(i))*(1-s(i));
            Nmat1(1,3) = 0.25*(1+r(i))*(1+s(i));
            Nmat1(1,4) = 0.25*(1-r(i))*(1+s(i));
            
            dNdr = 0.25.*[-(1-s(i)) (1-s(i)) (1+s(i)) -(1+s(i))];
            dNds = 0.25.*[-(1-r(i)) -(1+r(i)) (1+r(i)) (1-r(i))];
            
            dN = [dNdr; dNds];
            Nmat = [Nmat1; Nmat1];
            
%             J = dN*xl(:,1:nel);
%             inJ = inv(J);
%             DetJ = abs(det(J));
%            
% Find a better approach to compute Jacobian
            if elem == 1
                J = dN*xl(1:nel,:);
                inJ = inv(J);
                DetJ = abs(det(J));
            elseif elem == 2
                J = dN*xl(5:8,:);
                inJ = inv(J);
                DetJ = abs(det(J));
            elseif elem ==3
                J = dN*xl(9:12,:);
                inJ = inv(J);
                DetJ = abs(det(J));
            end
                 
            for j=1:nen
            
                 Bmat(1,(2*j-1)) = (inJ(1,:))*dN(:,j);
                 Bmat(2,(2*j)) = (inJ(2,:))*dN(:,j);
                 Bmat(3,(2*j-1)) = (inJ(2,:))*dN(:,j);
                 Bmat(3,(2*j)) = (inJ(1,:))*dN(:,j);
            end
                 
            c1 = w*DetJ*Patcht;
            ElemK = ElemK+c1*Bmat'*Dmat*Bmat;
        end

         ElemK
end