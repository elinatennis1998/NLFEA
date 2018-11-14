function [Nmat,Bmat,DetJ,c] = shf1d_Elina_M(lint,xl,nel,ndm,ndf,p)
%Elina Geut 11/3/2018
%   Function is used for assembling shape functions and it's derivatives in
%   terms of glabal coordinated by multiplying 

 
if lint == 1
    r(1) = -1;
    r(2) = 1;
    w = 2;

elseif lint == 2
    r(1) = -1/(sqrt(3));
    r(2) = 1/(sqrt(3));
    w = 1;

elseif lint == 3
    Nmat = zeros(ndm,ndf*nen);
    DerN = zeros(ndm,ndf*nen);
    r(1) = -(sqrt(6));
    r(2) = sqrt(6);
    r(3) = 0;
    w(1) = 5/9;
    w(2) = 8/9;
end

if p == 1
    Nmat1 = (1-r(1))/2;
    Nmat2 = (1+r(2))/2;
    DerN1 = -1/2;
    DerN2 =1/2;
    DerN = [DerN1 DerN2];
    J = DerN*xl(:,1:nel)';
    Bmat = (1/J)*DerN;
    Nmat = (1/J)*[Nmat1 Nmat2];
    DetJ = det(J);
    
    c = DetJ*w;
    
elseif p == 2
     Nmat = zeros(ndm,ndf*nel);
     DerN = zeros(ndm,ndf*nel);
     
    for i=1:lint
        Nmat(i,1) = (1-r(i))/2*(-r(i));
        Nmat(i,2) =  (1+r(i))*(1-r(i));
        Nmat(i,3) = (1+r(i))/2*(r(i));
        
        DerN(i,1) = -1/2+r(i);
        DerN(i,2) = -2*r(i);
        DerN(i,3) = 1/2+r(i);
        
        J = DerN*xl(:,1:nel)';
        Bmat = (1/J)*DerN;
        Nmat = (1/J)*Nmat;
        DetJ = det(J);
        c = w*DetJ;
    end
    
elseif p == 3

    for ii = 1:2
        for i=1:lint
         Nmat(i,1) =  1/16*(9*r(i)^2 - 1)*(1 - r(i));
         Nmat(i,2) =  9/16*(r(i)^2 - 1)*(3*r(i) - 1);
         Nmat(i,3) =  9/16*(1 - r(i)^2)*(3*r(i) + 1);
         Nmat(i,4) =  1/16*(9*r(i)^2 - 1)*(r(i) + 1);
            
         DerN(i,1) =  1/16*(1 + r(i)*(18 - 27*r(i)));
         DerN(i,2) =  9/16*(-3 - r(i)*(2 - 9*r(i)));
         DerN(i,3) =  9/16*(3 - r(i)*(2 + 9*r(i)));
         DerN(i,4) =  1/16*(-1 + r(i)*(18 + 27*r(i)));
         
         J = DerN*xl(:,1:nel)';
         Bmat = (1/J)*DerN;
         Nmat = (1/J)*Nmat;
         DetJ = det(J);
        
        end
        c = w(ii)*DetJ;
     end
end