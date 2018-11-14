function [sigma, D] = SigmaCmatNSCST3i(F,J,mateprop)
%
% Shear constitutive material model evaluation for NL elasticity
% 10/16/2011
% Output: sigma = J*sigma = tau = Kirchhoff stress
%         D = F_iI*F_jJ*F_kK*F_lL*C_IJKL

% 07/02/2013 - add a zero placeholder in mateprop(3) for matp to all input files 

one = [1; 1; 1; 0; 0; 0; 0; 0; 0];
mat1 = one*one';
matE = diag([2,2,2,1,1,1,0,0,0]);

mati = mateprop(1);

if mati <= 4
    if mati <= 2
        if mati == 1 % Standard Neo-Hookean Material

            PatchE = mateprop(4);
            Patchv = mateprop(5);
            mu = PatchE/(2*(1+Patchv));
            
            b = F*F';
            bv = [b(1,1); b(2,2); b(3,3); b(1,2); b(2,3); b(1,3); 0; 0; 0];

            sigma = mu*(bv - one);

            D = mu*matE;
            
        else % mati == 2 % Isochoric Neo-Hookean Material

            PatchE = mateprop(4);
            Patchv = mateprop(5);
            mu = PatchE/(2*(1+Patchv));%5.67;%1.61148;%
            
            b = F*F';
            bv = [b(1,1); b(2,2); b(3,3); b(1,2); b(2,3); b(1,3); 0; 0; 0];
            trb = trace(b);
            J23 = J^(-2/3);

            sigma = mu*J23*(bv - 1/3*trb*one);
            
            D = -2/3*(one*sigma' + sigma*one' + mu*J23*trb*(1/3*mat1 - 1/2*matE));

        end
    elseif mati == 3 % Mooney-Rivlin Material
        
%         PatchE = mateprop(4);
%         Patchv = mateprop(5);
        c1 = mateprop(6);
        c2 = mateprop(7);
% mu1/2=c1, mu2/2=c2, mu=mu1+mu2=2*(c1+c2) where mu is the small-strain shear modulus
        
        b = F*F';
        bv = [b(1,1); b(2,2); b(3,3); b(1,2); b(2,3); b(1,3); 0; 0; 0];
        b2 = b*b;
        b2 = [b2(1,1); b2(2,2); b2(3,3); b2(1,2); b2(2,3); b2(1,3); 0; 0; 0];
        b3 = (bv*bv');
        b4 = [2*bv(1)*bv(1) 2*bv(4)*bv(4) 2*bv(6)*bv(6) 2*bv(1)*bv(4) 2*bv(4)*bv(6) 2*bv(6)*bv(1) 0 0 0
              2*bv(4)*bv(4) 2*bv(2)*bv(2) 2*bv(5)*bv(5) 2*bv(2)*bv(4) 2*bv(5)*bv(2) 2*bv(4)*bv(5) 0 0 0
              2*bv(6)*bv(6) 2*bv(5)*bv(5) 2*bv(3)*bv(3) 2*bv(6)*bv(5) 2*bv(5)*bv(3) 2*bv(3)*bv(6) 0 0 0
              2*bv(1)*bv(4) 2*bv(2)*bv(4) 2*bv(6)*bv(5) (bv(1)*bv(2)+bv(4)*bv(4)) (bv(4)*bv(5)+bv(6)*bv(2)) (bv(6)*bv(4)+bv(1)*bv(5)) 0 0 0
              2*bv(4)*bv(6) 2*bv(5)*bv(2) 2*bv(5)*bv(3) (bv(4)*bv(5)+bv(6)*bv(2)) (bv(2)*bv(3)+bv(5)*bv(5)) (bv(5)*bv(6)+bv(4)*bv(3)) 0 0 0
              2*bv(6)*bv(1) 2*bv(4)*bv(5) 2*bv(3)*bv(6) (bv(6)*bv(4)+bv(1)*bv(5)) (bv(5)*bv(6)+bv(4)*bv(3)) (bv(3)*bv(1)+bv(6)*bv(6)) 0 0 0
              zeros(3,9)];
        I1 = trace(b);
        I2 = 1/2*(I1^2 - (trace(b*b)));
        J23 = J^(-2/3);
        J43 = J^(-4/3);

        sigma1 = 2*c1*J23*(bv - 1/3*I1*one);
        sigma2 = 2*c2*J43*(I1*bv - b2 - 2/3*I2*one);

        sigma = sigma1 + sigma2;
        
        D1 = -2/3*(one*sigma1' + sigma1*one' + 2*c1*J23*I1*(1/3*mat1 - 1/2*matE));
        D2 = -4/3*(one*sigma2' + sigma2*one') + 4*c2*J43*(b3 - 1/2*b4 - 4/9*I2*mat1 + 1/3*I2*matE);

        D = D1 + D2;
        
    else % mati == 4 % Ogden Material, Reese IJNME38
        
        n = 3;
        mu = mateprop(6:6+n-1);
        a = mateprop(6+n:6+2*n-1);

        b = F*F';
        [m,bp] = eig(b);
        bp = sqrt(diag(bp));
        Q = m';
        P = [Q(1,1)^2 Q(1,2)^2 Q(1,3)^2 2*Q(1,1)*Q(1,2) 2*Q(1,2)*Q(1,3) 2*Q(1,3)*Q(1,1)
             Q(2,1)^2 Q(2,2)^2 Q(2,3)^2 2*Q(2,1)*Q(2,2) 2*Q(2,2)*Q(2,3) 2*Q(2,3)*Q(2,1)
             Q(3,1)^2 Q(3,2)^2 Q(3,3)^2 2*Q(3,1)*Q(3,2) 2*Q(3,2)*Q(3,3) 2*Q(3,3)*Q(3,1)
             Q(1,1)*Q(2,1) Q(1,2)*Q(2,2) Q(1,3)*Q(2,3) Q(1,1)*Q(2,2)+Q(1,2)*Q(2,1) Q(1,2)*Q(2,3)+Q(1,3)*Q(2,2) Q(1,3)*Q(2,1)+Q(1,1)*Q(2,3)
             Q(2,1)*Q(3,1) Q(2,2)*Q(3,2) Q(2,3)*Q(3,3) Q(2,1)*Q(3,2)+Q(2,2)*Q(3,1) Q(2,2)*Q(3,3)+Q(2,3)*Q(3,2) Q(2,3)*Q(3,1)+Q(2,1)*Q(3,3)
             Q(3,1)*Q(1,1) Q(3,2)*Q(1,2) Q(3,3)*Q(1,3) Q(3,1)*Q(1,2)+Q(3,2)*Q(1,1) Q(3,2)*Q(1,3)+Q(3,3)*Q(1,2) Q(3,3)*Q(1,1)+Q(3,1)*Q(1,3)];
        sig = zeros(6,1);
        for i = 1:3
            si = 0;
            for r = 1:n
                si = si + mu(r)*(bp(i)^a(r) - 1);
            end
            sig(i) = si;
        end

        sigma = [P*sig; zeros(3,1)];

        C = zeros(3,1);
        for i = 1:3
            si = 0;
            for r = 1:n
                si = si +  mu(r)*(2+(a(r)-2)*bp(i)^a(r));
            end
            C(i) = si;
        end
        % C1 = mu(1)*(2+(a(1)-2)*bp(1)^a(1)) + mu(2)*(2+(a(2)-2)*bp(1)^a(2));
        % C2 = mu(1)*(2+(a(1)-2)*bp(2)^a(1)) + mu(2)*(2+(a(2)-2)*bp(2)^a(2));
        C1 = C(1);
        C2 = C(2);
        C3 = C(3);
        if bp(1) ~= bp(2)
        C4 = (bp(2)^2*sig(1) - bp(1)^2*sig(2))/(bp(1)^2-bp(2)^2);
        else
        C4 = 1/2*C1;
        end
        if bp(2) ~= bp(3)
        C5 = (bp(3)^2*sig(2) - bp(2)^2*sig(3))/(bp(2)^2-bp(3)^2);
        else
        C5 = 1/2*C2;
        end
        if bp(3) ~= bp(1)
        C6 = (bp(1)^2*sig(3) - bp(3)^2*sig(1))/(bp(3)^2-bp(1)^2);
        else
        C6 = 1/2*C3;
        end

        Diso = P*diag([C1,C2,C3,C4,C5,C6])*P';
        Diso = [Diso zeros(6,3); zeros(3,9)];

        D = Diso;
        
    end
elseif mati > 5 && mati <= 8
    if mati <= 6
        if mati == 5 % Yeoh Material, Stein CMAME190

            matp = 1/2*matE - 1/3*mat1;

            C1 = mateprop(6);
            C2 = mateprop(7);
            C3 = mateprop(8);
            
            b = (F*F');
            J23 = J^(-2/3);
            bv = J23*[b(1,1); b(2,2); b(3,3); 2*b(1,2); 2*b(2,3); 2*b(3,1); 0; 0; 0];
            I1 = bv'*one;
            C1b = C1 + 2*C2*(I1-3) + 3*C3*(I1-3)^2;

            sigma = 2*C1b*matp*bv;
            
            D = -2/3*(one*sigma' + sigma*one' - 2*C1b*I1*matp) ...
                + 8*(C2 + 3*C3*(I1-3))*matp*(bv*bv')*matp;
            
        else % mati == 6 % Anisotropic, Stein CMAME190

one = [1; 1; 1; 0; 0; 0];
mat1 = one*one';
matE = diag([2,2,2,1,1,1]);
matp = 1/2*matE - 1/3*mat1;

            C1 = mateprop(6);
            C2 = mateprop(7);
            C3 = mateprop(8);
            C4 = mateprop(9);
            alpha = degtorad(mateprop(10));
            A = [cos(alpha), sin(alpha), 0];
            
            b = F*F';
            J23 = J^(-2/3);
            C = J23*(F'*F);
            bv = J23*[b(1,1); b(2,2); 1; 2*b(1,2); 0; 0];
            I1 = bv'*one;
            C1b = C1 + 2*C2*(I1-3) + 3*C3*(I1-3)^2;
            I4 = A*C*A';
            a = F*A';
            aa = a*a';
            C4b = 4*C4*(I4-1);
            C6 = 8*C4;
            av = J23*[aa(1,1); aa(2,2); 0; 2*aa(1,2); 0; 0];
            tra = av'*one;

            sigma = 2*C1b*matp*bv + C4b*matp*av;
            
            D = -2/3*(one*sigma' + sigma*one' - (2*C1b*I1 + C4b*tra)*matp) ...
                + 8*(C2 + 3*C3*(I1-3))*matp*(bv*bv')*matp ...
                + C6*(matp*(av*av')*matp);
            
            sigma = [sigma; 0; 0; 0];
            D = [D zeros(6,3); zeros(3,9)];
            
        end
    elseif mati == 7 % Ogden Material, Simo CMAME85
        
        %%%%%% FIX THE bpt and bp DISCREPANCIES
        
        n = 3;
        mu = mateprop(6:6+n-1);
        c = mateprop(6+n:6+2*n-1);

        b = F*F';
        bv = [b(1,1); b(2,2); b(3,3); b(1,2); b(2,3); b(1,3); 0; 0; 0];
        Ib = [2*bv(1)*bv(1) 2*bv(4)*bv(4) 2*bv(6)*bv(6) 2*bv(1)*bv(4) 2*bv(4)*bv(6) 2*bv(6)*bv(1) 0 0 0
              2*bv(4)*bv(4) 2*bv(2)*bv(2) 2*bv(5)*bv(5) 2*bv(2)*bv(4) 2*bv(5)*bv(2) 2*bv(4)*bv(5) 0 0 0
              2*bv(6)*bv(6) 2*bv(5)*bv(5) 2*bv(3)*bv(3) 2*bv(6)*bv(5) 2*bv(5)*bv(3) 2*bv(3)*bv(6) 0 0 0
              2*bv(1)*bv(4) 2*bv(2)*bv(4) 2*bv(6)*bv(5) (bv(1)*bv(2)+bv(4)*bv(4)) (bv(4)*bv(5)+bv(6)*bv(2)) (bv(6)*bv(4)+bv(1)*bv(5)) 0 0 0
              2*bv(4)*bv(6) 2*bv(5)*bv(2) 2*bv(5)*bv(3) (bv(4)*bv(5)+bv(6)*bv(2)) (bv(2)*bv(3)+bv(5)*bv(5)) (bv(5)*bv(6)+bv(4)*bv(3)) 0 0 0
              2*bv(6)*bv(1) 2*bv(4)*bv(5) 2*bv(3)*bv(6) (bv(6)*bv(4)+bv(1)*bv(5)) (bv(5)*bv(6)+bv(4)*bv(3)) (bv(3)*bv(1)+bv(6)*bv(6)) 0 0 0
              zeros(3,9)];
        [m,bp] = eig(b);
        bp = sqrt(diag(bp));
        bpt = bp/J^(1/3);
        I1 = trace(b);
        I3 = J^2;
        DA = [2*bp(1)^4-I1*bp(1)^2+I3*bp(1)^-2
              2*bp(2)^4-I1*bp(2)^2+I3*bp(2)^-2
              2*bp(3)^4-I1*bp(3)^2+I3*bp(3)^-2];
        DAp = [8*bp(1)^3-2*I1*bp(1)-2*I3*bp(1)^-3
               8*bp(2)^3-2*I1*bp(2)-2*I3*bp(2)^-3
               8*bp(3)^3-2*I1*bp(3)-2*I3*bp(3)^-3];
        bA = zeros(3,1);
        for A = 1:3
            si = 0;
            for a = 1:n
                si = si + c(a)/mu(a)*(bpt(A)^mu(a));
                for B = 1:3
                    si = si - 1/3*c(a)/mu(a)*(bpt(B)^mu(a));
                end
            end
            bA(A) = si;
        end
        gAB = zeros(3,3);
        for A = 1:3
            for B = 1:3
                si = 0;
                if A == B
                    for a = 1:n
                        si = si + 1/3*c(a)*(bpt(A)^mu(a));
                        for C = 1:3
                            si = si + 1/9*c(a)*(bpt(C)^mu(a));
                        end
                    end
                else
                    for a = 1:n
                        si = si - 1/3*c(a)*(bpt(A)^mu(a) + bpt(B)^mu(a));
                        for C = 1:3
                            si = si + 1/9*c(a)*(bpt(C)^mu(a));
                        end
                    end
                end
                gAB(A,B) = si;
            end
        end
        
        mmat = [m(1,1)^2      m(1,2)^2      m(1,3)^2
                m(2,1)^2      m(2,2)^2      m(2,3)^2
                m(3,1)^2      m(3,2)^2      m(3,3)^2
                m(1,1)*m(2,1) m(1,2)*m(2,2) m(1,3)*m(2,3)
                m(2,1)*m(3,1) m(2,2)*m(3,2) m(2,3)*m(3,3)
                m(3,1)*m(1,1) m(3,2)*m(1,2) m(3,3)*m(1,3)
                zeros(3,3)];
         
        if (bp(1) == bp(2)) && (bpt(2) == bp(3))
            sigma = zeros(6,1);
            si = 0;
            for a = 1:n
                si = si + c(a)*(bpt(A)^mu(a));
            end
            D = si*(matE/2 - 1/3*mat1);
        elseif bp(1) == bpt(2)
        dgm3 = 1/DA(3)*(Ib/2-(bv*bv')+I3*bp(3)^-2*(mat1-matE/2)) ...
             + 1/DA(3)*(bp(3)^2*(bv*mmat(:,3)'+mmat(:,3)*bv')-1/2*DAp(3)*bp(3)*mmat(:,3)*mmat(:,3)') ...
             - 1/DA(3)*I3*bp(3)^-2*(one*mmat(:,3)'+mmat(:,3)*one');
            sigma = bA(1)*one + (bA(3) - bA(1))*mmat(:,3);
            D = -bA(1)*matE + (bA(3)-bA(1))*2*dgm3 + gAB(1,1)*(one-mmat(:,3))*(one-mmat(:,3))' ...
              + gAB(3,3)*mmat(:,3)*mmat(:,3)' + gAB(1,3)*(mmat(:,3)*(one-mmat(:,3))' + (one-mmat(:,3))*mmat(:,3)');
        elseif bp(2) == bpt(3)
        dgm1 = 1/DA(1)*(Ib/2-(bv*bv')+I3*bp(1)^-2*(mat1-matE/2)) ...
             + 1/DA(1)*(bp(1)^2*(bv*mmat(:,1)'+mmat(:,1)*bv')-1/2*DAp(1)*bp(1)*mmat(:,1)*mmat(:,1)') ...
             - 1/DA(1)*I3*bp(1)^-2*(one*mmat(:,1)'+mmat(:,1)*one');
            sigma = bA(2)*one + (bA(1) - bA(2))*mmat(:,1);
            D = -bA(2)*matE + (bA(1)-bA(2))*2*dgm1 + gAB(2,2)*(one-mmat(:,1))*(one-mmat(:,1))' ...
              + gAB(1,1)*mmat(:,1)*mmat(:,1)' + gAB(1,2)*(mmat(:,1)*(one-mmat(:,1))' + (one-mmat(:,1))*mmat(:,1)');
        elseif bp(3) == bpt(1)
        dgm2 = 1/DA(2)*(Ib/2-(bv*bv')+I3*bp(2)^-2*(mat1-matE/2)) ...
             + 1/DA(2)*(bp(2)^2*(bv*mmat(:,2)'+mmat(:,2)*bv')-1/2*DAp(2)*bp(2)*mmat(:,2)*mmat(:,2)') ...
             - 1/DA(2)*I3*bp(2)^-2*(one*mmat(:,2)'+mmat(:,2)*one');
            sigma = bA(2)*one + (bA(2) - bA(1))*mmat(:,2);
            D = -bA(1)*matE + (bA(2)-bA(1))*2*dgm2 + gAB(1,1)*(one-mmat(:,2))*(one-mmat(:,2))' ...
              + gAB(2,2)*mmat(:,2)*mmat(:,2)' + gAB(1,2)*(mmat(:,2)*(one-mmat(:,2))' + (one-mmat(:,2))*mmat(:,2)');
        else
        dgm1 = 1/DA(1)*(Ib/2-(bv*bv')+I3*bp(1)^-2*(mat1-matE/2)) ...
             + 1/DA(1)*(bp(1)^2*(bv*mmat(:,1)'+mmat(:,1)*bv')-1/2*DAp(1)*bp(1)*mmat(:,1)*mmat(:,1)') ...
             - 1/DA(1)*I3*bp(1)^-2*(one*mmat(:,1)'+mmat(:,1)*one');
        dgm2 = 1/DA(2)*(Ib/2-(bv*bv')+I3*bp(2)^-2*(mat1-matE/2)) ...
             + 1/DA(2)*(bp(2)^2*(bv*mmat(:,2)'+mmat(:,2)*bv')-1/2*DAp(2)*bp(2)*mmat(:,2)*mmat(:,2)') ...
             - 1/DA(2)*I3*bp(2)^-2*(one*mmat(:,2)'+mmat(:,2)*one');
        dgm3 = 1/DA(3)*(Ib/2-(bv*bv')+I3*bp(3)^-2*(mat1-matE/2)) ...
             + 1/DA(3)*(bp(3)^2*(bv*mmat(:,3)'+mmat(:,3)*bv')-1/2*DAp(3)*bp(3)*mmat(:,3)*mmat(:,3)') ...
             - 1/DA(3)*I3*bp(3)^-2*(one*mmat(:,3)'+mmat(:,3)*one');
            sigma = mmat*bA;
            C = zeros(9,9);
            for A = 1:3
                for B = 1:3
                    C = C + gAB(A,B)*mmat(:,A)*mmat(:,B)';
                end
            end
            D = 2*(bA(1)*dgm1+bA(2)*dgm2+bA(3)*dgm3) + C;
        end
        
    else % mati == 8 % Hencky model
        
    end
    
else
    
    if mati == 9 % anisotropic without volumetric projection

            PatchE = mateprop(4);
            Patchv = mateprop(5);
            mu = PatchE/(2*(1+Patchv));
            C4 = mateprop(6);
            alpha = degtorad(mateprop(7));
            A = [cos(alpha), sin(alpha), 0];
            
            b = F*F';
            C = (F'*F);
            bv = [b(1,1); b(2,2); b(3,3); b(1,2); b(2,3); b(1,3); 0; 0; 0];
            I4 = A*C*A';
            a = F*A';
            aa = a*a';
            C4b = 4*C4*(I4-1);
            C6 = 8*C4;
            av = [aa(1,1); aa(2,2); aa(3,3); aa(1,2); aa(2,3); aa(3,1); 0; 0; 0];

            sigma = mu*(bv - one) + C4b*av;

            D = mu*matE + C6*(av*av');
    
    elseif mati == 10 % Compressible double anisotropic

            PatchE = mateprop(4);
            Patchv = mateprop(5);
            mu = PatchE/(2*(1+Patchv));
            C4 = mateprop(6);
            A = [mateprop(9) mateprop(10) mateprop(11)]; % Unit vector one
            
            b = F*F';
            C = (F'*F);
            bv = [b(1,1); b(2,2); b(3,3); b(1,2); b(2,3); b(1,3); 0; 0; 0];
            I4 = A*C*A';
            a = F*A';
            aa = a*a';
            C4b = 4*C4*(I4-1);
            C6 = 8*C4;
            av = [aa(1,1); aa(2,2); aa(3,3); aa(1,2); aa(2,3); aa(3,1); 0; 0; 0];

            sigma = mu*(bv - one) + C4b*av;

            D = mu*matE + C6*(av*av');
            
            C4 = mateprop(7);
            A = [mateprop(12) mateprop(13) mateprop(14)]; % Unit vector two
            
            C = (F'*F);
            I4 = A*C*A';
            a = F*A';
            aa = a*a';
            C4b = 4*C4*(I4-1);
            C6 = 8*C4;
            av = [aa(1,1); aa(2,2); aa(3,3); aa(1,2); aa(2,3); aa(3,1); 0; 0; 0];

            sigma = sigma + C4b*av;

            D = D + C6*(av*av');
            
    end
end