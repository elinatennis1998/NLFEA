function D2 = J2Radial6thOrder(eps_n1,ep_n,beta_n,a_n,mu,bulk,K,H,sigy)
% 08/17/2014
% Tim Truster
% Function to compute 6th order curvature material tensor for J2 plasticity
% Verified against the numerical differentiation routine that both methods
% give a nearly identical value for the tensor. All symmetry properties are
% present.

D = zeros(3,3,3,3,3,3); % Tensor form

% e_n1 = zeros(6,1);
One = [1.0d0; 1.0d0; 1.0d0; 0.0d0; 0.0d0; 0.0d0];
I4 = diag([1.0d0 1.0d0 1.0d0 0.5d0 0.5d0 0.5d0]);
OneOne = (One*One');

% Constants
% R0 = 0.0d0;
% RP5 = 0.5d0;
R1 = 1.0d0;
R2 = 2.0d0;
R3 = 3.0d0;
sqr23 = sqrt(R2/R3);

% Set additional constants
Hb = K + H;
% Special case for perfect plasticity to avoid dividing by zero
if Hb == 0
bet = 0;
else
bet = K/Hb;
end
R2G = R2*mu;
R3G = R3*mu;

% Elastic predictors
eev = eps_n1(1) + eps_n1(2) + eps_n1(3); %volumetric strain
press = bulk*eev;
eevd3 = eev/R3;
e_n1 = eps_n1 - eevd3*One; %deviatoric strain
%convert engineering strain to physical component
str_n1 = R2G*I4*(e_n1 - ep_n); %trial elastic deviatoric stress
xtr_n1 = str_n1 - beta_n; %trial relative stress

% Check yield condition
normxtr = sqrt(   xtr_n1(1)*xtr_n1(1) +    xtr_n1(2)*xtr_n1(2) ...
             +    xtr_n1(3)*xtr_n1(3) + R2*xtr_n1(4)*xtr_n1(4) ...
             + R2*xtr_n1(5)*xtr_n1(5) + R2*xtr_n1(6)*xtr_n1(6));
ftr_n1 = normxtr - sqr23*(sigy + bet*Hb*a_n);

if ftr_n1 <= 1e-12
    
else
    
    % Plastic step
    n_n1 = xtr_n1/normxtr; %unit normal to yield surface
    dgama = ftr_n1/(R2G*(R1 + Hb/R3G)); %delta-gamma
    
    nij = [n_n1(1) n_n1(4) n_n1(6)
           n_n1(4) n_n1(2) n_n1(5)
           n_n1(6) n_n1(5) n_n1(3)];
       
    dnde = dnde(R2G,nij,normxtr);

    ddgde = ddgde(nij,Hb,R3G);

    d2ndede = d2ndede(dnde,normxtr,R2G,nij);

    d2dgdede = d2dgdede(dnde,Hb,R3G);
    
    D = -R2G*dgama*d2ndede;
    
    for I = 1:3
        for J = 1:3
            for K = 1:3
                for L = 1:3
                    for M = 1:3
                        for N = 1:3
          D(I,J,K,L,M,N) = D(I,J,K,L,M,N) + ...
                           -R2G*dnde(I,J,M,N)*ddgde(K,L) ...
                           -R2G*nij(I,J)*d2dgdede(K,L,M,N) ...
                           -R2G*dnde(I,J,K,L)*ddgde(M,N);
                        end
                    end
                end
            end
        end
    end
 
    
    
end

% Matrix form; diag will pull out the entries [11 22 33] from the 3x3
% submatrices in D.
D2 = [ diag(D(:,:,1,1,1,1)) diag(D(:,:,2,2,1,1)) diag(D(:,:,3,3,1,1))   diag(D(:,:,1,2,1,1)) diag(D(:,:,2,3,1,1)) diag(D(:,:,3,1,1,1))
      [diag(D(:,:,1,2,1,1)) diag(D(:,:,2,3,1,1)) diag(D(:,:,3,1,1,1))]' [D(1,2,1,2,1,1) D(1,2,2,3,1,1) D(1,2,3,1,1,1); D(2,3,1,2,1,1) D(2,3,2,3,1,1) D(2,3,3,1,1,1); D(3,1,1,2,1,1) D(3,1,2,3,1,1) D(3,1,3,1,1,1)]
       diag(D(:,:,1,1,2,2)) diag(D(:,:,2,2,2,2)) diag(D(:,:,3,3,2,2))   diag(D(:,:,1,2,2,2)) diag(D(:,:,2,3,2,2)) diag(D(:,:,3,1,2,2)) 
      [diag(D(:,:,1,2,2,2)) diag(D(:,:,2,3,2,2)) diag(D(:,:,3,1,2,2))]' [D(1,2,1,2,2,2) D(1,2,2,3,2,2) D(1,2,3,1,2,2); D(2,3,1,2,2,2) D(2,3,2,3,2,2) D(2,3,3,1,2,2); D(3,1,1,2,2,2) D(3,1,2,3,2,2) D(3,1,3,1,2,2)]
       diag(D(:,:,1,1,3,3)) diag(D(:,:,2,2,3,3)) diag(D(:,:,3,3,3,3))   diag(D(:,:,1,2,3,3)) diag(D(:,:,2,3,3,3)) diag(D(:,:,3,1,3,3)) 
      [diag(D(:,:,1,2,3,3)) diag(D(:,:,2,3,3,3)) diag(D(:,:,3,1,3,3))]' [D(1,2,1,2,3,3) D(1,2,2,3,3,3) D(1,2,3,1,3,3); D(2,3,1,2,3,3) D(2,3,2,3,3,3) D(2,3,3,1,3,3); D(3,1,1,2,3,3) D(3,1,2,3,3,3) D(3,1,3,1,3,3)]
       diag(D(:,:,1,1,1,2)) diag(D(:,:,2,2,1,2)) diag(D(:,:,3,3,1,2))   diag(D(:,:,1,2,1,2)) diag(D(:,:,2,3,1,2)) diag(D(:,:,3,1,1,2)) 
      [diag(D(:,:,1,2,1,2)) diag(D(:,:,2,3,1,2)) diag(D(:,:,3,1,1,2))]' [D(1,2,1,2,1,2) D(1,2,2,3,1,2) D(1,2,3,1,1,2); D(2,3,1,2,1,2) D(2,3,2,3,1,2) D(2,3,3,1,1,2); D(3,1,1,2,1,2) D(3,1,2,3,1,2) D(3,1,3,1,1,2)]
       diag(D(:,:,1,1,2,3)) diag(D(:,:,2,2,2,3)) diag(D(:,:,3,3,2,3))   diag(D(:,:,1,2,2,3)) diag(D(:,:,2,3,2,3)) diag(D(:,:,3,1,2,3)) 
      [diag(D(:,:,1,2,2,3)) diag(D(:,:,2,3,2,3)) diag(D(:,:,3,1,2,3))]' [D(1,2,1,2,2,3) D(1,2,2,3,2,3) D(1,2,3,1,2,3); D(2,3,1,2,2,3) D(2,3,2,3,2,3) D(2,3,3,1,2,3); D(3,1,1,2,2,3) D(3,1,2,3,2,3) D(3,1,3,1,2,3)]
       diag(D(:,:,1,1,3,1)) diag(D(:,:,2,2,3,1)) diag(D(:,:,3,3,3,1))   diag(D(:,:,1,2,3,1)) diag(D(:,:,2,3,3,1)) diag(D(:,:,3,1,3,1)) 
      [diag(D(:,:,1,2,3,1)) diag(D(:,:,2,3,3,1)) diag(D(:,:,3,1,3,1))]' [D(1,2,1,2,3,1) D(1,2,2,3,3,1) D(1,2,3,1,3,1); D(2,3,1,2,3,1) D(2,3,2,3,3,1) D(2,3,3,1,3,1); D(3,1,1,2,3,1) D(3,1,2,3,3,1) D(3,1,3,1,3,1)]];

end

function eye4 = eye4()
    %4th order symmetric identity tensor
    eye4 = zeros(3,3,3,3);
    eye4(:,:,1,1) = ...
     [1     0     0
     0     0     0
     0     0     0];
eye4(:,:,2,1) =...
     [                0   0.500000000000000                   0
   0.500000000000000                   0                   0
                   0                   0                   0];
eye4(:,:,3,1) =...
     [                0                   0   0.500000000000000
                   0                   0                   0
   0.500000000000000                   0                   0];
eye4(:,:,1,2) =...
     [                0   0.500000000000000                   0
   0.500000000000000                   0                   0
                   0                   0                   0];
eye4(:,:,2,2) =...
     [0     0     0
     0     1     0
     0     0     0];
eye4(:,:,3,2) =...
     [               0                   0                   0
                   0                   0   0.500000000000000
                   0   0.500000000000000                   0];
eye4(:,:,1,3) =...
     [                0                   0   0.500000000000000
                   0                   0                   0
   0.500000000000000                   0                   0];
eye4(:,:,2,3) =...
     [                0                   0                   0
                   0                   0   0.500000000000000
                   0   0.500000000000000                   0];
eye4(:,:,3,3) =...
     [0     0     0
     0     0     0
     0     0     1];
%     for I = 1:3
%         for J = 1:3
%             for K = 1:3
%                 for L = 1:3
%                     eye4(I,J,K,L) = CIJKL2(I,J,K,L);
%                 end
%             end
%         end
%     end
    
end

% function modul = CIJKL2(I,J,K,L)
% % Copied from NL Mixed Elasticity/Terms folder
% % dC_IJ^-1/dC_KL
% 
% C = eye(3);
% 
% modul = C(I,K)*C(J,L)+C(I,L)*C(J,K);
% modul = modul/2;
% 
% end

function ab = outer(a,b)
% Computes outer product of 2 second order tensors
    ab = zeros(3,3,3,3);
    for I = 1:3
        for J = 1:3
            ab(:,:,I,J) = a*b(I,J);
        end
    end
end

function dnde = dnde(R2G,n,normxtr)
% Derivative of flow rule normal tensor n_ij with respect to strain e_kl
    eye4 = eye4;
    dnde = R2G/normxtr*(eye4 - outer(n,n) - 1/3*outer(eye(3),eye(3)));
end

function ddgde = ddgde(n,Hb,R3G)
% Derivative of plastic multiplier dgamma w.r.t. strain e_kl
    ddgde = n/(1+Hb/R3G);
end

function d2ndede = d2ndede(dnde,normxtr,R2G,n)
% Second derivative of normal tensor n_ij w.r.t. strain e_kl, e_mn
    d2ndede = zeros(3,3,3,3,3,3);
    for I = 1:3
        for J = 1:3
            for K = 1:3
                for L = 1:3
                    for M = 1:3
                        for N = 1:3
    d2ndede(I,J,K,L,M,N) = -n(M,N)/normxtr*dnde(I,J,K,L) ...
                           -n(K,L)/normxtr*dnde(I,J,M,N) ...
                           -n(I,J)/normxtr*dnde(M,N,K,L);
                        end
                    end
                end
            end
        end
    end
    
    d2ndede = R2G*d2ndede;
end

function d2dgdede = d2dgdede(dnde,Hb,R3G)
% Second derivative of dgamma w.r.t. strain e_kl, e_mn
    d2dgdede = dnde/(1+Hb/R3G);
end