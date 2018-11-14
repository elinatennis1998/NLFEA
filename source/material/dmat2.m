function [dmat]=dmat2(J,mateprop,lam,p)
 %d_ijklmn term
        I1 = [1; 1; 0];
        cpmat1 = I1*I1';
        dpmat1 = [cpmat1; cpmat1; zeros(3,3)];
        dpmat2 =-[6 2 0 
                  2 2 0 
                  0 0 1 
                  2 2 0 
                  2 6 0 
                  0 0 1 
                  0 0 1 
                  0 0 1 
                  1 1 0 ];
        dpmat3 = [8 0 0 
                  0 0 0 
                  0 0 2 
                  0 0 0 
                  0 8 0 
                  0 0 2 
                  0 0 2 
                  0 0 2 
                  2 2 0 ];
        dpmat0 = dpmat1 + dpmat2 + dpmat3;
if mateprop(1) > 0
    [dmati] = dmat2i(mateprop);
    [theta1,theta2,theta3] = ThetaNS(J,mateprop);
    A = J*theta1 + 3*J^2*theta2 + J^3*theta3;
    B = J*theta1 + J^2*theta2;
    C = J*theta1;
    dmatp = p*J*dpmat0;
    dmat = dmatp-dmati;
end