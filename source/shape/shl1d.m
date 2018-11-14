function [basislist,basisder,basisder2,bubble,basisder3] = shl1d(ulist,unum,p)
%
% Tim Truster
% 5/2009
% UIUC
% copied to NLFEAver2 on 7/8/2013
%
% Function to evaluate 1-D Lagrange Basis Functions of degree p at 
% integration points (ulist). Basis functions are ordered for
% tensor-product formation of the shape functions in 2-D or 3-D.
%
% For polynomials of 4th or higher order, the functions are evaluated by
% applying the following formulas:
%
%         p+1     (u - xi_b)
% Na(u) = prod   -------------
%         b =  1 (xi_a - xi_b)
%         b /= a
%
%            p+1      ( p+1       (u - xi_c)   )
% Na_xi(u) = sum      ( prod     ------------- )
%            b = 1    ( c =  1   (xi_a - xi_c) )
%            b /= a   ( c /= a,b               )
%
% numertable(a,b) contains the prod. of the numerators (u - xi_c) c /= a,b

% SOMEDAY, CONVERT TO VECTOR OPERATIONS USING .
% Also need to fix p>3 derivative calculation

basislist = zeros(unum,(p+1));
basisder = basislist;
basisder2 = basislist;
basisder3 = basislist;
bubble = zeros(unum,3);

switch p
    case 0
        for i = 1:unum
            u = ulist(i);

            basislist(i,1) = 1;
            
            basisder(i,1) = 0;
            
            bubble(i,1) = 1-u^2;
            bubble(i,2) = -2*u;
            bubble(i,3) = -2;
        end %i
    case 1
        for i = 1:unum
            u = ulist(i);

            basislist(i,1) = (1 - u)/2;
            basislist(i,2) = (1 + u)/2;
            
            basisder(i,1) = -1/2;
            basisder(i,2) = 1/2;
            
            bubble(i,1) = 1-u^2;
            bubble(i,2) = -2*u;
            bubble(i,3) = -2;
        end %i
    case 2
        for i = 1:unum
            u = ulist(i);

            basislist(i,1) = u*(u - 1)/2;
            basislist(i,2) = 1-u^2;
            basislist(i,3) = u*(u + 1)/2;
            
            basisder(i,1) = u - 1/2;
            basisder(i,2) = -2*u;
            basisder(i,3) = u + 1/2;
            
            basisder2(i,1) = 1;
            basisder2(i,2) = -2;
            basisder2(i,3) = 1;
            
            
        end %i
    case 3
        for i = 1:unum
            u = ulist(i);

            basislist(i,1) =  1/16*(9*u^2 - 1)*(1 - u);
            basislist(i,2) =  9/16*(u^2 - 1)*(3*u - 1);
            basislist(i,3) =  9/16*(1 - u^2)*(3*u + 1);
            basislist(i,4) =  1/16*(9*u^2 - 1)*(u + 1);
            
            basisder(i,1) =  1/16*(1 + u*(18 - 27*u));
            basisder(i,2) =  9/16*(-3 - u*(2 - 9*u));
            basisder(i,3) =  9/16*(3 - u*(2 + 9*u));
            basisder(i,4) =  1/16*(-1 + u*(18 + 27*u));
            
            basisder2(i,1) =  1/16*((18 - 54*u));
            basisder2(i,2) =  9/16*(-(2 - 18*u));
            basisder2(i,3) =  9/16*(-(2 + 18*u));
            basisder2(i,4) =  1/16*((18 + 54*u));
            
            basisder3(i,1) =  1/16*((- 54));
            basisder3(i,2) =  9/16*(-( - 18));
            basisder3(i,3) =  9/16*(-( + 18));
            basisder3(i,4) =  1/16*(( + 54));
            
            
        end %i
    otherwise % Does NOT give correct values for derivatives
        %Set xi_a list
        xi = (-1:2/p:1);
        %Form numertable, denominators
        for i = 1:unum
            u = ulist(i);
            numertable = ones(p+1,p+1);
            denom = ones(p+1,1);
            %Loop over xi_a
            for a = 1:p+1
                numeraa = 1;
                denoma = 1;
                xia = xi(a);
                %Compute product (u - xi_b) b /= a
                for c = 1:a-1
                    xic = xi(c);
                    numeraa = numeraa*(u - xic);
                    denoma = denoma*(xia - xic);
                end
                for c = a+1:p+1
                    xic = xi(c);
                    numeraa = numeraa*(u - xic);
                    denoma = denoma*(xia - xic);
                end
                %Compute product (u - xi_c) c /= a,b
                for b = a+1:p+1
                    numerab = 1;
                    for c = 1:a-1
                        xic = xi(c);
                        numerab = numerab*(u - xic);
                    end
                    for c = a+1:b-1
                        xic = xi(c);
                        numerab = numerab*(u - xic);
                    end
                    for c = b+1:p+1
                        xic = xi(c);
                        numerab = numerab*(u - xic);
                    end
                    numertable(a,b) = numerab;
                    numertable(b,a) = numerab;
                end
                %Compute basis function Na
                basislist(i,a) = numeraa/denoma;
                denom(a) = denoma;
            end
            %Compute first derivatives Na_xi
            for a = 1:p+1
                %Compute sum [ prod (u - xi_c) ]_b, b /= a
                numer1 = 0;
                for b = 1:a-1
                    numer1 = numer1 + numertable(b,a);
                end
                for b = a+1:p+1
                    numer1 = numer1 + numertable(b,a);
                end
                basisder(i,a) = numer1/denoma;
            end
            %Compute second derivatives Na_xixi
        end
end