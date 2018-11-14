function [B,len] = PurgeList(A, SizeA, numcol)
    
    j = 1;
    epsilon = 1e-10;
    
    for i = 2:SizeA
        if abs(A(j, 1) - A(i, 1)) > epsilon
            j = j + 1;
            for k = 1:numcol
                A(j, k) = A(i, k);
            end
        end
    end
    len = j;
    B = zeros(len,numcol);
    for i = 1:len
        for j = 1:numcol
            B(i,j) = A(i,j);
        end
    end