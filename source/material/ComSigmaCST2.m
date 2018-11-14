function sigma = ComSigmaCST2(F,mu,p)

sigma = zeros(3,1);

sigma(1) = mu*(sum2(F,1,1) - 1)+p;
sigma(2) = mu*(sum2(F,2,2) - 1)+p;
sigma(3) = mu*sum2(F,1,2);

end

function sum = sum2(F,ind1,ind2)

    sum = 0;
    
    for I = 1:2
        sum = sum + F(ind1,I)*F(ind2,I);
    end

end