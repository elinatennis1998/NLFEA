function sigma = ComSigma09(J,F,mu,lam)

sigma = zeros(3,1);

sigma(1) = lam/J*log(J) + mu/J*(sum2(F,1,1) - 1);
sigma(2) = lam/J*log(J) + mu/J*(sum2(F,2,2) - 1);
sigma(3) = mu/J*sum2(F,1,2);

end

function sum = sum2(F,ind1,ind2)

    sum = 0;
    
    for I = 1:2
        sum = sum + F(ind1,I)*F(ind2,I);
    end

end