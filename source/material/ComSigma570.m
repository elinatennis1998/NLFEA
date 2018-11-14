function sigma = ComSigma570(F,mu,lam)

sigma = zeros(3,1);
epskk = F(1,1) + F(2,2);
sigma(1) = lam*epskk + mu*2*F(1,1);
sigma(2) = lam*epskk + mu*2*F(2,2);
sigma(3) = mu*(F(1,2)+F(2,1));

end