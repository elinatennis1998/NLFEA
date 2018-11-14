function D = ComCijklCST2(mu,p)

mu = (mu-p);
D = zeros(4,4);
D(1,1) = 2*mu;
% D(1,2) = lam;
% D(2,1) = D(1,2);
D(2,2) = 2*mu;
D(3,3) = mu;

end

% function sum = sum1(F,ind1,ind2,ind3,ind4)
% 
%     sum = 0;
%     
%     for I = 1:2
%         for J = 1:2
%             sum = sum + F(ind1,I)*F(ind2,I)*F(ind3,J)*F(ind4,J);
%         end
%     end
% 
% end