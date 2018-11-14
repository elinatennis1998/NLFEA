function [ nye ] = mm10_calc_nye_omar( ReT_grad93 )


% Levi-Civita Matrix:
lcMat = zeros(3,9);
lcMat([8 12 22]) = 1; lcMat([6 16 20]) = -1;
% Multply L-C Matrix by a 9x3 reshaped gradRe to get -curl
% (Nye tensor)
nye = -(lcMat * reshape(ReT_grad93([1:3:27,2:3:27,3:3:27]),9,3));


end

