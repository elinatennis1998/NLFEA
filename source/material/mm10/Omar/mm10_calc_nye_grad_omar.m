function [ nye_grad39 ] = mm10_calc_nye_grad_omar( gradgrad_ReT_96 )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

trnsfm = [1,28,46,4,31,49,7,34,52,...
          2,29,47,5,32,50,8,35,53,...
          3,30,48,6,33,51,9,36,54,...
          28,10,37,31,13,40,34,16,43,...
          29,11,38,32,14,41,35,17,44,...
          30,12,39,33,15,42,36,18,45,...
          46,37,19,49,40,22,52,43,25,...
          47,38,20,50,41,23,53,44,26,...
          48,39,21,51,42,24,54,45,27];

% gradgrad_ReT_96 = ReT_98 * shp2;
temp = reshape(gradgrad_ReT_96,54,1);
gradgrad_ReT_99 = temp(trnsfm);
gradgrad_ReT_99 = reshape(gradgrad_ReT_99,9,9);

% Levi-Civita Matrix:
lcMat = zeros(3,9);
lcMat([8 12 22]) = 1; lcMat([6 16 20]) = -1;

nye_grad39 = -(lcMat * gradgrad_ReT_99 );

end

