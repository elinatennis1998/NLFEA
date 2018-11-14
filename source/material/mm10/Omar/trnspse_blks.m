function [ B ] = trnspse_blks( A )
%Function to transpose 3x3 blocked matrix
%
% eg: [a11 a12 a13 b11 b12 b13          [a11 a12 a13
%      a21 a22 a23 b21 b22 b23    ==>    a21 a22 a23
%      a31 a32 a33 b31 b32 b33]          a31 a32 a33
%                                        b11 b12 b13
%                                        b21 b22 b23
%                                        b31 b32 b33]
% Variables:
% A: Blocked matrix
x = size(A,1);
y = size(A,2);
num_mtx_x = x/3;
num_mtx_y = y/3;
if x>y
    B = zeros(size(A))';
    for i =1:num_mtx_x
        for j = 1:num_mtx_y
            B((1:3)+3*(j-1),(1:3)+3*(i-1)) = A((1:3)+3*(i-1),(1:3)+3*(j-1));
        end
    end
else
    B = zeros(size(A))';
    for i =1:num_mtx_y
        for j = 1:num_mtx_x
            B((1:3)+3*(i-1),(1:3)+3*(j-1)) = A((1:3)+3*(j-1),(1:3)+3*(i-1));
        end
    end
end
end

