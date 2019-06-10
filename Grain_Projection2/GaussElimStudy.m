% Tim Truster
% 07/09/17
%
% Looking at inverse characteristics of tridiagonal matrices
n=8;
e = ones(n,1);
A = spdiags([e -2*e e], -1:1, n, n);
Ainv = full(inv(A));
[L,U] = lu(A);
Linv = full(inv(L));
Uinv = full(inv(U));
Ainv2 = Uinv*Linv;