% c
% c ****************************************************************************
% c *                                                                          *
% c *    mm10_invsym                                                                *
% c *                                                                          *
% c *         written by : mcm                                                 *
% c *         last modified : 3/22/12 mcm                                      *
% c *                                                                          *
% c *         Calculate the inversion of a symmetric matrix using LAPACK       *
% c *                                                                          *
% c ****************************************************************************
% c
function [A, n] = mm10_invsym(A, n)
%             implicit none
%             double precision, intent(inout), dimension(n,n) :: A
%             integer, intent(in) :: n
% 
%             integer :: i,j,info,lwork
%             integer, allocatable :: ipiv(:)
%             double precision, allocatable :: work(:)
% c
% c           Allocate storage
%             allocate(ipiv(n))
%              lwork = n*n
%             allocate(work(lwork))
% c           Factor
            ['U',~,A,n,ipiv,work,~,info] = DSYTRF('U',n,A,n,ipiv,...
                work,lwork,info); %DSYTRF computes the factorization 
% of a real symmetric matrix A using the Bunch-Kaufman diagonal pivoting method.  The form of the
% factorization is A = U*D*U**T  or  A = L*D*L**T
% c           Inverse
            [~,~,A,n,~,~,info] = DSYTRI(U,n,A,n,ipiv,...
                work,info); % DSYTRI computes the inverse of a real symmetric indefinite matrix
% A using the factorization A = U*D*U**T or A = L*D*L**T computed by DSYTRF.
% c           Sym -> Full
            for i = 1:n;
                  for j = 1:i - 1;
                        A(i,j) = A(j,i);
                  end
            end
% c           Free storage
%             deallocate(ipiv)
%             deallocate(work)
% 
%             return
end