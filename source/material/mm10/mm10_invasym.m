% c
% c ****************************************************************************
% c *                                                                          *
% c *    mm10_invasym                                                               *
% c *                                                                          *
% c *         written by : mcm                                                 *
% c *         last modified : 3/22/12 mcm                                      *
% c *                                                                          *
% c *         Calculate the inversion of a non-symmetric matrix using LAPACK   *
% c *                                                                          *
% c ****************************************************************************
% c
function [A,n] = mm10_invasym(A,n)
%             implicit none
%             double precision, intent(inout), dimension(n,n) :: A
%             integer, intent(in) :: n
% 
%             integer :: i,j,info,lwork
%             integer, allocatable :: ipivt(:)
%             double precision, allocatable :: work(:)
% c
% c           Allocate storage
%             allocate(ipivt(n))
            lwork = n*n;
%             allocate(work(lwork))
% c           Factor
[~,~,A,n,ipivt,info] = DGETRF(n,n,A,n,ipivt,info); % http://tinyurl.com/nku66ys
% c           Inverse
[~,A,n,~,~,~,info] = DGETRI(n,A,n,ipivt,work,lwork,info); % http://tinyurl.com/pbpe8q5
% c           Free storage
%             deallocate(ipivt)
%             deallocate(work)   
% 
%             return
end