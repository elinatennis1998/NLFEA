function [A, SizeA] = BubbleSortArray(A, SizeA, numcol, col)
%
%...  Written by Tim Truster (Fall 2009)
%     No modifications when copied to NLFEA ver2
%...  Program to numerically sort a small array
%
%**********************************************************************

%       implicit none
% 
% %     Input Variables
%       integer SizeA,numcol,col
% 	integer A(SizeA,numcol)
% 
% %     Output Variables
% 
% %     Local Variables
%       integer i,j,k,temp,PassNum
    
      PassNum = SizeA - 1;
      for i = 1:PassNum
        for j = 1:SizeA-i
            if (A(j,col)>A(j+1,col))
                for k = 1:numcol
                    temp = A(j,k);
                    A(j,k) = A(j+1,k);
                    A(j+1,k) = temp;
                end
            end
        end
      end
    
