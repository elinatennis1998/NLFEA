function [B,SizeB] = ShuffleIn(A,B,SizeA,SizeB)
%
%...  Written by Tim Truster (Fall 2009)
%     No modifications when copied to NLFEA ver2
%...  Program to combine two arrays, numerically sorted, and output in
%     array B
%
%**********************************************************************    

%       implicit none
% 
% %     Input Variables
%       integer SizeA,SizeB
% 	integer A(SizeA),B(SizeB)
% 
% %     Output Variables
% 
% %     Local Variables
%       integer i,j,k,l,ind(SizeA),elemA,elemB,indo,m
% 	logical dprt,imatch
    ind = zeros(SizeA,1);

% Determine index location of entries of A within array B
	k = 0;
	l = 1;
	for j = 1:SizeA
         elemA = A(j);
	
         [elemB,indo,imatch] = binsearch(elemA,B,l,SizeB);
	
         if(imatch)
            ind(j) = -1;
         elseif(indo==SizeB)
            ind(j) = indo + 1;
            k = k + 1;	        
         else
           ind(j) = indo;
           k = k + 1;
         end
	end

	i = SizeB + k;
	j = SizeA;
	l = SizeB;

	SizeB = i;

% Transfer A into B
      while(k>0)
%   Skip entries in A that were in B
         while(ind(j)==-1)
            j = j - 1;
         end
%   Shift entries of B
         while(l>ind(j))
	        B(i) = B(l);
	        i = i - 1;
	        l = l - 1;
         end
%   Transfer entry from A
	     B(i) = A(j);
	     i = i - 1;
	     j = j - 1;
	     k = k - 1;
      end