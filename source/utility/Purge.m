function [C,SizeC,A,SizeA] = Purge(A,B,SizeA,SizeB)
%
%...  Written by Tim Truster (Fall 2009)
%     No modifications when copied to NLFEA ver2
%...  Program to combine two arrays, numerically sorted, and output in
%     array B; entries common to A and B are removed from A
%
%**********************************************************************

%       implicit none
% 
% %     Input Variables
%       integer SizeA,SizeB
% 	integer A(SizeA),B(SizeB)
% 
% %     Output Variables
%       integer SizeC,C(*)
% 
% %     Local Variables
%       integer i,j,k,l,Ai
% 	logical dprt
    C = zeros(SizeA+SizeB,1);

	i = 1;
      j = 1;
	k = 1;
	l = 1;
    
      while((i<=SizeA)&&(j<=SizeB))
        if(B(j)<A(i))
            C(k) = B(j);
            j = j + 1;
            k = k + 1;
        elseif(A(i)<B(j))
	      Ai = A(i);
            C(k) = Ai;
	      A(l) = Ai;
            i = i + 1;
	      l = l + 1;
            k = k + 1;
	  else
            i = i + 1;
        end
      end

      if(i>SizeA)
        while(j<=SizeB)
            C(k) = B(j);
            j = j + 1;
	      k = k + 1;
        end
      else
       while(i<=SizeA)
	      Ai = A(i);
            C(k) = Ai;
	      A(l) = Ai;
            i = i + 1;
	      l = l + 1;
	      k = k + 1;
       end
      end

	SizeA = l - 1;
	SizeC = k - 1;