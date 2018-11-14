function t = tranr4(fl,fr)

% c      * * F E A P * * A Finite Element Analysis Program
% 
% c....  Copyright (c) 1984-2008: Regents of the University of California
% c                               All rights reserved
% 
% c-----[--.----+----.----+----.-----------------------------------------]
% c     Modification log                                Date (dd/mm/year)
% c       Original version                                    01/11/2006
% c-----[--.----+----.----+----.-----------------------------------------]
% c      Purpose:   Transformation array for 4th rank tensor
% c                 t(a,b) = fl(i,I)*fr(j,J) : a -> I,J ; b -> i,j
% c                   a,b  |  1    2    3    4    5    6
% c                  ------+-----------------------------
% c                  (I,J) | 1,1  2,2  3,3  1,2  2,3  3,1
% c               or (i,j) |                2,1  3,2  1,3
% 
% c     Input:
% c       fl(3,3) - left  deformation gradient
% c       fr(3,3) - right deformation gradient
% c     Output:
% c       t(6,6) - transformation array
% c-----[--.----+----.----+----.-----------------------------------------]
%       implicit  none

%       integer   i,j, i1(6),i2(6)
%       real*8    fl(3,3),fr(3,3), t(6,6)
% 
%       data      i1 /1,2,3,1,2,3/
%       data      i2 /1,2,3,2,3,1/
      t = zeros(6,6);

%       i1 = [1,2,3,1,2,3];
%       i2 = [1,2,3,2,3,1];
      i113 = [1,2,3];
      i213 = [1,2,3];
      i146 = [1,2,3];
      i246 = [2,3,1];

% c     Form transformation array for a 4th rank tensor in matrix form

%       for i = 1:3
%         for j = 1:3
%           t(i,j) =  fl(i1(j),i1(i))*fr(i2(j),i2(i));
%         end %do ! j
%         for j = 4:6
%           t(i,j) = (fl(i1(j),i1(i))*fr(i2(j),i2(i)) ...
%                  +  fl(i2(j),i2(i))*fr(i1(j),i1(i)))*0.5d0;
%         end %do ! j
%       end %do ! i
%         t(1:3,1:3) = fl(i1(1:3),i1(1:3))'.*fr(i2(1:3),i2(1:3))';
%         t(1:3,4:6) = (fl(i1(4:6),i1(1:3))'.*fr(i2(4:6),i2(1:3))' ...
%                    +  fl(i2(4:6),i2(1:3))'.*fr(i1(4:6),i1(1:3))')*0.5d0;
        t(1:3,1:3) = fl(i113,i113)'.*fr(i213,i213)';
        t(1:3,4:6) = (fl(i146,i113)'.*fr(i246,i213)' ...
                   +  fl(i246,i213)'.*fr(i146,i113)')*0.5d0;

%       for i = 4:6
%         for j = 1:3
%           t(i,j) =  fl(i1(j),i1(i))*fr(i2(j),i2(i)) ...
%                  +  fl(i2(j),i2(i))*fr(i1(j),i1(i));
%         end %do ! j
%         for j = 4:6
%           t(i,j) = (fl(i1(j),i1(i))*fr(i2(j),i2(i)) ...
%                  +  fl(i2(j),i1(i))*fr(i1(j),i2(i)) ...
%                  +  fl(i1(j),i2(i))*fr(i2(j),i1(i)) ...
%                  +  fl(i2(j),i2(i))*fr(i1(j),i1(i)))*0.5d0;
%         end %do ! j
%       end %do ! i
        t(4:6,1:3) = fl(i113,i146)'.*fr(i213,i246)' ...
                  +  fl(i213,i246)'.*fr(i113,i146)';
        t(4:6,4:6) = (fl(i146,i146)'.*fr(i246,i246)' ...
                   +  fl(i246,i146)'.*fr(i146,i246)' ...
                   +  fl(i146,i246)'.*fr(i246,i146)' ...
                   +  fl(i246,i246)'.*fr(i146,i146)')*0.5d0;

      end