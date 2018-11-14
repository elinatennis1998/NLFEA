function [NodeBC,BCLIndex] = getBCDG_2d(Coordinates)

% Function to assign BC to surface 5 of lap joint


NodeBC = [];
   j=1;
   for i =1:size(Coordinates,1)
      if Coordinates(i,2)==0
          NodeBC(j,:) = [i 1 0];
          j=j+1;
          NodeBC(j,:) = [i 2 0];
          j=j+1;

      end
   
   end
 numBC = size(NodeBC,1);
 BCLIndex = [numBC 0];   