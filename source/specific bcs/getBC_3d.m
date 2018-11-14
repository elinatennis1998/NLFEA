function [NodeBC] = getBC_3d(Coordinates)
j = 1;
k = 1;
for i =1:size(Coordinates,1)
   if Coordinates(i,1)==0
      if k==1
      NodeBC(j,:) = [i 1 0];
      j = j+1;
      NodeBC(j,:) = [i 2 0];
      j = j+1;
      NodeBC(j,:) = [i 3 0];
      j = j+1;
      elseif k ==2
      NodeBC(j,:) = [i 1 0];   
      j=j+1;
      NodeBC(j,:) = [i 3 0];   
      j=j+1;
      else
      NodeBC(j,:) = [i 1 0];   
      j=j+1;          
      end
   k = k+1; 
   end
end
