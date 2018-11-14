function [NodeDisp]=getDisp_3d(Coordinates,L,alpha,lamda)
   j=1;
for i =1:size(Coordinates,1)
   x1 = (L./(alpha*lamda)+Coordinates(i,3))*sin(alpha*lamda)-Coordinates(i,1);
   x3 = Coordinates(i,3) + (L./(alpha*lamda)+Coordinates(i,3))*(cos(alpha*lamda)-1)-Coordinates(i,3);
   if Coordinates(i,1)==L
      NodeDisp(j,:) = [i 1 x1];
      j=j+1;
      NodeDisp(j,:) = [i 3 x3];
      j=j+1;
   end
end