function NodeDisp = getDisp_2d(Coordinates,theta)
   NodeDisp = [];
   j=1;
   for i =1:size(Coordinates,1)
      if Coordinates(i,2)>=0.9999999
          u2 = theta*(Coordinates(i,1)-0.25)*(Coordinates(i,1)-0.75);
          NodeDisp(j,:) = [i 1 0];
           j=j+1;
          NodeDisp(j,:) = [i 2 u2];
          j=j+1;

      end
   
   end