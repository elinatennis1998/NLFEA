% plot the stability terms ep based on the interface
close all
         ep_11 = squeeze(ep_list(5,1,:));
         ep_12 = squeeze(ep_list(5,2,:));
         ep_21 = squeeze(ep_list(5,3,:));  
         ep_22 = squeeze(ep_list(5,4,:));    
         ep_11max= max(ep_11)
         ep_11min = min(ep_11)  
         ep_11maxmin = max(ep_11)-min(ep_11);
         ep_12max= max(ep_12)
         ep_12min = min(ep_12) 
         ep_12maxmin = max(ep_12)-min(ep_12);  
         ep_21max= max(ep_21)
         ep_21min = min(ep_21)  
         ep_21maxmin = max(ep_21)-min(ep_21);  
         ep_22max= max(ep_22)
         ep_22min = min(ep_22)  
         ep_22maxmin = max(ep_22)-min(ep_22);  
for ii_inter = 1:numCL
         ep_norm(ii_inter) = sqrt(ep_11(ii_inter)^2+ep_12(ii_inter)^2+ep_21(ii_inter)^2+ep_22(ii_inter)^2);

end
         ep_norm_max = max(ep_norm)
         ep_norm_min = min(ep_norm)
         ep_norm_maxmin = ep_norm_max-ep_norm_min;
for ii_inter=1:numCL
    ep_color_11 = (ep_11(ii_inter)-ep_11min)./ep_11maxmin;
    ep_color_12 = (ep_12(ii_inter)-ep_12min)./ep_12maxmin;    
    ep_color_21 = (ep_21(ii_inter)-ep_21min)./ep_21maxmin;    
    ep_color_22 = (ep_22(ii_inter)-ep_22min)./ep_22maxmin;    
    ep_color_norm = (ep_norm(ii_inter)-ep_norm_min)./ep_norm_maxmin;        
    figure(1) % ep_11
    hold on
    plot([Coordinates(SurfacesI(ii_inter,1),1),Coordinates(SurfacesI(ii_inter,2),1)],[Coordinates(SurfacesI(ii_inter,1),2),Coordinates(SurfacesI(ii_inter,2),2)],'Color',[ep_color_11,0,0],'LineWidth',4)
    hold off
    figure(2)  %ep_12
    hold on
    plot([Coordinates(SurfacesI(ii_inter,1),1),Coordinates(SurfacesI(ii_inter,2),1)],[Coordinates(SurfacesI(ii_inter,1),2),Coordinates(SurfacesI(ii_inter,2),2)],'Color',[ep_color_12,0,0],'LineWidth',4)
    hold off
    figure(3)   %ep_21
    hold on
    plot([Coordinates(SurfacesI(ii_inter,1),1),Coordinates(SurfacesI(ii_inter,2),1)],[Coordinates(SurfacesI(ii_inter,1),2),Coordinates(SurfacesI(ii_inter,2),2)],'Color',[ep_color_21,0,0],'LineWidth',4)
    hold off
    figure(4)   %ep_22
    hold on
    plot([Coordinates(SurfacesI(ii_inter,1),1),Coordinates(SurfacesI(ii_inter,2),1)],[Coordinates(SurfacesI(ii_inter,1),2),Coordinates(SurfacesI(ii_inter,2),2)],'Color',[ep_color_22,0,0],'LineWidth',4)
    hold off
    figure(5)   %ep_norm
    hold on
    plot([Coordinates(SurfacesI(ii_inter,1),1),Coordinates(SurfacesI(ii_inter,2),1)],[Coordinates(SurfacesI(ii_inter,1),2),Coordinates(SurfacesI(ii_inter,2),2)],'Color',[ep_color_norm,0,0],'LineWidth',4)
    hold off
end