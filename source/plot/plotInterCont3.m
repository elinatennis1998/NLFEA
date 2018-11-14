function plotInterCont3
% Placeholder for a 3d plotting function of interfaces
% 04/04/2015
% I could envision making a version just like the 2d one, which plots
% element-wise contours. Those could easily come from elemental stresses,
% isw=26, just like the 2d version.
% For better plots, a node-based method would be needed.
% From the EMI paper in CM journal, Interface4 folder, there are functions
% called plotdbondCont2.m and colorflip2beta.m which made the 3d cylinder
% plots for the composite interface. That would be a place to start looking
% for 3d ideas, namely to make a wire mesh by using the elements in
% SurfacesI and the proper faces.
% Currently, NL_FEA_Program has isw=60,61 reserved for making interface
% variables with FormI. A 2D example that works(?) is NL_Elem21_2d_2.m.
% NL_Elem3_3d_bubble has something but it is not current. NL_Elem21_3d_7.m
% has more in it and looks current, with InterQuant being used.
% So, the best ideas are to use isw=26 or isw=60 to get quantities around
% the interface(s) in order to plot.