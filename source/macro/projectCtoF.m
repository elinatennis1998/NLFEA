
% Tim Truster
% 11/7/2010
%     No modifications when copied to NLFEA ver2
%
% Routine to project coarse scale field onto fine scale field by evaluating
% the field at all the fine scale nodes; it loops over coarse elements and
% stores the coarse field at each of the fine nodes in the cells within
% that element.

ixe = zeros(nen,1);
ixc = zeros(nen,1);
xe = zeros(ndm,nen);
ue = zeros(ndf,nen);
uc = zeros(ndf,nen);
up = zeros(ndf,numnps);

proj36 = [0 1 0 .5 .5 0
          0 0 1 0 .5 .5];
proj49 = [-1 1 1 -1 0 1 0 -1 0
          -1 -1 1 1 -1 0 1 0 0];

NodesOnElement = NodesOnElement';
x = Coordinates';
      
for elem = 1:numel
    
    nele = cis(2,elem);
    if((nele==3)||(nele==6))
        neleB = 3;
    else
        neleB = 4;
    end
    
%     Extract element coordinates
    for j = 1:nele
       node = NodesOnElement(j,elem);
       ixe(j) = node;
       for i = 1:ndm
          xe(i,j) = x(i,node);
       end
       for i = 1:ndf
          ue(i,j) = u(i,node);
       end
    end  
    
%    -----------------------------------------------------
%     Loop over cells in element
%    -----------------------------------------------------
%     If Triangular Element
            if(neleB==3)

               celle = 0;

               for cellj = 1:minc
                  for celli = 1:2*(minc-cellj)+1
                 
                  celle = celle + 1;
                
%       Load cell coordinates, dof IDs

                  cell = cel*(elem-1) + celle;
                  ixc = zeros(nen,1);
   
                  for j = 1:nele
                     node = ixs(j,cell);
                     ixc(j) = node;
                  end

%       Set element-cell local coordinate map
                  evenodd = floor((celli+1)/2);
                  if(celli/2==evenodd) %even celli

	               %adjust order of cell flags for proper assembly
	               %due to reversal of element local coordinates
                   if(nele==3)
	                  i = ixc(1);
	                  ixc(1) = ixc(2);
	                  ixc(2) = ixc(3);
	                  ixc(3) = i;
	               else
	                  i = ixc(1);
	                  ixc(1) = ixc(2);
	                  ixc(2) = ixc(3);
	                  ixc(3) = i;
	                  i = ixc(4);
	                  ixc(4) = ixc(5);
	                  ixc(5) = ixc(6);
	                  ixc(6) = i;
                   end
                     MR = -1.d0/minc;
                     BR = celli/(2.d0*minc);
                     MS = -1.d0/minc;
                     BS = cellj/(minc);
                  else
                     MR = 1.d0/minc;
                     BR = ((celli+1)/2.d0-1.d0)/minc;
                     MS = 1.d0/minc;
                     BS = (cellj-1.d0)/minc;
                  end

%      -----------------------------------------------------
%       Project coarse quantities
%      -----------------------------------------------------
         for l=1:nele

            litr = proj36(1,l);
            lits = proj36(2,l);
            
            bigR = MR*litr + BR;
            bigS = MS*lits + BS;

%             if(nele==3||nele==6)
              shel = shlt(bigR,bigS,nele,nele,0,0);
%             else
%               shel = shlq(bigR,bigS,nele,0,0);
%             end

            un = zeros(ndf,1);

            for j = 1:nele
                for i = 1:ndf
                    un(i) = un(i) + shel(j)*ue(i,j);
                end
            end
            
            for i = 1:ndf
               uc(i,l) = un(i);
            end
            
         end
         
         
          for j = 1:nele
             node = ixc(j);
             for i = 1:ndf
                up(i,node) = uc(i,j);
             end
          end

%    -----------------------------------------------------
%     End loop over cells
%    -----------------------------------------------------
                  end
               end
%%
%     Else Quadrilateral Element
            else
              
               celle = 0;
               for cellj = 1:minc
                  for celli = 1:minc
                 
                  celle = celle + 1;

%       Load cell coordinates, dof IDs

                  cell = cel*(elem-1) + celle;
                  ixc = zeros(nen,1);
   
                  for j = 1:nele
                     node = ixs(j,cell);
                     ixc(j) = node;
                  end

%       Set element-cell local coordinate map
                  %MR = 1.d0/dble(m) 
                  MR = 1.d0/minc;
                  BR = -1.d0+(2.d0*celli-1.d0)/minc;
                  MS = 1.d0/minc;
                  BS = -1.d0+(2.d0*cellj-1.d0)/minc;

%      -----------------------------------------------------
%       Project coarse quantities
%      -----------------------------------------------------
         for l=1:nele

            litr = proj49(1,l);
            lits = proj49(2,l);
            
            bigR = MR*litr + BR;
            bigS = MS*lits + BS;

%             if(nele==3||nele==6)
%               shel = shlt(bigR,bigS,nele,0,0);
%             else
              shel = shlq(bigR,bigS,nele,nele,0,0);
%             end

            un = zeros(ndf,1);

            for j = 1:nele
                for i = 1:ndf
                    un(i) = un(i) + shel(j)*ue(i,j);
                end
            end
            
            for i = 1:ndf
               uc(i,l) = un(i);
            end
            
         end
         
         
          for j = 1:nele
             node = ixc(j);
             for i = 1:ndf
                up(i,node) = uc(i,j);
             end
          end

%    -----------------------------------------------------
%     End loop over cells
%    -----------------------------------------------------
                  end
               end

            end
    
end

NodesOnElement = NodesOnElement';