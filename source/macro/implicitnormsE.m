NodesOnElement = NodesOnElement';
ut = us(1:ndf,:)-up(1:ndf,:);


    
    ixc = zeros(nen,1);
	xc = zeros(ndm,celn);
    uc = zeros(ndf,nen);

    utl2 = zeros(ndf,1);
    utix = zeros(ndf,1);
    utiy = zeros(ndf,1);
	utieff = zeros(numel,1);
    
%....	set shape function flags
    ib = 0;
	der = 0;1;
	bf = 0;1;



%%
%-----------------------------------------------------
% Loop over elements in domain
%-----------------------------------------------------
   for elem = 1:numel

	   mat = RegionOnElement(elem);
       iel = MatTypeTable(2,mat);
       nonlin = MatTypeTable(3,mat);
       mateprop = MateT(mat,:);
         nele = cis(2,elem);
         if((nele==3)||(nele==6))
            neleB = 3;
         else
            neleB = 4;
         end
         lint = IntPoint(nele);
         
       ixe = zeros(nen,1);
       xe = zeros(ndm,nen);
%        ue = zeros(ndf,nen);
   
%     Extract element coordinates
            for j = 1:nele
               node = NodesOnElement(j,elem);
               ixe(j) = node;
               for i = 1:ndm
                  xe(i,j) = xs(i,node);
               end
%                for i = 1:ndf
%                   ue(i,j) = Node_U_V(node,i);
%                end
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
                     for i = 1:ndm
                        xc(i,j) = xs(i,node);
                     end
                     for i = 1:ndf
                        uc(i,j) = ut(i,node);
                     end
                  end

%       Set element-cell local coordinate map
                  evenodd = floor((celli+1)/2);
                  if(celli/2==evenodd) %even celli

	               %adjust order of cell flags for proper assembly
	               %due to reversal of element local coordinates
                   if(nele==3)
                      for i = 1:ndm
	                     temp = xc(i,1);
	                     xc(i,1) = xc(i,2);
	                     xc(i,2) = xc(i,3);
	                     xc(i,3) = temp;
                      end
                      for i = 1:ndfs
	                     temp = uc(i,1);
	                     uc(i,1) = uc(i,2);
	                     uc(i,2) = uc(i,3);
	                     uc(i,3) = temp;
                      end
                   else
                      for i = 1:ndm
	                     temp = xc(i,1);
	                     xc(i,1) = xc(i,2);
	                     xc(i,2) = xc(i,3);
	                     xc(i,3) = temp;
	                     temp = xc(i,4);
	                     xc(i,4) = xc(i,5);
	                     xc(i,5) = xc(i,6);
	                     xc(i,6) = temp;
                      end
                      for i = 1:ndfs
	                     temp = uc(i,1);
	                     uc(i,1) = uc(i,2);
	                     uc(i,2) = uc(i,3);
	                     uc(i,3) = temp;
	                     temp = uc(i,4);
	                     uc(i,4) = uc(i,5);
	                     uc(i,5) = uc(i,6);
	                     uc(i,6) = temp;
                      end
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

                    starcellerrorE

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
                     for i = 1:ndm
                        xc(i,j) = xs(i,node);
                     end
                     for i = 1:ndf
                        uc(i,j) = ut(i,node);
                     end
                  end

%       Set element-cell local coordinate map
                  %MR = 1.d0/dble(m) 
                  MR = 1.d0/minc;
                  BR = -1.d0+(2.d0*celli-1.d0)/minc;
                  MS = 1.d0/minc;
                  BS = -1.d0+(2.d0*cellj-1.d0)/minc;

                    starcellerrorE

%    -----------------------------------------------------
%     End loop over cells
%    -----------------------------------------------------
                  end
               end
            end

%-----------------------------------------------------
% End loop over elements in domain
%-----------------------------------------------------
    end

	ul2 = sqrt(utl2(1)+utl2(2));
	uh1 = sqrt(utix(1)+utix(2)+utiy(1)+utiy(2));

%
%....	calculate the log
%
      dl10  = log(10.d0);
      ul2  = log(ul2) / dl10;
      uh1  = log(uh1) / dl10;
        fprintf('Total u  %i  %1.7e  %1.7e\n',numel,ul2,uh1)
        
    
NodesOnElement = NodesOnElement';
