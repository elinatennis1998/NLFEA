function [NodeBC,NodesOnElement,RegionOnElement,nummat,MateT,MatTypeTable,numBC,numel,ccell] = Meso_Locking_patch2(NodesOnElement,num_locked_g,NodeBC,locked_g,....
    meso_nen,grainG,nen_bulk,numelemg,GrainA,numel,numnpMicro,numnpMeso,MateT,numgrain,nummat,MatTypeTable,RegionOnElement,RegionsOnInterface)
%Elina Geut
%Created 3/4/2019
%Last Modified 4/16/2019
%Reassigning BC with locked grains

if num_locked_g == 1 %Checks how many grains are being locked
     %Calculations id only one grain is locked
    
     micro = unique(reshape(NodesOnElement(grainG(locked_g,:),1:nen_bulk),numelemg*nen_bulk,1));
     l = length(micro);
     r_g = numgrain-locked_g;
     NodeBC = [NodeBC
         [(numnpMicro+1:numnpMicro+(locked_g-1)*meso_nen)' ones((locked_g-1)*meso_nen,1) zeros((locked_g-1)*meso_nen,1)
         (numnpMicro+1:numnpMicro+(locked_g-1)*meso_nen)' 2*ones((locked_g-1)*meso_nen,1) zeros((locked_g-1)*meso_nen,1)]
         [(numnpMicro+locked_g*meso_nen+1:numnpMicro+numnpMeso)'   ones(r_g*meso_nen,1) zeros(r_g*meso_nen,1)
         (numnpMicro+locked_g*meso_nen+1:numnpMicro+numnpMeso)' 2*ones(r_g*meso_nen,1) zeros(r_g*meso_nen,1)]
         [micro ones(l,1) zeros(l,1)
         micro 2*ones(l,1) zeros(l,1)]];
     
     % Correction for grain (i) of mesoscale
     index = find(MateT(:,1) == locked_g);
     MateT(index(1),5:8) = [1 0 0 1];
     MateT(index(2),5:8) = [1 0 0 1];
     index2 = find(MateT(:,2) == locked_g);
     MateT(index2(1),5:8) = [0 1 1 0];
     MateT(index2(2),5:8) = [0 1 1 0];
     l_MT = length(MateT)+1;
     MateT(l_MT,1:8) = [MateT(locked_g,1:3) GrainA 0 0 0 0];
     %Modification for CF subroutine
     %4/16/2019
     MatTypeTable(1:3,l_MT) = [l_MT 11 0]';
     %End of modification 4/16/2019
     numel = numel + 1;
     %Recised up until this point
     NodesOnElement(numel,1:12) = [numnpMicro+(locked_g-1)*meso_nen+1 numnpMicro+(locked_g-1)*meso_nen+2 numnpMicro+(locked_g-1)*meso_nen+3 zeros(1,12-3)];
     RegionOnElement(numel) = l_MT;
     nummat = nummat + 1;
     ccell = numel; %List of meso elements
else  %Calculations for a case where more than 1 grain is locked
    
    %Data initialization/Node,Material reassignment/Zeroing out micro scale
    %and replacing by meso scale
        for i = 1:num_locked_g
%             micro = unique(reshape(NodesOnElement(grainG(locked_g(i),:),1:nen_bulk),numelemg*nen_bulk,1)); %Corrected
           micro = unique(reshape(NodesOnElement(grainG(locked_g(i),:),1:nen_bulk),numelemg*nen_bulk,1)); %Corrected
            l(i) = length(micro);
            r_g = numgrain-locked_g(i);
            NodeBC = [NodeBC
                [micro ones(l(i),1) zeros(l(i),1)
                micro 2*ones(l(i),1) zeros(l(i),1)]]; %Zero out elements ina grain
            
            index1 = find(MateT(:,1) == locked_g(i));
            %The MateT arrangement is [cR cL fR fL]
            %Value of 1st column 
            
            for k = 1:2
                if MateT(index1(k),1) == locked_g(i)
                    flag_meso_top = 1;
                    flag_micro_top = 0;
                else
                    flag_meso_top = 0;
                    flag_micro_top = 1;
                end
                
                if ismember(MateT(index1(k),2),locked_g) == 1
                    flag_meso_bot = 1;
                    flag_micro_bot = 0;
                else
                    flag_meso_bot = 0;
                    flag_micro_bot = 1;
                end
                    
                 MateT(index1(k),5:8) = [flag_meso_top flag_meso_bot flag_micro_top flag_micro_bot];
            end
            
            %Value of 2nd column
            
            index2 = find(MateT(:,2) == locked_g(i));
            for j = 1:2
                
                if MateT(index2(j),2) == locked_g(i)
                    flag_meso_l = 1;
                    flag_micro_l = 0;
                else
                    flag_meso_l = 0;
                    flag_micro_l = 1;
                end
                
                if ismember(MateT(index2(j),1),locked_g) == 1
                    flag_meso_r = 1;
                    flag_micro_r = 0;
                else
                    flag_meso_r = 0;
                    flag_micro_r = 1;
                end
                MateT(index2(j),5:8) = [flag_meso_r flag_meso_l flag_micro_r flag_micro_l];
                MateT(index2(j),5:8) = [0 1 1 0];
            end
            
            l_MT = length(MateT)+1;
            MateT(l_MT,1:8) = [MateT(locked_g(i),1:3) GrainA 0 0 0 0];
            MatTypeTable(1:3,l_MT) = [l_MT 11 0]';
            BoundGrain = find(sum(MateT(:,1) == MateT(:,2),2) >= 1);
            for count = 1:length(BoundGrain)
                DBCposition = BoundGrain(count);
                MateT(DBCposition,5:8) = [1 1 1 0];
            end
            numel = numel + 1;
            NodesOnElement(numel,1:12) = [numnpMicro+(locked_g(i)-1)*meso_nen+1 numnpMicro+(locked_g(i)-1)*meso_nen+2 ...
                numnpMicro+(locked_g(i)-1)*meso_nen+3 zeros(1,12-3)];
            RegionOnElement(numel) = l_MT;
            nummat = nummat + 1;
            ccell(i) = numel; %List of meso elements
        end
%         Assigning NodeBC for the first grain 
% The set up procedure is written by hand in a notebook
            NodeBC = [NodeBC 
            [(numnpMicro+1:numnpMicro+(locked_g(1)-1)*meso_nen)' ones((locked_g(1)-1)*meso_nen,1) zeros((locked_g(1)-1)*meso_nen,1)
            (numnpMicro+1:numnpMicro+(locked_g(1)-1)*meso_nen)' 2*ones((locked_g(1)-1)*meso_nen,1) zeros((locked_g(1)-1)*meso_nen,1)]];
%          NodeBC for all grains other than first and last 
% r_g number of grains between the two locked grains 

            for n = 1:num_locked_g-1
                r_g(n) = locked_g(n+1)-(locked_g(n)+1);
                NodeBC = [NodeBC
                    [(numnpMicro+locked_g(n)*meso_nen+1:numnpMicro+(locked_g(n+1)-1)*meso_nen)'   ones(r_g(n)*meso_nen,1) zeros(r_g(n)*meso_nen,1)
                    (numnpMicro+locked_g(n)*meso_nen+1:numnpMicro+(locked_g(n+1)-1)*meso_nen)' 2*ones(r_g(n)*meso_nen,1) zeros(r_g(n)*meso_nen,1)]];
            end
%         Assigning NodeBC for the last grain

        r_g = numgrain-locked_g(num_locked_g);
        NodeBC = [NodeBC
            [(numnpMicro+locked_g(num_locked_g)*meso_nen+1:numnpMicro+numnpMeso)'  ones(r_g*meso_nen,1) zeros(r_g*meso_nen,1)
            (numnpMicro+locked_g(num_locked_g)*meso_nen+1:numnpMicro+numnpMeso)' 2*ones(r_g*meso_nen,1) zeros(r_g*meso_nen,1)]];
        
end
   numBC = length(NodeBC); %Reassign number of boudary conditions
   
   %4/25/2019
   %Loop added to material property for CS flag 
   if exist('FSon','var') && FSon == 1
        for n = 1:num_locked_g
           MateT(locked_g,9) = 2; %Flag = 2 to turn on both CS and FS
        end
   else
       for n = 1:num_locked_g
           MateT(locked_g,9) = 1; %Flag = 1 to turn on CS only
       end
   end