% c$Id:$
%       subroutine sethis(ie,ix,rben,nie,nen,nen1,
%      &                  numel,nummat,prt)
% 
% c      * * F E A P * * A Finite Element Analysis Program
% 
% c....  Copyright (c) 1984-2008: Regents of the University of California
% c                               All rights reserved
% 
% c-----[--.----+----.----+----.-----------------------------------------]
% c     Modification log                                Date (dd/mm/year)
% c       Original version                                    01/11/2006
% c       1.                Conversion to Matlab by TJT       19/04/2013
% c-----[--.----+----.----+----.-----------------------------------------]
% c      Purpose: Set up history addresses in ix array
% 
% c      Inputs:
% c         ie(nie,*) - Material set assembly information
% c         rben(*)   - List of rigid body to nodes
% c         nie       - Dimension of ie array
% c         nen       - Number of nodes/element
% c         nen1      - Dimension of ix array
% c         numel     - Number of elements in mesh
% c         nummat    - Number of material sets in mesh
% c         prt       - Flag, output results if true
% 
% c      Outputs:
% c         ix(nen1,*)- History data pointers added to positions nen+1,
% c                     nen+2 and nen+3
% c-----[--.----+----.----+----.-----------------------------------------]
%       implicit  none
% 
%       include  'hdatam.h'
%       include  'iofile.h'
% 
%       logical   flag,prt,setvar,palloc
%       integer   i,n,nh0,nhf,nie,nen,nen1,numel,nummat,ma
%       integer   ie(nie,*),ix(nen1,*),rben(*)
% 
%       save

%     Compute maximum length necessary to store history variables

      nhmax  = 0;
      nh3max = 0;
      for mg = 1:nummat
        nh0 = 0;
        nhf = 0;
        for ma = 1:nummat
          if(ieFEAP(nie-2,ma) == mg) %then
            ieFEAP(nie-3,ma) = nh0; %#ok<*SAGROW>
            ieFEAP(nie-4,ma) = nhf;
            nh0         = nh0 + ieFEAP(nie,  ma);
            nhf         = nhf + ieFEAP(nie-5,ma);
          end
        end % ma
        nhmax  = max(nhmax, ieFEAP(nie,mg));
        nh3max = max(nh3max,ieFEAP(nie-5,mg));
      end % mg

%     Set pointers for history variables into ix-array

      nh0 = 0;
      for elem = 1:numel

%       Set number of history terms for rigid elements to zero

%         if(rben(elem) <= 0) %then

          ma  = RegionOnElement(elem);
          ixFEAP(7,elem) = ma;

%         Variable storage history

          nhf = 0;
          for mg = 1:nummat
            if(ieFEAP(nie-2,mg) == ma) 
              nhf = nhf + ieFEAP(nie,ma);
            end
          end % mg
          if(nhf > 0) %then
            ixFEAP(1,elem) = nh0;
            nh0 = nh0 + nhf;
            ixFEAP(2,elem) = nh0;
            nh0 = nh0 + nhf;
            nhf = 0;
          end

%         Fixed storage history

          for mg = 1:nummat
            if(ieFEAP(nie-2,mg) == ma) 
              nhf = nhf + ieFEAP(nie-5,ma);
            end
          end % mg
          if(nhf > 0) %then
            ixFEAP(3,elem) = nh0;
            nh0 = nh0 + nhf;
          end
%         end
      end % elem
      nhf = nh0;
      hrvec = zeros(nhf,1);
      hr = zeros(2*nhmax+nh3max,1);
      if hrstore == 1
      hrmat = zeros(nhf,datastep+1); % use only for small meshes
      end
