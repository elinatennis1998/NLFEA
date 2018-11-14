% c$Id:$
%       subroutine plink(id,ntyp,ndf,numnp,neq,prt)
% 
% c      * * F E A P * * A Finite Element Analysis Program
% 
% c....  Copyright (c) 1984-2008: Regents of the University of California
% c                               All rights reserved
% 
% c-----[--.----+----.----+----.-----------------------------------------]
% c     Modification log                                Date (dd/mm/year)
% c       Original version                                    01/11/2006
% c       1. Move 'go to 40' from line 103 to line 101        21/02/2007
% c          Modify line with 'sign' to check for error combinations.
% c       2  Allow up to 48 dof/node for links                26/06/2007
% c       3. Remove the go to 40 statement after endif.       29/06/2007
% c          Change 'mode1.lt.node1' and 'mode2.lt.node2' to 'le' in
% c          the 'go to 20' expression.
% c       4. Return checks in 3. to 'lt'; add 'exit' so that  08/07/2007
% c          last non-zero input record is processed.
% c       5.                Conversion to Matlab by TJT       29/04/2013
% c       6.                Optimized, switching for -> :     30/04/2013
% c-----[--.----+----.----+----.-----------------------------------------]
% c      Purpose: Link degrees of freedom to have same solution value
% 
% c      Inputs:
% c         id(ndf,*) - Equation numbers before link
% c         ntyp(*)   - Node type: >= 0 exist; < 0 not active
% c         ndf       - Number dof/node
% c         numnp     - Number of nodes in mesh
% c         prt       - Output links performed if true
% 
% c      Outputs:
% c         id(ndf,*) - Equation numbers for each dof after link
% c         neq       - Number of equations active after link
% c-----[--.----+----.----+----.-----------------------------------------]
%       implicit  none
% 
%       include  'cdat2.h'
%       include  'comfil.h'
%       include  'conval.h'
%       include  'iodata.h'
%       include  'iofile.h'
%       include  'ioincl.h'
%       include  'part0.h'
%       include  'trdata.h'
% 
%       logical   prt,lsave,errck, pinput, oprt,prth, exit
%       character fnamr*132,fext*4, type*4
%       integer   ndf, numnp, neq, iosfil
%       integer   i, ii, i1,i2, j, j1,j2
%       integer   mode1,mode2, node1,node2, nmax, id(ndf,*),idl(6),jdl(12),ntyp(*)
%       real*8    td(48)
% 
%       save

%     Routine to link degrees of freedom together

      mode1 = 0;
      ncom = 0;
      dontexit = 1;
      while dontexit
        ncom = ncom + 1;
        if ncom <= numComp
        node1 = NodeComp(ncom,1);
        node2 = NodeComp(ncom,2);
        i1 = NodeComp(ncom,3);
        i2 = NodeComp(ncom,4);
%         for i = 1:ndf
%           idl(i) = NodeComp(ncom,i+4);
%         end % i
        idl = NodeComp(ncom,5:4+ndf);
        if((node1 == 0) || (node1 > numnp) || (node2 == 0) || (node2 > numnp)) %then
          ncom,mode1,mode2 %#ok<NOPTS>
          error('linking nodes failed because of values out of range')
        end
        else
        node1 = 0;
        node2 = 0;
        i1 = 0;
        i2 = 0;
        dontexit = 0;
        end
        if(mode1 > 0) %then
          linking = 1;
          while linking
%           if(ntyp(mode1).lt.0 .or. ntyp(mode2).lt.0) then
%             if(prt) then
%               write(iow,3001) mode1,mode2
%               if(iosfil.lt.0) write(*,3001) mode1,mode2
%             endif
%           elseif(mode1.eq.mode2) then
%             if(prt) then
%               write(iow,3002) mode1,mode2
%               if(iosfil.lt.0) write(*,3002) mode1,mode2
%             endif
%           else
%             if(prt) then
%               write(iow,2001) mode1,mode2,(jdl(i),i=1,ndf)
%               if(iosfil.lt.0) then
%                 write(*,2001) mode1,mode2,(jdl(i),i=1,ndf)
%               endif
%             endif

%           Check that node pair has not already linked d.o.f.

            for j = 1:ndf
%               if((ndfp(j) == npart) && (jdl(j) == 0)) %then
              if((jdl(j) == 0)) %then
                if((idFEAP(j,mode1) >0) && (idFEAP(j,mode2) > 0)) %then

%                 Select node to renumber dof

                  if(idFEAP(j,mode1) == idFEAP(j,mode2)) %then
                    disp('nodes already tied')
                    mode1,mode2,j
                    continue
                  elseif(idFEAP(j,mode1) < idFEAP(j,mode2)) %then
                    nmax     = idFEAP(j,mode2);
                    idFEAP(j,mode2) = idFEAP(j,mode1); %#ok<*SAGROW>
                  else
                    nmax     = idFEAP(j,mode1);
                    idFEAP(j,mode1) = idFEAP(j,mode2);
                  end
%                   for ii = 1:numnp
%                     if(idFEAP(j,ii) == nmax) %then
%                       idFEAP(j,ii) = idFEAP(j,mode1);
%                     end
%                   end % ii
                  idFEAP(j,(idFEAP(j,ii) == nmax)) = idFEAP(j,mode1);

%                 Loop through all nodes to reduce equation numbers

%                   errck = 0;
%                   for i = 1:ndf
% %                     if(ndfp(i) == npart) %then
%                       for ii = 1:numnp
%                         if(idFEAP(i,ii) > nmax) %then
%                           idFEAP(i,ii) = idFEAP(i,ii) - 1;
%                           errck    = 1;
%                         end
%                       end % ii
% %                     end
%                   end % i
                  dofs = find(idFEAP>nmax);
                  errck = ~isempty(dofs);
                  idFEAP(dofs) = idFEAP(dofs) - 1;
%                   idFEAP(idFEAP>nmax) = idFEAP(idFEAP>nmax) - 1;
                  if(errck)
                    neq = neq - 1;
                  end
                else

%                 Attempt to link restrained dof

%                   write(ilg,3003) mode1,mode2,j
%                   write(iow,3003) mode1,mode2,j
%                   if(ior.lt.0) then
%                     write(*,3003) mode1,mode2,j
%                     go to 40
%                   end
                  mode1,mode2,j %#ok<NOPTS>
                  error('Attempt to link restrained dof')

                end
              end

            end % j
%             if(exit) 
%               go to 40
%             end
          mode1 = mode1 + j1;
          mode2 = mode2 + j2;
          if ~( (j1 > 0 && mode1 < node1) || (j2 > 0 && mode2 < node2) )
              linking = 0;
          end
          end

%         Check for error combinations

          if(node1 > 0) %then
            j1 =  1;
          else
            j1 = -1;
          end
          if(node2 > 0) %then
            j2 =  1;
          else
            j2 = -1;
          end
          if(j1*j2 <= 0) 
              dontexit = 0;
          end
        end
        mode1 = node1;
        mode2 = node2;
        j1 = i1;
        j2 = i2;
%         for i = 1:ndf
%           jdl(i) = idl(i);
%         end % i
        jdl = idl;

      end

%     Check that the number of equations is correct

%       neq = 0;
%       for i = 1:ndf
% %         if(ndfp(i) == npart) %then
%           for ii = 1:numnp
%             neq = max(neq,idFEAP(i,ii));
%           end % ii
% %         end
%       end % i
      neq = max(max(idFEAP));

