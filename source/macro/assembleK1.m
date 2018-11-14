% %**********************************************************************
% %
%       subroutine assembleK(nel,nen,ndf,edoft,eK,Kdd11)
% %
% %...  Written by Tim Truster (Spring 2010)
%     No modifications when copied to NLFEA ver2
% %...  Program to assmeble stiffness matrix for cell into global matrix
% %
% %**********************************************************************
% 
%       Implicit None
% 
% 	include 'starprob.h'
% 
%       logical debug
%       common /debugs/ debug
% 
% %     Input Variables
%       integer nel,nen,ndf,edoft(ndf*nen)
%       real*8  eK(nen*ndf,nen*ndf),Kdd11(neqp,neqp)
% 
% %     Output Variables
% 
% %     Local Variables
% 	integer	locind1,locind2,grow,gcol

    for locind1 = 1:nele*ndfs

         grow = edoft(locind1);

         if(grow<=neqp)

            %Fd1(grow) = Fd1(grow) + ElemF(locind1)

            for locind2 = 1:nele*ndfs

               gcol = edoft(locind2);

               if(gcol<=neqp)
                  KDD11(grow,gcol) = KDD11(grow,gcol) + eK(locind1,locind2); %#ok<AGROW>
               end %gcol

            end

         end

    end

