% %**********************************************************************
% %
%       subroutine assembleKF(nel,nen,ndf,edoft,eK,eF,Kdd11,Fd1)
% %
% %...  Written by Tim Truster (Fall 2009)
%     No modifications when copied to NLFEA ver2
% %...  Program to assemble stiffness matrix and force vector for cell
% %     into the global matrix and vector
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
% 	real*8  eF(nen*ndf),Fd1(neqp)
% 
% %     Output Variables
% 
% %     Local Variables
% 	integer	locind1,locind2,grow,gcol

    for locind1 = 1:nele*ndfs

         grow = edoft(locind1);

         if(grow<=neqp)

            FD1(grow) = FD1(grow) + eF(locind1); %#ok<AGROW>

            for locind2 = 1:nele*ndfs

               gcol = edoft(locind2);

               if(gcol<=neqp)
                  KDD11(grow,gcol) = KDD11(grow,gcol) + eK(locind1,locind2); %#ok<AGROW>
               end %gcol

            end

         end

    end

