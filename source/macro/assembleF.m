% %**********************************************************************
% %
%       subroutine assembleF(nel,nen,ndf,edoft,eF,Fd1)
% %
% %...  Written by Tim Truster (Fall 2009)
%     No modifications when copied to NLFEA ver2
% %...  Program to assemble force vector for cell into global vector
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
% 	real*8  eF(nen*ndf),Fd1(neqp)
% 
% %     Output Variables
% 
% %     Local Variables
% 	integer	locind1,grow

    for locind1 = 1:nele*ndfs

         grow = edoft(locind1);

         if(grow<=neqp)

            FD1(grow) = FD1(grow) + eF(locind1); %#ok<AGROW>

         end

    end

