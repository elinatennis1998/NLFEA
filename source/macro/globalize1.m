% %**********************************************************************
% %
%       subroutine globalize(usi,globid,ldoft,us)
% %
% %...  Written by Tim Truster (Spring 2010)
%     No modifications when copied to NLFEA ver2
% %...  Program to get DOF IDs for cell of submesh
% %
% %**********************************************************************
% 
%       Implicit None
% 
% 	integer         numnp,numel,nummat,nen,neq,ipr
%       common /cdata/  numnp,numel,nummat,nen,neq,ipr
% 
%       logical debug
%       common /debugs/ debug
% 
% 	include 'smparam.h'
% 	include 'stardata.h'
% 	include 'starprob.h'
% 
% %     Input Variables
%       integer globid(numnpp),ldoft(ndfs,numnpp)
% 	real*8  usi(neqp),us(ndfs,numel*celn-numnp)
% 
% %     Output Variables
% 
% %     Local Variables
% 	integer	i,j,node,dof

    for i = 1:numnpp

	   node = globid(i);

       for j = 1:ndfs

	      dof = ldoft(j,i);

          if(dof<=neqp)
	
	         us(j,node) = us(j,node) + usi(dof); %#ok<AGROW>

          end

       end

    end