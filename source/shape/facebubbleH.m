function [b,db] = facebubbleH(ss)
% Tim Truster
% 03/05/2014
% face bubble function for stabilized DG Hexahedral elements
   one = 1.d0;
   two = 2.d0;
   db=zeros(3,1);
   r=ss(1);
   s=ss(2);
   t=ss(3);
   onemrsq = one-r*r;
   onemssq = one-s*s;
   bubble(1) =-two*r*onemssq*(1-t);
   bubble(2) =-two*s*onemrsq*(1-t);
   bubble(3) = -onemrsq*onemssq;
   b = onemrsq*onemssq*(1-t);
%    bubble(1) = 32.d0*r*s*s*(-1.d0+r*r)*(-1.d0+s*s)*(1-t)+32.d0*r*r*r*s*s*(-1.d0+s*s)*(1-t);
%    bubble(2) = 32.d0*r*r*s*(-1.d0+r*r)*(-1.d0+s*s)*(1-t)+32.d0*r*r*s*s*s*(-1.d0+r*r)*(1-t);
%    bubble(3) = -16.d0*r*r*s*s*(-1.d0+r*r)*(-1.d0+s*s);
%    b = 16.d0*r*r*s*s*(-1.d0+r*r)*(-1.d0+s*s)*(1-t);
   db(1)=bubble(1);
   db(2)=bubble(2);
   db(3)=bubble(3);