function f = getf(b,c,K0,deltad,deltal,ADOFN)

%HW2 Q1 Arc-Length - Modified Newton-Ralphson
%Tim Truster
%10/08/2008
%CEE 576 UIUC
%Prof. Masud

%Modified 06/29/2009

%This routine computes f for this problem

% f = 0;
% for i = 1:ADOFN
%     f = f + c*deltad(i)*K0(i)*deltad(i);
% end
f = c*sum(deltad(1:ADOFN).*K0(1:ADOFN).*deltad(1:ADOFN));
f = f + b*deltal^2;
f = sqrt(f);