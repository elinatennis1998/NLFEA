function [b,bd] = bubbleCSTT(ss,xs)

bd = zeros(1,3);
u=1.d0-ss(1)-ss(2)-ss(3);

% % First - works poorly for T4
% bd(1) = 256.d0*ss(2)*ss(3)*(u-ss(1));
% bd(2) = 256.d0*ss(1)*ss(3)*(u-ss(2));
% bd(3) = 256.d0*ss(1)*ss(2)*(u-ss(3));
% b = 256.d0*ss(1)*ss(2)*ss(3)*u;

% % Second - works worst
% bd(1) = 2*ss(1)*ss(2)^2*ss(3)^2*u^2 - 2*ss(1)^2*ss(2)^2*ss(3)^2*u;
% bd(2) = 2*ss(1)^2*ss(2)*ss(3)^2*u^2 - 2*ss(1)^2*ss(2)^2*ss(3)^2*u;
% bd(3) = 2*ss(1)^2*ss(2)^2*ss(3)*u^2 - 2*ss(1)^2*ss(2)^2*ss(3)^2*u;
% b = ss(1)^2*ss(2)^2*ss(3)^2*u^2;

% % Third - works best but is not zero on faces
%    shld = zeros(10,3);
%    shl = zeros(10,1);
% 
%     shld(5,1)=4.d0*ss(2);
%     shld(5,2)=4.d0*ss(1);
%     shld(5,3)=0.d0;
%     shl(5)=4.d0*ss(1)*ss(2);
% 
%     shld(9,1)=0.d0;
%     shld(9,2)=4.d0*ss(3);
%     shld(9,3)=4.d0*ss(2);
%     shl(9)=4.d0*ss(2)*ss(3);
% 
%      shld(8,1)=4.d0*ss(3);
%     shld(8,2)=0.d0;
%     shld(8,3)=4.d0*ss(1);
%     shl(8)=4.d0*ss(1)*ss(3);
% 
%     shld(7,1)=4.d0*(u-ss(1));
%     shld(7,2)=-4.d0*ss(1);
%     shld(7,3)=-4.d0*ss(1);
%     shl(7)=4.d0*ss(1)*u;
% 
%     shld(6,1)=-4.d0*ss(2);
%     shld(6,2)=4.d0*(u-ss(2));
%     shld(6,3)=-4.d0*ss(2);
%     shl(6)=4.d0*ss(2)*u;
% 
%     shld(10,1)=-4.d0*ss(3);
%     shld(10,2)=-4.d0*ss(3);
%     shld(10,3)=4.d0*(u-ss(3));
%     shl(10)=4.d0*ss(3)*u;
%     
% %     b = sum(shl);
% %     bd = sum(shld);
% 
% % Fourth - works poorly
%     ba = sum(shl);
%     bda = sum(shld);
% 
%     bdb = zeros(1,3);
%     bdb(1) = ss(2)*ss(3)*(u-ss(1));
%     bdb(2) = ss(1)*ss(3)*(u-ss(2));
%     bdb(3) = ss(1)*ss(2)*(u-ss(3));
%     bb = ss(1)*ss(2)*ss(3)*u;
% 
%     b = ba*bb + bb;
%     bd = ba*bdb + bb*bda + bdb;

% % Fifth - equivalent to first
%     bd(1) = 2*ss(1)*ss(2)*ss(3)*u-ss(1)^2*ss(2)*ss(3)+ss(2)^2*ss(3)*(u-ss(1))+ss(2)*ss(3)^2*(u-ss(1))+ss(2)*ss(3)*u^2-2*ss(1)*ss(2)*ss(3)*u;
%     bd(2) = ss(1)^2*ss(3)*(u-ss(2))+2*ss(1)*ss(2)*ss(3)*u-ss(1)*ss(2)^2*ss(3)+ss(1)*ss(3)^2*(u-ss(2))+ss(1)*ss(3)*u^2-2*ss(1)*ss(2)*ss(3)*u;
%     bd(3) = ss(1)^2*ss(2)*(u-ss(3))+ss(1)*ss(2)^2*(u-ss(3))+2*ss(1)*ss(2)*ss(3)*u-ss(1)*ss(2)*ss(3)^2+ss(1)*ss(2)*u^2-2*ss(1)*ss(2)*ss(3)*u;
%     b = ss(1)^2*ss(2)*ss(3)*u+ss(1)*ss(2)^2*ss(3)*u+ss(1)*ss(2)*ss(3)^2*u+ss(1)*ss(2)*ss(3)*u^2;

% Sixth - works better than first
    bd(1) = 3*ss(1)^2*ss(2)*ss(3)*u-ss(1)^3*ss(2)*ss(3)+ss(2)^3*ss(3)*(u-ss(1))+ss(2)*ss(3)^3*(u-ss(1))+ss(2)*ss(3)*u^3-3*ss(1)*ss(2)*ss(3)*u^2;
    bd(2) = ss(1)^3*ss(3)*(u-ss(2))+3*ss(1)*ss(2)^2*ss(3)*u-ss(1)*ss(2)^3*ss(3)+ss(1)*ss(3)^3*(u-ss(2))+ss(1)*ss(3)*u^3-3*ss(1)*ss(2)*ss(3)*u^2;
    bd(3) = ss(1)^3*ss(2)*(u-ss(3))+ss(1)*ss(2)^3*(u-ss(3))+3*ss(1)*ss(2)*ss(3)^2*u-ss(1)*ss(2)*ss(3)^3+ss(1)*ss(2)*u^3-3*ss(1)*ss(2)*ss(3)*u^2;
    b = ss(1)^3*ss(2)*ss(3)*u+ss(1)*ss(2)^3*ss(3)*u+ss(1)*ss(2)*ss(3)^3*u+ss(1)*ss(2)*ss(3)*u^3;

bd = bd*xs;