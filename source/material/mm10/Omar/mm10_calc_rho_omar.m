function [ rho_grad18_24, rho_M18_8 ] = mm10_calc_rho_omar( g, ReT_grad93,...
    ReT_M3_24, ngp, nye, B, b, nye_grad)

ReT_grad39_blks = reshape(ReT_grad93,3,9);

%% Get Q = g(ReT) and Q transpose 
Q.grad93_blks = trnspse_blks(g*ReT_grad39_blks);
QT.grad39_blks = trnspse_blks(trnspse_blks(...
    reshape(ReT_grad93([1,4,7,2,5,8,3,6,9],:),3,9))*g');

Q.M24_3 = trnspse_blks(g*ReT_M3_24);
Q.M24_24 = zeros(3*ngp,3*ngp);
for i =1:ngp
    Q.M24_24(1+(i-1)*3:3+(i-1)*3,1+(i-1)*3:3+(i-1)*3)= Q.M24_3(1+(i-1)*3:3+(i-1)*3,1:3);
end
QT.M3_24 = (Q.M24_3)';
    
RHS = ((Q.grad93_blks * nye * QT.M3_24) + ...
    trnspse_blks(Q.M24_3 * nye * QT.grad39_blks));
if (nargin == 8)
    RHS = RHS + (trnspse_blks ( Q.M24_3 * nye_grad) * QT.M3_24);
end
RHS = reshape(RHS,27,ngp);
%% Compute grad(rho)
rho_grad18_24 = ((B * reshape(RHS ([1:3,10:12,19:21, 4:6,13:15,22:24,...
    7:9,16:18,25:27],:), 9,3*ngp)) ./ b)* Q.M24_24;

%% rho computations
    rho.temp = zeros(3,3*ngp);
    for i =1:ngp
        rho.temp(:,(1:3)+3*(i-1)) = Q.M24_3((1:3)+3*(i-1),:) *...
            nye * QT.M3_24(:,(1:3)+3*(i-1));
    end
    rho.temp = reshape (rho.temp,9,ngp);
    
    rho_M18_8 = (B * rho.temp) ./ b;
end

