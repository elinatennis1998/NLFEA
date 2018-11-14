function [ basis ] = mm10_basis_omarf( g0, ReT_3_3ngp, basis )
ngp = size(basis.si,3);
for i =1:ngp
    ReT = ReT_3_3ngp(1:3,1+(i-1)*3:3+(i-1)*3);
    g = g0*ReT;
    f = 1/sqrt(2.0D0);
    
    basis.si(1, :, i)=g'*(f.*[-1  1  0]');
    basis.si(2, :, i)=g'*(f.*[ 1  0 -1]');
    basis.si(3, :, i)=g'*(f.*[ 0 -1  1]');
    basis.si(4, :, i)=g'*(f.*[-1 -1  0]');
    basis.si(5, :, i)=g'*(f.*[ 1  0  1]');
    basis.si(6, :, i)=g'*(f.*[ 0  1 -1]');
    basis.si(7, :, i)=g'*(f.*[ 1  1  0]');
    basis.si(8, :, i)=g'*(f.*[-1  0  1]');
    basis.si(9, :, i)=g'*(f.*[ 0 -1 -1]');
    basis.si(10,:, i)=g'*(f.*[ 1 -1  0]');
    basis.si(11,:, i)=g'*(f.*[-1  0 -1]');
    basis.si(12,:, i)=g'*(f.*[ 0  1  1]');
    basis.si(13,:, i)=g'*(f.*[ 1  1  0]');
    basis.si(14,:, i)=g'*(f.*[ 1  0  1]');
    basis.si(15,:, i)=g'*(f.*[ 0  1  1]');
    basis.si(16,:, i)=g'*(f.*[-1  1  0]');
    basis.si(17,:, i)=g'*(f.*[ 1  0 -1]');
    basis.si(18,:, i)=g'*(f.*[ 0 -1  1]');
    %     basis.si(19,:, i)=f.*[ 1  1  0];
    %     basis.si(20,:, i)=f.*[ 1  0  1];
    %     basis.si(21,:, i)=f.*[ 0  1  1];
    %     basis.si(22,:, i)=f.*[-1  1  0];
    %     basis.si(23,:, i)=f.*[ 1  0 -1];
    %     basis.si(24,:, i)=f.*[ 0 -1  1];
    
    f = 1/sqrt(3.0D0);
    
    basis.ni(1, :, i)=g'*(f.*[ 1  1  1]');
    basis.ni(2, :, i)=g'*(f.*[ 1  1  1]');
    basis.ni(3, :, i)=g'*(f.*[ 1  1  1]');
    basis.ni(4, :, i)=g'*(f.*[ 1 -1 -1]');
    basis.ni(5, :, i)=g'*(f.*[ 1 -1 -1]');
    basis.ni(6, :, i)=g'*(f.*[ 1 -1 -1]');
    basis.ni(7, :, i)=g'*(f.*[-1  1 -1]');
    basis.ni(8, :, i)=g'*(f.*[-1  1 -1]');
    basis.ni(9, :, i)=g'*(f.*[-1  1 -1]');
    basis.ni(10,:, i)=g'*(f.*[-1 -1  1]');
    basis.ni(11,:, i)=g'*(f.*[-1 -1  1]');
    basis.ni(12,:, i)=g'*(f.*[-1 -1  1]');
    basis.ni(13,:, i)=g'*(f.*[ 1 -1 -1]');
    basis.ni(14,:, i)=g'*(f.*[ 1 -1 -1]');
    basis.ni(15,:, i)=g'*(f.*[-1  1 -1]');
    basis.ni(16,:, i)=g'*(f.*[ 1  1  1]');
    basis.ni(17,:, i)=g'*(f.*[ 1  1  1]');
    basis.ni(18,:, i)=g'*(f.*[ 1  1  1]');
    %     basis.ni(19,:, i)=f.*[-1  1 -1];
    %     basis.ni(20,:, i)=f.*[-1 -1  1];
    %     basis.ni(21,:, i)=f.*[-1 -1  1];
    %     basis.ni(22,:, i)=f.*[-1 -1  1];
    %     basis.ni(23,:, i)=f.*[-1  1 -1];
    %     basis.ni(24,:, i)=f.*[ 1 -1 -1];
    
    for j =1:12
        basis.pi(j,:, i) = cross(basis.si(j,:, i),basis.ni(j,:, i));
    end
    %     basis.pi(13:24,:, i) = basis.si(13:24,:, i);
    basis.pi(13:18,:, i) = basis.si(13:18,:, i);
end
end

