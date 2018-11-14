if ~(isfield(omar,'basis'))
    omar.basis.nslip = 18;
    g = cc_props.g;
    f = 1/sqrt(2.0D0);
    
    omar.basis.si(1, :)=g*(f.*[-1  1  0]');
    omar.basis.si(2, :)=g*(f.*[ 1  0 -1]');
    omar.basis.si(3, :)=g*(f.*[ 0 -1  1]');
    omar.basis.si(4, :)=g*(f.*[-1 -1  0]');
    omar.basis.si(5, :)=g*(f.*[ 1  0  1]');
    omar.basis.si(6, :)=g*(f.*[ 0  1 -1]');
    omar.basis.si(7, :)=g*(f.*[ 1  1  0]');
    omar.basis.si(8, :)=g*(f.*[-1  0  1]');
    omar.basis.si(9, :)=g*(f.*[ 0 -1 -1]');
    omar.basis.si(10,:)=g*(f.*[ 1 -1  0]');
    omar.basis.si(11,:)=g*(f.*[-1  0 -1]');
    omar.basis.si(12,:)=g*(f.*[ 0  1  1]');
    omar.basis.si(13,:)=g*(f.*[ 1  1  0]');
    omar.basis.si(14,:)=g*(f.*[ 1  0  1]');
    omar.basis.si(15,:)=g*(f.*[ 0  1  1]');
    omar.basis.si(16,:)=g*(f.*[-1  1  0]');
    omar.basis.si(17,:)=g*(f.*[ 1  0 -1]');
    omar.basis.si(18,:)=g*(f.*[ 0 -1  1]');
%     omar.basis.si(19,:)=f.*[ 1  1  0];
%     omar.basis.si(20,:)=f.*[ 1  0  1];
%     omar.basis.si(21,:)=f.*[ 0  1  1];
%     omar.basis.si(22,:)=f.*[-1  1  0];
%     omar.basis.si(23,:)=f.*[ 1  0 -1];
%     omar.basis.si(24,:)=f.*[ 0 -1  1];
    
    f = 1/sqrt(3.0D0);
    
    omar.basis.ni(1, :)=g*(f.*[ 1  1  1]');
    omar.basis.ni(2, :)=g*(f.*[ 1  1  1]');
    omar.basis.ni(3, :)=g*(f.*[ 1  1  1]');
    omar.basis.ni(4, :)=g*(f.*[ 1 -1 -1]');
    omar.basis.ni(5, :)=g*(f.*[ 1 -1 -1]');
    omar.basis.ni(6, :)=g*(f.*[ 1 -1 -1]');
    omar.basis.ni(7, :)=g*(f.*[-1  1 -1]');
    omar.basis.ni(8, :)=g*(f.*[-1  1 -1]');
    omar.basis.ni(9, :)=g*(f.*[-1  1 -1]');
    omar.basis.ni(10,:)=g*(f.*[-1 -1  1]');
    omar.basis.ni(11,:)=g*(f.*[-1 -1  1]');
    omar.basis.ni(12,:)=g*(f.*[-1 -1  1]');
    omar.basis.ni(13,:)=g*(f.*[ 1 -1 -1]');
    omar.basis.ni(14,:)=g*(f.*[ 1 -1 -1]');
    omar.basis.ni(15,:)=g*(f.*[-1  1 -1]');
    omar.basis.ni(16,:)=g*(f.*[ 1  1  1]');
    omar.basis.ni(17,:)=g*(f.*[ 1  1  1]');
    omar.basis.ni(18,:)=g*(f.*[ 1  1  1]');
%     omar.basis.ni(19,:)=f.*[-1  1 -1];
%     omar.basis.ni(20,:)=f.*[-1 -1  1];
%     omar.basis.ni(21,:)=f.*[-1 -1  1];
%     omar.basis.ni(22,:)=f.*[-1 -1  1];
%     omar.basis.ni(23,:)=f.*[-1  1 -1];
%     omar.basis.ni(24,:)=f.*[ 1 -1 -1];
    
    for i =1:12
        omar.basis.pi(i,:) = cross(omar.basis.si(i,:),omar.basis.ni(i,:));
    end
%     omar.basis.pi(13:24,:) = omar.basis.si(13:24,:);
    omar.basis.pi(13:18,:) = omar.basis.si(13:18,:);
end