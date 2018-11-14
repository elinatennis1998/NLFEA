% 04/30/2013

%Get BC
for i = 1:numBC
    node = NodeBC(i,1);
    dir = NodeBC(i,2);
    if dir > ndf
        errmsg = ['dof ID exceeds ndf for BC=' num2str(i)];
        error(errmsg)
    end
    displacement = NodeBC(i,3);
    idFEAP(dir,node) = -1; %#ok<*SAGROW>
    idFEAPBC(dir,node) = displacement;
end
for i = 1:numBCnp
    node = NodeBCnp(i,1);
    dir = NodeBCnp(i,2);
    if dir > ndf
        errmsg = ['dof ID exceeds ndf for BCnp=' num2str(i)];
        error(errmsg)
    end
    displacement = NodeBCnp(i,3);
    idFEAP(dir,node) = -1; %#ok<*SAGROW>
    idFEAPBCnp(dir,node) = displacement;
end

% Include proportional multiplier BCs
for i = 1:numBCmt
    len = size(NodeBCmt{2,i},1);
    for j = 1:len
        node = NodeBCmt{2,i}(j,1);
        dir = NodeBCmt{2,i}(j,2);
        if dir > ndf
            errmsg = ['dof ID exceeds ndf for BCmt=' num2str(i) 'entry' num2str(j)];
            error(errmsg)
        end
        idFEAP(dir,node) = -1; %#ok<*SAGROW>
    end
end