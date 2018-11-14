% 04/30/2013

%Get BC
for i = 1:numBC2
    node = NodeBC2(i,1);
    dir = NodeBC2(i,2);
    displacement = NodeBC2(i,3);
    idFEAP2(dir,node) = -1; %#ok<*SAGROW>
    idFEAPBC2(dir,node) = displacement;
end
for i = 1:numBCnp2
    node = NodeBCnp2(i,1);
    dir = NodeBCnp2(i,2);
    displacement = NodeBCnp2(i,3);
    idFEAP2(dir,node) = -1; %#ok<*SAGROW>
    idFEAPBCnp2(dir,node) = displacement;
end