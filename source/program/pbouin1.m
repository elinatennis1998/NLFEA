% 04/30/2013

%Get BC
for i = 1:numBC1
    node = NodeBC1(i,1);
    dir = NodeBC1(i,2);
    displacement = NodeBC1(i,3);
    idFEAP1(dir,node) = -1; %#ok<*SAGROW>
    idFEAPBC1(dir,node) = displacement;
end
for i = 1:numBCnp1
    node = NodeBCnp1(i,1);
    dir = NodeBCnp1(i,2);
    displacement = NodeBCnp1(i,3);
    idFEAP1(dir,node) = -1; %#ok<*SAGROW>
    idFEAPBCnp1(dir,node) = displacement;
end