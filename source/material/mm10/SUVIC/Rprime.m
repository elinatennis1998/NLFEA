function Rp = Rprime(eps,props)

R0 = props.R0;
m = props.m;
spds0 = sigmaprime(eps,props);
Rp = R0*(spds0)^m;