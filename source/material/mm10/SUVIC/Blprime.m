function Blp = Blprime(eps,props)

B0l = props.B0l;
m = props.m;
spds0 = sigmaprime(eps,props);
Blp = B0l*(spds0)^m;