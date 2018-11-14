function [ n ] = normalizer( v )
%NORMALIZER Return 1/norm(v) or 0
%   The quantity to multiply a vector by to normalize it respecting the
%   zero vector
tol = 1.0e-15;
nv = norm(v);

if abs(nv) < tol
    n = 0.0;
else
    n = 1.0/nv;
end


end

