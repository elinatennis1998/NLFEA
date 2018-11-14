if ~(isfield(omar,'B'))
    o.a = sqrt(3)/9;
    o.c = sqrt(3)/84;
    o.d = 1/18;
    o.e = 3/14;
    
    omar.B = [  o.a,  -7*o.c,     o.c,   7*o.c,    -o.a,    -o.c, -13*o.c,  13*o.c,     0
               -o.a,    -o.c,   7*o.c,  13*o.c,       0, -13*o.c,  -7*o.c,     o.c,   o.a
                  0, -13*o.c,  13*o.c,     o.c,     o.a,  -7*o.c,    -o.c,   7*o.c,  -o.a 
                o.a,   7*o.c,    -o.c,  -7*o.c,    -o.a,    -o.c,  13*o.c,  13*o.c,     0
               -o.a,     o.c,  -7*o.c, -13*o.c,       0, -13*o.c,   7*o.c,     o.c,   o.a
                  0,  13*o.c, -13*o.c,    -o.c,     o.a,  -7*o.c,     o.c,   7*o.c,  -o.a 
                o.a,   7*o.c,     o.c,  -7*o.c,    -o.a,     o.c, -13*o.c, -13*o.c,     0
               -o.a,     o.c,   7*o.c, -13*o.c,       0,  13*o.c,  -7*o.c,    -o.c,   o.a
                  0,  13*o.c,  13*o.c,    -o.c,     o.a,   7*o.c,    -o.c,  -7*o.c,  -o.a 
                o.a,  -7*o.c,    -o.c,   7*o.c,    -o.a,     o.c,  13*o.c, -13*o.c,     0
               -o.a,    -o.c,  -7*o.c,  13*o.c,       0,  13*o.c,   7*o.c,    -o.c,   o.a
                  0, -13*o.c, -13*o.c,     o.c,     o.a,   7*o.c,     o.c,  -7*o.c,  -o.a
              5*o.d,     o.e,       0,     o.e,   5*o.d,       0,       0,       0,  -o.d
              5*o.d,       0,     o.e,       0,    -o.d,       0,     o.e,       0, 5*o.d
               -o.d,       0,       0,       0,   5*o.d,     o.e,       0,     o.e, 5*o.d
              5*o.d,    -o.e,       0,    -o.e,   5*o.d,       0,       0,       0,  -o.d
              5*o.d,       0,    -o.e,       0,    -o.d,       0,    -o.e,       0, 5*o.d
               -o.d,       0,       0,       0,   5*o.d,    -o.e,       0,    -o.e, 5*o.d ];
           
    clear('o')
end