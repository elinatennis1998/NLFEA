% formR for pressure-precipitation model
function work_vec1 = mm10_halite_formRpp( props, work_vec1, stress, tinc )

I_dev = [0.666666666666667,-0.333333333333333,-0.333333333333333,0,0,0;
    -0.333333333333333,0.666666666666667,-0.333333333333333,0,0,0;
    -0.333333333333333,-0.333333333333333,0.666666666666667,0,0,0;
    0,0,0,1,0,0;
    0,0,0,0,1,0;
    0,0,0,0,0,1];
temp = I_dev * stress';
eps_norton = 3.0/2.0 * props.cp_033 * tinc * temp;
eps_norton(4:6) = eps_norton(4:6) * 2.0;
work_vec1 = work_vec1 - eps_norton;

end