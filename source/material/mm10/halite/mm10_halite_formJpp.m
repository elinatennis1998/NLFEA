% formJ for pressure-precipitation model
function J11 = mm10_halite_formJpp( props, J11 , tinc )

I_dev = [0.666666666666667,-0.333333333333333,-0.333333333333333,0,0,0;
    -0.333333333333333,0.666666666666667,-0.333333333333333,0,0,0;
    -0.333333333333333,-0.333333333333333,0.666666666666667,0,0,0;
    0,0,0,1,0,0;
    0,0,0,0,1,0;
    0,0,0,0,0,1];
temp = props.stiffness * I_dev;
temp(4:6, 1:6) = temp(4:6, 1:6) * 2.0;
% J11 = J11 - 3.0/2.0*props.cp_033*temp*tinc;
J11 = J11 + 3.0/2.0*props.cp_033*temp*tinc;

end