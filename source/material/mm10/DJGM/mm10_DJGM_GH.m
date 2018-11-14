function [h_0, gamma_dot_tilde, g_tilde, r, n, m, g_0, q] = mm10_DJGM_GH(props)

% material parameters for titanium alpha phase slipping
% taken from Deka, Joseph, 2006 paper
    
if props.k_0 == -100 % User supplied
    
    if props.s_type == 9
        h_0 = [ones(3,1)*props.c1;ones(3,1)*props.c2];
        gamma_dot_tilde = [ones(3,1)*props.c4;ones(3,1)*props.c5];
        g_tilde = [ones(3,1)*props.c7;ones(3,1)*props.c8];
        r = [ones(3,1)*props.Qbulk;ones(3,1)*props.p_y];
        n = [ones(3,1)*props.q_y;ones(3,1)*props.q_v];
        m = [ones(3,1)*props.burgers;ones(3,1)*props.tau_a];
        g_0 = [ones(3,1)*props.G_0_y;ones(3,1)*props.eps_dot_0_y];
        q = [1 1 1 1 1 1;
             1 1 1 1 1 1;
             1 1 1 1 1 1;
             1 1 1 1 1 1;
             1 1 1 1 1 1;
             1 1 1 1 1 1;];
    elseif props.s_type == 10
        h_0 = [ones(3,1)*props.c1;ones(3,1)*props.c2;ones(12,1)*props.c3];
        gamma_dot_tilde = [ones(3,1)*props.c4;ones(3,1)*props.c5;ones(12,1)*props.c6];
        g_tilde = [ones(3,1)*props.c7;ones(3,1)*props.c8;ones(12,1)*props.Qslip];
        r = [ones(3,1)*props.Qbulk;ones(3,1)*props.p_y;ones(12,1)*props.p_v];
        n = [ones(3,1)*props.q_y;ones(3,1)*props.q_v;ones(12,1)*props.boltz];
        m = [ones(3,1)*props.b;ones(3,1)*props.tau_a;ones(12,1)*props.tau_hat_v];
        g_0 = [ones(3,1)*props.G_0_y;ones(3,1)*props.eps_dot_0_y;ones(12,1)*props.G_0_v];

        q = ones(18,18);
    else
        error('Unvalid slip type! Please try again!');
    end
else
if props.s_type == 9
    if props.theta_0 == -1
    h_0 = -1e6*[2000;2000;2000;4153;4153;4153];
    gamma_dot_tilde = [0.0023;0.0023;0.0023;0.0023;0.0023;0.0023];
    g_tilde = 1e6*[504;504;504;504;504;504]/2;
    r = [0.3;0.3;0.3;0.29;0.29;0.29];
    n = [0.14;0.14;0.14;0.15;0.15;0.15];
    m = [0.02;0.02;0.02;0.02;0.02;0.02];
    g_0 = 1e6*[300;300;300;240;240;240];
    q = [1 1 1 1 1 1;
         1 1 1 1 1 1;
         1 1 1 1 1 1;
         1 1 1 1 1 1;
         1 1 1 1 1 1;
         1 1 1 1 1 1;];
    else
    h_0 = 1e6*[2000;2000;2000;4153;4153;4153];
    gamma_dot_tilde = [0.0023;0.0023;0.0023;0.0023;0.0023;0.0023];
    g_tilde = 1e6*[504;504;504;504;504;504];
    r = [0.3;0.3;0.3;0.29;0.29;0.29];
    n = [0.14;0.14;0.14;0.15;0.15;0.15];
    m = [0.02;0.02;0.02;0.02;0.02;0.02];
    g_0 = 1e6*[300;300;300;240;240;240];
    q = [1 1 1 1 1 1;
         1 1 1 1 1 1;
         1 1 1 1 1 1;
         1 1 1 1 1 1;
         1 1 1 1 1 1;
         1 1 1 1 1 1;];
    end
elseif props.s_type == 10
    h_0 = 1e6*[2000;2000;2000;4153;4153;4153;ones(12,1)*4153*3];
    gamma_dot_tilde = [0.0023;0.0023;0.0023;0.0023;0.0023;0.0023;ones(12,1)*0.0023*3];
    g_tilde = 1e6*[504;504;504;504;504;504;ones(12,1)*504*3];
    r = [0.3;0.3;0.3;0.29;0.29;0.29;ones(12,1)*0.29];
    n = [0.14;0.14;0.14;0.15;0.15;0.15;ones(12,1)*0.15];
    m = [0.02;0.02;0.02;0.02;0.02;0.02;ones(12,1)*0.02];
    g_0 = 1e6*[300;300;300;240;240;240;ones(12,1)*240*3];
    q = ones(18,18);
else
    error('Unvalid slip type! Please try again!');
end
end

end