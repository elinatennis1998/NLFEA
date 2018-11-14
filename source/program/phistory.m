nhmax = 13*8*2; % maximum number of internal variables/element
switch iel
    case 1
        nha = (6+6+1)*8; % total number of time dependent hist vars/element
        nhb = nha; % total number of t_n+1 variables
        nhc = 0; % total number of time indep hist vars/element (e.g. tau)
    case 2
        nha = (6+6+1)*8; % total number of time dependent hist vars/element
        nhb = nha; % total number of t_n+1 variables
        nhc = 0; % total number of time indep hist vars/element (e.g. tau)
    case 3
        nha = (7)*6; % total number of time dependent hist vars/element
        nhb = nha; % total number of t_n+1 variables
        nhc = 0; % total number of time indep hist vars/element (e.g. tau)
    case 4
        nha = (7)*6; % total number of time dependent hist vars/element
        nhb = nha; % total number of t_n+1 variables
        nhc = 0; % total number of time indep hist vars/element (e.g. tau)
    case 5
        nha = (7)*6; % total number of time dependent hist vars/element
        nhb = nha; % total number of t_n+1 variables
        nhc = 0; % total number of time indep hist vars/element (e.g. tau)
    case 6
%         nha = (3+3+1)*4; % total number of time dependent hist vars/element
        nha = (12)*4; % total number of time dependent hist vars/element
        nhb = nha; % total number of t_n+1 variables
        nhc = 0; % total number of time indep hist vars/element (e.g. tau)
    case 7
        nha = (6+5+1)*8; % total number of time dependent hist vars/element
%         nha = (5+5+1)*8; % total number of time dependent hist vars/element
        nhb = nha; % total number of t_n+1 variables
        nhc = 0; % total number of time indep hist vars/element (e.g. tau)
    case 8
        nha = (3+3+1)*9; % total number of time dependent hist vars/element
        nhb = nha; % total number of t_n+1 variables
        nhc = 0; % total number of time indep hist vars/element (e.g. tau)
    case 9
        nha = (3+3+1)*4; % total number of time dependent hist vars/element
        nhb = nha; % total number of t_n+1 variables
        nhc = 3; % total number of time indep hist vars/element (e.g. tau)
    otherwise
        nha = 0;
        nhb = 0;
        nhc = 0;
end