
D = 75000;
rho_c = 2900;
rho_m = 3750;
g = 3.72;
R = 3396000;

function [phi] = flexural_inf (n)

    phi = ( 1+ (D/((rho_m-rho_c)*g))*((2*n+1)/(2*R))^4) ^ -1;

end
